# Single-cell tRNA quantification (Clover-Seq extension) — OCM

Extends [GDSC-Clover-Seq](https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq)
(bulk tRNA-seq) to 10x On-Chip Multiplexing (OCM) single-cell data. Standard
scRNA-seq pipelines (Cell Ranger) systematically undercount tRNA reads,
mostly because tRNA genes exist as near-identical copies at multiple genomic
loci and Cell Ranger's generic aligner treats that as low-confidence
multi-mapping and drops the reads. Clover-Seq's `choosemappings.py` handles
exactly this case for bulk data; this workflow adds cell barcode/UMI
awareness so the same logic works per-cell, per-OB-sample.

**Scope**: this Snakefile pipeline only recovers the **tRNA side** of the
data (isoacceptor/isodecoder counts per cell). Gene expression / protein_coding
counts come from the lab's own `cellranger multi` runs (already produced
upstream, outside this repo, at `.../outs/per_sample_outs/{sample_id}/`) —
there is no "pass 2 cellranger" step in the current pipeline. See
`grant_data/code/extract_proteincoding_cellranger.py` for how the
protein_coding side is pulled in for combined analyses.

## How it works

1. Each **well** (a physical 10x GEM lane, e.g. `236860-1`) is described by a
   `cellranger multi` CSV (`OCM1.csv` / `OCM2.csv`), which lists both the
   physical sequencing libraries (`[libraries]`, GEX + FBC fastq paths) and
   the OCM sample tags (`[samples]`, `OB1`-`OB4`, each with its own
   already-called-cells barcode list from the `cellranger multi` output).
   The Snakefile parses this CSV directly — nothing about OB1-OB4 identity is
   hardcoded here, it all comes from the CSV.
2. Selected libraries (`library_filter`: `"all"` = GEX+FBC pooled, matching
   Cell Ranger's own input scope; or a substring like `"FBC1"` to isolate
   just the dedicated small-RNA library) are concatenated per well.
3. A union cell whitelist is built per well from all 4 OBs' own
   `cellranger multi`-called barcodes (so we only ever count UMIs for cells
   Cell Ranger itself already called real).
4. Cell barcode + UMI are extracted from R1 (`umi_tools extract`, whitelist-
   restricted) and appended to each R2 read name — `{original_name}_{CB}_{UMI}`.
5. R2 is trimmed (`cutadapt`: 5' TSO, 3' adapter + poly-A), then aligned with
   `bowtie2 -k <maxMaps> --local --very-sensitive` against the same tRNA+genome
   index bulk Clover-Seq uses.
6. Clover-Seq's **unmodified** `choosemappings.py` resolves multi-mapping and
   picks the best tRNA (or non-tRNA genomic) alignment per read. It only
   groups alignments by read name and compares alignment scores — it never
   inspects read-name content — so the barcode/UMI suffix rides through
   completely unmodified.
7. `count_singlecell_tRNA.py` parses the resolved, sorted BAM, pulls the
   barcode/UMI back out of each read name, splits reads per OB (via each
   OB's own whitelist), and UMI-deduplicates within each (cell, tRNA gene)
   pair to build one count matrix per well per OB.

## Setup

1. **Conda env**: `env_config/clover-seq.yaml` has everything the pipeline
   itself needs (bowtie2, cutadapt, `umi_tools`, pysam, samtools). Note: the
   *named* env (`conda activate clover-seq`) is currently broken on this
   machine (scipy/umi_tools import errors) — Snakemake's own
   `--use-conda`-built env under `.snakemake/conda/` works and is what
   actually ran successfully; see `CONTEXT.md` §4 for the exact binary path
   and a known-working standalone alternative
   (`/dartfs-hpc/scratch/f00873z/envs/trna_env`).
2. **`config.yaml`** — driven entirely by per-well `cellranger multi` CSVs,
   see the `wells:` block. Key params:
   - `maxMaps` (bowtie2 `-k`, max alignments reported per read before giving
     up on a repetitive locus) — currently `50`. **Meeting note 2026-08-10:
     re-run with `maxMaps=100`, see below.**
   - `library_filter` / `library_exclude` / `results_dir` — these three at
     the bottom of the checked-in `config.yaml` are leftovers from the last
     manual run, **not** the actual settings any given job uses — every
     `sbatch/run_*.sbatch` script overrides all three via
     `snakemake --config library_filter=... results_dir=...` on the command
     line. Always check the specific `.sbatch` script being run, not just
     `config.yaml`, to know what a given job actually targets.
   - `trna_db` / `bt2_index` — prebuilt mm10 tRNA database, only change for a
     different genome.
   - `smrna_gtf` / `fbc_biotype_obs` — used by the (separate, non-Snakemake)
     full biotype-classification scripts under `grant_data/code/` and
     `code/` (import-locked scripts), not by the Snakefile rules themselves.

## Running it

```bash
snakemake --configfile config.yaml --cores <N> --use-conda \
  --config library_filter=FBC1 library_exclude=__none__ results_dir=fbc_only
```

In practice, every real run goes through one of the scripts in `sbatch/`
(each already `cd`s to the repo's absolute path and sets its own
`--config` overrides) — submit with `sbatch sbatch/<name>.sbatch`. See
`sbatch/run_fbc_align.sbatch` for the FBC-only tRNA recovery run and
`sbatch/run_gex_align.sbatch` for the GEX-library cross-check.

**`--config key=` with an empty value crashes** Snakemake's config parser —
use a placeholder like `library_exclude=__none__` instead of
`library_exclude=`.

## Pipeline steps (from the actual `Snakefile`)

| # | Rule | What it does |
|---|---|---|
| 0 | `pool_runs` | Concatenates R1/R2 fastqs across every selected library for a well into one pooled fastq pair. |
| 0b | `cell_whitelist` | Builds the union cell whitelist for a well from all 4 OBs' `cellranger multi` barcode lists (`-1` suffix stripped). |
| 1 | `barcode_extract` | `umi_tools extract` pulls the cell barcode + UMI off R1 (whitelist-restricted) and appends them to each R2 read name. |
| 2 | `sc_trimming` | `cutadapt` trims the 5' TSO and 3' adapter/poly-A off the barcode-tagged R2, same parameters as bulk Clover-Seq. |
| 3 | `sc_tRNA_bowtie2` | Aligns trimmed R2 to the tRNA+genome bowtie2 index (`db-tRNAgenome`), `-k maxMaps --local --very-sensitive`. |
| 4 | `sc_tRNA_choosemappings` | Runs Clover-Seq's unmodified `choosemappings.py` to resolve multi-mapping and pick the best alignment per read; sorts + indexes the result. |
| 5 | `sc_tRNA_count` | Per OB: splits the resolved BAM by that OB's whitelist, UMI-deduplicates per (cell, gene), writes `{well}/{OB}/matrix.mtx` + `barcodes.tsv` + `features.tsv` (isoacceptor level by default; `--level isodecoder --isodecoder-count {primary,unique}` for finer resolution, used by the `run_isodecoder_*.sbatch` scripts). |

## Known limitations / extension points

- UMI deduplication is exact-match only (no 1-mismatch UMI correction).
- Reads tied across genuinely different isoacceptors are assigned to
  whichever locus `choosemappings.py` designated primary, not split
  fractionally.
- Multi-lane fastqs must already be concatenated per library before this
  pipeline runs; Cell Ranger's own lane-handling isn't reproduced here.
- `count_biotype_umi_dedup.py`'s own internal tRNA bucket is known
  inaccurate (pools all isodecoders before UMI-clustering) — always use
  `count_singlecell_tRNA.py`'s isodecoder_unique output for tRNA numbers,
  see `CONTEXT.md` §1/§5.

## Repo layout

- **`code/`** — the Snakemake pipeline itself (`choosemappings.py`,
  `count_singlecell_tRNA.py`) plus `trnasequtils.py` (shared dependency) and
  a handful of other scripts that are import-locked to this directory (see
  `CONTEXT.md` §0 for the full list and why).
- **`code/diagnostics/`** — validation/sensitivity/comparison scripts for
  the grant analyses (not needed to reproduce a final reported number).
- **`code/archive/`** — superseded early script versions.
- **`grant_data/code/`** — scripts that generate the specific tables/figures
  sent to Esteban for the grant (biotype composition, tRNA UMAPs, isodecoder
  tables), kept separate from the Snakemake pipeline.
- **`grant_data/`** (root) — the actual output tables/figures/pptx sent to
  Esteban.
- **`sbatch/`** — all cluster-submission scripts for pipeline runs.
- **`diagnostics/`** (top-level, distinct from `code/diagnostics/`) — the
  original standalone toolkit used to diagnose *why* tRNA counts were low in
  Cell Ranger in the first place (`extract_trna_beds.py`,
  `run_gtf_overlap_analysis.sh`, `analyze_bam_tags.py`,
  `analyze_antisense.py`), plus later one-off investigation subfolders. This
  showed, for the Orellana `236860_3` run: ~41% of tRNA loci overlap another
  annotated gene (mostly intronic), only ~5% of reads at tRNA loci survive to
  a counted UMI, the dominant cause is raw genomic multi-mapping (`NH>1`)
  rather than gene-annotation conflicts (`GX>1`), and antisense strand
  orientation is not a contributing factor.

See `CONTEXT.md` for the full working history: every grant item's status,
every design decision and why, every error hit and its fix, and exact
commands to reproduce or resume any piece of this.

## Meeting notes

### 2026-08-10 — action items
- **Re-run the pipeline with `maxMaps=100`** (currently `50` in
  `config.yaml`) to check whether the tRNA recovery rate changes when
  bowtie2 is allowed to report more alignments per read before giving up on
  a repetitive locus.
- **Redo the tRNA UMAP for OB2 vs OB3** (`grant_data/code/umap_ob2_ob3_trna_local.R`)
  on the `maxMaps=100` output, to check whether the population-level
  dominant-isodecoder signal (grant item 4) is sensitive to this parameter.
