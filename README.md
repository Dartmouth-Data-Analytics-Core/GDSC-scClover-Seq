# Single-cell tRNA quantification (Clover-Seq extension, OCM)

The Single-cell Clover-Seq extension recovers per-cell tRNA expression from
10x On-Chip Multiplexing (OCM) sequencing data, implemented through
Snakemake for use on the Dartmouth `dartfs-hpc` cluster. It extends
[GDSC-Clover-Seq](https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq)
(bulk tRNA-seq) by adding cell-barcode/UMI awareness to its
`choosemappings.py` multi-mapping resolution step, recovering tRNA reads
that a standard scRNA-seq pipeline (Cell Ranger) systematically drops,
because tRNA genes exist as near-identical copies at multiple genomic loci
and Cell Ranger's generic aligner treats that as low-confidence
multi-mapping. Software dependencies are installed per rule by Snakemake
from the conda environment file in `env_config/`.

This pipeline assumes `cellranger multi` has already been run for every
pool; it reads that run's own CSV and `outs/` directory directly (cell
calling, barcode whitelists, and the GEX/protein-coding side all come
from there) and never re-runs or replaces it.

## Documentation

- [Summary](#summary)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Optional Features](#optional-features)
- [Understanding the Outputs](#understanding-the-outputs)
- [Contact](#contact)

## Summary

Currently the pipeline performs the following, per configured well:

- Pooling of the selected sequencing library or libraries (the dedicated
  tRNA-enriched small-RNA library, the gene-expression library, or both
  combined) across sequencing runs
- Cell barcode + UMI extraction from R1, restricted to the union of
  barcodes `cellranger multi` already called as real cells
- Adapter/TSO trimming of R2 (`cutadapt`)
- Alignment to a combined tRNA+genome Bowtie2 reference (`--local
  --very-sensitive`), with reads that never align anywhere in that
  reference captured separately rather than discarded
- Multi-mapping resolution using Clover-Seq's unmodified `choosemappings.py`
- Per-cell, per-sample UMI-deduplicated tRNA counting, at both isoacceptor
  and isodecoder level
- A flat, externally-numbered sample view spanning every well, plus a
  manifest mapping each external sample back to its well/internal origin
- Per-sample biotype composition (tRNA, pre-tRNA, every Ensembl
  `gene_biotype`, then "other"), reusing Clover-Seq's unmodified
  `count_all_smRNA.py`
- Pull-through of `cellranger multi`'s own already-computed, splice-aware
  gene expression matrix, filtered to protein-coding genes, for combined
  tRNA + gene-expression analyses

## Quick Start

> **Important**
> `cellranger multi` must already have been run for every well before
> starting here. This pipeline does not run it for you.

1. Point `config.yaml`'s `wells:` block at your `cellranger multi` CSV(s)
   and their `outs/` directories.
2. Set `library_filter` / `results_dir` in `config.yaml`, or use the
   FBC1/`fbc_only` defaults already checked in (see
   [Configuration](#configuration)).
3. Submit one of the scripts in `sbatch/` to SLURM.

```bash
sbatch sbatch/run_fbc_align.sbatch
```

## Installation

Clone the repository:

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/GDSC-scClover-Seq
cd GDSC-scClover-Seq
```

Then get the shared Snakemake environment on your `PATH` (shared Snakemake environment):

```bash
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake
```

## Configuration

### 1. Well configuration

There is no separate sample sheet. Each **well** (e.g. `236860-1`) is
described by its own `cellranger multi` CSV, which Snakemake parses
directly. `config.yaml`'s `wells:` block just points at where that CSV and
its `outs/` directory live.

> **Important**
> What this pipeline calls a "well" is the final sequencing pool
> (`236860-1`, `236860-2`), not one of the OCM chip's own physical loading
> wells. The chip has 8 physical wells, but it multiplexes 4 samples into
> each final pool, so 2 pools come out of one chip run. Cell Ranger only
> accepts those 4 samples per pool being named `OB1`-`OB4`, which is why
> every `multi_csv` below only ever declares `OB1`-`OB4`, never more,
> regardless of how many pools exist. Any name beyond that (e.g. this
> pipeline's own external `OB5`-`OB8` numbering for a second pool) is
> assigned outside Cell Ranger.

| Field | Description |
|---|---|
| well name (key) | Short well identifier, e.g. `236860-1`. Any number of wells is supported. |
| `multi_csv` | Path to that well's `cellranger multi` CSV (e.g. `OCM1.csv`), listing both the sequencing libraries (`[libraries]`) and the OCM sample tags (`[samples]`, each with its own already-called-cells barcode list). |
| `multi_out` | Path to that well's `cellranger multi` `outs/` directory. |

An example `multi_csv` (well `236860-1`) is shown below. The Snakefile
only reads the `[libraries]` and `[samples]` sections; any other section
(here, `[gene-expression]`, Cell Ranger's own reference/create-bam
settings) is ignored. `[libraries]` can list the same library sequenced
across multiple runs (here, `run1` and `run2`); the pipeline pools them.

```csv
[gene-expression],,
reference,/dartfs-hpc/rc/lab/G/GSR_Active/scripts/refs/Mouse/mm10_smRNA,
create-bam,TRUE,
[libraries],,
fastq_id,fastqs,feature_types
236860-1_GEX1_DIL1,/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Orellana/260608_10X/run1,gene expression
236860-1_GEX1,/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Orellana/260608_10X/run2,gene expression
236860-1_GEX1_FBC1,/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Orellana/260608_10X/run1,gene expression
236860-1_GEX1_FBC1,/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Orellana/260608_10X/run2,gene expression
[samples],,
sample_id,ocm_barcode_ids,
1,OB1,
2,OB2,
3,OB3,
4,OB4,
```

### 2. Pipeline parameters

All settings live in `config.yaml`. At minimum, check the following
before a first run:

| Parameter | Description |
|---|---|
| `library_filter` | Which library or libraries to pool per well: `"FBC1"` (the dedicated tRNA-enriched library), `"all"` (GEX+FBC pooled, matching Cell Ranger's own input scope), or  `"GEX1"`. |
| `library_exclude` | A `fastq_id` token to exclude even if `library_filter` would otherwise include it. Use the placeholder `"__none__"` to exclude nothing (an empty string crashes Snakemake's config parser). |
| `results_dir` | Output root for a given run. Keep this distinct per `library_filter`/parameter combination so runs never overwrite each other (checked-in default: `fbc_only`). |
| `maxMaps` | Bowtie2 `-k`: max alignments reported per read before giving up on a repetitive locus (currently `100`). |
| `trna_db` / `bt2_index` | Prebuilt mm10 tRNA database and its Bowtie2 index. Only change for a different genome. |
| `bc_pattern` / `umi_separator` | Cell barcode + UMI structure (10x 3' v4: 16 bp CB + 12 bp UMI) and the separator used when appending them to read names. |
| `tso` / `adapter_1` / `minlength` | 5' TSO and 3' adapter/poly-A trimmed from R2 before Bowtie2 alignment, since Bowtie2 has no trimming of its own, unlike Cell Ranger's internal STAR step. `tso` is the 30 bp Template Switch Oligo 10x's own chemistry flanks every cDNA construct with (Cell Ranger trims this same sequence internally before its own alignment). |
| `smrna_gtf` / `reclass_gtf` / `fbc_biotype_obs` | Gene annotation used for biotype classification (rule `sc_biotype_by_sample`) and protein-coding extraction (rule `sc_protein_coding_matrix`). |

### 3. Job submission scripts

Every run goes through one of the scripts in `sbatch/`.

> **Important**
> Check the `--account`/`--partition` and email fields at the top of the
> `.sbatch` script you're submitting before running it on your own
> allocation.

### 4. Submitting the job

The pipeline's main mode, FBC1-only tRNA recovery, is what `config.yaml`'s
own defaults already match:

```yaml
library_filter: "FBC1"
library_exclude: "__none__"
results_dir: fbc_only
```

```bash
snakemake --cores 16 --use-conda --conda-frontend conda --rerun-incomplete \
  --config library_filter=FBC1 library_exclude=__none__ results_dir=fbc_only \
  --snakefile Snakefile
```

`library_filter` matching is exact-token, not substring, so `"FBC1"` never
accidentally matches `"FBC2"` or a `"GEX1"`-only row.

| `library_filter` value | What gets pooled | Typical `results_dir` |
|---|---|---|
| `"FBC1"` | Dedicated small-RNA library only | `fbc_only` |
| `"all"` | GEX + FBC pooled, matching Cell Ranger's own input | `combined_gex_fbc` |
| `"GEX1"` | Gene-expression library only | `gex_only` |

Or just submit the already-configured script directly:

```bash
sbatch sbatch/run_fbc_align.sbatch
```

## Optional Features

### Running just one rule

Ask Snakemake for a specific output file instead of the whole `rule all`;
it resolves the dependency graph for just that target, skipping any
upstream step whose output is already current.

```bash
# rebuild only the isodecoder-level tRNA matrix for one well/OB (rule sc_tRNA_count_isodecoder)
snakemake --cores 4 --use-conda --conda-frontend conda \
  --config results_dir=fbc_only \
  --snakefile Snakefile \
  fbc_only/03_tRNA_matrix/236860-1/OB2_isodecoder_unique/matrix.mtx

# rebuild only the per-sample biotype composition table (rule sc_biotype_by_sample)
snakemake --cores 8 --use-conda --conda-frontend conda \
  --config results_dir=fbc_only \
  --snakefile Snakefile \
  fbc_only/03_biotype_by_sample/biotype_by_sample_raw.txt
```


### Unmapped-read capture

Reads that never align anywhere in the tRNA+genome reference are written
to `{results_dir}/02_sc_alignment/{well}.unmapped.fastq.gz` alongside the
main alignment. This happens automatically as part of the normal
alignment rule; no extra target is needed to produce it.

> **Note**
> Reads that DID align somewhere (a tRNA locus or not, ambiguous or not)
> are already fully characterized elsewhere: `choosemappings.py` doesn't
> drop them, it just picks a best reference, and the per-sample biotype
> table (`sc_biotype_by_sample`) already classifies every one of them.
> `unmapped.fastq.gz` is specifically the reads that reached neither.

## Understanding the Outputs

| Path | Contents |
|---|---|
| `{results_dir}/03_tRNA_matrix_by_sample/{OB}[_isodecoder_{unique,primary}]/` | Flat, externally-numbered per-sample tRNA count matrices (`matrix.mtx` + `barcodes.tsv` + `features.tsv`, `Read10X()`-compatible). This is what most downstream analysis reads from. |
| `{results_dir}/sample_manifest.tsv` | Maps each external sample (`ext_ob`) back to its well and internal OB label. |
| `{results_dir}/03_biotype_by_sample/biotype_by_sample_{raw,norm}.txt` | Per-sample biotype composition. The `protein_coding` row here should not be read as real mRNA capture: Bowtie2 isn't splice-aware, so it aligns equally well to a gene's introns as to its exons, and the classifier counts a read as `protein_coding` anywhere in that gene's full genomic span, not just the exons. |
| `{results_dir}/03_protein_coding_matrix_by_sample/{OB}/` | Same flat per-sample layout as the tRNA matrix, for Cell Ranger's own splice-aware protein-coding counts. |
| `{results_dir}/02_sc_alignment/{well}.unmapped.fastq.gz` | Reads that never aligned anywhere in the tRNA+genome reference; see [Unmapped-read capture](#unmapped-read-capture). |

## Contact

**Contact and questions**: Please address questions to [DataAnalyticsCore@groups.dartmouth.edu](mailto:DataAnalyticsCore@groups.dartmouth.edu) or submit an issue in the GitHub repository.

This pipeline extends [GDSC-Clover-Seq](https://github.com/Dartmouth-Data-Analytics-Core/GDSC-clover-Seq),
itself adapted from the [tRAX tool](https://github.com/UCSC-LoweLab/tRAX)
(GPL v3.0). If you use this pipeline, please cite:
[Holmes AD, Howard JM, Chan PP, and Lowe TM.](https://www.biorxiv.org/content/10.1101/2022.07.02.498565v1)
