#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Single-cell tRNA quantification (Clover-Seq extension) - OCM
# Multi-well; driven by the cellranger multi CSVs.
#
# LIBRARY SELECTION (config: library_filter)
#   "all"   -> pool ALL libraries in the CSV [libraries] block (GEX + FBC).
#              This matches cellranger multi's own input, so the tRNA
#              comparison against Cell Ranger is apples-to-apples, and
#              maximizes tRNA molecules per cell.
#   "FBC1"  -> pool only libraries whose fastq_id contains that token
#              (the dedicated tRNA-enriched small-RNA library only).
#
# Outputs go under config: results_dir (default "results"), so FBC-only and
# GEX+FBC runs coexist without clobbering each other.
#
# OCM sample identity is in the 16 bp cell barcode (validated on 236860-2):
# cellranger multi's per-OB called-cell lists demultiplex for us, so we align
# ONCE per well and count once per OB in that well (any number of OBs).
# bulk Clover-Seq logic (bowtie2 -k + choosemappings.py) is UNCHANGED.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import os
import glob

configfile: "config.yaml"

RES = config.get("results_dir", "results").rstrip("/")
LIBRARY_FILTER = config.get("library_filter", "all")


def parse_multi_csv(csv_path, multi_out, library_filter):
    """Parse a cellranger multi CSV into {fbc_r1, fbc_r2, ob_barcodes}."""
    section = None
    libs = []          # (fastq_id, fastqs_dir)
    samples = []       # (sample_id, ob_id)
    with open(csv_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                section = line.split(",", 1)[0].strip().strip("[]").strip().lower()
                continue
            fields = [f.strip() for f in line.split(",")]
            if section == "libraries":
                if fields[0] == "fastq_id":
                    continue
                libs.append((fields[0], fields[1]))
            elif section == "samples":
                if fields[0] == "sample_id":
                    continue
                samples.append((fields[0], fields[1]))

    # select libraries: "all" -> every row; otherwise EXACT token match on
    # fastq_id split by "_" (e.g. "236860-1_GEX1_FBC1" -> tokens
    # ["236860-1", "GEX1", "FBC1"]). Exact match, not substring -- substring
    # matching let "FBC" match "FBC1" and "GEX" match "GEX1", a real bug
    # (see CONTEXT.md) that forced every FBC1-only job to also override
    # library_exclude to a neutral placeholder just to dodge it. Token
    # match makes that unnecessary for filter/exclude values that are
    # themselves exact tokens (GEX1, FBC1, ...).
    exclude = config.get("library_exclude", None)
    r1, r2 = [], []
    for fastq_id, d in libs:
        tokens = fastq_id.split("_")
        if library_filter != "all" and library_filter not in tokens:
            continue
        if exclude and exclude in tokens:
            continue
        r1 += sorted(glob.glob(os.path.join(d, f"{fastq_id}_S*_R1_001.fastq.gz")))
        r2 += sorted(glob.glob(os.path.join(d, f"{fastq_id}_S*_R2_001.fastq.gz")))

    if not r1 or len(r1) != len(r2):
        raise ValueError(
            f"{csv_path}: found {len(r1)} R1 and {len(r2)} R2 fastqs for "
            f"library_filter '{library_filter}'. Check the [libraries] dirs and fastq_ids.")

    ob_barcodes = {}
    for sid, ob in samples:
        ob_barcodes[ob] = os.path.join(
            multi_out, "per_sample_outs", sid,
            "sample_filtered_feature_bc_matrix", "barcodes.tsv.gz")

    return {"fbc_r1": r1, "fbc_r2": r2, "ob_barcodes": ob_barcodes}


WELLS = {
    well: parse_multi_csv(w["multi_csv"], w["multi_out"], LIBRARY_FILTER)
    for well, w in config["wells"].items()
}

wildcard_constraints:
    well    = r"236860-\d+",
    ob      = r"OB\d+",
    ext_ob  = r"OB\d+",
    mode    = r"unique|primary",
    variant = r"|_isodecoder_unique|_isodecoder_primary"

# "unique" (discards reads ambiguous across >1 isodecoder) is the
# convention every actual analysis in this project uses -- it's the fair,
# apples-to-apples analog of how Cell Ranger counts every other biotype
# (see README). "primary" (keeps ambiguous reads, tie-broken) is kept as
# an OPTIONAL target only -- the rule still supports it (request it
# explicitly, e.g. `snakemake .../OB2_isodecoder_primary/matrix.mtx`), it
# just isn't built by a plain `snakemake`/`rule all` run.
TARGETS = [
    f"{RES}/03_tRNA_matrix/{well}/{ob}/matrix.mtx"
    for well, wd in WELLS.items()
    for ob in wd["ob_barcodes"]
] + [
    f"{RES}/03_tRNA_matrix/{well}/{ob}_isodecoder_unique/matrix.mtx"
    for well, wd in WELLS.items()
    for ob in wd["ob_barcodes"]
]

#----- OB1-N external sample naming (per Esteban's request, 2026-09-02):
# each well's own internal OB labels get relabeled onto one flat,
# sequential OB1-N range spanning ALL wells, so anyone consuming
# {RES}/*_by_sample/ doesn't need to know about the well/OB nesting. A
# well that's a TECHNICAL duplicate of another (same physical samples run
# through OCM twice, not independent biological replicates -- confirmed
# with Elizabeth Sergison for this project's well2, see CONTEXT.md) just
# reuses the same internal OB labels as the well it duplicates; nothing
# here assumes any fixed number of wells or samples per well.
#
# Cumulative offset per well, driven by each well's ACTUAL sample count
# (WELLS[well]["ob_barcodes"], already parsed from that well's own
# multi_csv [samples] section) -- NOT a hardcoded block size. Well order
# is fixed by sorted well name (not dict/config order), so this is
# deterministic regardless of how config.yaml lists its wells. Previously
# `i * 4`, which silently broke as soon as a well's sample count differed
# from 4: fewer samples left a gap in the OB1-N numbering (harmless but
# confusing), more samples caused a genuine collision (EXT_OB_SOURCE is a
# plain dict with no collision check, so two different physical samples
# could end up sharing one external OB name and one clobbers the other).
WELL_OFFSET = {}
_cumulative_ob = 0
for _well in sorted(config["wells"].keys()):
    WELL_OFFSET[_well] = _cumulative_ob
    _cumulative_ob += len(WELLS[_well]["ob_barcodes"])

EXT_OB_SOURCE = {}   # "OB<n>" (n = 1..total samples across all wells) -> (well, internal "OB<n>")
for well, wd in WELLS.items():
    for ob in wd["ob_barcodes"]:
        ext_ob = f"OB{int(ob[2:]) + WELL_OFFSET[well]}"
        EXT_OB_SOURCE[ext_ob] = (well, ob)

BY_SAMPLE_TARGETS = [
    f"{RES}/03_tRNA_matrix_by_sample/{ext_ob}{variant}/matrix.mtx"
    for ext_ob in EXT_OB_SOURCE
    for variant in ["", "_isodecoder_unique"]   # "_isodecoder_primary" optional, see above
]

BIOTYPE_TARGETS = [
    f"{RES}/03_biotype_by_sample/biotype_by_sample_raw.txt",
    f"{RES}/03_biotype_by_sample/biotype_by_sample_norm.txt",
]

PC_TARGETS = [
    f"{RES}/03_protein_coding_matrix_cellranger/{well}_{ob}/matrix.mtx"
    for well, wd in WELLS.items()
    for ob in wd["ob_barcodes"]
] + [
    f"{RES}/03_protein_coding_matrix_by_sample/{ext_ob}/matrix.mtx"
    for ext_ob in EXT_OB_SOURCE
]

MANIFEST_TARGETS = [f"{RES}/sample_manifest.tsv"]

TARGETS = TARGETS + BY_SAMPLE_TARGETS + BIOTYPE_TARGETS + PC_TARGETS + MANIFEST_TARGETS


rule all:
    input:
        TARGETS


#----- Rule 0: pool the selected libraries across sequencing runs (per well).
rule pool_runs:
    input:
        r1 = lambda w: WELLS[w.well]["fbc_r1"],
        r2 = lambda w: WELLS[w.well]["fbc_r2"]
    output:
        r1 = temp(f"{RES}/00_pool/{{well}}.R1.fastq.gz"),
        r2 = temp(f"{RES}/00_pool/{{well}}.R2.fastq.gz")
    message: "Pooling reads across runs/libraries: {wildcards.well}"
    threads: 1
    resources: maxtime="4:00:00", mem_mb="4gb"
    shell: """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
    """


#----- Rule 0b: union cell whitelist for the well (one list per OB), "-1" stripped.
rule cell_whitelist:
    input:
        bcs = lambda w: list(WELLS[w.well]["ob_barcodes"].values())
    output:
        wl = f"{RES}/00_pool/{{well}}.cell_whitelist.txt"
    message: "Building union cell whitelist: {wildcards.well}"
    threads: 1
    resources: maxtime="0:30:00", mem_mb="4gb"
    shell: """
        zcat {input.bcs} | sed 's/-1$//' | sort -u > {output.wl}
    """


#----- Rule 1: extract CB + UMI from R1, embed in R2 read names.
rule barcode_extract:
    input:
        r1 = f"{RES}/00_pool/{{well}}.R1.fastq.gz",
        r2 = f"{RES}/00_pool/{{well}}.R2.fastq.gz",
        wl = f"{RES}/00_pool/{{well}}.cell_whitelist.txt"
    output:
        r2_tagged  = f"{RES}/01_barcode_tagged/{{well}}.R2.tagged.fastq.gz",
        r1_discard = temp(f"{RES}/01_barcode_tagged/{{well}}.R1.discard.fastq")
    log: f"{RES}/01_barcode_tagged/logs/{{well}}.umi_tools_extract.log"
    message: "Extracting cell barcode + UMI: {wildcards.well}"
    conda: "env_config/clover-seq.yaml"
    threads: 2
    resources: maxtime="12:00:00", mem_mb="16gb"
    params:
        bc_pattern = config["bc_pattern"]
    shell: """
        umi_tools extract \
            -I {input.r1} \
            --read2-in={input.r2} \
            --bc-pattern={params.bc_pattern} \
            --whitelist={input.wl} \
            -S {output.r1_discard} \
            --read2-out={output.r2_tagged} \
            2> {log}
    """


#----- Rule 2: trim tagged R2, Cell Ranger style (TSO 5', poly-A + adapter 3').
# On long GEX reads the TSO simply won't be present at the 5' end, so -g is a
# no-op there; it only fires on short (tRNA-length) molecules. Safe for both.
rule sc_trimming:
    input:
        r2_tagged = f"{RES}/01_barcode_tagged/{{well}}.R2.tagged.fastq.gz"
    output:
        trimmed = f"{RES}/02_sc_alignment/{{well}}.R2.trim.fastq.gz"
    log: f"{RES}/02_sc_alignment/logs/{{well}}.cutadapt.log"
    message: "Trimming TSO (5') + poly-A / adapter (3'): {wildcards.well}"
    conda: "env_config/clover-seq.yaml"
    threads: 8
    resources: maxtime="4:00:00", mem_mb="60gb"
    params:
        tso       = config["tso"],
        adapter_1 = config["adapter_1"],
        minlength = config["minlength"]
    shell: """
        cutadapt \
            -g {params.tso} \
            -a {params.adapter_1} \
            --poly-a \
            -m {params.minlength} \
            -j {threads} \
            -o {output.trimmed} \
            {input.r2_tagged} > {log}
    """


#----- Rule 3: bowtie2 to the tRNA+genome index. Core alignment logic
# UNCHANGED from bulk -- the only addition (2026-09-09) is `--un-gz`,
# bowtie2's own native flag to ALSO write out every read that failed to
# align anywhere in the tRNA+genome index (not even a non-tRNA genomic
# locus), alongside the main `--no-unal`-filtered BAM. Reads that DID
# align somewhere -- tRNA locus or not, ambiguous or not -- are already
# fully characterized downstream: choosemappings.py (rule 4) doesn't drop
# them, it just picks a best reference, and sc_biotype_by_sample (rule 9)
# already classifies every read in the resulting per-sample BAM by
# biotype (that's what most of biotype_by_sample_raw.txt's non-tRNA rows
# already are). This flag exists to capture the one genuinely
# uncharacterized bucket: reads that never matched the combined
# tRNA+genome reference at all. See CONTEXT.md for the fuller reasoning
# (why a separate Cell-Ranger-pass-2-style script was NOT the right tool
# for this question).
rule sc_tRNA_bowtie2:
    input:
        trimmed = f"{RES}/02_sc_alignment/{{well}}.R2.trim.fastq.gz"
    output:
        rawBam   = temp(f"{RES}/02_sc_alignment/{{well}}.raw.bam"),
        unmapped = f"{RES}/02_sc_alignment/{{well}}.unmapped.fastq.gz"
    log: f"{RES}/02_sc_alignment/logs/{{well}}.bowtie2.log"
    message: "Bowtie2 alignment: {wildcards.well}"
    conda: "env_config/clover-seq.yaml"
    threads: 16
    resources: maxtime="24:00:00", mem_mb="60gb"
    params:
        bt2_index = config["bt2_index"],
        maxMaps   = config["maxMaps"],
        nPenalty  = config["nPenalty"]
    shell: """
        bowtie2 \
            --local \
            -x {params.bt2_index} \
            -U {input.trimmed} \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --ignore-quals \
            --no-unal \
            --un-gz {output.unmapped} \
            --reorder \
            -p {threads} \
            2> {log} \
        | samtools view -b -@ 4 -o {output.rawBam} -
    """


#----- Rule 4: choosemappings.py (unmodified).
rule sc_tRNA_choosemappings:
    input:
        rawBam = f"{RES}/02_sc_alignment/{{well}}.raw.bam"
    output:
        srtBam = f"{RES}/02_sc_alignment/{{well}}.srt.bam",
        bai    = f"{RES}/02_sc_alignment/{{well}}.srt.bam.bai"
    log: f"{RES}/02_sc_alignment/logs/{{well}}.choosemappings.log"
    message: "Selecting best tRNA mappings: {wildcards.well}"
    conda: "env_config/clover-seq.yaml"
    threads: 10
    resources: maxtime="12:00:00", mem_mb="60gb"
    params:
        tRNA_db = config["trna_db"]
    shell: """
        python code/choosemappings.py {params.tRNA_db}/db-trnatable.txt \
            --input {input.rawBam} \
            --progname sc-clover-seq \
            --fqname {wildcards.well} \
            --expname {wildcards.well} \
            --minnontrnasize 20 \
            2> {log} \
        | samtools sort -@ {threads} -o {output.srtBam} -

        samtools index {output.srtBam}
    """


#----- Rule 5: per-cell UMI dedup -> tRNA isoacceptor x cell matrix, per OB.
rule sc_tRNA_count:
    input:
        bam = f"{RES}/02_sc_alignment/{{well}}.srt.bam",
        wl  = lambda w: WELLS[w.well]["ob_barcodes"][w.ob]
    output:
        matrix   = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}/matrix.mtx",
        barcodes = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}/barcodes.tsv",
        features = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}/features.tsv"
    log: f"{RES}/03_tRNA_matrix/logs/{{well}}.{{ob}}.count.log"
    message: "Counting single-cell tRNA UMIs: {wildcards.well} {wildcards.ob}"
    conda: "env_config/clover-seq.yaml"
    threads: 1
    resources: maxtime="2:00:00", mem_mb="32gb"
    params:
        tRNA_db       = config["trna_db"],
        umi_separator = config["umi_separator"]
    shell: """
        python code/count_singlecell_tRNA.py \
            --input {input.bam} \
            --trnatable {params.tRNA_db}/db-trnatable.txt \
            --whitelist {input.wl} \
            --umi-separator {params.umi_separator} \
            --outdir $(dirname {output.matrix}) \
            2> {log}
    """


#----- Rule 5b: same as rule 5, but at isodecoder level (individual tRNA
# gene, not amino-acid+anticodon) -- this is the resolution every OB2/OB3
# analysis actually consumes, previously only produced by standalone
# run_isodecoder_*.sbatch scripts outside any Snakemake tracking. `mode`
# selects count_singlecell_tRNA.py's --isodecoder-count: "unique" discards
# reads ambiguous across >1 isodecoder, "primary" keeps them tie-broken.
rule sc_tRNA_count_isodecoder:
    input:
        bam = f"{RES}/02_sc_alignment/{{well}}.srt.bam",
        wl  = lambda w: WELLS[w.well]["ob_barcodes"][w.ob]
    output:
        matrix   = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}_isodecoder_{{mode}}/matrix.mtx",
        barcodes = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}_isodecoder_{{mode}}/barcodes.tsv",
        features = f"{RES}/03_tRNA_matrix/{{well}}/{{ob}}_isodecoder_{{mode}}/features.tsv"
    log: f"{RES}/03_tRNA_matrix/logs/{{well}}.{{ob}}.isodecoder_{{mode}}.count.log"
    message: "Counting single-cell tRNA UMIs (isodecoder, {wildcards.mode}): {wildcards.well} {wildcards.ob}"
    conda: "env_config/clover-seq.yaml"
    threads: 1
    resources: maxtime="2:00:00", mem_mb="32gb"
    params:
        tRNA_db       = config["trna_db"],
        umi_separator = config["umi_separator"]
    shell: """
        python code/count_singlecell_tRNA.py \
            --input {input.bam} \
            --trnatable {params.tRNA_db}/db-trnatable.txt \
            --whitelist {input.wl} \
            --umi-separator {params.umi_separator} \
            --level isodecoder \
            --isodecoder-count {wildcards.mode} \
            --outdir $(dirname {output.matrix}) \
            2> {log}
    """


#----- Rule 6: OB1-N external-naming view (see EXT_OB_SOURCE above).
# Symlinks ONLY -- never a second copy of the data. 03_tRNA_matrix/{well}/
# stays the single source of truth; every existing script that already
# reads it (UMAP scripts, isodecoder tables, the objects already shared
# with Esteban) is completely untouched by this rule.
rule sc_tRNA_by_sample_view:
    input:
        matrix   = lambda w: f"{RES}/03_tRNA_matrix/{EXT_OB_SOURCE[w.ext_ob][0]}/{EXT_OB_SOURCE[w.ext_ob][1]}{w.variant}/matrix.mtx",
        barcodes = lambda w: f"{RES}/03_tRNA_matrix/{EXT_OB_SOURCE[w.ext_ob][0]}/{EXT_OB_SOURCE[w.ext_ob][1]}{w.variant}/barcodes.tsv",
        features = lambda w: f"{RES}/03_tRNA_matrix/{EXT_OB_SOURCE[w.ext_ob][0]}/{EXT_OB_SOURCE[w.ext_ob][1]}{w.variant}/features.tsv"
    output:
        matrix   = f"{RES}/03_tRNA_matrix_by_sample/{{ext_ob}}{{variant}}/matrix.mtx",
        barcodes = f"{RES}/03_tRNA_matrix_by_sample/{{ext_ob}}{{variant}}/barcodes.tsv",
        features = f"{RES}/03_tRNA_matrix_by_sample/{{ext_ob}}{{variant}}/features.tsv"
    message: "OB1-N sample view: {wildcards.ext_ob}{wildcards.variant}"
    threads: 1
    resources: maxtime="0:10:00", mem_mb="1gb"
    run:
        import os
        for src, dst in zip(
            [input.matrix, input.barcodes, input.features],
            [output.matrix, output.barcodes, output.features],
        ):
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)


#----- Rule 7: split each well's resolved BAM into one BAM per external
# sample (OB1-N, see EXT_OB_SOURCE), so count_all_smRNA.py (rule 9) --
# which has no --whitelist option, unlike count_singlecell_tRNA.py -- can
# classify each sample separately. code/split_bam_by_whitelist.py
# already supports passing just one --ob label per call.
rule sc_biotype_split_bam:
    input:
        bam = lambda w: f"{RES}/02_sc_alignment/{EXT_OB_SOURCE[w.ext_ob][0]}.srt.bam",
        wl  = lambda w: WELLS[EXT_OB_SOURCE[w.ext_ob][0]]["ob_barcodes"][EXT_OB_SOURCE[w.ext_ob][1]]
    output:
        bam = f"{RES}/02_sc_alignment/by_sample/biosample_{{ext_ob}}.srt.bam",
        bai = f"{RES}/02_sc_alignment/by_sample/biosample_{{ext_ob}}.srt.bam.bai"
    log: f"{RES}/02_sc_alignment/logs/biosample_{{ext_ob}}.split.log"
    message: "Splitting BAM for external sample: {wildcards.ext_ob}"
    conda: "env_config/clover-seq.yaml"
    threads: 1
    resources: maxtime="2:00:00", mem_mb="16gb"
    shell: """
        python code/split_bam_by_whitelist.py \
            --input {input.bam} \
            --ob {wildcards.ext_ob}={input.wl} \
            --outdir {RES}/02_sc_alignment/by_sample \
            --prefix biosample \
            2> {log}
    """


#----- Rule 8: samplefile for count_all_smRNA.py, one row per external
# sample (OB1-N) -- each gets its own `replicate` value so nothing gets
# pooled across wells (unlike run_fbc_biotype.sbatch's OB1/OB2 tables,
# which deliberately DO pool well1+well2 under one replicate value).
rule sc_biotype_samplefile:
    output:
        tsv = f"{RES}/03_biotype_by_sample/samples.txt"
    message: "Writing count_all_smRNA.py samplefile for OB1-N"
    threads: 1
    resources: maxtime="0:10:00", mem_mb="1gb"
    run:
        import os
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        # NO header line -- trnasequtils.samplefile's parser (code/trnasequtils.py)
        # treats every line as data, doesn't skip a header row. A header
        # here gets parsed as a 9th bogus sample (sample_id="sample_id",
        # bamdir="bamdir"), causing "No such file or directory:
        # bamdir/sample_id.srt.bam". Matches the format of the
        # already-working fbc_results/fbc_ob_samples.txt, which also has
        # no header.
        with open(output.tsv, "w") as fh:
            for ext_ob in EXT_OB_SOURCE:
                fh.write(f"biosample_{ext_ob}\t{ext_ob}\t{RES}/02_sc_alignment/by_sample\n")


#----- Rule 9: full biotype classification per sample (OB1-N), general
# pipeline output -- NOT tied to any one comparison (Esteban's request,
# 2026-09-02). Reuses code/count_all_smRNA.py UNMODIFIED, the original
# bulk Clover-Seq classifier (mature tRNA -> pre-tRNA loci -> Ensembl GTF
# -> other, see its own docstring) -- categories are not hardcoded here,
# they come from whatever gene_biotype values appear in smrna_gtf for
# each sample's reads. Deliberately does NOT report protein_coding: this
# rule's input BAM is bowtie2 (non-splice-aware) + choosemappings.py,
# which has a structural bias toward resolving ambiguous reads as tRNA
# over anything else, so a protein_coding number out of THIS classifier
# would be systematically biased low, not a fair measurement -- that
# question is already answered better by the separate, splice-aware
# cellranger-native path (code/extract_proteincoding_cellranger.py).
# See CONTEXT.md.
rule sc_biotype_by_sample:
    input:
        samplefile = f"{RES}/03_biotype_by_sample/samples.txt",
        bams = [f"{RES}/02_sc_alignment/by_sample/biosample_{ext_ob}.srt.bam" for ext_ob in EXT_OB_SOURCE]
    output:
        norm = f"{RES}/03_biotype_by_sample/biotype_by_sample_norm.txt",
        raw  = f"{RES}/03_biotype_by_sample/biotype_by_sample_raw.txt"
    log: f"{RES}/03_biotype_by_sample/logs/count_all_smRNA.log"
    message: "Classifying biotype composition per sample (OB1-N)"
    conda: "env_config/clover-seq.yaml"
    threads: 8
    resources: maxtime="6:00:00", mem_mb="32gb"
    params:
        ensemblgtf  = config["smrna_gtf"],
        maturetrnas = f"{config['trna_db']}/db-maturetRNAs.bed",
        trnaloci    = f"{config['trna_db']}/db-trnaloci.bed",
        trnatable   = f"{config['trna_db']}/db-trnatable.txt"
    shell: """
        python code/count_all_smRNA.py \
            --samplefile {input.samplefile} \
            --ensemblgtf {params.ensemblgtf} \
            --maturetrnas {params.maturetrnas} \
            --trnaloci {params.trnaloci} \
            --trnatable {params.trnatable} \
            --countfile {output.norm} \
            --realcountfile {output.raw} \
            --cores {threads} \
            2> {log}
    """


#----- Rule 10: protein_coding matrix, per well/OB -- brings
# code/extract_proteincoding_cellranger.py into the tracked
# pipeline (previously a standalone manual script). This
# is a pure filter, not a recount: takes cellranger multi's own
# already-computed sample_filtered_feature_bc_matrix (GEX+FBC1-pooled,
# splice-aware STAR counts, already-called cells, already UMI-deduplicated
# -- nothing here re-aligns or re-counts anything), keeps only genes tagged
# protein_coding in smrna_gtf, strips cellranger's "-1" barcode suffix to
# match this pipeline's own raw-barcode convention. Deliberately filtered,
# not the raw all-biotypes matrix: mixing rRNA/tRNA capture differences
# into one shared normalization would reintroduce a compositional confound
# already established and fixed elsewhere in this project -- see CONTEXT.md.
rule sc_protein_coding_matrix:
    input:
        matrix_dir = lambda w: os.path.dirname(WELLS[w.well]["ob_barcodes"][w.ob])
    output:
        matrix   = f"{RES}/03_protein_coding_matrix_cellranger/{{well}}_{{ob}}/matrix.mtx",
        barcodes = f"{RES}/03_protein_coding_matrix_cellranger/{{well}}_{{ob}}/barcodes.tsv",
        features = f"{RES}/03_protein_coding_matrix_cellranger/{{well}}_{{ob}}/features.tsv"
    log: f"{RES}/03_protein_coding_matrix_cellranger/logs/{{well}}.{{ob}}.extract.log"
    message: "Extracting protein_coding matrix (cellranger-native): {wildcards.well} {wildcards.ob}"
    # Uses the shared clover-seq.yaml (scipy added there 2026-09-09) --
    # previously an isolated env_config/protein_coding.yaml specifically to
    # avoid marking every other clover-seq.yaml rule "environment changed"
    # (see CONTEXT.md for the full-rebuild scare that caused). Merged back
    # now because the maxMaps=100 rerun already requires a from-scratch
    # build under its own results_dir (fbc_only_mm100), so that risk no
    # longer applies for this run -- env_config/protein_coding.yaml is kept
    # on disk but unused/superseded. Re-isolate it again before ever
    # touching an already-built results_dir (e.g. fbc_only/) with this env.
    conda: "env_config/clover-seq.yaml"
    threads: 1
    resources: maxtime="1:00:00", mem_mb="16gb"
    params:
        gtf = config["smrna_gtf"]
    shell: """
        python code/extract_proteincoding_cellranger.py \
            --matrix {input.matrix_dir} \
            --gtf {params.gtf} \
            --outdir $(dirname {output.matrix}) \
            2> {log}
    """


#----- Rule 11: OB1-N external-naming view for the protein_coding
# matrix, mirroring rule 6 (sc_tRNA_by_sample_view) exactly -- same
# EXT_OB_SOURCE mapping, symlinks only. Goal: Esteban (or anyone) can
# Read10X() a clean OB1-N-numbered matrix for EITHER assay without
# knowing the well1/well2 internal layout, then build their own Seurat
# object (see grant_data/code/umap_ob2_ob3_trna_local.R for a worked
# example that also adds the tRNA assay on the same object).
rule sc_protein_coding_by_sample_view:
    input:
        matrix   = lambda w: f"{RES}/03_protein_coding_matrix_cellranger/{EXT_OB_SOURCE[w.ext_ob][0]}_{EXT_OB_SOURCE[w.ext_ob][1]}/matrix.mtx",
        barcodes = lambda w: f"{RES}/03_protein_coding_matrix_cellranger/{EXT_OB_SOURCE[w.ext_ob][0]}_{EXT_OB_SOURCE[w.ext_ob][1]}/barcodes.tsv",
        features = lambda w: f"{RES}/03_protein_coding_matrix_cellranger/{EXT_OB_SOURCE[w.ext_ob][0]}_{EXT_OB_SOURCE[w.ext_ob][1]}/features.tsv"
    output:
        matrix   = f"{RES}/03_protein_coding_matrix_by_sample/{{ext_ob}}/matrix.mtx",
        barcodes = f"{RES}/03_protein_coding_matrix_by_sample/{{ext_ob}}/barcodes.tsv",
        features = f"{RES}/03_protein_coding_matrix_by_sample/{{ext_ob}}/features.tsv"
    message: "OB1-N sample view (protein_coding): {wildcards.ext_ob}"
    threads: 1
    resources: maxtime="0:10:00", mem_mb="1gb"
    run:
        import os
        for src, dst in zip(
            [input.matrix, input.barcodes, input.features],
            [output.matrix, output.barcodes, output.features],
        ):
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)


#----- Rule 12: sample manifest -- writes EXT_OB_SOURCE (the external
# OB-name -> well/internal-OB mapping computed above) to disk, so
# downstream scripts (e.g. grant_data/code/build_seurat_all_samples_local.R)
# read the mapping instead of re-deriving or hardcoding it themselves.
# Single source of truth stays here in the Snakefile: if the number of
# wells or samples per well ever changes, this file is the only thing
# that needs to be regenerated, not every consumer script.
rule sample_manifest:
    output:
        tsv = f"{RES}/sample_manifest.tsv"
    message: "Writing OB -> well/internal-OB sample manifest"
    threads: 1
    resources: maxtime="0:10:00", mem_mb="1gb"
    run:
        ordered_ext_obs = sorted(EXT_OB_SOURCE, key=lambda x: int(x[2:]))
        with open(output.tsv, "w") as fh:
            fh.write("ext_ob\twell\tinternal_ob\n")
            for ext_ob in ordered_ext_obs:
                well, ob = EXT_OB_SOURCE[ext_ob]
                fh.write(f"{ext_ob}\t{well}\t{ob}\n")
