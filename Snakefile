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
# ONCE per well and count 4x per well. bulk Clover-Seq logic (bowtie2 -k +
# choosemappings.py) is UNCHANGED.
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

    # select libraries: "all" -> every row; otherwise substring match on fastq_id
    exclude = config.get("library_exclude", None)
    r1, r2 = [], []
    for fastq_id, d in libs:
        if library_filter != "all" and library_filter not in fastq_id:
            continue
        if exclude and exclude in fastq_id:
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
    well = r"236860-\d+",
    ob   = r"OB\d+"

TARGETS = [
    f"{RES}/03_tRNA_matrix/{well}/{ob}/matrix.mtx"
    for well, wd in WELLS.items()
    for ob in wd["ob_barcodes"]
]


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
    resources: cpus="1", maxtime="4:00:00", mem_mb="4gb"
    shell: """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
    """


#----- Rule 0b: union cell whitelist for the well (4 OB lists), "-1" stripped.
rule cell_whitelist:
    input:
        bcs = lambda w: list(WELLS[w.well]["ob_barcodes"].values())
    output:
        wl = f"{RES}/00_pool/{{well}}.cell_whitelist.txt"
    message: "Building union cell whitelist: {wildcards.well}"
    resources: cpus="1", maxtime="0:30:00", mem_mb="4gb"
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
    resources: cpus="2", maxtime="12:00:00", mem_mb="16gb"
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
    resources: cpus="8", maxtime="4:00:00", mem_mb="60gb"
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
            -j {resources.cpus} \
            -o {output.trimmed} \
            {input.r2_tagged} > {log}
    """


#----- Rule 3: bowtie2 to the tRNA+genome index. UNCHANGED from bulk.
rule sc_tRNA_bowtie2:
    input:
        trimmed = f"{RES}/02_sc_alignment/{{well}}.R2.trim.fastq.gz"
    output:
        rawBam = temp(f"{RES}/02_sc_alignment/{{well}}.raw.bam")
    log: f"{RES}/02_sc_alignment/logs/{{well}}.bowtie2.log"
    message: "Bowtie2 alignment: {wildcards.well}"
    conda: "env_config/clover-seq.yaml"
    resources: cpus="16", maxtime="24:00:00", mem_mb="60gb"
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
            --reorder \
            -p {resources.cpus} \
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
    resources: cpus="10", maxtime="12:00:00", mem_mb="60gb"
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
        | samtools sort -@ {resources.cpus} -o {output.srtBam} -

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
    resources: cpus="1", maxtime="2:00:00", mem_mb="32gb"
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