#!/usr/bin/env python3
"""
extract_proteincoding_cellranger.py — protein_coding-only cell x gene matrix,
subset from cellranger's own per-OB, per-well GEX+FBC1 combined count matrix.

Replaces count_protein_coding_umi_dedup.py (custom BAM-based re-counting) for
the boss's hypothesis test (does dominant tRNA isodecoder track protein_coding
expression clusters?). Reasoning, decided in discussion before writing this:

  - cellranger multi's config.csv (see e.g.
    /Volumes/GSR_Active/.../236860_1/outs/config.csv) declares BOTH the pure
    GEX1 fastqs AND the GEX1_FBC1 fastqs as feature_types "gene expression" --
    confirmed by qc_library_metrics.csv, where "Physical library ID GEX_1"
    read count == exact sum of all 4 fastq_id read counts. cellranger already
    pools GEX + FBC1 reads into one combined, splice-aware (STAR), per-cell
    gene count -- this is the same source already used for grant items 1/2
    (GEX biotype table + GEX UMAP), just not yet demuxed to protein_coding-only
    and OB2/OB3-only for this specific test.
  - Using this instead of a custom re-implementation removes an entire class
    of bug (we already found and fixed one real one: intron-embedded ncRNA
    getting mis-credited to its protein_coding host gene, because the custom
    script's biotype-priority check was missing). cellranger's own gene
    assignment is already splice-aware and doesn't need that heuristic at all.
  - tRNA identity (dominant_trna, joined in downstream by the R script) stays
    on the FBC1-only custom pipeline (fbc_only/03_tRNA_matrix), where the
    multi-mapping recovery that pipeline exists for actually matters. Each
    half of the test uses the library source it's actually good at:
    GEX(+FBC1 pooled) for protein_coding, FBC1-only custom pipeline for tRNA.

Usage:
    python extract_proteincoding_cellranger.py \\
        --matrix /Volumes/GSR_Active/Labs/Orellana/260608_10X/analysis/236860_1/outs/per_sample_outs/2/sample_filtered_feature_bc_matrix \\
        --gtf grant_data/ref/genes.reclass.gtf.gz \\
        --outdir fbc_only/03_protein_coding_matrix_cellranger/236860-1_OB2
"""
import argparse
import gzip
import os
import sys

import scipy.io
import scipy.sparse


def load_gene_biotypes(gtf_path):
    """gene_id -> gene_biotype, from the GTF's 'gene'-type rows."""
    opener = gzip.open if gtf_path.endswith(".gz") else open
    biotype = {}
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attr = cols[8]
            if 'gene_id "' not in attr or 'gene_biotype "' not in attr:
                continue
            gid = attr.split('gene_id "', 1)[1].split('"', 1)[0]
            bt = attr.split('gene_biotype "', 1)[1].split('"', 1)[0]
            biotype[gid] = bt
    return biotype


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True,
                     help="cellranger sample_filtered_feature_bc_matrix dir (per OB, per well)")
    ap.add_argument("--gtf", required=True, help="genes.reclass.gtf.gz (or genes.gtf.gz)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    print(f"loading biotypes from {a.gtf} ...", file=sys.stderr)
    gene_biotype = load_gene_biotypes(a.gtf)
    print(f"  {len(gene_biotype)} genes annotated", file=sys.stderr)

    print(f"loading cellranger matrix from {a.matrix} ...", file=sys.stderr)
    with gzip.open(os.path.join(a.matrix, "features.tsv.gz"), "rt") as f:
        features = [line.rstrip("\n").split("\t") for line in f]  # [gene_id, gene_name, feature_type]
    with gzip.open(os.path.join(a.matrix, "barcodes.tsv.gz"), "rt") as f:
        barcodes = [line.strip() for line in f]
    mat = scipy.io.mmread(os.path.join(a.matrix, "matrix.mtx.gz")).tocsr()
    if mat.shape[0] != len(features):
        sys.exit(f"matrix has {mat.shape[0]} rows but {len(features)} features")
    if mat.shape[1] != len(barcodes):
        sys.exit(f"matrix has {mat.shape[1]} columns but {len(barcodes)} barcodes")

    # protein_coding only, unannotated gene_ids dropped (not expected, but
    # fail loud rather than silently keep them if the reference ever drifts)
    unannotated = [gid for gid, _, _ in features if gid not in gene_biotype]
    if unannotated:
        print(f"  WARNING: {len(unannotated)} gene_ids in matrix not found in GTF "
              f"(dropped, e.g. {unannotated[:3]})", file=sys.stderr)
    keep_idx = [i for i, (gid, _, _) in enumerate(features)
                if gene_biotype.get(gid) == "protein_coding"]
    print(f"  {len(keep_idx)} / {len(features)} features are protein_coding", file=sys.stderr)

    sub = mat[keep_idx, :]
    kept_features = [features[i] for i in keep_idx]

    # strip cellranger's "-1" GEM-well suffix so barcodes match the raw
    # 16bp convention used by the FBC-only custom tRNA pipeline (same
    # convention split_bam_by_whitelist.py already relies on)
    stripped_barcodes = [bc.rsplit("-", 1)[0] for bc in barcodes]

    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "barcodes.tsv"), "w") as f:
        for b in stripped_barcodes:
            f.write(b + "\n")
    with open(os.path.join(a.outdir, "features.tsv"), "w") as f:
        for gid, gname, _ in kept_features:
            f.write(f"{gid}\t{gname}\tGene Expression\n")

    sub_coo = sub.tocoo()
    nnz = sub_coo.nnz
    with open(os.path.join(a.outdir, "matrix.mtx"), "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{sub.shape[0]} {sub.shape[1]} {nnz}\n")
        for r, c, v in zip(sub_coo.row, sub_coo.col, sub_coo.data):
            f.write(f"{r + 1} {c + 1} {int(v)}\n")

    print(f"{sub.shape[0]} genes, {sub.shape[1]} cells, {int(sub.sum()):,} total UMIs", file=sys.stderr)
    print(f"wrote {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
