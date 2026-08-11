#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clover-Seq single-cell extension
# Builds a cell x tRNA-isoacceptor UMI count matrix from a
# tRNA_choosemappings-resolved BAM whose read names carry a
# cell barcode + UMI suffix (umi_tools extract convention:
# "{original_read_name}_{CELL_BARCODE}_{UMI}").
#
# Feature resolution is selectable (--level):
#   isoacceptor (default) -> amino acid + anticodon (e.g. Gly-GCC).
#       tRNA gene copies are near-identical across loci, so a read
#       that multi-maps across copies of the SAME isoacceptor is not
#       a real ambiguity; collapsing here recovers those multimappers.
#   isodecoder            -> the individual tRNA transcript (e.g.
#       tRNA-Gly-GCC-1). This RE-INTRODUCES the cross-copy ambiguity
#       that isoacceptor avoids, so --isodecoder-count controls it:
#         primary (default) -> count the primary alignment's isodecoder
#             for every read (choosemappings' tie-break; some noise).
#         unique            -> count a read only if it maps to exactly
#             one isodecoder (YR tag == 1); conservative, drops the
#             multi-isodecoder reads (mirrors the bulk unique counts).
# Only the PRIMARY alignment per read is used to name its feature.
#
# UMI counting mirrors Cell Ranger's two steps, at isoacceptor
# resolution:
#   (1) UMI error correction: within each (cell, isoacceptor), UMIs
#       one Hamming apart are collapsed (directional; umi_tools
#       UMIClusterer, the same algorithm Cell Ranger uses).
#   (2) collision resolution: if one (cell, UMI) molecule is seen on
#       more than one isoacceptor, it is assigned to the isoacceptor
#       with the most supporting reads; on a tie it is discarded.
# Each surviving (cell, UMI, isoacceptor) is one molecule = one count.
# This replaces the previous len(set(UMIs)) per (cell, isoacceptor),
# which skipped correction (1) and double-counted molecules across
# isoacceptors (2).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import argparse
import gzip
import os
import sys
from collections import defaultdict, Counter

import pysam
from trnasequtils import transcriptfile

# umi_tools' reference directional clusterer (same method as Cell Ranger).
_clusterer = None
try:
    try:
        from umi_tools import UMIClusterer
    except ImportError:
        from umi_tools.network import UMIClusterer
    _clusterer = UMIClusterer(cluster_method="directional")
except Exception:
    _clusterer = None


def parse_barcode_umi(qname, separator="_"):
    parts = qname.rsplit(separator, 2)
    if len(parts) != 3:
        return None, None
    _, cb, umi = parts
    return cb, umi


def load_whitelist(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        # strip GEM-well suffix (e.g. "-1") so raw sequenced barcodes match
        return {line.strip().split("-")[0] for line in f if line.strip()}


def cluster_umis(umi_read_counts):
    """Collapse UMIs 1 Hamming apart (directional). Input/return: {umi: reads}.

    Returns {representative_umi: summed_reads}. Falls back to the raw dict if
    umi_tools is unavailable or there is nothing to collapse.
    """
    if _clusterer is None or len(umi_read_counts) <= 1:
        return dict(umi_read_counts)
    bc = {u.encode(): c for u, c in umi_read_counts.items()}
    try:
        clusters = _clusterer(bc, 1)          # umi_tools >= 1.0: (counts, threshold)
    except TypeError:
        clusters = _clusterer(bc, threshold=1)
    out = {}
    for clust in clusters:
        rep = clust[0].decode()
        out[rep] = sum(bc[u] for u in clust)
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Build a cell x tRNA-isoacceptor UMI count matrix (Cell Ranger-style UMI counting).")
    parser.add_argument("--input", required=True,
                        help="Coordinate-sorted BAM from the sc_tRNA_choosemappings rule")
    parser.add_argument("--trnatable", required=True,
                        help="db-trnatable.txt from the tRNA reference database")
    parser.add_argument("--whitelist",
                        help="Optional barcodes.tsv(.gz) of called cells to restrict counting to real cells")
    parser.add_argument("--umi-separator", default="_",
                        help="Separator umi_tools used to embed CB/UMI in read names (default: _)")
    parser.add_argument("--level", choices=["isoacceptor", "isodecoder"], default="isoacceptor",
                        help="Feature resolution: isoacceptor (amino-anticodon, default) or "
                             "isodecoder (individual tRNA transcript)")
    parser.add_argument("--isodecoder-count", choices=["primary", "unique"], default="primary",
                        help="Only for --level isodecoder: 'primary' counts the primary alignment's "
                             "isodecoder for every read; 'unique' counts a read only if it maps to a "
                             "single isodecoder (YR tag == 1). Ignored for isoacceptor.")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for matrix.mtx / barcodes.tsv / features.tsv")
    args = parser.parse_args()

    if _clusterer is None:
        print("WARNING: umi_tools UMIClusterer not importable; UMI error correction "
              "is DISABLED (counts will be slightly inflated). Collision resolution "
              "still applies.", file=sys.stderr)

    trnadata = transcriptfile(args.trnatable)
    trnatranscripts = trnadata.gettranscripts()
    valid_barcodes = load_whitelist(args.whitelist) if args.whitelist else None

    # (cb, feature) -> Counter(umi -> supporting reads)
    reads = defaultdict(Counter)
    kept = skipped_not_trna = skipped_no_bcumi = skipped_not_whitelisted = 0
    skipped_multi_isodecoder = 0

    bam = pysam.AlignmentFile(args.input, "rb")
    refnames = bam.references
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue

        refname = refnames[read.reference_id]
        if refname not in trnatranscripts:
            skipped_not_trna += 1
            continue

        cb, umi = parse_barcode_umi(read.query_name, args.umi_separator)
        if cb is None:
            skipped_no_bcumi += 1
            continue

        if valid_barcodes is not None and cb not in valid_barcodes:
            skipped_not_whitelisted += 1
            continue

        if args.level == "isoacceptor":
            amino = trnadata.amino.get(refname, "Und")
            anticodon = trnadata.anticodon.get(refname, "NNN")
            feature = f"{amino}-{anticodon}"
        else:  # isodecoder: the individual transcript
            if args.isodecoder_count == "unique":
                # YR = # of tRNA transcripts (isodecoders) this read maps to,
                # written by choosemappings on every tRNA alignment. Absent -> assume 1.
                yr = int(read.get_tag("YR")) if read.has_tag("YR") else 1
                if yr > 1:
                    skipped_multi_isodecoder += 1
                    continue
            feature = refname
        reads[(cb, feature)][umi] += 1
        kept += 1
    bam.close()

    # --- Step 1: UMI error correction WITHIN each (cell, isoacceptor) ---
    # (cb, umi) -> {feature: supporting reads}  (after within-group collapse)
    molecule = defaultdict(dict)
    for (cb, feature), umi_read_counts in reads.items():
        for umi, nreads in cluster_umis(umi_read_counts).items():
            molecule[(cb, umi)][feature] = nreads

    # --- Step 2: resolve (cell, UMI) seen on >1 isoacceptor: max reads wins,
    # tie discards. Each surviving molecule = one count. ---
    counts = defaultdict(int)   # (cb, feature) -> molecule count
    n_molecules = collisions = ties_discarded = 0
    for (cb, umi), feats in molecule.items():
        n_molecules += 1
        if len(feats) == 1:
            counts[(cb, next(iter(feats)))] += 1
            continue
        collisions += 1
        top = max(feats.values())
        winners = [f for f, c in feats.items() if c == top]
        if len(winners) == 1:
            counts[(cb, winners[0])] += 1
        else:
            ties_discarded += 1  # ambiguous molecule, cannot assign -> drop

    # --- write matrix (features x barcodes, Matrix Market) ---
    barcodes = sorted({cb for cb, _ in counts})
    features = sorted({feat for _, feat in counts})
    barcode_idx = {b: i for i, b in enumerate(barcodes)}
    feature_idx = {f: i for i, f in enumerate(features)}

    os.makedirs(args.outdir, exist_ok=True)

    with open(os.path.join(args.outdir, "barcodes.tsv"), "w") as f:
        for b in barcodes:
            f.write(b + "\n")

    feature_type = "tRNA_isoacceptor" if args.level == "isoacceptor" else "tRNA_isodecoder"
    with open(os.path.join(args.outdir, "features.tsv"), "w") as f:
        for feat in features:
            f.write(f"{feat}\t{feat}\t{feature_type}\n")

    nonzero = [(feature_idx[feat] + 1, barcode_idx[cb] + 1, c)
               for (cb, feat), c in counts.items() if c > 0]
    with open(os.path.join(args.outdir, "matrix.mtx"), "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{len(features)} {len(barcodes)} {len(nonzero)}\n")
        for feat_i, bc_i, c in nonzero:
            f.write(f"{feat_i} {bc_i} {c}\n")

    total_umis = sum(counts.values())
    print(f"Reads counted (tRNA, barcode-tagged, primary): {kept}", file=sys.stderr)
    print(f"Reads skipped (not a tRNA transcript alignment): {skipped_not_trna}", file=sys.stderr)
    print(f"Reads skipped (could not parse CB/UMI from name): {skipped_no_bcumi}", file=sys.stderr)
    print(f"Reads skipped (barcode not in whitelist): {skipped_not_whitelisted}", file=sys.stderr)
    print(f"Feature level: {args.level}"
          + (f" ({args.isodecoder_count})" if args.level == "isodecoder" else ""), file=sys.stderr)
    if args.level == "isodecoder" and args.isodecoder_count == "unique":
        print(f"Reads skipped (multi-isodecoder, YR>1): {skipped_multi_isodecoder}", file=sys.stderr)
    print(f"UMI error correction: {'ON (umi_tools directional)' if _clusterer else 'OFF'}", file=sys.stderr)
    print(f"Molecules (post-correction): {n_molecules}", file=sys.stderr)
    print(f"  multi-isoacceptor molecules (collisions): {collisions}", file=sys.stderr)
    print(f"  of those discarded on a tie: {ties_discarded}", file=sys.stderr)
    print(f"UMI counts written: {total_umis}", file=sys.stderr)
    print(f"Cells: {len(barcodes)}, tRNA {args.level} features: {len(features)}", file=sys.stderr)


if __name__ == "__main__":
    main()