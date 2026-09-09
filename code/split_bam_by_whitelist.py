#!/usr/bin/env python3
"""
split_bam_by_whitelist.py — Split a barcode-tagged well BAM into one BAM per OB.

Input BAM has cell barcode + UMI embedded in the read name (umi_tools
extract convention: "{original_read_name}_{CELL_BARCODE}_{UMI}"), same
format count_singlecell_tRNA.py already parses. This partitions reads by
which OB's whitelist their CB falls in, so a whole-BAM classifier
(count_all_smRNA.py) can be run per OB instead of pooled across OBs.

Usage:
    python split_bam_by_whitelist.py \\
        --input fbc_only/02_sc_alignment/236860-1.srt.bam \\
        --ob OB1=/path/236860_1/outs/per_sample_outs/1/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \\
        --ob OB2=/path/236860_1/outs/per_sample_outs/2/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \\
        --outdir fbc_only/02_sc_alignment \\
        --prefix 236860-1

Writes fbc_only/02_sc_alignment/236860-1_OB1.srt.bam,
       fbc_only/02_sc_alignment/236860-1_OB2.srt.bam
(coordinate order preserved from the input, indexed).
"""
import argparse
import gzip
import os
import sys

import pysam


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_whitelist(path):
    # strip GEM-well suffix ("-1") so it matches the raw CB embedded in read names
    with opener(path) as f:
        return {line.strip().split("-")[0] for line in f if line.strip()}


def parse_barcode(qname, separator="_"):
    parts = qname.rsplit(separator, 2)
    if len(parts) != 3:
        return None
    return parts[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="coordinate-sorted, barcode-tagged BAM")
    ap.add_argument("--ob", action="append", required=True,
                     help="OB_LABEL=path/to/barcodes.tsv.gz, repeatable")
    ap.add_argument("--umi-separator", default="_")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True, help="well id, e.g. 236860-1")
    a = ap.parse_args()

    ob_whitelists = {}
    for spec in a.ob:
        label, path = spec.split("=", 1)
        wl = load_whitelist(path)
        ob_whitelists[label] = wl
        print(f"{label}: {len(wl)} whitelisted barcodes ({path})", file=sys.stderr)

    # barcode -> OB label (barcodes should be disjoint across OBs; warn if not)
    cb_to_ob = {}
    for label, wl in ob_whitelists.items():
        for cb in wl:
            if cb in cb_to_ob:
                print(f"WARNING: barcode {cb} in both {cb_to_ob[cb]} and {label}, "
                      f"keeping first assignment", file=sys.stderr)
                continue
            cb_to_ob[cb] = label

    os.makedirs(a.outdir, exist_ok=True)
    bam_in = pysam.AlignmentFile(a.input, "rb")
    outs = {
        label: pysam.AlignmentFile(
            os.path.join(a.outdir, f"{a.prefix}_{label}.srt.bam"), "wb", template=bam_in)
        for label in ob_whitelists
    }

    counts = {label: 0 for label in ob_whitelists}
    unmatched = 0
    for read in bam_in.fetch(until_eof=True):
        cb = parse_barcode(read.query_name, a.umi_separator)
        label = cb_to_ob.get(cb) if cb is not None else None
        if label is None:
            unmatched += 1
            continue
        outs[label].write(read)
        counts[label] += 1
    bam_in.close()
    for fh in outs.values():
        fh.close()

    for label in ob_whitelists:
        out_path = os.path.join(a.outdir, f"{a.prefix}_{label}.srt.bam")
        pysam.index(out_path)
        print(f"{label}: {counts[label]} reads written -> {out_path}", file=sys.stderr)
    print(f"reads with CB not in any OB whitelist: {unmatched}", file=sys.stderr)


if __name__ == "__main__":
    main()
