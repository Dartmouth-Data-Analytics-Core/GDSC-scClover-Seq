#!/usr/bin/env python3
"""
rRNA summary for one well: how much of each sample's (OB's) library is rRNA,
split into reads that MAPPED (rRNA / Mt_rRNA rows of
03_biotype_by_sample/biotype_by_sample_raw.txt) and reads that did NOT map
(the well's unmapped.fastq.gz, aligned here against the mouse rDNA repeat
by rule unmapped_vs_rdna; whatever fails to align lands in --remainder).

Each unmapped read is attributed to its OB through the cell barcode
umi_tools extract embedded in the read name ("{name}_{16bp CB}_{12bp UMI}")
matched against that OB's Cell Ranger barcode whitelist ("-1" GEM suffix
stripped, since embedded barcodes never carry it).

Outputs (per well):
  --out-tsv    one row per OB + ALL: mapped / unmapped / total rRNA counts and
               percentages, reads <= --short-max nt, and a profile of the
               remainder (unmapped reads that did not align to the rRNA reference).
  --out-hist   remainder read-length histogram per OB (long format).
  --out-fasta  top --top-n unique remainder sequences longer than --short-max,
               with counts, ready for BLAST.

The bowtie2 log of the rDNA alignment is cross-checked against the fastqs
(total reads and reads that aligned 0 times must match) so a wrong file
pairing fails loudly instead of producing a plausible-looking table.

Definitions worth knowing before reading the table:
  - mapped_reads   = sum of every biotype row for that OB's column in
                     biotype_by_sample_raw.txt (as count_all_smRNA.py counted
                     them), not an independent BAM read count.
  - unmapped_rRNA  = unmapped_reads - remainder_reads.
  - total_rRNA     = mapped rRNA + mapped Mt_rRNA + unmapped_rRNA.
  - "short"        = length <= --short-max (default 20). With bowtie2's default
                     local-mode --score-min (G,20,8) a 20 nt read cannot align
                     even with zero mismatches, so every short read ends up in
                     the remainder (checked, and the script aborts otherwise; that
                     is why there is no separate remainder_le column).
"""

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_whitelist(path):
    """barcodes.tsv[.gz] -> set of 16bp barcodes, "-1" suffix stripped."""
    bcs = set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            bc = line.strip().split("-")[0]
            if bc:
                bcs.add(bc)
    return bcs


def ob_of(header, barcode_to_ob, sep):
    parts = header[1:].split()[0].rsplit(sep, 2)  # [..., CB, UMI]
    if len(parts) != 3:
        return "unassigned"
    return barcode_to_ob.get(parts[1], "unassigned")


def read_records(path):
    """Yield (header, sequence) for each 4-line FASTQ record."""
    with open_maybe_gz(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield header, seq


def scan_unmapped(path, barcode_to_ob, sep, short_max):
    total, short = Counter(), Counter()
    for header, seq in read_records(path):
        ob = ob_of(header, barcode_to_ob, sep)
        total[ob] += 1
        if len(seq) <= short_max:
            short[ob] += 1
    return total, short


def scan_remainder(path, barcode_to_ob, sep, short_max):
    total, short, long_ = Counter(), Counter(), Counter()
    hist = Counter()                       # (ob, length) -> reads
    uniq = defaultdict(set)                # ob -> unique long sequences
    seq_counts = Counter()                 # long sequence -> reads (whole well)
    for header, seq in read_records(path):
        ob = ob_of(header, barcode_to_ob, sep)
        n = len(seq)
        total[ob] += 1
        hist[(ob, n)] += 1
        if n <= short_max:
            short[ob] += 1
        else:
            long_[ob] += 1
            uniq[ob].add(seq)
            seq_counts[seq] += 1
    return total, short, long_, hist, uniq, seq_counts


def read_biotype(path, ext_obs):
    """Per external OB: total mapped reads (all rows), rRNA row, Mt_rRNA row."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        sys.exit(f"{path}: no data rows")
    offset = len(rows[0]) - len(header)   # 1 when data rows carry a leading label column
    col = {name: i + offset for i, name in enumerate(header)}
    out = {}
    for ext in ext_obs:
        key = f"biosample_{ext}"
        if key not in col:
            sys.exit(f"{path}: no column '{key}' (found: {', '.join(header)})")
        j = col[key]
        mapped = rrna = mt = 0
        for r in rows:
            try:
                v = int(r[j])
            except (ValueError, IndexError):
                continue
            mapped += v
            if r[0] == "rRNA":
                rrna += v
            elif r[0] == "Mt_rRNA":
                mt += v
        out[ext] = {"mapped": mapped, "rRNA": rrna, "MtrRNA": mt}
    return out


def parse_bowtie2_log(path):
    total = unaligned = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*(\d+) reads; of these:", line)
            if m:
                total = int(m.group(1))
            m = re.match(r"\s*(\d+) \([\d.]+%\) aligned 0 times", line)
            if m:
                unaligned = int(m.group(1))
    if total is None or unaligned is None:
        sys.exit(f"{path}: could not find the bowtie2 alignment summary")
    return total, unaligned


def pct(part, whole):
    return f"{100 * part / whole:.2f}" if whole else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well", required=True)
    ap.add_argument("--unmapped", required=True, help="{well}.unmapped.fastq.gz")
    ap.add_argument("--remainder", required=True, help="unmapped reads that did NOT align to the rDNA reference")
    ap.add_argument("--bowtie2-log", required=True, help="summary written by the rDNA bowtie2 run")
    ap.add_argument("--biotype-raw", required=True, help="03_biotype_by_sample/biotype_by_sample_raw.txt")
    ap.add_argument("--ob", action="append", required=True,
                    help="EXTERNAL_OB=path/to/barcodes.tsv.gz, repeatable (this well's OBs)")
    ap.add_argument("--short-max", type=int, default=20)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--umi-separator", default="_")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-hist", required=True)
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    S = args.short_max
    barcode_to_ob, ext_obs = {}, []
    for spec in args.ob:
        label, path = spec.split("=", 1)
        ext_obs.append(label)
        for bc in load_whitelist(path):
            barcode_to_ob[bc] = label

    mapped = read_biotype(args.biotype_raw, ext_obs)
    un_total, un_short = scan_unmapped(args.unmapped, barcode_to_ob, args.umi_separator, S)
    (rem_total, rem_short, rem_long,
     hist, uniq, seq_counts) = scan_remainder(args.remainder, barcode_to_ob, args.umi_separator, S)

    log_total, log_unaligned = parse_bowtie2_log(args.bowtie2_log)
    if sum(un_total.values()) != log_total or sum(rem_total.values()) != log_unaligned:
        sys.exit(
            f"Inconsistent inputs for {args.well}: unmapped fastq has {sum(un_total.values()):,} reads and "
            f"remainder fastq {sum(rem_total.values()):,}, but the bowtie2 log reports {log_total:,} reads / "
            f"{log_unaligned:,} aligned 0 times. These files are not from the same run.")

    labels = list(ext_obs)
    if un_total.get("unassigned"):
        labels.append("unassigned")

    header = ["well", "ob",
              "mapped_reads", "mapped_rRNA", "mapped_rRNA_pct", "mapped_MtrRNA", "mapped_MtrRNA_pct",
              "unmapped_reads", "unmapped_rRNA", "unmapped_rRNA_pct", f"unmapped_le{S}nt", f"unmapped_le{S}nt_pct",
              "total_reads", "total_rRNA", "total_rRNA_pct", "unmapped_pct_of_total",
              "remainder_reads", f"remainder_gt{S}nt",
              f"remainder_gt{S}nt_unique_seqs"]

    def row_for(label, m_reads, m_rrna, m_mt, un, un_le, rem, rem_le, rem_gt, rem_uniq):
        if rem_le != un_le:   # a read this short cannot align to rDNA, so every short unmapped read must be in the remainder
            sys.exit(f"{args.well} {label}: {un_le:,} unmapped reads <= {S} nt but {rem_le:,} in the remainder")
        un_rrna = un - rem
        tot = m_reads + un
        tot_rrna = m_rrna + m_mt + un_rrna
        return [args.well, label,
                m_reads, m_rrna, pct(m_rrna, m_reads), m_mt, pct(m_mt, m_reads),
                un, un_rrna, pct(un_rrna, un), un_le, pct(un_le, un),
                tot, tot_rrna, pct(tot_rrna, tot), pct(un, tot),
                rem, rem_gt, rem_uniq]

    rows, sums = [], Counter()
    for lab in labels:
        m = mapped.get(lab, {"mapped": 0, "rRNA": 0, "MtrRNA": 0})
        vals = dict(m_reads=m["mapped"], m_rrna=m["rRNA"], m_mt=m["MtrRNA"], un=un_total[lab], un_le=un_short[lab],
                    rem=rem_total[lab], rem_le=rem_short[lab], rem_gt=rem_long[lab])
        rows.append(row_for(lab, rem_uniq=len(uniq[lab]), **vals))
        sums.update(vals)
    rows.append(row_for("ALL", rem_uniq=len(seq_counts), **sums))

    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    all_hist = Counter()
    for (lab, n), c in hist.items():
        all_hist[n] += c
    with open(args.out_hist, "w") as fh:
        fh.write("well\tob\tlength\treads\n")
        for lab in labels:
            for n in sorted(n for (l, n) in hist if l == lab):
                fh.write(f"{args.well}\t{lab}\t{n}\t{hist[(lab, n)]}\n")
        for n in sorted(all_hist):
            fh.write(f"{args.well}\tALL\t{n}\t{all_hist[n]}\n")

    with open(args.out_fasta, "w") as fh:
        for i, (seq, c) in enumerate(seq_counts.most_common(args.top_n), start=1):
            fh.write(f">{args.well}_seq{i}_count{c}\n{seq}\n")


if __name__ == "__main__":
    main()
