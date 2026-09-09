#!/usr/bin/env python3
"""
count_all_smRNA.py — Classify and count reads across tRNA and small RNA biotypes.

For each sample BAM, reads are assigned to feature categories in priority order:
  1. Mature tRNA   — O(1) dict lookup by BAM chromosome name
  2. Pre-tRNA loci — RangeBin spatial search on genomic loci
  3. Ensembl GTF   — highest-priority biotype wins per read
  4. BED features  — arbitrary extra feature files
  5. Extra sequences
  6. Unclassified  — counted as 'other'

Outputs normalized biotype counts grouped by replicate, raw per-sample counts,
amino-acid and anticodon tRNA counts, mismatch distributions, and read-length
distributions broken down by RNA class. Parallelises across samples via
multiprocessing.Pool when --cores > 1.
"""

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
from multiprocessing import Pool


class counttypes:
    """Per-sample read classification counters for tRNA and smRNA features."""

    def __init__(self, samplename, bamfile, trnas=list(), trnaloci=list(), emblgenes=list(), otherfeats=list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats
        self.emblbiotypes = set()
        self.aminos = set()
        self.bedtypes = set()
        self.extraseqtypes = set()
        self.anticodons = set()

        self.embltypecounts = defaultdict(int)
        self.bedtypecounts = defaultdict(int)
        self.trnafragtypes = defaultdict(int)
        self.totalreads = 0
        self.trnareads = 0
        self.otherreads = 0
        self.readlengths = defaultdict(int)
        self.trnareadlengths = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.aminocounts = defaultdict(int)
        self.anticodoncounts = defaultdict(int)
        self.indelreads = defaultdict(int)
        self.trnaanticounts = defaultdict(int)
        self.trnaemblcounts = defaultdict(int)
        self.pretrnareadlengths = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.partiallocuscounts = defaultdict(int)
        self.fulllocuscounts = defaultdict(int)
        self.extraseqcounts = defaultdict(int)
        self.trnaantilocuscounts = defaultdict(int)
        self.mismatchcounts = defaultdict(int)
        self.trnamismatchcounts = defaultdict(int)

    def addsamplecounts(self):
        self.totalreads += 1
    def addreadlengths(self, length):
        self.readlengths[length] += 1
    def addtrnareadlengths(self, length):
        self.trnareadlengths[length] += 1
    def addpretrnareadlengths(self, length):
        self.pretrnareadlengths[length] += 1
    def addpartiallocuscounts(self, currbed):
        self.partiallocuscounts[currbed] += 1
    def addtrnasamplecounts(self):
        self.trnareads += 1
    def addtrnacounts(self, currbed):
        self.trnacounts[currbed] += 1
    def addtrnaantilocuscounts(self, currbed):
        self.trnaantilocuscounts[currbed] += 1
    def addtrnalocuscounts(self, currbed):
        self.trnalocuscounts[currbed] += 1
    def addfulllocuscounts(self, currbed):
        self.fulllocuscounts[currbed] += 1
    def addaminocounts(self, curramino):
        self.aminos.add(curramino)
        self.aminocounts[curramino] += 1
    def addanticodoncounts(self, curranticodon):
        self.anticodons.add(curranticodon)
        self.anticodoncounts[curranticodon] += 1
    def addindelreads(self, curramino):
        self.indelreads[curramino] += 1
    def addmismatchcounts(self, mismatchcounts):
        self.mismatchcounts[mismatchcounts] += 1
    def addtrnamismatchcounts(self, mismatchcounts):
        self.trnamismatchcounts[mismatchcounts] += 1
    def addotherreads(self):
        self.otherreads += 1
    def addtrnaantisense(self, currbed):
        self.trnaanticounts[currbed] += 1
    def addemblcounts(self, currtype):
        self.emblbiotypes.add(currtype)
        self.embltypecounts[currtype] += 1
    def addbedcounts(self, genetype):
        self.bedtypes.add(genetype)
        self.bedtypecounts[genetype] += 1
    def addextracounts(self, genetype):
        self.extraseqtypes.add(genetype)
        self.extraseqcounts[genetype] += 1


def counttypereads(bamfile, samplename, trnainfo, trnaloci, trnalist, maturenames,
                   featurelist=dict(), otherseqlist=list(), embllist=None,
                   maxmismatches=None, bamnofeature=False):
    """Classify every primary-alignment read in bamfile into a read-type category.

    Mature tRNA reads are resolved first via a pre-built dict keyed on chromosome
    name, which covers the majority of reads in O(1). Remaining reads fall through
    to the pre-tRNA loci, Ensembl GTF, BED, and extra-sequence layers in order.

    Returns a counttypes object populated with per-category read tallies.
    """
    bedlist = list(featurelist.keys())
    readtypecounts = counttypes(samplename, bamfile, trnas=trnalist, trnaloci=trnaloci,
                                emblgenes=embllist, otherfeats=bedlist)
    mitochrom = None
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    ncrnaorder = defaultdict(int)
    currbam = bamfile

    for i, curr in enumerate(reversed(["snoRNA", "miRNA", "rRNA", "snRNA", "misc_RNA", "lincRNA", "protein_coding"])):
        ncrnaorder[curr] = i + 1

    # Pre-build flat chrom -> (currbed, feature) mapping for O(1) mature tRNA lookup.
    # Each mature tRNA has its own BAM chromosome, so a dict hit replaces RangeBin
    # search for the large majority of reads.
    chrom_to_mature = {}
    for currbed in trnalist:
        for tname, feat in maturenames[currbed].items():
            chrom_to_mature[tname] = (currbed, feat)

    try:
        if not os.path.isfile(currbam + ".bai") or os.path.getmtime(currbam + ".bai") < os.path.getmtime(currbam):
            pysam.index("" + currbam)
        bamfile = pysam.Samfile("" + currbam, "rb")
        if bamnofeature:
            outname = os.path.splitext(currbam)[0] + "_nofeat.bam"
            outbamnofeature = pysam.Samfile(outname, "wb", template=bamfile)
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()

    for currread in getbam(bamfile, primaryonly=True):
        readlength = currread.getlength()
        gotread = False
        readtypecounts.totalreads += 1
        readtypecounts.mismatchcounts[currread.getmismatches()] += 1
        if currread.hasindel():
            readtypecounts.indelreads[readlength] += 1
        readtypecounts.readlengths[readlength] += 1

        # --- Layer 1: mature tRNA (O(1) dict lookup) ---
        if currread.chrom in chrom_to_mature:
            currbed, currfeat = chrom_to_mature[currread.chrom]
            if currread.strand == "+":
                readtypecounts.trnareadlengths[readlength] += 1
                readtypecounts.trnareads += 1
                readtypecounts.trnacounts[currbed] += 1
                amino = trnainfo.getamino(currfeat.name)
                readtypecounts.aminos.add(amino)
                readtypecounts.aminocounts[amino] += 1
                anti = trnainfo.getanticodon(currfeat.name)
                readtypecounts.anticodons.add(anti)
                readtypecounts.anticodoncounts[anti] += 1
                readtypecounts.trnamismatchcounts[currread.getmismatches()] += 1
            elif currfeat.antisense().coverage(currread) > 10:
                readtypecounts.trnaanticounts[currbed] += 1
            continue

        # --- Layer 2: pre-tRNA genomic loci ---
        for currbed in trnaloci:
            for currfeat in trnaloci[currbed].getbin(currread):
                expandfeat = currfeat.addmargin(30)
                if currfeat.coverage(currread) > 10:
                    if currfeat.strand != currread.strand:
                        readtypecounts.trnaantilocuscounts[currbed] += 1
                        gotread = True
                        break
                    if (currread.start + minpretrnaextend <= currfeat.start or
                            currread.end - minpretrnaextend >= currfeat.end):
                        pass
                    readtypecounts.pretrnareadlengths[readlength] += 1
                    readtypecounts.trnalocuscounts[currbed] += 1
                    if (currread.start + fullpretrnathreshold < currfeat.start and
                            currread.end - fullpretrnathreshold + 3 > currfeat.end):
                        readtypecounts.fulllocuscounts[currbed] += 1
                    else:
                        readtypecounts.partiallocuscounts[currbed] += 1
                    gotread = True
                    break
                if currfeat.getdownstream(30).coverage(currread) > 10:
                    readtypecounts.trnaanticounts[currbed] += 1
                    gotread = True
                    break
                elif expandfeat.antisense().coverage(currread) > 5:
                    readtypecounts.trnaanticounts[currbed] += 1
                    gotread = True
                    break
            if gotread:
                break
        if gotread:
            continue

        # --- Layer 3: Ensembl GTF biotypes (highest-priority biotype wins) ---
        if embllist is not None:
            currtype = None
            for currfeat in embllist.getbin(currread):
                if currfeat.coverage(currread) > 10:
                    if currtype is None or ncrnaorder[currfeat.data["biotype"]] > ncrnaorder[currtype]:
                        currtype = currfeat.data["biotype"]
            if currtype is not None:
                readtypecounts.emblbiotypes.add(currtype)
                readtypecounts.embltypecounts[currtype] += 1
                gotread = True
        if gotread:
            continue

        # --- Layer 4: BED features ---
        for currbed in bedlist:
            for currfeat in featurelist[currbed].getbin(currread):
                if currfeat.coverage(currread) > 10:
                    readtypecounts.bedtypes.add(currbed)
                    readtypecounts.bedtypecounts[currbed] += 1
                    gotread = True
                    break
            if gotread:
                break
        if gotread:
            continue

        # --- Layer 5: extra sequences ---
        for currbed in otherseqlist:
            for currfeat in otherseqlist[currbed][currread.chrom]:
                if currfeat.coverage(currread) > 10:
                    readtypecounts.extraseqtypes.add(currbed)
                    readtypecounts.extraseqcounts[currbed] += 1
                    gotread = True
                    break
            if gotread:
                break
        if gotread:
            continue

        # --- Layer 6: unclassified ---
        readtypecounts.otherreads += 1
        if embllist is not None and mitochrom == currread.chrom:
            readtypecounts.emblbiotypes.add("Mitochondrial_other")
            readtypecounts.embltypecounts["Mitochondrial_other"] += 1
        if not gotread and bamnofeature:
            outbamnofeature.write(currread.bamline)

    return readtypecounts


def printtypefile(countfile, samples, sampledata, allcounts, trnalist, trnaloci,
                  bedtypes, emblbiotypes, sizefactor, extraseqtypes=set()):
    """Write normalized biotype counts per replicate group to countfile.

    Biotypes are written in a fixed order: snoRNA/snRNA/scaRNA/sRNA/miRNA first,
    then other biotypes, then Mt_rRNA/Mt_tRNA/rRNA last.
    """
    replicates = list(sampledata.allreplicates())
    print("\t".join(replicates), file=countfile)

    print("other\t" + "\t".join(
        str(sum(allcounts[s].otherreads / sizefactor[s]
                for s in sampledata.getrepsamples(rep)))
        for rep in replicates), file=countfile)

    for currbed in bedtypes:
        print(os.path.basename(currbed).split(".")[0] + "\t" + "\t".join(
            str(sum(allcounts[s].bedtypecounts[currbed] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)

    for currname in extraseqtypes:
        print(currname + "_seq\t" + "\t".join(
            str(sum(allcounts[s].extraseqcounts[currname] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)

    biotypefirst = ["snoRNA", "snRNA", "scaRNA", "sRNA", "miRNA"]
    biotypelast = ["Mt_rRNA", "Mt_tRNA", "rRNA"]
    otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
    biotypeorder = biotypefirst + otherbiotypes + biotypelast

    for currbiotype in biotypeorder:
        print(currbiotype + "\t" + "\t".join(
            str(sum(allcounts[s].embltypecounts[currbiotype] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)

    for currbed in trnaloci:
        print("pretRNA_antisense\t" + "\t".join(
            str(sum(allcounts[s].trnaantilocuscounts[currbed] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)
        print("pretRNA\t" + "\t".join(
            str(sum(allcounts[s].trnalocuscounts[currbed] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)

    for currbed in trnalist:
        print("tRNA_antisense\t" + "\t".join(
            str(sum(allcounts[s].trnaanticounts[currbed] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)
        print("tRNA\t" + "\t".join(
            str(sum(allcounts[s].trnacounts[currbed] / sizefactor[s]
                    for s in sampledata.getrepsamples(rep)))
            for rep in replicates), file=countfile)


def printrealcounts(countfile, samples, sampledata, allcounts, trnalist, trnaloci,
                    bedtypes, emblbiotypes, extraseqtypes=set()):
    """Write raw (un-normalized) biotype counts per sample to countfile."""
    biotypefirst = ["snoRNA", "snRNA", "scaRNA", "sRNA", "miRNA"]
    biotypelast = ["Mt_rRNA", "Mt_tRNA", "rRNA"]
    otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
    biotypeorder = biotypefirst + otherbiotypes + biotypelast

    print("\t".join(samples), file=countfile)
    print("other\t" + "\t".join(str(allcounts[s].otherreads) for s in samples), file=countfile)
    for currbed in bedtypes:
        print(os.path.basename(currbed) + "\t" +
              "\t".join(str(allcounts[s].bedtypecounts[currbed]) for s in samples), file=countfile)
    for currname in extraseqtypes:
        print(currname + "_seq\t" +
              "\t".join(str(allcounts[s].extraseqcounts[currname]) for s in samples), file=countfile)
    for currbiotype in reversed(biotypeorder):
        print(currbiotype + "\t" +
              "\t".join(str(allcounts[s].embltypecounts[currbiotype]) for s in samples), file=countfile)
    for currbed in trnaloci:
        print("pretRNA\t" +
              "\t".join(str(allcounts[s].trnalocuscounts[currbed]) for s in samples), file=countfile)
    for currbed in trnalist:
        print("tRNA_antisense\t" +
              "\t".join(str(allcounts[s].trnaanticounts[currbed]) for s in samples), file=countfile)
        print("tRNA\t" +
              "\t".join(str(allcounts[s].trnacounts[currbed]) for s in samples), file=countfile)


def printaminocounts(trnaaminofilename, sampledata, trnainfo, allcounts, sizefactor):
    """Write normalized amino-acid-level tRNA counts per replicate group."""
    replicates = list(sampledata.allreplicates())
    aminos = trnainfo.allaminos()
    with open(trnaaminofilename, "w") as f:
        print("\t".join(replicates), file=f)
        for curramino in aminos:
            print(curramino + "\t" + "\t".join(
                str(sum(allcounts[s].aminocounts[curramino] / sizefactor[s]
                        for s in sampledata.getrepsamples(rep)))
                for rep in replicates), file=f)


def printanticodoncounts(trnaanticodonfilename, sampledata, trnainfo, allcounts, sizefactor):
    """Write normalized anticodon-level tRNA counts per sample."""
    anticodons = trnainfo.allanticodons()
    allsamples = list(sampledata.getsamples())
    with open(trnaanticodonfilename, "w") as f:
        print("\t".join(allsamples), file=f)
        for curranticodon in anticodons:
            print(curranticodon + "\t" + "\t".join(
                str(allcounts[s].anticodoncounts[curranticodon] / sizefactor[s])
                for s in allsamples), file=f)


def printmismatchcounts(trnamismatchname, sampledata, trnainfo, allcounts, sizefactor):
    """Write per-sample mismatch count distributions split into tRNA vs. non-tRNA reads."""
    allsamples = list(sampledata.getsamples())
    mismatchrange = list(range(10))
    with open(trnamismatchname, "w") as f:
        print("count\ttype\t" + "\t".join(allsamples), file=f)
        for currmismatch in mismatchrange:
            print(str(currmismatch) + "\ttrna\t" + "\t".join(
                str(allcounts[s].trnamismatchcounts[currmismatch] / sizefactor[s])
                for s in allsamples), file=f)
            print(str(currmismatch) + "\tnontrna\t" + "\t".join(
                str(allcounts[s].mismatchcounts[currmismatch] / sizefactor[s]
                    - allcounts[s].trnamismatchcounts[currmismatch] / sizefactor[s])
                for s in allsamples), file=f)


def printlengthfile(readlengthfile, samples, allcounts):
    """Write per-sample read-length distributions broken down by RNA type."""
    with open(readlengthfile, "w") as f:
        print("Length\tSample\tother\ttrnas\tpretrnas", file=f)
        for currsample in samples:
            maxlen = max(allcounts[currsample].readlengths.keys())
            for curr in range(0, maxlen + 1):
                trna_pre = (allcounts[currsample].trnareadlengths[curr]
                            + allcounts[currsample].pretrnareadlengths[curr])
                other = allcounts[currsample].readlengths[curr] - trna_pre
                print("\t".join([str(curr), currsample, str(other),
                                 str(allcounts[currsample].trnareadlengths[curr]),
                                 str(allcounts[currsample].pretrnareadlengths[curr])]), file=f)


def counttypereadspool(args):
    return counttypereads(*args[0], **args[1])


def compressargs(*args, **kwargs):
    return tuple([args, kwargs])


def main(**argdict):
    """Load reference files, classify reads in each sample BAM, write output tables."""
    argdict = defaultdict(lambda: None, argdict)
    ensemblgtf = argdict["ensemblgtf"]
    bamnofeature = argdict["bamnofeature"]
    trnatable = argdict["trnatable"]
    trnaaminofilename = argdict["trnaaminofile"]
    trnaanticodonfilename = argdict["trnaanticodonfile"]

    bamdir = argdict["bamdir"] if argdict["bamdir"] is not None else "./"
    sampledata = samplefile(argdict["samplefile"], bamdir=bamdir)

    cores = argdict["cores"] if argdict["cores"] is not None else 1
    threadmode = cores > 1

    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"])
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print("Size factor file " + argdict["sizefactors"] + " missing " + currsample,
                      file=sys.stderr)
                sys.exit(1)

    bedfiles = argdict["bedfile"] if argdict["bedfile"] is not None else list()
    locifiles = argdict["trnaloci"] if argdict["trnaloci"] is not None else list()
    maturetrnafiles = argdict["maturetrnas"] if argdict["maturetrnas"] is not None else list()

    readlengthfile = argdict["readlengthfile"]
    mismatchfilename = argdict["mismatchfile"]

    if argdict["realcountfile"] is None or argdict["realcountfile"] == "stdout":
        realcountfile = sys.stdout
    else:
        realcountfile = open(argdict["realcountfile"], "w")

    if argdict["countfile"] is None or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")

    maturenames = dict()
    trnainfo = transcriptfile(trnatable)
    samples = list(sampledata.getsamples())

    try:
        featurelist = dict()
        trnaloci = dict()
        trnalist = dict()
        for currfile in bedfiles:
            featurelist[currfile] = RangeBin(readfeatures(currfile))
        for currfile in locifiles:
            trnaloci[currfile] = RangeBin(readbed(currfile), binfactor=10000)
        for currfile in maturetrnafiles:
            matlist = list(readbed(currfile))
            trnalist[currfile] = list(matlist)
            maturenames[currfile] = {curr.name: curr for curr in matlist}
        embllist = RangeBin(readgtf(ensemblgtf, filtertypes=set())) if ensemblgtf is not None else None
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()

    maxmismatches = None
    allcounts = dict()
    otherseqlist = dict()

    if threadmode:
        countpool = Pool(processes=cores)
        arglist = []
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            arglist.append(compressargs(currbam, currsample, trnainfo, trnaloci, trnalist, maturenames,
                                        otherseqlist=otherseqlist, embllist=embllist,
                                        featurelist=featurelist, maxmismatches=maxmismatches,
                                        bamnofeature=bamnofeature))
        results = countpool.map(counttypereadspool, arglist)
        for i, curr in enumerate(samples):
            allcounts[curr] = results[i]
    else:
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = counttypereads(currbam, currsample, trnainfo, trnaloci, trnalist,
                                                   maturenames, otherseqlist=otherseqlist,
                                                   embllist=embllist, featurelist=featurelist,
                                                   maxmismatches=maxmismatches,
                                                   bamnofeature=bamnofeature)

    emblbiotypes = set(itertools.chain.from_iterable(c.emblbiotypes for c in allcounts.values()))
    bedtypes = set(itertools.chain.from_iterable(c.bedtypes for c in allcounts.values()))
    extraseqtypes = set(itertools.chain.from_iterable(c.extraseqtypes for c in allcounts.values()))

    printtypefile(countfile, samples, sampledata, allcounts, trnalist, trnaloci,
                  bedtypes, emblbiotypes, sizefactor, extraseqtypes=extraseqtypes)
    printrealcounts(realcountfile, samples, sampledata, allcounts, trnalist, trnaloci,
                    bedtypes, emblbiotypes, extraseqtypes=extraseqtypes)

    if readlengthfile is not None:
        printlengthfile(readlengthfile, samples, allcounts)
    if trnaaminofilename is not None:
        printaminocounts(trnaaminofilename, sampledata, trnainfo, allcounts, sizefactor)
    if trnaanticodonfilename is not None:
        printanticodoncounts(trnaanticodonfilename, sampledata, trnainfo, allcounts, sizefactor)
    if mismatchfilename is not None:
        printmismatchcounts(mismatchfilename, sampledata, trnainfo, allcounts, sizefactor)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads mapping to tRNA and other small RNA features.")
    parser.add_argument("--samplefile",
                        help="Sample file listing sample IDs, groups, and BAM directories")
    parser.add_argument("--sizefactors",
                        help="Optional tab-delimited size factors file for normalization")
    parser.add_argument("--bedfile", nargs="*", default=list(),
                        help="BED file(s) with non-tRNA features")
    parser.add_argument("--ensemblgtf",
                        help="Ensembl GTF with biotype annotations")
    parser.add_argument("--trnaloci", nargs="+", default=list(),
                        help="BED file(s) with pre-tRNA genomic loci")
    parser.add_argument("--maturetrnas", nargs="+", default=list(),
                        help="BED file(s) with mature tRNA features")
    parser.add_argument("--trnatable",
                        help="tRNA transcript table (db-trnatable.txt)")
    parser.add_argument("--trnaaminofile",
                        help="Output file for amino-acid-level tRNA counts")
    parser.add_argument("--readlengthfile",
                        help="Output file for read-length distributions by RNA class")
    parser.add_argument("--realcountfile",
                        help="Output file for raw (un-normalized) biotype counts per sample")
    parser.add_argument("--countfile",
                        help="Output file for normalized biotype counts by replicate group")
    parser.add_argument("--mismatchfile",
                        help="Output file for mismatch count distributions")
    parser.add_argument("--trnaanticodonfile",
                        help="Output file for anticodon-level tRNA counts per sample")
    parser.add_argument("--bamnofeature", action="store_true", default=False,
                        help="Write unclassified reads to a separate BAM file")
    parser.add_argument("--mitochrom",
                        help="Mitochondrial chromosome name in the database")
    parser.add_argument("--cores", type=int, default=1,
                        help="Number of parallel worker processes (default: 1)")

    args = parser.parse_args()
    main(**vars(args))
