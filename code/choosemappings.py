#!/usr/bin/env python3

import sys
import os.path
import itertools
import pysam
import argparse

from trnasequtils import *

from collections import defaultdict

defminnontrnasize = 20
maxmaps = 50


def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)


def getbesttrnamappings(trnafile, inputbam="-", bamout=True, logfile=sys.stderr,
                        progname=None, fqname=None, libname=None,
                        setcountfile=None, extraseqfilename=None,
                        minnontrnasize=defminnontrnasize):

    trnadata = transcriptfile(trnafile)
    trnatranscripts = set(trnadata.gettranscripts())
    extraseqs = set()
    seqdata = None
    if extraseqfilename is not None:
        seqdata = extraseqfile(extraseqfilename)
        extraseqs = set(trnadata.gettranscripts())

    totalreads = 0
    multimaps = 0
    duperemove = 0
    shortened = 0
    mapsremoved = 0
    totalmaps = 0

    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0])) + "/"
    gitversion, githash = getgithash(scriptdir)

    trnareads = 0
    maxreads = 0
    diffreads = 0
    uniquenontrnas = 0
    nonuniquenontrnas = 0

    # Open BAM from file path or stdin
    mode = "rb" if inputbam != "-" else "r"
    bamfile = pysam.AlignmentFile(inputbam, mode)

    # Cache reference names upfront — avoids a C-boundary crossing per call
    refnames = bamfile.references

    ambanticodon = 0
    ambamino = 0
    ambtrna = 0
    acsets = defaultdict(int)
    aminosets = defaultdict(int)
    trnasetcounts = defaultdict(int)
    newheader = bamfile.header.to_dict()
    newheader["RG"] = [{}]
    imperfect = 0
    if "PG" not in newheader:
        newheader["PG"] = []
    if progname is not None:
        newheader["PG"].append({"PN": progname, "ID": progname, "VN": gitversion})
    if fqname is not None:
        newheader["RG"][0]["ID"] = fqname
    if libname is not None:
        newheader["RG"][0]["LB"] = libname

    outmode = "wb" if bamout else "w"
    outfile = pysam.AlignmentFile("-", outmode, header=newheader)

    for pairedname, allmaps in itertools.groupby(bamfile, lambda x: x.query_name):
        allmaps = list(allmaps)
        # any() short-circuits; sum() would iterate the whole list
        if any(curr.flag & 0x004 for curr in allmaps):
            continue
        totalreads += 1

        mappings = 0
        currscore = None
        newset = set()
        readlength = None

        for currmap in allmaps:
            totalmaps += 1
            if currmap.tid == -1:
                continue
            readlength = len(currmap.query_sequence)
            mappings += 1
            # get_tag is a single pysam call; avoids building a tag dict per alignment
            score = currmap.get_tag("AS")
            if currscore is None or currscore < score:
                newset = {currmap}
                currscore = score
            elif currscore == score:
                newset.add(currmap)

        if mappings > 1:
            multimaps += 1
        if len(newset) < mappings:
            shortened += 1
        if len(newset) >= 50:
            maxreads += 1

        finalset = []

        # Single pass over newset — refnames[] is a tuple index, no C call
        trnamappings = [curr for curr in newset if refnames[curr.tid] in trnatranscripts]

        if trnamappings:
            trnareads += 1
            diff = len(newset) - len(trnamappings)
            anticodons = frozenset(trnadata.getanticodon(refnames[curr.tid]) for curr in trnamappings)
            aminos     = frozenset(trnadata.getamino(refnames[curr.tid])     for curr in trnamappings)
            locusmaps  = list(itertools.chain.from_iterable(
                trnadata.transcriptdict[refnames[curr.tid]] for curr in trnamappings))

            if trnamappings[0].get_tag("XM") + trnamappings[0].get_tag("XO") > 0:
                imperfect += 1
            if len(anticodons - frozenset(['NNN'])) > 1:
                ambanticodon += 1
                acsets[anticodons] += 1
            if len(aminos - frozenset(['Und'])) > 1:
                ambamino += 1
                aminosets[aminos] += 1
            if diff > 0:
                diffreads += 1
            if len(trnamappings) > 1:
                ambtrna += 1
            if setcountfile is not None:
                trnasetcounts[frozenset(refnames[curr.tid] for curr in trnamappings)] += 1
            for currtrnamap in trnamappings:
                currtrnamap.tags = (currtrnamap.tags +
                    [("YA", len(anticodons)), ("YM", len(aminos)),
                     ("YR", len(trnamappings)), ("YL", len(locusmaps))])
            finalset = trnamappings

        elif extraseqs and any(refnames[curr.tid] in extraseqs for curr in newset):
            for currseqset in seqdata.seqlist():
                if any(refnames[curr.tid] in extraseqdict[currseqset] for curr in newset):
                    finalset = [curr for curr in newset
                                if refnames[curr.tid] in extraseqdict[currseqset]]
                    break

        else:
            if readlength is None or readlength < minnontrnasize:
                continue
            finalset = list(newset)
            if len(finalset) > maxmaps:
                duperemove += 1
                continue
            if len(newset) > 1:
                nonuniquenontrnas += 1
            else:
                uniquenontrnas += 1

        mapsremoved += mappings - len(finalset)
        if not any(isprimarymapping(curr) for curr in finalset):
            for i, curr in enumerate(finalset):
                if i == 0:
                    curr.flag &= ~0x0100
                outfile.write(curr)
        else:
            for curr in finalset:
                outfile.write(curr)

    if setcountfile is not None:
        with open(setcountfile, "w") as setcounts:
            for currset in trnasetcounts:
                print(",".join(currset) + "\t" + str(trnasetcounts[currset]), file=setcounts)

    print("tRNA Reads with multiple transcripts:" + str(ambtrna),        file=logfile)
    print("tRNA Reads with multiple anticodons:"  + str(ambanticodon),   file=logfile)
    print("tRNA Reads with multiple aminos:"      + str(ambamino),       file=logfile)
    print("Total tRNA Reads:"                     + str(trnareads),      file=logfile)
    print("Single mapped non-tRNAs:"              + str(uniquenontrnas), file=logfile)
    print("Multiply mapped non-tRNAs:"            + str(nonuniquenontrnas), file=logfile)
    print("Imperfect matches:"                    + str(imperfect) + "/" + str(trnareads), file=logfile)
    outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select best tRNA mappings from a bowtie2 BAM')
    parser.add_argument('trnaname',
                        help='tRNA table file (db-trnatable.txt)')
    parser.add_argument('--input', default="-",
                        help='Input BAM file (default: stdin)')
    parser.add_argument('--progname',
                        help='Program name written into BAM @PG header')
    parser.add_argument('--fqname',
                        help='Fastq file name written into BAM @RG header')
    parser.add_argument('--expname',
                        help='Library/experiment name')
    parser.add_argument('--trnasetcounts',
                        help='Output file for per-tRNA-set read counts')
    parser.add_argument('--minnontrnasize', type=int, default=20,
                        help='Minimum read length to retain non-tRNA alignments (default: 20)')

    args = parser.parse_args()
    getbesttrnamappings(args.trnaname, inputbam=args.input,
                        progname=args.progname, fqname=args.fqname, libname=args.expname,
                        setcountfile=args.trnasetcounts, minnontrnasize=args.minnontrnasize)
