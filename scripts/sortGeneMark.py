#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Sort genemark.gtf by coordinates; change gene and transcript IDs to appear in
# ascending order.
# ==============================================================


import argparse
import csv
import re
import sys
import tempfile
import os
import subprocess


def extractFeatureGtf(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def replaceFeatureGtf(text, feature, oldValue, newValue):
    old = feature + " \"" + oldValue + "\""
    new = feature + " \"" + newValue + "\""
    return text.replace(old, new)


def sort(gtf, sortedOut):
    cmd = "sort -k1,1 -k4,4n -k5,5n -o " + sortedOut + " " + gtf

    if subprocess.call(cmd, shell=True) != 0:
        sys.exit("error in " + __file__ + " while executing the following " +
                 "command: " + cmd)


def updateIdsAndPrint(input, sortedInput):
    geneIdCounter = 1
    trIdCounter = 1
    geneIdMap = {}
    trIdMap = {}
    # This map exists just for checking purposes
    gene2tr = {}

    for row in csv.reader(open(sortedInput), delimiter='\t'):
        if row[0][0] == "#":
            continue

        geneId = extractFeatureGtf(row[8], "gene_id")
        trId = extractFeatureGtf(row[8], "transcript_id")

        if geneId not in geneIdMap:
            geneIdMap[geneId] = str(geneIdCounter) + "_g"
            geneIdCounter += 1
            gene2tr[geneId] = trId

        if trId not in trIdMap:
            trIdMap[trId] = str(trIdCounter) + "_t"
            trIdCounter += 1

        if gene2tr[geneId] != trId:
            sys.exit("error in " + __file__ + ": The GeneMark file " + input +
                     " is expected to contain one transcript per gene. More " +
                     "than 1 transcript were detected for gene " + geneId + ".")

        row[8] = replaceFeatureGtf(row[8], "gene_id", geneId,
                                   geneIdMap[geneId])
        row[8] = replaceFeatureGtf(row[8], "transcript_id", trId,
                                   trIdMap[trId])

        print("\t".join(row))


def main():
    args = parseCmd()

    inputDir = os.path.dirname(os.path.realpath(args.input))
    sortedInput = tempfile.NamedTemporaryFile(prefix="genemark", dir=inputDir,
                                              delete=False).name

    sort(args.input, sortedInput)
    updateIdsAndPrint(args.input, sortedInput)
    os.remove(sortedInput)


def parseCmd():

    parser = argparse.ArgumentParser(description='Sort genemark.gtf by\
        coordinates; change gene and transcript IDs to appear in ascending\
        order.')

    parser.add_argument('input', metavar='genemark.gtf', type=str,
                        help='Input GeneMark file in gtf format.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
