#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Print longest coding isoforms in a gtf file
# ==============================================================


import csv
import sys
import re
import argparse


def extractFeature(text, feature):
    regex = feature + ' "([^"]+)"'
    result = re.search(regex, text)
    if result:
        return result.groups()[0]
    else:
        return None


def computeLengths(input):
    transcriptLengths = dict()
    for row in csv.reader(open(input), delimiter='\t'):
        if (row[2] == 'CDS'):
            gene = extractFeature(row[8], 'gene_id')
            transcript = extractFeature(row[8], 'transcript_id')
            if not gene or not transcript:
                continue
            if gene not in transcriptLengths:
                transcriptLengths[gene] = dict()
            if transcript not in transcriptLengths[gene]:
                transcriptLengths[gene][transcript] = 0
            transcriptLengths[gene][transcript] += int(row[4]) - int(row[3])
    return transcriptLengths


def getLongestTranscript(transcriptLengths):
    longestTranscripts = dict()
    for gene in transcriptLengths:
        max = 0
        longestTranscript = ""
        for transcript in transcriptLengths[gene]:
            length = transcriptLengths[gene][transcript]
            if (length > max):
                max = length
                longestTranscript = transcript
        longestTranscripts[gene] = longestTranscript
    return longestTranscripts


def printLongest(input, longestTranscripts):
    for row in csv.reader(open(input), delimiter='\t'):
        gene = extractFeature(row[8], 'gene_id')
        transcript = extractFeature(row[8], 'transcript_id')
        if not gene or not transcript:
            continue
        if (longestTranscripts[gene] == transcript):
            print('\t'.join(row))


def main():
    args = parseCmd()
    transcriptLengths = computeLengths(args.input)
    longestTranscripts = getLongestTranscript(transcriptLengths)
    printLongest(sys.argv[1], longestTranscripts)


def parseCmd():

    parser = argparse.ArgumentParser(description='Print longest isoforms in a \
                                     gtf file')

    parser.add_argument('input', type=str,
                        help='Input gtf file')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
