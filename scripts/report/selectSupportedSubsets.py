#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Split supplied gene set into different subsets based on the amount of support
# by external evidence.
# ==============================================================


import argparse
import predictionAnalysis as analysis


def main():
    args = parseCmd()

    prediction = analysis.PredictionAnalysis(args.prediction, args.hints)
    prediction.saveSupportedSubsets(args.fullSupport, args.anySupport,
                                    args.noSupport)


def parseCmd():

    parser = argparse.ArgumentParser(description='Split supplied gene set \
        into different subsets based on the amount of support by external \
        evidence.')

    parser.add_argument('prediction', metavar='prediction.gtf', type=str,
                        help='Prediction file.')

    parser.add_argument('hints', metavar='hints.gff', type=str,
                        help='File with external hints.')

    parser.add_argument('--fullSupport', type=str, required=True,
                        help='Output transcripts fully supported by external \
        evidence to this file. All introns in a transcript must be supported \
        by external evidence. In case of single-exon transcripts, both start \
        and stop codon must be supported by external evidence. On top of \
        these criteria, transcripts must be complete.')

    parser.add_argument('--anySupport', type=str, required=True,
                        help='Output transcripts with any external support to \
        this file. At least one intron, start or stop codon of a predicted \
        transcript must be supported.')

    parser.add_argument('--noSupport', type=str, required=True,
                        help='Output transcripts with no external support to \
        this file')

    return parser.parse_args()


if __name__ == '__main__':
    main()
