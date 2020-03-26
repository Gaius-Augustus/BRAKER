#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Ensure that the file with "good" GeneMark training genes contains at least N
# genes. If the file contains < N genes, add random genes from the file with
# "bad" genes.
# ==============================================================


import argparse
import csv
import re
import sys
import random
from collections import OrderedDict


def extractFeature(text, feature):
    regex = feature + ' "([^"]+)"'
    return re.search(regex, text).groups()[0]


def loadGenes(genesFile):
    genes = OrderedDict()
    for row in csv.reader(open(genesFile), delimiter='\t'):
        geneId = extractFeature(row[8], "gene_id")
        if geneId not in genes:
            genes[geneId] = [row]
        else:
            genes[geneId].append(row)
    return genes


def ensureNGoodGenes(goodGenes, badGenes, N):
    goodN = len(goodGenes)
    if goodN >= N:
        return

    toAdd = N - goodN
    if toAdd > len(badGenes):
        scriptName = sys.argv[0][(sys.argv[0].rfind("/") + 1):]
        sys.stderr.write(scriptName + ': warning: File with "bad" genes does' +
                         ' not have enough genes to ensure ' + str(N) +
                         ' genes. Adding all "bad" genes to the "good" gene' +
                         ' set.\n')
        toAdd = len(badGenes)

    for i in range(toAdd):
        geneId = random.choice(list(badGenes.keys()))
        goodGenes[geneId] = badGenes.pop(geneId)


def printGenes(genes, output):
    f = open(output, 'w')
    for key in genes:
        for row in genes[key]:
            f.write("\t".join(row) + "\n")
    f.close()


def main():
    args = parseCmd()

    if args.randomSeed:
        random.seed(args.randomSeed)

    goodGenes = loadGenes(args.goodGenes)
    badGenes = loadGenes(args.badGenes)

    ensureNGoodGenes(goodGenes, badGenes, args.N)

    printGenes(goodGenes, args.goodGenes)
    printGenes(badGenes, args.badGenes)


def parseCmd():

    parser = argparse.ArgumentParser(description='Ensure that the file with\
                                     "good" GeneMark training genes contains\
                                     at least N genes. If the file contains\
                                     < N genes, add random genes from the file\
                                     with "bad" genes.')

    parser.add_argument('--goodGenes', type=str, required=True,
                        help='"Good" genemark training genes')

    parser.add_argument('--badGenes', type=str, required=True,
                        help='"Bad" genemark training genes')

    parser.add_argument('--N', type=int, required=True,
                        help='Minimum required number of good genes.')

    parser.add_argument('--randomSeed', type=int, required=False,
                        help='Use this random seed for adding genes from\
                        "bad" file.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
