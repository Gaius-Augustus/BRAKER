#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
#
# Ensure that the file with "good" GeneMark training genes contains at least N
# genes. If the file contains < N genes, add additional genes from the file with
# "bad" genes. The bad genes are added based on their support by hints: Genes
# which are supported by most hints (normalized by the number of exons) are
# added first.
# ==============================================================


import argparse
import csv
import re
import sys
from collections import OrderedDict


def extractFeature(text, feature):
    regex = feature + ' "([^"]+)"'
    search = re.search(regex, text)
    if search:
        return search.groups()[0]
    else:
        return None


class Gene:

    def __init__(self, row, geneId):
        self.rows = [row]
        self.anchoredCount = 0
        self.anchoredSum = 0
        self.geneId = geneId
        self.countSupport(row)

    def addRow(self, row):
        self.countSupport(row)
        self.rows.append(row)

    def countSupport(self, row):
        if row[2] != "CDS":
            return

        evidence = extractFeature(row[8], "anchored")
        if evidence == "1_1":
            self.anchoredSum += 2
        elif evidence is not None:
            self.anchoredSum += 1
        self.anchoredCount += 1

    def getSupportRatio(self):
        return float(self.anchoredSum) / self.anchoredCount

    def __lt__(self, other):
        return self.getSupportRatio() > other.getSupportRatio()


def loadGenes(genesFile):
    genes = OrderedDict()
    for row in csv.reader(open(genesFile), delimiter='\t'):
        geneId = extractFeature(row[8], "gene_id")
        if geneId not in genes:
            genes[geneId] = Gene(row, geneId)
        else:
            genes[geneId].addRow(row)
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

    sortedBad = sorted(badGenes.values())
    for i in range(toAdd):
        goodGenes[sortedBad[i].geneId] = badGenes.pop(sortedBad[i].geneId)


def printGenes(genes, output):
    f = open(output, 'w')
    for key in genes:
        for row in genes[key].rows:
            f.write("\t".join(row) + "\n")
    f.close()


def main():
    args = parseCmd()

    goodGenes = loadGenes(args.goodGenes)
    badGenes = loadGenes(args.badGenes)

    ensureNGoodGenes(goodGenes, badGenes, args.N)

    printGenes(goodGenes, args.goodGenes)
    printGenes(badGenes, args.badGenes)


def parseCmd():

    parser = argparse.ArgumentParser(description='Ensure that the file with\
                                     "good" GeneMark training genes contains\
                                     at least N genes. If the file contains\
                                     < N genes, add additional genes from the\
                                     file with "bad" genes. The bad genes are\
                                     added based on their support by hints:\
                                     Genes which are supported by most hints\
                                     (normalized by the number of exons) are \
                                     added first.')

    parser.add_argument('--goodGenes', type=str, required=True,
                        help='"Good" genemark training genes')

    parser.add_argument('--badGenes', type=str, required=True,
                        help='"Bad" genemark training genes')

    parser.add_argument('--N', type=int, required=True,
                        help='Minimum required number of good genes.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
