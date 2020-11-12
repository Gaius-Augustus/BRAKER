#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Generate BRAKER annotation report
#
# The code for outlier aware histogram comes from
# https://github.com/bdoughty/outlier-aware-histogram
# ==============================================================


import argparse
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import predictionAnalysis as analysis
from matplotlib.backends.backend_pdf import PdfPages
import tempfile
import os
import subprocess


def printBasicStatistics(prediction, report):
    fig = plt.figure()
    fig.text(0.5, 0.9, 'BRAKER Prediction Report', ha='center', va='center',
             size=20)

    fig.text(0.1, 0.8, 'All genes in this report are represented by the '
             'longest coding transcript.', ha='left', va='top', size=9)

    fig.text(0.1, 0.72, 'General Statistics', ha='left', va='top', size=15)

    geneCount = prediction.getTranscriptCount()
    anySupport = prediction.getAnySupportedTranscriptCount()
    fullSupport = prediction.getFullySupportedTranscriptCount()
    noSupport = prediction.getUnsupportedTranscriptCount()
    complete = prediction.getCompleteCount()
    partial = prediction.getIncompleteCount()

    text = 'Gene count: ' + str(prediction.getTranscriptCount()) + "\n"
    text += '    Single-exon genes: ' + \
            str(prediction.getSingleTranscriptCount()) + "\n"
    text += '    Multi-exon genes: ' + \
            str(prediction.getMultiTranscriptCount()) + "\n"
    text += "\n"
    text += 'Introns per gene: ' + \
            str(round(prediction.getIntronsPerTranscript(), 2)) + "\n"
    text += 'Introns per multi-exon gene: ' + \
            str(round(prediction.getIntronsPerMultiTranscript(), 2)) + "\n"
    text += "\n"
    text += "Genes fully supported by external evidence: " + str(fullSupport) \
            + " (" + str(round(100 * fullSupport / geneCount, 2)) + "%)\n"
    text += "Genes partially supported by external evidence: " + \
            str(anySupport) + " (" + \
            str(round(100 * anySupport / geneCount, 2)) + "%)\n"
    text += "Genes unsupported by any external evidence: " + str(noSupport) + \
            " (" + str(round(100 * noSupport / geneCount, 2)) + "%)\n"
    text += "\n"
    text += "Complete genes: " + str(complete) + \
            " (" + str(round(100 * complete / geneCount, 2)) + "%)\n"
    text += "Partial genes: " + str(partial) + \
            " (" + str(round(100 * partial / geneCount, 2)) + "%)\n"

    fig.text(0.1, 0.64, text, ha='left', va='top', size=9, linespacing=1.5)

    report.savefig(fig)


def getBounds(data, zScore):
    std = np.std(data)
    median = np.median(data)
    return (median + zScore * std)


def histogram(data, report, title, xlabel, zScore,
              minimum=-1, maximum=-1, bins=-1):
    data = np.asarray(data)

    if maximum == -1:
        upper = getBounds(data, zScore)
    else:
        upper = maximum

    if upper > data.max():
        upper = data.max()
        upper_outliers = False
    else:
        upper_outliers = True

    fig = plt.figure()

    color = "tab:blue"

    if bins == -1:
        bins = 'auto'
        histData = data
    elif bins == "single":
        upper = int(upper)

        # The two lines below deal with the inclusivity of the last boundary
        # (https://matplotlib.org/3.3.2/api/_as_gen/matplotlib.pyplot.hist.html)
        bins = range(1, upper + 2)
        histData = data[data != upper + 1]

        plt.xticks(bins[:-1])
        color = "forestgreen"
    else:
        sys.error("Error: Unexpected bins argument: " + bins)

    # Note that range has no effect if "bins" is a sequence
    if minimum != -1:
        n, bins, patches = plt.hist(histData, range=(minimum, maximum),
                                    bins=bins, align='left')
    else:
        n, bins, patches = plt.hist(histData, range=(data.min(), upper),
                                    bins=bins, align='left', color=color)

    if upper_outliers:
        outliersSum = (data > upper).sum()
        lastColumn = patches[-1].get_height() + outliersSum

        if plt.ylim()[1] < lastColumn + lastColumn / 20:
            plt.ylim(0, lastColumn + lastColumn / 20)

        patches[-1].set_height(lastColumn)
        patches[-1].set_facecolor('m')
        patches[-1].set_label('Range of upper outliers: (' + str(int(upper))
                              + ', ' + str(data.max()) + ')')

    if upper_outliers:
        plt.legend()

    ax = plt.gca()
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    # plt.text(10, 20, '* Genes are represented by\n' +
    #          '   the longest coding transcript',
    #          ha='left', va='top', size=7, transform=None)
    report.savefig(fig)
    return data.min(), upper


def printGeneHistogram(prediction, report):
    data = prediction.getTranscriptLengths()
    histogram(data, report, "Histogram of gene lengths", "Gene length", 3)


def printExonHistograms(prediction, report):
    minimum, maximum = histogram(prediction.getExonLengths("all"), report,
                                 "Histogram of exon lengths", "Exon length", 3)
    histogram(prediction.getExonLengths("initial"), report,
              "Histogram of intial exon lengths", "Initial exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengths("internal"), report,
              "Histogram of internal exon lengths", "Internal exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengths("terminal"), report,
              "Histogram of terminal exon lengths", "Terminal exon length", 3,
              minimum, maximum)
    histogram(prediction.getExonLengths("single"), report,
              "Histogram of single exon lengths", "Single exon length", 3,
              minimum, maximum)


def printIntronHistogram(prediction, report):
    data = np.asarray(prediction.getIntronLengths())
    z = 3
    if data.max() > 10000:
        z = 1
    histogram(data, report, "Histogram of intron lengths",
              "Intron length", z)


def printExonsPerGene(prediction, report):
    data = prediction.getExonsPerTranscript()
    histogram(data, report, "Exons per gene", "Exon number", 5,
              bins='single')


def main():
    args = parseCmd()
    report = PdfPages(args.output)

    longestIsoforms = tempfile.NamedTemporaryFile(prefix="longest",
                                                  delete=False).name
    returnVal = subprocess.call(os.path.dirname(os.path.realpath(__file__)) +
                                "/printLongestIsoforms.py " + args.prediction +
                                " > " + longestIsoforms, shell=True)
    if returnVal != 0:
        os.remove(longestIsoforms)
        sys.exit("Error: Something went wrong during the selection of " +
                 "longest coding isoforms.")

    prediction = analysis.PredictionAnalysis(longestIsoforms, args.hints)
    printBasicStatistics(prediction, report)
    printGeneHistogram(prediction, report)
    printExonHistograms(prediction, report)
    printExonsPerGene(prediction, report)
    printIntronHistogram(prediction, report)

    os.remove(longestIsoforms)
    report.close()


def parseCmd():

    parser = argparse.ArgumentParser(description='Generate BRAKER \
        annotation report.')

    parser.add_argument('prediction', metavar='prediction.gtf', type=str,
                        help='BRAKER prediction file.')

    parser.add_argument('hints', metavar='hints.gff', type=str,
                        help='File with external hints.')

    parser.add_argument('output', metavar='output.pdf', type=str,
                        help='Name of the output pdf report.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
