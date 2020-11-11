#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Collect various statistics about the prediction in gtf format
# ==============================================================

import csv
import re


def extractFeatureGtf(row, feature):
    regex = feature + ' "([^"]+)"'
    result = re.search(regex, row[8])
    if result:
        return result.groups()[0]
    else:
        return None


def getSignature(row):
    return row[0] + "_" + row[3] + "_" + row[4] + "_" + row[6]


class PredictionAnalysis():

    def __init__(self, prediction, hints):
        self.loadHints(hints)
        self.loadPrediction(prediction)

    def loadPrediction(self, prediction):
        self.transcripts = {}
        self.prediction = prediction

        for row in csv.reader(open(prediction), delimiter='\t'):
            transcriptID = extractFeatureGtf(row, "transcript_id")
            if not transcriptID:
                continue
            if transcriptID not in self.transcripts:
                self.transcripts[transcriptID] = Transcript(self)
            self.transcripts[transcriptID].addFeature(row)

        self.collectOverallStatistics()

    def loadHints(self, hints):
        self.hints = hints
        self.intronHints = set()
        self.startHints = set()
        self.stopHints = set()
        for row in csv.reader(open(hints), delimiter='\t'):
            if row[2].lower() == "intron":
                self.intronHints.add(getSignature(row))
            elif row[2].lower() == "start_codon" or row[2].lower() == "start":
                self.startHints.add(getSignature(row))
            elif row[2].lower() == "stop_codon" or row[2].lower() == "stop":
                self.stopHints.add(getSignature(row))

    def collectOverallStatistics(self):
        self.singleTranscriptCount = 0
        self.intronCount = 0
        self.completeCount = 0
        self.fullSupportCount = 0
        self.anySupportCount = 0
        self.transcriptLengths = []
        self.intronLengths = []

        exonTypes = ["single", "initial", "internal", "terminal",
                     "unknown", "all"]
        self.exonLengths = {}
        for exonType in exonTypes:
            self.exonLengths[exonType] = []
        self.exonCounts = []

        for transcript in self.transcripts.values():
            if len(transcript.exons) == 1:
                self.singleTranscriptCount += 1

            self.intronCount += len(transcript.introns)

            if transcript.isComplete():
                self.completeCount += 1

            if transcript.fullSupport():
                self.fullSupportCount += 1

            if transcript.anySupport:
                self.anySupportCount += 1

            transcript.categorizeExons()

            self.transcriptLengths.append(transcript.length)

            for exon in transcript.exons:
                self.exonLengths["all"].append(exon.length)
                self.exonLengths[exon.type].append(exon.length)

            for intron in transcript.introns:
                self.intronLengths.append(intron.length)

            self.exonCounts.append(len(transcript.exons))

    def getTranscriptCount(self):
        return len(self.transcripts)

    def getSingleTranscriptCount(self):
        return self.singleTranscriptCount

    def getMultiTranscriptCount(self):
        return self.getTranscriptCount() - self.getSingleTranscriptCount()

    def getIntronCount(self):
        return self.intronCount

    def getIntronsPerTranscript(self):
        return self.getIntronCount() / self.getTranscriptCount()

    def getIntronsPerMultiTranscript(self):
        return self.getIntronCount() / self.getMultiTranscriptCount()

    def getCompleteCount(self):
        return self.completeCount

    def getIncompleteCount(self):
        return self.getTranscriptCount() - self.getCompleteCount()

    def getTranscriptLengths(self): 
        return self.transcriptLengths

    def getExonLengths(self, exonType):
        return self.exonLengths[exonType]

    def getIntronLengths(self):
        return self.intronLengths

    def getExonsPerTranscript(self):
        return self.exonCounts

    def getFullySupportedTranscriptCount(self):
        return self.fullSupportCount

    def getAnySupportedTranscriptCount(self):
        return self.anySupportCount

    def getUnsupportedTranscriptCount(self):
        return self.getTranscriptCount() - \
            self.getAnySupportedTranscriptCount()


class Feature():
    def __init__(self, row):
        self.support = False
        self.length = int(row[4]) - int(row[3]) + 1
        self.signature = getSignature(row)
        self.chr = row[0]
        self.beginning = int(row[3])
        self.end = int(row[4])
        self.strand = row[6]
        self.row = row
        self.type = None

    def __gt__(self, other):
        return self.beginning > other.beginning


class Transcript():

    def __init__(self, analysis):
        self.predictionAnalysis = analysis
        self.exons = []
        self.introns = []
        self.start = None
        self.stop = None
        self.length = 0
        self.startFound = False
        self.stopFound = False
        self.fullIntronSupport = True
        self.anySupport = False
        self.exonsCategorized = False

    def addFeature(self, row):
        if row[2] == "CDS":
            self.addExon(row)
        elif row[2] == "intron":
            self.addIntron(row)
        elif row[2] == "start_codon":
            self.addStart(row)
        elif row[2] == "stop_codon":
            self.addStop(row)

    def addExon(self, row):
        self.exonsCategorized = False
        exon = Feature(row)
        self.length += exon.length
        self.exons.append(exon)

    def addIntron(self, row):
        intron = Feature(row)
        if intron.signature in self.predictionAnalysis.intronHints:
            self.support = True
            self.anySupport = True
        else:
            self.fullIntronSupport = False
        self.introns.append(intron)

    def addStart(self, row):
        self.start = Feature(row)
        start = Feature(row)
        if start.signature in self.predictionAnalysis.startHints:
            self.anySupport = True
            self.support = True

    def addStop(self, row):
        self.stop = Feature(row)
        stop = Feature(row)
        if stop.signature in self.predictionAnalysis.stopHints:
            self.anySupport = True
            self.support = True

    def isComplete(self):
        return self.start and self.stop

    def fullSupport(self):
        # This is debatable, needs discussion...
        if not self.isComplete():
            return False

        if len(self.introns) == 0:
            if self.start and self.start.support \
               and self.stop and self.stop.support:
                return True
        else:
            return self.fullIntronSupport

        return False

    def categorizeExons(self):
        self.exonsCategorized = True

        if len(self.exons) == 1:
            if self.isComplete():
                self.exons[0].type = "single"
            else:
                self.exons[0].type = "unknown"
            return

        self.exons.sort()

        for i in range(len(self.exons)):
            exon = self.exons[i]

            if i == 0:
                if exon.strand == "+" and self.start:
                    exon.type = "initial"
                elif exon.strand == "-" and self.stop:
                    exon.type = "terminal"
                else:
                    exon.type = "unknown"

            elif i != len(self.exons) - 1:
                exon.type = "internal"

            else:
                if exon.strand == "+" and self.stop:
                    exon.type = "terminal"
                elif exon.strand == "-" and self.start:
                    exon.type = "initial"
                else:
                    exon.type = "unknown"

    def print(self):
        first = self.start
        last = self. stop
        if self.exons[0].strand == "-":
            first = self.stop
            last = self. start

        if first:
            print("\t".join(first.row))

        self.exons.sort()
        for exon in self.exons:
            print("\t".join(exon.row) + ' cds_type "' + exon.type + '";')

        self.introns.sort()
        for intron in self.introns:
            print("\t".join(intron.row))

        if last:
            print("\t".join(last.row))
