#!/usr/bin/env python3

####################################################################################################
#                                                                                                  #
# gmst2globalCoords.pl                                                                             #
# Script to transform transcript level coordinates of gene predictions made by GeneMarkS-T on      #
# transcripts to genome level coordinates and print the genes to a file in GTF-format              #
#                                                                                                  #
# Authors: Hannah Thierfeldt, Lars Gabriel                                                         #
#                                                                                                  #
# Credits : Mario Stanke, Katharina Hoff                                                           #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# Version: 1.0.0                                                                                   #
#                                                                                                  #
# Copyright 2019-2021. All rights reserved.                                                        #
####################################################################################################

import argparse
import re
import sys
import csv

# author attributes
__author__ = "Hannah Thierfeldt, Lars Gabriel"
__copyright__ = "Copyright 2019-2021. All rights reserved."
__credits__ = "Mario Stanke, Katharina Hoff"
__license__ = "Artistic License"
__version__ = "1.0.0"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

parser = argparse.ArgumentParser(description='Script to transform transcript level coordinates of gene predictions' +
                    'made by GeneMarkS-T on transcripts to genome level coordinates and print the genes to a ' +
                    'file in GTF-format')
parser.add_argument('-t', '--transcripts', required=True, type=str,
                    help='GTF-file of (StringTie) transcripts the GeneMarkS-T predictions were made on')
parser.add_argument('-p', '--predictions', required=True, type=str,
                    help='GFF-file with GeneMarkS-T predictions made on (StringTie) transcripts')
parser.add_argument('-r', '--remove_alt_tx', required=False, action='store_true',
                    help='When predictions were made on several alternative transcripts of a transcript cluster, ' +
                    'print only the gene with the most introns')
parser.add_argument('-g', '--genome', required=False, type=str,
                    help='Fasta file of the genome used for the GeneMarkS-T predictions. ' +
                    'If a genome file is provided, only transcripts whose coding sequence ' +
                    'starts with a start-codon are printed.')
parser.add_argument('-o', '--outfile', required=False, type=str, default='traingenes.gtf',
                    help='Output file in GTF format')
args = parser.parse_args()

# read transcripts
# dict with gene id as key and a 2nd dict as value
# 2nd dict has transcript ids as key and a 3rd dict as value
tx = {}
try:
    with open(args.transcripts, 'r') as file_handle:
        for line in file_handle:
            match = re.match(r'(\S+)\s+\S+\s+exon\s+(\d+)\s+(\d+)\s+\S+\s+(\S)\s+(\S)\s+(.*)',
                    line)
            if match:
                seq, start, stop, strand, frame, attributes = match.groups()
                if not ('gene_id' in attributes or 'transcript_id' in attributes):
                    continue
                gene, = re.match(r'.*gene_id\s+"(\S+)";.*', attributes).groups()
                tx_id, = re.match(r'.*transcript_id\s+"(\S+)";.*', attributes).groups()
                if gene not in tx:  # add new gene
                    tx[gene] = {'seq': seq}
                if tx_id not in tx[gene]:  # add new transcript
                    tx[gene][tx_id] = {'strand': strand, 'frame': frame, 'exons': []}
                tx[gene][tx_id]['exons'].append([int(start), int(stop)])  # add new exon
except IOError:
    sys.stderr.write('Error: Failed to open ' + args.transcripts + ' for reading!\n')
    quit(1)

# read predictions
# dict with gene id as key and a 2nd dict as value
# 2nd dict has transcript ids as key and a 3rd dict as value
pred_gene = {}
try:
    with open(args.predictions, 'r') as file_handle:
        for line in file_handle:
            if re.match(r'(\S+)\.(\S+)\s+\S+\s+CDS\s+(\d+)\s+(\d+)\s+(\S+)\s+\S\s+(\S)', line):
                # only one prediction per transcript
                gene, tx_id, start, stop, score, frame = re.match(
                    r'(\S+)\.(\S+)\s+\S+\s+CDS\s+(\d+)\s+(\d+)\s+(\S+)\s+\S\s+(\S)', line).groups()
                if gene not in pred_gene:  # add new gene
                    pred_gene[gene] = {}
                # add prediction on transcript
                pred_gene[gene][gene+'.'+tx_id] = {'start': int(start), 'stop': int(stop), 'score': score,
                                                   'frame': frame}
except IOError:
    sys.stderr.write('Error: Failed to open ' + args.predictions + ' for reading!\n')
    quit(1)

# find genome level coordinates
for gene in pred_gene:
    for tx_id in pred_gene[gene]:
        # strand on genome level
        if tx[gene][tx_id]['strand'] == '-':
            g_st = '-'
            # start = stop coord of transcript - stop coord of prediction + 1
            start = tx[gene][tx_id]['exons'][-1][1] - pred_gene[gene][tx_id]['stop'] + 1
            # stop = stop coord of transcript - start coord of prediction + 1
            stop = tx[gene][tx_id]['exons'][-1][1] - pred_gene[gene][tx_id]['start'] + 1

            intron_len = 0
            # start with last exon and go to first
            for t_ind in range(len(tx[gene][tx_id]['exons']) - 1, -1, -1):
                if tx[gene][tx_id]['exons'][t_ind][0] <= stop + intron_len <= tx[gene][tx_id]['exons'][t_ind][1]:
                    g_stop = stop + intron_len
                    break
                else:  # add intron length to bridge missing stretch
                    intron_len += tx[gene][tx_id]['exons'][t_ind - 1][1] - tx[gene][tx_id]['exons'][t_ind][0] + 1

            # compute stop coord
            t_start = t_ind  # start loop with last exon of prev loop
            for t_ind in range(t_start, -1, -1):
                if tx[gene][tx_id]['exons'][t_ind][0] <= start + intron_len <= tx[gene][tx_id]['exons'][t_ind][1]:
                    g_start = start + intron_len
                    break
                else:  # add intron length to bridge missing stretch
                    intron_len += tx[gene][tx_id]['exons'][t_ind - 1][1] - tx[gene][tx_id]['exons'][t_ind][0] + 1

            cds_f = t_ind
            cds_l = t_start

        else:
            g_st = '+'
            start = pred_gene[gene][tx_id]['start'] + tx[gene][tx_id]['exons'][0][0] - 1
            stop = pred_gene[gene][tx_id]['stop'] + tx[gene][tx_id]['exons'][0][0] - 1

            # compute start coord
            # start with first exon
            intron_len = 0
            for t_ind in range(len(tx[gene][tx_id]['exons'])):
                if tx[gene][tx_id]['exons'][t_ind][0] <= start + intron_len <= tx[gene][tx_id]['exons'][t_ind][1]:
                    g_start = start + intron_len
                    break
                else:  # add intron length to bridge missing stretch
                    intron_len += tx[gene][tx_id]['exons'][t_ind + 1][0] - tx[gene][tx_id]['exons'][t_ind][1] - 1

            # compute stop coord
            t_start = t_ind  # start loop with last exon of prev loop
            for t_ind in range(t_start, len(tx[gene][tx_id]['exons'])):
                if tx[gene][tx_id]['exons'][t_ind][0] <= stop + intron_len <= tx[gene][tx_id]['exons'][t_ind][1]:
                    g_stop = stop + intron_len
                    break
                else:  # add intron length to bridge missing stretch
                    intron_len += tx[gene][tx_id]['exons'][t_ind + 1][0] - tx[gene][tx_id]['exons'][t_ind][1] - 1
            cds_f = t_start
            cds_l = t_ind

        # add genome level strand, start and stop coord, as well as indices of exons in CDS
        pred_gene[gene][tx_id]['genome_start'] = g_start
        pred_gene[gene][tx_id]['genome_stop'] = g_stop
        pred_gene[gene][tx_id]['1stCDS'] = cds_f
        pred_gene[gene][tx_id]['lastCDS'] = cds_l
        pred_gene[gene][tx_id]['genome_strand'] = g_st

# when predictions were made on several alternative transcripts of a transcript cluster
# only keep the prediction that includes the most introns
if args.remove_alt_tx:
    for gene in pred_gene:
        if len(pred_gene[gene]) > 1:  # predictions on alternative transcripts
            nr_intr = {}  # save number of introns
            for tx_id in pred_gene[gene]:
                nr_intr[tx_id] = pred_gene[gene][tx_id]['lastCDS'] - pred_gene[gene][tx_id]['1stCDS']

            max_intr = max(nr_intr, key=nr_intr.get)  # get gene with most introns

            remove = []
            for tx_id in pred_gene[gene]:  # get entry keys to remove
                if tx_id != max_intr:
                    remove.append(tx_id)
            for tx_id in remove:  # remove entries
                pred_gene[gene].pop(tx_id)

# remove all transcripts whose first CDS does not start with a start-codon

# ids of transcripts that are fragmented (1st CDS doesn't start with a start-codon)
fragmented_transcripts = []
if args.genome:
    def make_complement(string):
        complement_dict = {'a' : 't', 't' : 'a', 'u' : 'a', 'g': 'c', 'c' : 'g'}
        return ''.join([complement_dict[s] for s in string[::-1]])
    # dict with seq names as keys and the sequence as values
    start_codon = ['atg', 'aug']

    # prepare a dict with first codon of each transcript for iterating through genome file
    transcript_start_codons = {}
    for gene in pred_gene:
        if tx[gene]['seq'] not in transcript_start_codons.keys():
            transcript_start_codons.update({ tx[gene]['seq'] : [] })
        for tx_id in pred_gene[gene]:
            if pred_gene[gene][tx_id]['genome_strand'] == '+':
                start = pred_gene[gene][tx_id]['genome_start'] - 1
            else:
                start = pred_gene[gene][tx_id]['genome_stop'] - 3
            transcript_start_codons[tx[gene]['seq']].append([start, gene, tx_id, pred_gene[gene][tx_id]['genome_strand']])
    # read genome
    with open(args.genome, 'r') as file_handle:
        seq_name = ''
        seq = ''
        line = file_handle.readline()
        while True:
            while line and not line[0] == '>':
                seq += line.strip('\n')
                line = file_handle.readline()
            if seq_name in transcript_start_codons.keys():
                for tx_start in transcript_start_codons[seq_name]:
                    if tx_start[0] + 3 > len(seq):
                        sys.stderr.write('### Genome sequence {} is shorter than the '.format(seq_name)
                            + 'start position of transcript {},{},{}.').format(*tx_start)
                    codon = seq[tx_start[0]:tx_start[0] + 3].lower()
                    if tx_start[3] == '-':
                        codon = make_complement(codon)
                    if codon not in start_codon:
                        fragmented_transcripts.append(tx_start[1:3])
            if not line:
                break
            seq_name = line.split('>')[1].split('\n')[0].strip()
            seq = ''
            line = file_handle.readline()
else :
    sys.stderr.write('### It is assumed that coding sequences of all transcripts \n'
        + '### that were predicted by GeneMarkS-T start with a start-codon.\n'
        + '### The start of the coding sequences will not be check as no genome file is provided.\n')

# write genes to outfile
try:
    with open(args.outfile, 'w') as file_handle:
        out_writer = csv.writer(file_handle, delimiter='\t', quotechar = "|", lineterminator = '\n')
        for gene in pred_gene:
            for tx_id in pred_gene[gene]:
                if [gene, tx_id] in fragmented_transcripts:
                    continue
                # exons with CDS
                cds_f = pred_gene[gene][tx_id]['1stCDS']
                cds_l = pred_gene[gene][tx_id]['lastCDS']
                iteration_list = range(len(tx[gene][tx_id]['exons']))
                current_frame = int(pred_gene[gene][tx_id]['frame'])
                tx_lines = []
                # UTR directions for strand
                if pred_gene[gene][tx_id]['genome_strand'] == '+':
                    start_utr = "5'-UTR"
                    stop_utr = "3'-UTR"
                else:
                    start_utr = "3'-UTR"
                    stop_utr = "5'-UTR"
                    iteration_list = reversed(iteration_list)

                for ex in iteration_list:
                    if ex == cds_f:
                        if cds_f == cds_l:  # only one CDS
                            # if CDS does not begin at start coord of exon
                            if pred_gene[gene][tx_id]['genome_start'] > tx[gene][tx_id]['exons'][ex][0]:
                                # UTR on exon before CDS
                                # start of exon to start of CDS-1
                                tx_lines.append([tx[gene]['seq'], 'manual', start_utr,
                                        tx[gene][tx_id]['exons'][ex][0],
                                        pred_gene[gene][tx_id]['genome_start']-1, '.',
                                        pred_gene[gene][tx_id]['genome_strand'], '.',
                                        'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                            # CDS
                            tx_lines.append([tx[gene]['seq'], 'GeneMark.hmm', 'CDS',
                                pred_gene[gene][tx_id]['genome_start'],
                                pred_gene[gene][tx_id]['genome_stop'], '.',
                                pred_gene[gene][tx_id]['genome_strand'], current_frame,
                                'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                            # if CDS does not stop at stop coord of exon
                            if pred_gene[gene][tx_id]['genome_stop'] < tx[gene][tx_id]['exons'][ex][1]:
                                # UTR on exon after CDS
                                # end of CDS+1 to end of exon
                                tx_lines.append([tx[gene]['seq'], 'manual', stop_utr,
                                        pred_gene[gene][tx_id]['genome_stop']+1,
                                        tx[gene][tx_id]['exons'][ex][1], '.',
                                        pred_gene[gene][tx_id]['genome_strand'], '.',
                                        'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])

                        else:  # first CDS of multiple
                            # if CDS does not start at first coord of exon
                            if pred_gene[gene][tx_id]['genome_start'] > tx[gene][tx_id]['exons'][ex][0]:
                                # UTR on exon before CDS
                                # start of exon to start of CDS-1
                                tx_lines.append([tx[gene]['seq'], 'manual', start_utr,
                                        tx[gene][tx_id]['exons'][ex][0],
                                        pred_gene[gene][tx_id]['genome_start']-1, '.',
                                        pred_gene[gene][tx_id]['genome_strand'], '.',
                                        'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                            # CDS
                            tx_lines.append([tx[gene]['seq'], 'GeneMark.hmm', 'CDS',
                                    pred_gene[gene][tx_id]['genome_start'],
                                    tx[gene][tx_id]['exons'][ex][1], '.',
                                    pred_gene[gene][tx_id]['genome_strand'], current_frame,
                                    'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                            current_frame = (3 - (tx[gene][tx_id]['exons'][ex][1] - \
                                pred_gene[gene][tx_id]['genome_start'] + 1 - current_frame) % 3) %3

                    elif ex == cds_l:  # last CDS of multiple
                        # CDS
                        tx_lines.append([tx[gene]['seq'], 'GeneMark.hmm', 'CDS',
                                tx[gene][tx_id]['exons'][ex][0],
                                pred_gene[gene][tx_id]['genome_stop'], '.',
                                pred_gene[gene][tx_id]['genome_strand'], current_frame,
                                'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                        current_frame = (3 - (pred_gene[gene][tx_id]['genome_stop'] - \
                            tx[gene][tx_id]['exons'][ex][0] + 1 - current_frame) % 3) %3

                        # if CDS does not stop at stop coord of exon
                        if pred_gene[gene][tx_id]['genome_stop'] < tx[gene][tx_id]['exons'][ex][1]:
                            # UTR on exon after CDS
                            # end of CDS+1 to end of exon
                            tx_lines.append([tx[gene]['seq'], 'manual', stop_utr,
                                    pred_gene[gene][tx_id]['genome_stop'] + 1,
                                    tx[gene][tx_id]['exons'][ex][1], '.',
                                    pred_gene[gene][tx_id]['genome_strand'], '.',
                                    'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])

                    elif cds_f < ex < cds_l:  # internal CDS
                        tx_lines.append([tx[gene]['seq'], 'GeneMark.hmm', 'CDS',
                                tx[gene][tx_id]['exons'][ex][0], tx[gene][tx_id]['exons'][ex][1],
                                '.', pred_gene[gene][tx_id]['genome_strand'],
                                current_frame, 'gene_id "' + gene + '"; transcript_id "' +
                                tx_id + '";'])
                        current_frame = (3 - (tx[gene][tx_id]['exons'][ex][1] - \
                            tx[gene][tx_id]['exons'][ex][0] + 1 - current_frame) % 3) %3

                    elif ex < cds_f:  # UTR before CDS
                        tx_lines.append([tx[gene]['seq'], 'manual', start_utr,
                                tx[gene][tx_id]['exons'][ex][0],
                                tx[gene][tx_id]['exons'][ex][1], '.',
                                pred_gene[gene][tx_id]['genome_strand'], '.',
                                'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])

                    else:  # UTR after CDS
                        tx_lines.append([tx[gene]['seq'], 'manual', stop_utr,
                                tx[gene][tx_id]['exons'][ex][0],
                                tx[gene][tx_id]['exons'][ex][1], '.',
                                pred_gene[gene][tx_id]['genome_strand'], '.',
                                'gene_id "' + gene + '"; transcript_id "' + tx_id + '";'])
                tx_lines.sort(key=lambda t:(t[3],t[4]))
                for line in tx_lines:
                    out_writer.writerow(line)

except IOError:
    print('Error: Failed to open ' + args.outfile + ' for writing!')
