#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Get range of GC content for a file of DNA sequences.
# ==============================================================

import argparse
import os
import sys
import re
import importlib.util

class FileDoesNotExist(Exception):
    pass

class SequenceFileError(Exception):
    pass

class InputError(Exception):
    pass

config = {
    'percentiles' : [5,95],
    'window_size' : 10000,
    'draw' : False,
    'draw_size' : 80 ,
    'draw_height' : 10 ,
    'delimiter' : ' ',
    'print_sequence_length' : False,
    'mem_size' : 4e+9
}

numpy_imported = False
if importlib.util.find_spec('numpy') is not None:
    import numpy as np
    numpy_imported = True

def main():
    args = parseCmd()
    seq_file_path = make_path_abs(args.sequences)
    set_config(args)

    # dict with seq name as key and gc content of sequence slices as values
    gc_content = {}
    seq_name = ''
    len_last_slice = 0
    last_slice = 0
    total_genome_length = 0

    with open(seq_file_path, 'r') as seq_file:
        text = seq_file.read(int(config['mem_size']))
        while text:
            text = text.split('>')
            for j, seq in enumerate(text):
                if j > 0:
                    if len_last_slice > 0.2*config['window_size']:
                        gc_content[seq_name].append(last_slice/len_last_slice)
                    total_genome_length += len_last_slice
                    len_last_slice = 0
                    last_slice = 0
                    seq_name_index = seq.index('\n')
                    seq_name = seq[:seq_name_index]
                    seq = seq[seq_name_index:]
                    if seq_name in gc_content:
                        raise SequenceFileError(f'Seq name {seq_name} is not unique!')
                    else:
                        gc_content.update({seq_name : []})
                seq = seq.replace('\n', '').lower().replace('n', '')
                intendation = config['window_size'] - int(last_slice)
                slices = [(last_slice + len(re.findall('[gc]', \
                    seq[0:intendation]))) / config['window_size']] \
                    + [len(re.findall('[gc]', seq[i:i+config['window_size']])) \
                    / config['window_size'] for i in range(intendation, \
                    len(seq), config['window_size'])]
                last_slice = 0

                if len(slices) > 1:
                    len_last_slice = len(seq) % config['window_size']
                    if len_last_slice > 0:
                        last_slice = slices.pop() * config['window_size']
                    total_genome_length += config['window_size'] * len(slices)
                    gc_content[seq_name] += slices
            text = seq_file.read(int(config['mem_size']))
    print(f'sequence_length:{config["delimiter"]}{total_genome_length}')
    print_gc_range(gc_content)

    if config['draw']:
        draw_seq_content(gc_content)

def draw_seq_content(content_dic):
    draw_size = min([config['draw_size']-3] + [len(c) for c in content_dic.values()])
    percentiles = {}
    min_axis = 101.0
    max_axis = -1.0
    for seq in content_dic:
        iter_range = [round(i*len(content_dic[seq])/(draw_size)) for i in range(0, draw_size)]
        percentiles.update({seq : []})
        for i in range(1, len(iter_range)):
            percentiles[seq].append(sum(content_dic[seq][iter_range[i-1]:iter_range[i]])\
                /(iter_range[i]-iter_range[i-1]))
        min_axis = min([min_axis] + percentiles[seq])
        max_axis = max([max_axis] + percentiles[seq])
    min_axis = int(min_axis*100)/100
    max_axis = int(max_axis*100)/100
    step_size = (max_axis-min_axis)/(config['draw_height']-1)
    if numpy_imported:
        iter_array = np.linspace(max_axis, min_axis, config['draw_height'])
    else:
        iter_array =[i*step_size+min_axis for i in range(config['draw_height']-1,-1,-1)]
    for seq in content_dic:
        print(f'GC-content distribution for {seq}')
        for k, i in enumerate(iter_array):
            if k == 0 or k == config['draw_height']-1:
                out_str = '{:<3}|'.format(int(i*100))
            else:
                out_str = '{:<3}|'.format(' ')
            for p in percentiles[seq]:
                if p>=i:
                    out_str += '-'
                else:
                    if k==config['draw_height']-1:
                        print(p)
                    out_str += ' '
            print(out_str)

def print_gc_range(content_dic):
    all_slices = []
    for g in content_dic.values():
        all_slices += g
    if not numpy_imported:
        all_slices.sort()
    if all_slices:
        gc_range = []
        for p in config['percentiles']:
            gc_range.append(get_percentile(all_slices, p))
            print(f'{p}-percentile:{config["delimiter"]}{gc_range[-1]}')

def get_percentile(array, q):
    if numpy_imported:
        return np.percentile(array, q)
    else:
        return array[int(q*(len(array)-1)/100)]


def check_percentile_range(q):
    if not (q>=0 and q<=100):
        raise InputError(f'The percentile {q} hast to be between 0 and 100!')

def set_config(args):
    global config
    if not args.percentiles == None:
        config['percentiles'] = []
        for p in args.percentiles.split(','):
            p = int(p)
            check_percentile_range(p)
            config['percentiles'].append(p)
    if args.window_size:
        config['window_size'] = args.window_size
    if args.print_sequence_length:
        config['print_sequence_length'] = args.print_sequence_length
    if args.draw_gc_content:
        config['draw'] = args.draw_gc_content
    if args.drawing_length:
        config['draw_size'] = args.drawing_length
    if args.drawing_height:
        config['draw_height'] = args.drawing_height
    if args.delimiter:
        config['delimiter'] = args.delimiter
    if args.max_mem_size:
        config['mem_size'] = args.max_mem_size * 1e+9

def make_path_abs(path_to_file):
    abs_path = os.path.abspath(path_to_file)
    if not os.path.exists(abs_path):
        raise FileDoesNotExist(f'Path does not exist: {abs_path} .')
    return abs_path

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='Compute the range of GC content '\
        'of the sequences in a FASTA file. The fractions of Guanin and Cytosin in '\
        'sequence slices of window_size are computed from the input file. The output '\
        'is the lower_percentile (default 5th precentile) and higher_percentile (default '\
        '95th quantile) of the computed GC content fractions.')
    parser.add_argument('--sequences', type=str, required=True,
        help='FASTA file of genomic sequences.')
    parser.add_argument('--window_size', type=int,
        help='Size of the sequence slices. (default: 10000)')
    parser.add_argument('--percentiles', type=str,
        help='List of percentiles that are printed to stdout.' \
            'A value k of a percentile has to be 0<=k<=100 and the values are '\
            'separated by a comma. (default: 5,95)')
    parser.add_argument('--print_sequence_length', action='store_true',
        help='Prints the total number of nucleotides in any sequence. (default False)')
    parser.add_argument('--draw_gc_content', action='store_true',
        help='Draws the distribution of the GC content of each input sequence to stdout. (default False)')
    parser.add_argument('--drawing_length', type=int,
        help='Number of characters in one row of the drawing, if --draw_gc_content is used. (default 80)')
    parser.add_argument('--drawing_height', type=int,
        help='Height (in rows) of the drawing for one sequence, if --draw_gc_content is used. (default 10)')
    parser.add_argument('--delimiter', type=str,
        help='Character used as separator between descriptions and values in the output. (default " ")')
    parser.add_argument('--max_mem_size', type=float,
        help='Maximal memory size in Gigabyte that can be used by this script. (default 4)')


    return parser.parse_args()

if __name__ == '__main__':
    main()
