#!/usr/bin/env python3
# ==============================================================
# Lars Gabriel
#
# Create hintsfile from the hints of a GeneMark-ETP+ run
# ==============================================================
import argparse
import os
import sys
import subprocess as sp
import csv

class GeneMarkDirNotFound(Exception):
    pass
class ScriptNotFound(Exception):
    pass

def main():
    args = parseCmd()

    # figure the path to the correct directory out
    etp_sub_dir = ''
    if os.path.exists(f'{args.etp_wdir}/proteins.fa'):
        etp_sub_dir = f'proteins.fa'
    else:
        for dir in [d for d in os.listdir(args.etp_wdir) if os.path.isdir(d)]:
            if set(['hc', 'model', 'nonhc', 'penalty']).\
                issubset(os.listdir(f'{args.etp_wdir}/{d}/')):
                etp_sub_dir = f'{d}'
                break
    if not etp_sub_dir:
        raise GeneMarkDirNotFound(
            f'ERROR: {args.etp_wdir} doesn\'t seem to be the working directory '\
            + f'of a GeneMarkETP+ run.')

    # getting nonHC protein hints
    if not os.path.exists(f'{args.genemark_scripts}/format_back.pl'):
        raise ScriptNotFound(f'The scripts {args.genemark_scripts}/formal_back.pl '\
            f'was not found in {args.genemark_scripts}.')
    cmd = f'{args.genemark_scripts}/format_back.pl {args.etp_wdir}/{etp_sub_dir}/'\
        + f'nonhc/prothint/prothint_augustus.gff '\
        + f'{args.etp_wdir}/{etp_sub_dir}/nonhc/nonhc.trace '\
        + f' >> {args.out}'
    sp.call(cmd, shell=True)

    # getting HC protein hints and RNA-Seq hints
    cmd = f'cat {args.etp_wdir}/rnaseq/hints/hintsfile_merged.gff '\
        + f'{args.etp_wdir}/rnaseq/hints/{etp_sub_dir}/prothint/prothint_augustus.gff '\
        + f'{args.etp_wdir}/rnaseq/hints/hintsfile_merged.gff '\
        + f'>> {args.out}'
    sp.call(cmd, shell=True)

def csv_write(tab, path, style='w+'):
    with open(path, style) as f:
        writer = csv.writer(f, delimiter='\t')
        for line in tab:
            writer.writerow(line)

def csv_read(file_path, de='\t'):
    result = []
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            lines = csv.reader(file, delimiter=de, quotechar = "|",\
                lineterminator = '\n')
            for line in lines:
                if line:
                    result.append(line)
    return result

def make_path_abs(path_to_file):
    abs_path = os.path.abspath(path_to_file)
    if not os.path.exists(abs_path):
        sys.stderr.write('Path does not exist: {}\n'.format(abs_path))
    return abs_path

def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description=\
        'Fetch extrinscic evidence hints for AUGUSTUS from a GeneMark-ETP+ run.')
    parser.add_argument('--genemark_scripts', type=str,
        help='Path to GeneMark-ETP+ scripts')
    parser.add_argument('--etp_wdir', type=str,
        help='Location of a GeneMark-ETP+ run')
    parser.add_argument('--out', type=str,
        help='Output hintsfile in AUGUSTUS format.')
    return parser.parse_args()

if __name__ == '__main__':
    main()
