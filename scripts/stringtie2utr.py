#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description="Updates BRAKER gene models with UTRs from a StringTie assembly.")
    
    # Mandatory input arguments
    parser.add_argument("-b", "--braker", required=True, help="File with BRAKER gene models in GTF format")
    parser.add_argument("-s", "--stringtie", required=True, help="File with StringTie transcript models in GFF format")

    args = parser.parse_args()

    # Now you can access the input file paths using args.braker and args.stringtie
    braker_file = args.braker
    stringtie_file = args.stringtie

    # Add your code here to process the input files and perform the desired actions

if __name__ == "__main__":
    main()
