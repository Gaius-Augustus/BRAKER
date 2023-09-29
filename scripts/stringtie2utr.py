#!/usr/bin/env python3

import argparse
import re


def read_gtf(gtf_file):
    """
    Reads a GTF file and extracts gene and non-gene features.

    Args:
    - gtf_file (str): Path to the GTF file.

    Returns:
    tuple: (non_gene_dict, gene_dict, tx_to_gene_dict, tx_dict) where:
    - non_gene_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.
    """
    non_gene_dict = {}
    gene_dict = {}
    tx_to_gene_dict = {}
    tx_dict = {}
    transcript_id_pattern = re.compile(r'transcript_id "([^"]+)"')
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')

    with open(gtf_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                feature_type = fields[2]
                last_field = fields[-1]

                if feature_type == "gene":
                    gene_dict[last_field] = line.strip()  # Use the entire last field as gene_id
                elif feature_type == "transcript":
                    tx_dict[last_field] = line.strip()  # Use the entire last field as transcript_id
                else:
                    transcript_id_match = transcript_id_pattern.search(last_field)
                    gene_id_match = gene_id_pattern.search(last_field)
                    gene_id = None
                    transcript_id = None
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)
                    if gene_id_match:
                        gene_id = gene_id_match.group(1)
                    if transcript_id:  # Ensure we have a transcript ID
                        non_gene_dict.setdefault(transcript_id, []).append(line.strip())
                        if transcript_id not in tx_to_gene_dict:
                            tx_to_gene_dict[transcript_id] = gene_id

    return non_gene_dict, gene_dict, tx_to_gene_dict, tx_dict




def add_intron_features(gtf_dict):
    """
    Adds intron feature lines based on the exon feature lines.
    
    Args:
    - gtf_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Updated GTF dictionary with intron feature lines added.
    """
    
    for transcript_id, entries in gtf_dict.items():
        sorted_exons = sorted([entry for entry in entries if entry.split('\t')[2] == 'exon'], key=lambda x: int(x.split('\t')[3]))
        
        new_entries = []
        for i in range(len(sorted_exons) - 1):
            curr_exon = sorted_exons[i]
            next_exon = sorted_exons[i + 1]
            
            curr_end = int(curr_exon.split('\t')[4])
            next_start = int(next_exon.split('\t')[3])
            
            if next_start - curr_end > 1:
                # Construct intron feature line based on exon line, but replace the feature name with "intron"
                intron_line = curr_exon.split('\t')
                intron_line[2] = 'intron'
                intron_line[3] = str(curr_end + 1)
                intron_line[4] = str(next_start - 1)
                
                # Remove exon_number attribute
                intron_line[8] = re.sub(r'exon_number "\d+";', '', intron_line[8]).strip()

                # Add the intron line to new entries
                new_entries.append('\t'.join(intron_line))

        # Extend the original entries with new intron feature lines
        gtf_dict[transcript_id].extend(new_entries)

    return gtf_dict



def create_introns_hash(non_gene_dict):
    """
    Creates a dictionary with intron strings as keys and transcript IDs as values.
    
    Args:
    - non_gene_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Dictionary with intron strings as keys and transcript IDs as values.
    """
    
    introns_hash = {}
    
    for transcript_id, entries in non_gene_dict.items():
        for entry in entries:
            # Check if the feature is an intron
            if entry.split('\t')[2] == 'intron':
                seqname = entry.split('\t')[0]
                start = entry.split('\t')[3]
                end = entry.split('\t')[4]
                strand = entry.split('\t')[6]
                intron_key = f"{seqname}_{start}_{end}_{strand}"

                # Store the transcript ID associated with the intron in the hash
                introns_hash[intron_key] = transcript_id

    return introns_hash


def find_matching_transcripts(intron_hash1, intron_hash2):
    """
    Find matching transcript IDs based on intron patterns.
    
    Args:
    - intron_hash1 (dict): Dictionary with intron strings as keys and transcript IDs from the first dataset as values.
    - intron_hash2 (dict): Dictionary with intron strings as keys and transcript IDs from the second dataset as values.

    Returns:
    dict: Dictionary with transcript IDs from intron_hash1 as keys and lists of matching transcript IDs from intron_hash2 as values.
    """
    
    # Reverse the hashes for easy lookup of intron patterns for each transcript
    reverse_hash1 = {}
    for intron, transcript in intron_hash1.items():
        if transcript not in reverse_hash1:
            reverse_hash1[transcript] = []
        reverse_hash1[transcript].append(intron)

    reverse_hash2 = {}
    for intron, transcript in intron_hash2.items():
        if transcript not in reverse_hash2:
            reverse_hash2[transcript] = []
        reverse_hash2[transcript].append(intron)

    matching_transcripts = {}

    for transcript1, introns1 in reverse_hash1.items():
        matches = set()
        
        for transcript2, introns2 in reverse_hash2.items():
            if all(intron in introns2 for intron in introns1):
                matches.add(transcript2)

        if matches:
            matching_transcripts[transcript1] = list(matches)

    return matching_transcripts


def tx_len(non_gene_dict):
    """
    Calculate the length of each transcript.

    Args:
    - non_gene_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Dictionary with transcript IDs as keys and their lengths as values.
    """
    
    transcript_lengths = {}

    for transcript_id, entries in non_gene_dict.items():
        total_length = 0
        
        for entry in entries:
            fields = entry.split('\t')
            start = int(fields[3])
            end = int(fields[4])
            
            total_length += (end - start + 1)  # +1 because both start and end are inclusive

        transcript_lengths[transcript_id] = total_length

    return transcript_lengths


def select_longest_tx(matched_transcripts, tx_lens_stringtie):
    """
    For each key in matched_transcripts, select the longest transcript based on lengths provided in tx_lens_stringtie.

    Args:
    - matched_transcripts (dict): Dictionary with keys as transcript IDs and values as lists of matching transcript IDs.
    - tx_lens_stringtie (dict): Dictionary with transcript IDs as keys and their lengths as values.

    Returns:
    dict: Dictionary with keys from matched_transcripts and values as the IDs of the longest transcripts.
    """
    
    selected_transcripts = {}

    for original_tx, matched_tx_list in matched_transcripts.items():
        # Find the longest transcript in the matched_tx_list
        longest_tx = max(matched_tx_list, key=lambda tx: tx_lens_stringtie.get(tx, 0))

        selected_transcripts[original_tx] = longest_tx

    return selected_transcripts


def overlap(exon_start, exon_end, cds_start, cds_end):
    """Check if exon overlaps with CDS."""
    return exon_start <= cds_end and exon_end >= cds_start


def merge_features(braker_gtf, stringtie_gtf, selected_transcripts):
    for braker_tx, stringtie_tx in selected_transcripts.items():
        # Retrieve features for current transcripts
        braker_features = braker_gtf[braker_tx]
        stringtie_exons = [f for f in stringtie_gtf[stringtie_tx] if f.split('\t')[2] == 'exon']
        braker_cds_features = [f for f in braker_features if f.split('\t')[2] == 'CDS']
        
        for exon in stringtie_exons:
            exon_parts = exon.split('\t')
            exon_start, exon_end = int(exon_parts[3]), int(exon_parts[4])
            
            # Flag to determine if current exon should be merged
            merge_exon = True
            
            for cds in braker_cds_features:
                cds_parts = cds.split('\t')
                cds_start, cds_end = int(cds_parts[3]), int(cds_parts[4])
                
                if overlap(exon_start, exon_end, cds_start, cds_end):
                    # Exon overlaps with CDS. Check the length condition.
                    if exon_end - exon_start + 1 < cds_end - cds_start + 1:
                        merge_exon = False
                        break
            
            # Add exon to braker_features if it passed the checks
            if merge_exon:
                braker_features.append(exon)
        
        # Update the braker_gtf dictionary with new features
        braker_gtf[braker_tx] = braker_features

    return braker_gtf


def compute_utr_features(braker_gtf):
    """
    Compute the UTR features for each transcript in braker_gtf based on strand information.

    Args:
    - braker_gtf (dict): Dictionary with transcript IDs as keys and lists of GTF feature lines as values.

    Returns:
    dict: Updated braker_gtf dictionary with UTR features added.
    """
    
    for transcript_id, features in braker_gtf.items():
        # Sort features by start position
        features.sort(key=lambda x: int(x.split('\t')[3]))
        
        utr_features = []
        for feature in features:
            fields = feature.split('\t')
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # If feature is an exon, check if there are overlapping CDS features
            if feature_type == "exon":
                overlapping_cds = [f for f in features if f.split('\t')[2] == "CDS" and int(f.split('\t')[3]) <= end and int(f.split('\t')[4]) >= start]
                
                if overlapping_cds:
                    cds_start = int(overlapping_cds[0].split('\t')[3])
                    cds_end = int(overlapping_cds[0].split('\t')[4])

                    # Check for UTR based on strand
                    if strand == "+":
                        if start < cds_start:
                            utr5 = "\t".join(fields[:2] + ["five_prime_UTR"] + fields[3:])
                            utr5 = utr5.replace(str(end), str(cds_start - 1))
                            utr_features.append(utr5)

                        if end > cds_end:
                            utr3 = "\t".join(fields[:2] + ["three_prime_UTR"] + fields[3:])
                            utr3 = utr3.replace(str(start), str(cds_end + 1))
                            utr_features.append(utr3)
                    
                    elif strand == "-":
                        if end > cds_end:
                            utr5 = "\t".join(fields[:2] + ["five_prime_UTR"] + fields[3:])
                            utr5 = utr5.replace(str(start), str(cds_end + 1))
                            utr_features.append(utr5)

                        if start < cds_start:
                            utr3 = "\t".join(fields[:2] + ["three_prime_UTR"] + fields[3:])
                            utr3 = utr3.replace(str(end), str(cds_start - 1))
                            utr_features.append(utr3)

        # Add the computed UTR features to the list of features for this transcript
        features.extend(utr_features)

    return braker_gtf


def fix_feature_coordinates(gtf_dict, gene_dict, tx_to_gene_dict, tx_dict):
    """
    The list of features in gtf_dict has been expanded compared to the original version. We need to identify the left most and right most coordinate
     for each transcript. Then, update the tx_dict lines with these new coordinates. In the end, loop over all tx_dict entries and identify the left most and right
      most coordinate per gene, then update the gene lines in gene_dict.
    
    Args:
    - gtf_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.

    Returns:
    gene_dict, tx_dict with updated coordinates
    """
    # 1. Update transcript coordinates based on feature coordinates in gtf_dict
    for tx_id, features in gtf_dict.items():
        starts = []
        ends = []
        for feature in features:
            fields = feature.split('\t')
            starts.append(int(fields[3]))
            ends.append(int(fields[4]))
        
        # Find the min and max coordinates for the transcript
        min_coord = min(starts)
        max_coord = max(ends)
        
        tx_fields = tx_dict[tx_id].split('\t')
        tx_fields[3] = str(min_coord)
        tx_fields[4] = str(max_coord)
        tx_dict[tx_id] = '\t'.join(tx_fields)
    
    # 2. Update gene coordinates based on updated transcript coordinates in tx_dict
    for gene_id in gene_dict:
        associated_transcripts = [tx for tx, gid in tx_to_gene_dict.items() if gid == gene_id]
        
        starts = []
        ends = []
        for tx_id in associated_transcripts:
            tx_fields = tx_dict[tx_id].split('\t')
            starts.append(int(tx_fields[3]))
            ends.append(int(tx_fields[4]))
        
        # Find the min and max coordinates for the gene
        min_coord = min(starts)
        max_coord = max(ends)
        
        gene_fields = gene_dict[gene_id].split('\t')
        gene_fields[3] = str(min_coord)
        gene_fields[4] = str(max_coord)
        gene_dict[gene_id] = '\t'.join(gene_fields)
    
    return gene_dict, tx_dict


def print_gtf(gtf_dict, gene_dict, tx_to_gene_dict, tx_dict):
    """
    Print GTF lines based on gene_dict and gtf_dict.
    
    Args:
    - gtf_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.

    Returns:
    None: Prints the GTF lines to stdout.
    """
    printed_gene = {}
    # Iterate over gene_dict entries
    for tx_id, tx_line in tx_dict.items():
        if not tx_id in printed_gene:
            print(gene_dict[tx_to_gene_dict[tx_id]])
            gene_printed = True
        print(tx_line)
        sorted_features = sorted(gtf_dict.get(tx_id, []), key=lambda x: int(x.split('\t')[3]))
        for feature in sorted_features:
            # If the feature is a UTR line, remove the exon_number
            if "UTR" in feature:
                # split line into fields
                fields = feature.split("\t")
                # build new gtf line
                print(fields[0] + "\tstringtie2utr\t", "\t".join(fields[2:8]), "\ttranscript_id \"" + tx_id + "\"; gene_id \"" + tx_to_gene_dict[tx_id] + "\";")
            elif "StringTie" not in feature:
                print(feature)


    
def main():
    parser = argparse.ArgumentParser(description="Updates BRAKER gene models with UTRs from a StringTie assembly.")
    
    # Mandatory input arguments
    parser.add_argument("-b", "--braker", required=True, help="File with BRAKER gene models in GTF format")
    parser.add_argument("-s", "--stringtie", required=True, help="File with StringTie transcript models in GFF format")

    args = parser.parse_args()

    # Now you can access the input file paths using args.braker and args.stringtie
    braker_file = args.braker
    stringtie_file = args.stringtie

    # read the braker_file
    braker_non_gene_dict, braker_gene_line_dict, braker_tx_to_gene_dict, braker_tx_dict = read_gtf(braker_file)

    # read the stringtie_file
    stringtie_non_gene_dict, stringtie_gene_line_dict, stringie_tx_to_gene_dict, stringtie_tx_dict = read_gtf(stringtie_file)
    # add intron features to the stringtie_non_gene_dict
    stringtie_non_gene_dict = add_intron_features(stringtie_non_gene_dict)

    # create introns hash for braker and stringtie_non_gene_dict)
    braker_introns_hash = create_introns_hash(braker_non_gene_dict)
    stringtie_introns_hash = create_introns_hash(stringtie_non_gene_dict)
    # find matching transcripts
    matched_transcripts = find_matching_transcripts(braker_introns_hash, stringtie_introns_hash)

    # calculate the length of each transcript
    tx_lens_stringtie = tx_len(stringtie_non_gene_dict)

    # select the longest transcript for each key in matched_transcripts
    final_matching_tx = select_longest_tx(matched_transcripts, tx_lens_stringtie)

    # merge features from stringtie_gtf into braker_gtf based on the selected transcripts
    braker_gtf = merge_features(braker_non_gene_dict, stringtie_non_gene_dict, final_matching_tx)
    # compute UTR features
    braker_gtf = compute_utr_features(braker_gtf)
    # fix gene and transcript coordinates
    braker_gene_line_dict, braker_tx_dict = fix_feature_coordinates(braker_gtf, braker_gene_line_dict, braker_tx_to_gene_dict, braker_tx_dict)
    # print the updated braker_gtf
    print_gtf(braker_gtf, braker_gene_line_dict, braker_tx_to_gene_dict, braker_tx_dict)

if __name__ == "__main__":
    main()
