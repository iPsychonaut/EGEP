#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_aa_sequences.py

Created on Mon Mar 24 2025

This script extracts amino acid sequences from an Exonerate GTF file and saves them
as a single FASTA file.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1


def is_valid_dna_sequence(s):
    """Check if a string is a valid DNA sequence.

    Parameters:
        s (str): String to check.

    Returns:
        bool: True if the string contains only valid DNA characters, False otherwise.
    """
    valid_chars = set("ATCGatcg{}.-")
    return all(char in valid_chars for char in s)


def extract_aa_sequences(gtf_file, output_fasta):
    """Extract amino acid sequences from a GTF file and save as a single FASTA file.

    Parses the GTF file to extract amino acid sequences from alignment sections,
    converts them to one-letter codes, and writes them to a single FASTA file.

    Parameters:
        gtf_file (str): Path to the input GTF file.
        output_fasta (str): Path to the output FASTA file.

    Returns:
        None: Writes the amino acid sequences to the specified FASTA file.

    Raises:
        FileNotFoundError: If the GTF file does not exist.
    """
    if not os.path.exists(gtf_file):
        print(f"Error: GTF file not found: {gtf_file}")
        sys.exit(1)

    sequences = []
    with open(gtf_file, 'r') as file:
        extraction_start = False
        base_name = ""
        temp_name_data = ""
        temp_seq_data = ""
        for line in file:
            if "C4 Alignment" in line:
                extraction_start = True
                base_name = ""
                temp_name_data = ""
                temp_seq_data = ""
            if extraction_start:
                if "START OF GFF DUMP" in line:
                    extraction_start = False
                    if temp_seq_data:
                        seq = Seq(seq1(temp_seq_data))  #Convert to one-letter code
                        seq_record = SeqRecord(seq, id=f"{base_name} {temp_name_data}", description="")
                        sequences.append(seq_record)
                else:
                    if "Query:" in line:
                        temp_query = line.split(": ")[-1].replace(" \n", "")
                        temp_name_data = temp_query
                    if "Target:" in line:
                        temp_target = line.split(": ")[-1].replace("\n", "").replace('.', '-')
                        base_name = os.path.basename(gtf_file).replace(".gtf", "")
                    if "Target range:" in line:
                        temp_range = line.split(": ")[-1]
                        if " [revcomp]" in temp_target:
                            temp_start = temp_range.split(" -> ")[-1]
                            temp_end = temp_range.split(" -> ")[0]
                        else:
                            temp_end = temp_range.split(" -> ")[-1]
                            temp_start = temp_range.split(" -> ")[0]
                        start = temp_start.replace('\n', '')
                        end = temp_end.replace('\n', '')
                        temp_name_data = f"{temp_name_data} {start}:{end}"
                    if " : " in line:
                        mRNA_code_raw = line.split(" : ")[1]
                        mRNA_code = ''.join([char for char in mRNA_code_raw if not char.isdigit()])
                        remove_list = [">", ".", "Target Intron", " ", "{", "}", "<-", "!", "|", ":", "-"]
                        for item in remove_list:
                            mRNA_code = mRNA_code.replace(item, "")
                        if not is_valid_dna_sequence(mRNA_code):
                            temp_seq_data += mRNA_code

    #Write all sequences to a single FASTA file
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"PASS:\tExtracted {len(sequences)} amino acid sequences to {output_fasta}")


if __name__ == "__main__":
    """Command-line entry point for extracting amino acid sequences.

    Expects two arguments: input GTF file and output FASTA file.
    """
    if len(sys.argv) != 3:
        print("Usage: python extract_aa_sequences.py <gtf_file> <output_fasta>")
        sys.exit(1)
    gtf_file = sys.argv[1]
    output_fasta = sys.argv[2]
    extract_aa_sequences(gtf_file, output_fasta)