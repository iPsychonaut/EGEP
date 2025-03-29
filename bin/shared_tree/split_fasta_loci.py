#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_fasta_loci.py

Updated on Wed Mar 19 2025

This script splits a concatenated FASTA file into individual FASTA files per sequence (locus),
based on the full sequence ID in the header. It is part of the Entheome Genome Extraction
Pipeline (EGEP) for preparing locus-specific files for gene tree inference in both BUSCO
and FASTA modes.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import os
from Bio import SeqIO

def split_fasta_loci(input_fasta, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    locus_dict = {}
    with open(input_fasta, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            # Extract locus from header (e.g., ITS from /mnt/d/...|ITS.AY281023.1.PSICU)
            parts = record.id.split("|")
            if len(parts) > 1:
                locus_id = parts[1].split('.')[0]  # e.g., ITS from ITS.AY281023.1.PSICU
            else:
                locus_id = "unknown"
            if locus_id not in locus_dict:
                locus_dict[locus_id] = []
            locus_dict[locus_id].append(record)
            print(f"DEBUG: Wrote sequence {record.id} to locus_files/{locus_id}.fasta", file=sys.stderr)
    
    for locus_id, records in locus_dict.items():
        output_file = os.path.join(output_dir, f"{locus_id}.fasta")
        with open(output_file, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: split_fasta_loci.py <input_fasta> <output_dir>", file=sys.stderr)
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    split_fasta_loci(input_fasta, output_dir)