#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_busco_loci.py

Updated on Mon Mar 17 2025

This script splits a concatenated FASTA file containing BUSCO sequences into individual
FASTA files per locus, based on the BUSCO ID in the sequence header. It is part of the
Entheome Genome Extraction Pipeline (EGEP) for preparing locus-specific files for gene
tree inference.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
from Bio import SeqIO
import os

def split_busco_loci(input_fasta, output_dir="locus_files"):
    """
    Split a concatenated BUSCO FASTA file into separate files per locus.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory to store individual locus FASTA files (default: 'locus_files').
    """
    # Read all sequences from the input FASTA
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print(f"ERROR: No sequences found in {input_fasta}", file=sys.stderr)
        sys.exit(1)
    
    # Group sequences by locus (assuming BUSCO ID is the first part of the header)
    locus_dict = {}
    for rec in records:
        # Extract locus ID (e.g., 'BUSCO123' from 'BUSCO123_Psilocybe_...')
        locus = rec.id.split('_')[0]
        if locus not in locus_dict:
            locus_dict[locus] = []
        locus_dict[locus].append(rec)
    
    # Create output directory if it doesn’t exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write each locus to a separate FASTA file
    for locus, seqs in locus_dict.items():
        output_file = os.path.join(output_dir, f"{locus}.fasta")
        SeqIO.write(seqs, output_file, "fasta")
        print(f"DEBUG: Wrote {len(seqs)} sequences to {output_file}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python split_busco_loci.py <input_fasta> [output_dir]", file=sys.stderr)
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) == 3 else "locus_files"
    split_busco_loci(input_fasta, output_dir)