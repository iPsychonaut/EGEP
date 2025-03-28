#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_cluster.py

Created on Mon Mar 24 2025

This script extracts a cluster subsequence from a FASTA file based on GFF3 annotations,
applying a 5000 bp buffer around the min/max positions.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import os
import pandas as pd
from Bio import SeqIO


def extract_cluster(gff3_file, fasta_file, output_fasta):
    """Extract a cluster subsequence from a FASTA file using GFF3 annotations.

    Determines the cluster range from GFF3 (min start to max end ± 5000 bp buffer)
    and extracts the corresponding subsequence from the FASTA file.

    Parameters:
        gff3_file (str): Path to the input GFF3 file.
        fasta_file (str): Path to the input FASTA file.
        output_fasta (str): Path to the output cluster FASTA file.

    Returns:
        None: Writes the cluster subsequence to the specified FASTA file.

    Raises:
        FileNotFoundError: If the GFF3 or FASTA file does not exist.
    """
    if not os.path.exists(gff3_file):
        print(f"Error: GFF3 file not found: {gff3_file}")
        sys.exit(1)
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        sys.exit(1)

    #Read GFF3 into DataFrame
    build_df = pd.read_csv(gff3_file, sep='\t', header=None, comment='#',
                           names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

    #Generate cluster range with 5000 bp buffer
    buffer_bp = 5000
    cluster_start = build_df['start'].min() - buffer_bp
    cluster_end = build_df['end'].max() + buffer_bp

    if cluster_start < 1:
        cluster_start = 1  #Ensure start is not negative

    #Extract subsequence from FASTA
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        if cluster_end > len(record):
            cluster_end = len(record)  #Adjust end if beyond sequence length
        record.seq = record.seq[cluster_start - 1:cluster_end]  #0-based indexing
        record.id = f"{record.id}_{cluster_start}->{cluster_end}"
        SeqIO.write(records, output_fasta, "fasta")
        print(f"PASS:\tExtracted cluster {cluster_start}-{cluster_end} to {output_fasta}")
        return

    print(f"Error: No sequences found in {fasta_file}")
    sys.exit(1)


if __name__ == "__main__":
    """Command-line entry point for extracting a cluster subsequence.

    Expects three arguments: GFF3 file, FASTA file, and output FASTA file.
    """
    if len(sys.argv) != 4:
        print("Usage: python extract_cluster.py <gff3_file> <fasta_file> <output_fasta>")
        sys.exit(1)
    gff3_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_fasta = sys.argv[3]
    extract_cluster(gff3_file, fasta_file, output_fasta)