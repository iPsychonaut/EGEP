#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_and_extract_best_contig.py

Created on Wed Mar 26 2025

This script analyzes BLAST output to find the best contig containing primary and secondary genes,
calculates a trimmed region with buffer (enforcing a max length), and extracts the sequence to a FASTA file.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import pandas as pd
import sys
import json
import argparse
from Bio import SeqIO
import numpy as np
import os  # Needed for working with paths

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze BLAST output and extract trimmed contig region")
    parser.add_argument("assembly_fasta", help="Assembly FASTA file containing contigs")
    parser.add_argument("blast_output", help="BLAST output file in tabular format")
    parser.add_argument("primary_genes", help="Space-separated primary gene IDs")
    parser.add_argument("secondary_genes", help="Space-separated secondary gene IDs")
    parser.add_argument("output_json", help="Output JSON file with analysis results")
    parser.add_argument("output_fasta", help="Output FASTA file with trimmed sequence")
    parser.add_argument("--buffer", type=int, default=5000, help="Buffer size in bp (default: 5000)")
    parser.add_argument("--max_length", type=int, default=5000, help="Max trimmed length in bp (default: 5000)")
    return parser.parse_args()

def convert_to_serializable(obj):
    """Convert numpy types to JSON-serializable Python types."""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    return obj

args = parse_args()
assembly_fasta = args.assembly_fasta
blast_output = args.blast_output
primary_genes = args.primary_genes.split()
secondary_genes = args.secondary_genes.split()
output_json = args.output_json
output_fasta = args.output_fasta
buffer_size = args.buffer
max_length = args.max_length

# Log initial inputs
print(f"Processing assembly: {assembly_fasta}", file=sys.stderr)
print(f"Primary genes: {primary_genes}", file=sys.stderr)
print(f"Secondary genes: {secondary_genes}", file=sys.stderr)
print(f"Buffer size: {buffer_size}bp, Max length: {max_length}bp", file=sys.stderr)

# Get contig lengths
contig_lengths = {record.id: len(record.seq) for record in SeqIO.parse(assembly_fasta, "fasta")}
print(f"Loaded {len(contig_lengths)} contigs", file=sys.stderr)

# Load BLAST output
df = pd.read_csv(blast_output, sep='\t', header=None,
                 names=['Query Seq-ID', 'Subject Seq-ID', '% identical matches', 
                        'Alignm len', '# mismatches', '# Gap openings', 'SoA-Q', 
                        'EoA-Q', 'SoA-Sub', 'EoA-Sub', 'Expect value', 'Bit score'])
print(f"Loaded {len(df)} BLAST hits", file=sys.stderr)

# Filter for all genes
all_genes = primary_genes + secondary_genes
filtered_df = df[df['Query Seq-ID'].str.contains('|'.join(all_genes))]

# Get best hit per gene per contig
best_hits = {}
contig_gene_counts = {}
for gene in all_genes:
    gene_hits = filtered_df[filtered_df['Query Seq-ID'].str.contains(gene)]
    if len(gene_hits) > 0:
        for contig, group in gene_hits.groupby('Subject Seq-ID'):
            best_hit = group.loc[group['Bit score'].idxmax()]
            start = min(best_hit['SoA-Sub'], best_hit['EoA-Sub'])
            end = max(best_hit['SoA-Sub'], best_hit['EoA-Sub'])
            if contig not in best_hits:
                best_hits[contig] = {}
            best_hits[contig][gene] = {"start": start, "end": end, "bit_score": best_hit['Bit score']}
            contig_gene_counts[contig] = contig_gene_counts.get(contig, set()).union({gene})

# Find best contig (must contain ALL primary genes)
best_contig = None
all_primary_set = set(primary_genes)  # Set of all primary genes to check against
for contig, genes in contig_gene_counts.items():
    primary_hits = {g for g in genes if g in primary_genes}
    if primary_hits == all_primary_set and contig in contig_lengths:  # Must have ALL primary genes
        best_contig = contig
        break  # Take the first valid contig; could refine further if needed

# Prepare output and extract sequence
output_dict = {"best_contig": best_contig if best_contig else "None"}
if best_contig:
    contig_hits = best_hits[best_contig]
    contig_length = contig_lengths[best_contig]
    # Use only primary gene hits for region calculation
    primary_hits = {gene: hit for gene, hit in contig_hits.items() if gene in primary_genes}
    starts = [hit["start"] for hit in primary_hits.values()]
    ends = [hit["end"] for hit in primary_hits.values()]
    
    raw_start = min(starts)
    raw_end = max(ends)
    raw_length = raw_end - raw_start + 1
    
    # Calculate buffered region
    buffered_start = max(1, raw_start - buffer_size)
    buffered_end = min(contig_length, raw_end + buffer_size)
    region_length = buffered_end - buffered_start + 1
        
    output_dict[best_contig] = {
        "genes": contig_hits,  # Still include all hits for reference
        "contig_length": contig_length,
        "region": {
            "raw_start": raw_start,
            "raw_end": raw_end,
            "buffered_start": buffered_start,
            "buffered_end": buffered_end,
            "buffer_size": buffer_size,
            "final_length": buffered_end - buffered_start + 1
        }
    }
    
    # Extract sequence using the updated logic
    for record in SeqIO.parse(assembly_fasta, "fasta"):
        if record.id == best_contig:
            trimmed_seq = record.seq[buffered_start - 1:buffered_end]
            # Extract the basename of the assembly FASTA file and remove the extension.
            fasta_basename = os.path.basename(assembly_fasta)
            fasta_basename_no_ext = fasta_basename.replace(".fasta", "")
            trimmed_record_id = f"{fasta_basename_no_ext}_trimmed"
            # Write the trimmed sequence to the output FASTA file.
            with open(output_fasta, "w") as f:
                f.write(f">{trimmed_record_id}_{buffered_start}-{buffered_end}\n{trimmed_seq}\n")
            print(f"Extracted {best_contig}: {buffered_start}-{buffered_end} ({buffered_end - buffered_start + 1}bp) with all primary genes", file=sys.stderr)
            break
    
    print(f"Best contig {best_contig} contains all {len(primary_genes)} primary genes", file=sys.stderr)
else:
    print(f"No contig found with all {len(primary_genes)} primary genes", file=sys.stderr)

# Write JSON output, ensuring serializability
with open(output_json, "w") as f:
    json.dump(convert_to_serializable(output_dict), f, indent=2)
print(f"Output written to {output_json} and {output_fasta}", file=sys.stderr)
