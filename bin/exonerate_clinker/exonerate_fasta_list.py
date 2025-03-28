#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exonerate_fasta_list.py

Created on Wed Mar 26 14:22:33 2025

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import multiprocessing
import os
import subprocess
import sys
import pandas as pd

# Function to Exonerate a set of sequences from a cassette contig
def exonerate_fasta(input_fasta, sequence_fasta, target_percent, organism_kingdom, maxintron):
    """Run Exonerate to align query sequences against a target FASTA file and generate GTF output.

    Args:
        input_fasta (str): Path to the target FASTA file (genomic sequence).
        sequence_fasta (str): Path to the query FASTA file (nucleotide or protein sequences).
        target_percent (float): Percent identity threshold for Exonerate alignment.
        organism_kingdom (str): Kingdom of the organism ('Funga', 'Flora', 'Fauna', or other).
        maxintron (str or float): Maximum intron length for alignment; if 'nan', set based on kingdom.

    Returns:
        str or None: Path to the generated GTF file if successful, None if an error occurs.
    """
    print(f"Exonerate Input Fasta: {input_fasta}")
    # Generate output GTF filename from input FASTA
    if '.fasta' in input_fasta:
        exonerate_output = input_fasta.replace('.fasta', '.gtf')
    elif '.fna' in input_fasta:
        exonerate_output = input_fasta.replace('.fna', '.gtf')
    else:
        exonerate_output = input_fasta + '.gtf'  # Fallback if no extension match
    
    # Skip if output exists and is sufficiently large (>1KB)
    if os.path.exists(exonerate_output) and os.path.getsize(exonerate_output) > 1000:
        print(f"NOTE:\tSkipping Exonerate: {exonerate_output} already exists and is larger than 1KB")
        return exonerate_output
    
    print(f'Running Exonerate on {input_fasta} targeting {sequence_fasta}...')
    
    # Set maxintron based on kingdom if not provided
    if not pd.isna(maxintron):
        maxintron = maxintron
    elif organism_kingdom == "Funga":
        maxintron = 200
    elif organism_kingdom == "Flora":
        maxintron = 2000
    elif organism_kingdom == "Fauna":
        maxintron = 2000
    else:
        maxintron = 200        
        print(f"UNKNOWN KINGDOM (organism_kingdom) to set exonerate maxintron, using default = {maxintron}")
    
    # Determine model based on query sequence type and fix maxintron in command
    if "_nuc" in sequence_fasta:
        exonerate_cmd = (
            f'exonerate --model cdna2genome --percent {target_percent} '
            f'--showalignment yes --showtargetgff yes --showvulgar no '
            f'-E --query {sequence_fasta} --maxintron {maxintron} '
            f'--target {input_fasta} > {exonerate_output}'
        )
    else:
        exonerate_cmd = (
            f'exonerate --model protein2genome:bestfit --percent {target_percent} '
            f'--showalignment yes --showtargetgff yes --showvulgar no '
            f'-E --query {sequence_fasta} --maxintron {maxintron} '
            f'--target {input_fasta} > {exonerate_output}'
        )
    
    print(f"CMD:\t{exonerate_cmd}")
    
    # Execute Exonerate command
    process = subprocess.run(exonerate_cmd, shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"ERROR:\t{process.stderr}")
        return None
    else:
        print("PASS:\tSuccessfully Exonerated genes of interest and Generated GTF file")
        return exonerate_output

def exonerate_pipeline(input_fasta_list, sequence_fasta, exonerate_pid, organism_kingdom, maxintron, cpu_threads):
    """Execute Exonerate on a list of FASTA files in parallel.

    Args:
        input_fasta_list (str): Comma-separated string of FASTA file paths.
        sequence_fasta (str): Path to the query FASTA file.
        exonerate_pid (float): Percent identity threshold for Exonerate.
        organism_kingdom (str): Kingdom of the organism.
        maxintron (str or float): Maximum intron length; can be 'nan' to use kingdom-based default.
        cpu_threads (int): Number of CPU threads for parallel processing.
    """
    # Split the comma-separated string into a list
    fasta_files = input_fasta_list.split(',')
    
    # Run exonerate_fasta in parallel
    with multiprocessing.Pool(cpu_threads) as pool:
        results = [
            pool.apply_async(exonerate_fasta, args=(
                fasta_file, sequence_fasta, exonerate_pid, organism_kingdom, maxintron))
            for fasta_file in fasta_files
        ]
        for result in results:
            result.get()  # Wait for all tasks to complete and raise any exceptions

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python exonerate_fasta_list.py <input_fasta_list> <sequence_fasta> "
              "<exonerate_pid> <organism_kingdom> <maxintron> <cpu_threads>", file=sys.stderr)
        sys.exit(1)
    
    input_fasta_list = sys.argv[1]  # Comma-separated list of FASTA files
    sequence_fasta = sys.argv[2]
    exonerate_pid = float(sys.argv[3])
    organism_kingdom = sys.argv[4]
    maxintron = sys.argv[5] if sys.argv[5] != 'nan' else float('nan')
    cpu_threads = int(sys.argv[6])
    
    exonerate_pipeline(input_fasta_list, sequence_fasta, exonerate_pid, organism_kingdom, maxintron, cpu_threads)