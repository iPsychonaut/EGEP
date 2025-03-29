#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_sequences_from_fasta.py

Updated on Wed Mar 19 2025

Extracts and concatenates sequences from an assembly FASTA file based on BLAST hits to an input FASTA,
producing a single concatenated sequence per assembly.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq  # Import Seq from Bio.Seq
import os

def extract_and_concatenate_sequences(assembly, input_fasta, output_fasta):
    """
    Extract and concatenate sequences from the assembly based on BLAST hits to the input FASTA.

    Args:
        assembly (str): Path to the assembly FASTA file.
        input_fasta (str): Path to the input FASTA file with query sequences.
        output_fasta (str): Path to the output FASTA file for the concatenated sequence.
    """
    print(f"DEBUG: Starting with assembly {assembly}", file=sys.stderr)
    
    # Get assembly basename for unique ID
    assembly_name = os.path.basename(assembly).replace('.fasta', '')
    
    # Create BLAST database from the assembly
    try:
        subprocess.run(["makeblastdb", "-in", assembly, "-dbtype", "nucl"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("DEBUG: BLAST database created", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to create BLAST database: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Parse input sequences
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print(f"ERROR: No sequences in {input_fasta}", file=sys.stderr)
        sys.exit(1)
    print(f"DEBUG: Found {len(records)} sequences in input FASTA", file=sys.stderr)
    
    # Parse assembly sequences into a dictionary for quick lookup
    assembly_seqs = SeqIO.to_dict(SeqIO.parse(assembly, "fasta"))
    print(f"DEBUG: Loaded {len(assembly_seqs)} contigs from assembly", file=sys.stderr)
    
    # Temporary BLAST output file
    temp_blast = "temp_blast.tsv"
    
    # Run BLAST to find hits for all queries at once
    try:
        subprocess.run(
            ["blastn", "-query", input_fasta, "-db", assembly, "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-max_target_seqs", "1", "-out", temp_blast],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print("DEBUG: BLAST completed", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BLAST failed: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process BLAST results to find best hits
    best_hits = {}
    with open(temp_blast, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = fields
            bitscore = float(bitscore)
            if qseqid not in best_hits or bitscore > best_hits[qseqid][1]:
                best_hits[qseqid] = (sseqid, bitscore, int(sstart), int(send))
    print(f"DEBUG: Found {len(best_hits)} BLAST hits", file=sys.stderr)
    
    # Concatenate sequences
    concat_seq = ""
    hit_ids = []
    for rec in records:
        if rec.id in best_hits:
            sseqid, bitscore, sstart, send = best_hits[rec.id]
            if sseqid in assembly_seqs:
                seq = assembly_seqs[sseqid].seq
                start = min(sstart, send) - 1  # Convert to 0-based
                end = max(sstart, send)        # End is inclusive in BLAST
                hit_seq = seq[start:end]
                concat_seq += str(hit_seq)
                hit_ids.append(f"{rec.id}_{sseqid}")
            else:
                print(f"WARNING: Contig {sseqid} not in assembly", file=sys.stderr)
        else:
            print(f"WARNING: No hit for {rec.id}", file=sys.stderr)
    
    # Write the concatenated sequence
    if concat_seq:
        output_record = SeqIO.SeqRecord(
            Seq(concat_seq),  # Use Bio.Seq.Seq
            id=f"{assembly_name}-Concatenated_Genes",
            description=f"Concatenated hits: {','.join(hit_ids)}"
        )
        SeqIO.write([output_record], output_fasta, "fasta")
        if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
            print(f"DEBUG: Wrote concatenated sequence to {output_fasta}", file=sys.stderr)
        else:
            print(f"ERROR: Failed to write {output_fasta}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"ERROR: No sequences extracted for concatenation", file=sys.stderr)
        sys.exit(1)
    
    # Clean up temporary BLAST file
    os.remove(temp_blast)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_sequences_from_fasta.py <assembly> <input_fasta> <output_fasta>", file=sys.stderr)
        sys.exit(1)
    
    assembly = sys.argv[1]
    input_fasta = sys.argv[2]
    output_fasta = sys.argv[3]
    extract_and_concatenate_sequences(assembly, input_fasta, output_fasta)