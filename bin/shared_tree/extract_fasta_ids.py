#!/usr/bin/env python3
"""
extract_fasta_ids.py

Created on Fri Mar 21 2025

Extracts FASTA loci shared across ALL assemblies, selects the best representative sequence per locus per assembly,
and concatenates them into a single sequence per assembly for the Entheome Genome Extraction Pipeline (EGEP).

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import pandas as pd
from Bio import SeqIO

def extract_fasta_ids(assemblies_file, blast_hits, input_fasta, output_fasta, overlap_file, threshold):
    print(f"DEBUG: Starting with threshold={threshold}", file=sys.stderr)
    
    with open(assemblies_file) as f:
        assemblies = [line.strip() for line in f]
    print(f"DEBUG: Found {len(assemblies)} assemblies", file=sys.stderr)

    # Load BLAST hits into a dictionary by assembly
    blast_dfs = {}
    for hit_file in blast_hits.split():
        asm_name = hit_file.split('/')[-1].replace('_blast_hits.tsv', '')
        df = pd.read_csv(hit_file, sep="\t", header=None, names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
        df = df[df['pident'] > 70]  # Filter by percent identity
        blast_dfs[asm_name] = df
        print(f"DEBUG: Loaded {hit_file} with {len(df)} hits", file=sys.stderr)
    print(f"DEBUG: Loaded {len(blast_dfs)} BLAST hit files", file=sys.stderr)

    # Identify loci present in ALL assemblies
    locus_presence = {}
    for asm, df in blast_dfs.items():
        loci = set(df["sseqid"].str.split('.').str[0])  # e.g., ITS from ITS.AY281023.1.PSICU
        for locus in loci:
            locus_presence[locus] = locus_presence.get(locus, set()) | {asm}

    # Only keep loci present in every assembly
    all_assemblies = set(assemblies)
    shared_loci = [locus for locus, asm_set in locus_presence.items() if asm_set == all_assemblies]
    if not shared_loci:
        shared_loci = sorted(locus_presence.keys(), key=lambda x: len(locus_presence[x]), reverse=True)[:2]  # Fallback to top 2 most common
        print(f"DEBUG: No loci in all assemblies, using top {len(shared_loci)} loci: {shared_loci}", file=sys.stderr)
    print(f"DEBUG: Shared loci (in all {len(assemblies)} assemblies): {shared_loci}", file=sys.stderr)

    # Load assembly sequences
    asm_seqs = {}
    for asm in assemblies:
        asm_seqs[asm] = {rec.id: str(rec.seq) for rec in SeqIO.parse(asm, "fasta")}
    print(f"DEBUG: Loaded sequences from {len(asm_seqs)} assemblies", file=sys.stderr)

    # Select best representative sequence per shared locus per assembly and concatenate
    concatenated_seqs = {}
    for asm in assemblies:
        asm_key = asm.split('/')[-1].replace('.fasta', '')
        if asm_key in blast_dfs:
            locus_best_hits = {}
            missing_locus = False
            for locus in shared_loci:
                locus_hits = blast_dfs[asm_key][blast_dfs[asm_key]["sseqid"].str.startswith(locus)]
                if not locus_hits.empty:
                    # Select hit with highest bitscore
                    best_hit = locus_hits.loc[locus_hits["bitscore"].idxmax()]
                    locus_best_hits[locus] = best_hit
                    print(f"DEBUG: Best hit for {locus} in {asm_key}: {best_hit['sseqid']} (bitscore={best_hit['bitscore']})", file=sys.stderr)
                else:
                    print(f"DEBUG: No hit for {locus} in {asm_key}, skipping assembly", file=sys.stderr)
                    missing_locus = True
                    break

            if missing_locus:
                continue  # Skip this assembly if any shared locus is missing

            # Concatenate best hits into one sequence per assembly
            concat_seq = ""
            concat_desc = []
            for locus in shared_loci:  # Ensure order consistency
                hit = locus_best_hits[locus]
                qseqid = hit["qseqid"]
                start, end = int(hit["qstart"]) - 1, int(hit["qend"])
                if qseqid in asm_seqs[asm]:
                    seq = asm_seqs[asm][qseqid][start:end]
                    concat_seq += seq
                    concat_desc.append(f"{hit['sseqid']}:{start}-{end}")
                else:
                    print(f"DEBUG: Sequence {qseqid} not found in {asm}, skipping assembly", file=sys.stderr)
                    concat_seq = ""
                    break

            if concat_seq:
                concatenated_seqs[asm] = (concat_seq, ";".join(concat_desc))
                print(f"DEBUG: Concatenated sequence for {asm_key} with {len(concat_seq)} bp", file=sys.stderr)

    # Write concatenated sequences
    with open(output_fasta, "w") as f:
        for asm, (seq, desc) in concatenated_seqs.items():
            asm_key = asm.split('/')[-1].replace('.fasta', '')
            header = f"{asm_key}|concatenated"
            f.write(f">{header} {desc}\n{seq}\n")

    overlap_str = f"Shared FASTA loci used: {len(shared_loci)}"
    with open(overlap_file, "w") as f:
        f.write(overlap_str)
    print(f"DEBUG: Wrote overlap: {overlap_str}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python extract_fasta_ids.py <assemblies_file> <blast_hits> <input_fasta> <output_fasta> <overlap_file> <threshold>")
        sys.exit(1)
    extract_fasta_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])