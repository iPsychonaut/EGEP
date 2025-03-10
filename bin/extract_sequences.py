#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_sequences.py

Created on Sun Mar 9 18:44:50 2025

This script extracts protein sequences from GFF files for shared BUSCO IDs in the
Entheome Genome Extraction Pipeline (EGEP). It generates individual FASTA files for
each assembly and a single combined FASTA file for all shared BUSCO sequences.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def extract_sequences(assemblies_file, shared_busco_ids_file, compleasm_db, combined_fasta_output=None):
    """Extract protein sequences from GFF files and generate FASTA files.

    Reads assembly paths and shared BUSCO IDs, parses corresponding GFF files to
    extract protein sequences, writes them to species-specific FASTA files, and
    optionally combines them into a single FASTA file.

    Args:
        assemblies_file (str): Path to the file listing filtered assembly paths.
        shared_busco_ids_file (str): Path to the file listing shared BUSCO IDs.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        combined_fasta_output (str, optional): Path to save the combined FASTA file.
                                               If None, no combined file is generated.

    Returns:
        None: Outputs individual FASTA files named "<species_id>_shared_busco.fasta"
              and, if specified, a combined FASTA file.

    Raises:
        FileNotFoundError: If an assembly, GFF, or input file cannot be found.
    """
    try:
        with open(assemblies_file, "r") as f:
            assemblies = [line.strip() for line in f]
        with open(shared_busco_ids_file, "r") as f:
            shared_buscos = [line.strip() for line in f]
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {e.filename}")
        sys.exit(1)

    combined_records = []
    for asm in assemblies:
        species_id = "_".join(os.path.basename(asm).split("_")[:2])
        gff = asm.replace(
            ".fasta",
            f"_{compleasm_db}_odb10_compleasm/{compleasm_db}_odb10/miniprot_output.gff"
        )
        out_faa = f"{species_id}_shared_busco.fasta"

        try:
            with open(gff, "r") as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: GFF file not found for {asm}: {gff}")
            continue

        records = {}
        i = 0
        while i < len(lines):
            if lines[i].startswith("##PAF"):
                paf_id = lines[i].split("\t")[1].strip().split("_")[0]
                matched_busco = next((b for b in shared_buscos if b in paf_id), None)
                i += 1
                if i < len(lines) and lines[i].startswith("##STA"):
                    prot_seq = lines[i][5:].strip()
                    if matched_busco:
                        record = SeqRecord(
                            Seq(prot_seq),
                            id=paf_id,
                            description=f"Shared_{compleasm_db}_BUSCO:{species_id}|{matched_busco}"
                        )
                        records[matched_busco] = record
                i += 1
            else:
                i += 1

        ordered_records = [records[b] for b in shared_buscos if b in records]
        if ordered_records:
            SeqIO.write(ordered_records, out_faa, "fasta")
            combined_records.extend(ordered_records)
            print(f"PASS:\tSequences extracted to {out_faa}")
        else:
            print(f"Warning: No shared BUSCO sequences found for {asm}")

    # Generate combined FASTA file if output path is provided
    if combined_fasta_output and combined_records:
        SeqIO.write(combined_records, combined_fasta_output, "fasta")
        print(f"PASS:\tCombined FASTA file created: {combined_fasta_output}")


if __name__ == "__main__":
    """Command-line entry point for extracting sequences from GFF files.

    Expects three or four arguments: assemblies file, shared BUSCO IDs file,
    Compleasm database identifier, and optionally a combined FASTA output path.
    """
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(
            "Usage: python extract_sequences.py <assemblies_file> "
            "<shared_busco_ids_file> <compleasm_db> [<combined_fasta_output>]"
        )
        sys.exit(1)
    assemblies_file = sys.argv[1]
    shared_busco_ids_file = sys.argv[2]
    compleasm_db = sys.argv[3]
    combined_fasta_output = sys.argv[4] if len(sys.argv) == 5 else None
    extract_sequences(assemblies_file, shared_busco_ids_file, compleasm_db, combined_fasta_output)