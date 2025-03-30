#!/usr/bin/env python3
"""
extract_sequences_from_busco.py

Updated on Wed Mar 30 2025

Extracts sequences from Compleasm BUSCO output based on shared BUSCO IDs and concatenates them per assembly.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_sequences(csv_file, compleasm_db):
    fasta_output = f"{compleasm_db}_busco_shared_sequences.fasta"
    print(f"Starting: CSV={csv_file}, Output={fasta_output}, DB={compleasm_db}")

    # Read CSV
    with open(csv_file, "r") as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        print("ERROR: CSV is empty or has no data lines", file=sys.stderr)
        sys.exit(1)

    assemblies = []
    all_busco_ids = []
    for line in lines[1:]:
        parts = line.split(",")
        assembly = parts[0]
        busco_ids = parts[1:]
        assemblies.append(assembly)
        all_busco_ids.extend(busco_ids)
    shared_busco_ids = list(set(all_busco_ids))
    print(f"Found {len(assemblies)} assemblies, {len(shared_busco_ids)} unique BUSCO IDs")
    print(f"DEBUG: Sample CSV BUSCO IDs: {shared_busco_ids[:5]}", file=sys.stderr)

    records = []
    gff_found = False
    for assembly in assemblies:
        base_name = os.path.basename(assembly)
        organism_id = base_name.replace(".fasta", "")
        gff_path = os.path.join(
            os.path.dirname(assembly),
            f"{os.path.basename(assembly).replace('.fasta','')}_{compleasm_db}_odb10_compleasm",
            f"{compleasm_db}_odb10",
            "miniprot_output.gff"
        )
        print(f"DEBUG: Looking for GFF file at {gff_path}", file=sys.stderr)

        if not os.path.isfile(gff_path):
            print(f"ERROR: GFF not found at {gff_path} for {assembly}.", file=sys.stderr)
            continue
        gff_found = True

        with open(gff_path, "r") as f:
            lines = f.read().splitlines()
        print(f"DEBUG: GFF file {gff_path} has {len(lines)} lines", file=sys.stderr)
        if lines:
            print(f"DEBUG: First few lines of GFF: {lines[:5]}", file=sys.stderr)

        seq_dict = {}
        paf_ids = set()
        i = 0
        while i < len(lines):
            if lines[i].startswith("##PAF"):
                fields = lines[i].split("\t")
                if len(fields) < 2:
                    print(f"WARNING: Malformed ##PAF line: {lines[i]}", file=sys.stderr)
                    i += 1
                    continue
                full_paf_id = fields[1].strip()
                paf_id = full_paf_id.split("_")[0]  # e.g., "10022at5204_1036808_1:000f6b" -> "10022at5204"
                paf_ids.add(paf_id)
                i += 1
                if i < len(lines) and lines[i].startswith("##STA"):
                    prot_seq = lines[i][5:].strip()
                    # Try matching with and without 'at5204' suffix
                    base_id = paf_id.replace("at5204", "") if "at5204" in paf_id else paf_id
                    if paf_id in shared_busco_ids or base_id in shared_busco_ids:
                        matched_id = paf_id if paf_id in shared_busco_ids else base_id
                        seq_dict[matched_id] = prot_seq
                        print(f"DEBUG: Matched BUSCO ID {matched_id} in GFF for {assembly}", file=sys.stderr)
                    i += 1
                    while i < len(lines) and not lines[i].startswith("##"):
                        i += 1
                else:
                    i += 1
            else:
                i += 1

        print(f"DEBUG: Found {len(paf_ids)} unique core BUSCO IDs in GFF for {assembly}", file=sys.stderr)
        print(f"DEBUG: Sample GFF BUSCO IDs: {list(paf_ids)[:5]}", file=sys.stderr)
        print(f"DEBUG: Found {len(seq_dict)} matching BUSCOs in GFF for {assembly}", file=sys.stderr)
        if not seq_dict:
            print(f"WARNING: No matching BUSCOs in GFF for {assembly}.", file=sys.stderr)

        # Concatenate sequences
        concat_seq = "".join(seq_dict.get(busco, "") for busco in shared_busco_ids)
        if concat_seq:
            record = SeqRecord(Seq(concat_seq), id=f"{organism_id}-Shared_BUSCO_SEQ", description="")
            records.append(record)
            print(f"Added sequence for {organism_id} with {len(seq_dict)} BUSCOs")

    if not gff_found:
        print(f"ERROR: No GFF files found for any assembly in {csv_file}", file=sys.stderr)
        sys.exit(1)

    # Write output, even if empty, to allow pipeline to proceed
    SeqIO.write(records, fasta_output, "fasta")
    if not records:
        print(f"WARNING: No sequences extracted for {compleasm_db}, writing empty FASTA", file=sys.stderr)
    else:
        print(f"Success: Wrote {fasta_output} with {len(records)} records")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_sequences_from_busco.py <csv_file> <compleasm_db>", file=sys.stderr)
        sys.exit(1)
    csv_file = sys.argv[1]
    compleasm_db = sys.argv[2]
    extract_sequences(csv_file, compleasm_db)
