#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_sequences(csv_file, compleasm_db):
    # Define output file directly (no subdirectory in script)
    fasta_output = f"{compleasm_db}_busco_shared_sequences.fasta"

    print(f"Starting: CSV={csv_file}, Output={fasta_output}, DB={compleasm_db}")

    # Read CSV
    with open(csv_file, "r") as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        print("Error: CSV is empty or has no data lines", file=sys.stderr)
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

    records = []
    for assembly in assemblies:
        base_name = os.path.basename(assembly)
        if "_" not in base_name:
            print(f"Skipping {base_name}: no underscore in filename")
            continue
        organism_id = base_name.replace(".fasta", "")
        gff_path = assembly.replace(".fasta", f"_{compleasm_db}_odb10_compleasm/{compleasm_db}_odb10/miniprot_output.gff")
        if not os.path.isfile(gff_path):
            print(f"Error: GFF not found at {gff_path} for {assembly}", file=sys.stderr)
            sys.exit(1)

        with open(gff_path, "r") as f:
            lines = f.read().splitlines()

        seq_dict = {}
        i = 0
        while i < len(lines):
            if lines[i].startswith("##PAF"):
                fields = lines[i].split("\t")
                paf_id = fields[1].split("_")[0]
                i += 1
                if i < len(lines) and lines[i].startswith("##STA"):
                    prot_seq = lines[i][5:]
                    if paf_id in shared_busco_ids:
                        seq_dict[paf_id] = prot_seq
                    i += 1
                    while i < len(lines) and not lines[i].startswith("##"):
                        i += 1
                else:
                    i += 1
            else:
                i += 1

        if not seq_dict:
            print(f"Error: No matching BUSCOs in GFF for {assembly}", file=sys.stderr)
            sys.exit(1)

        concat_seq = "".join(seq_dict.get(busco, "") for busco in shared_busco_ids)
        if concat_seq:
            record = SeqRecord(Seq(concat_seq), id=f"{organism_id}-Shared_BUSCO_SEQ", description="")
            records.append(record)
            print(f"Added sequence for {organism_id} with {len(seq_dict)} BUSCOs")

    if not records:
        print("Error: No sequences extracted from any GFF files", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(records, fasta_output, "fasta")
    if not os.path.isfile(fasta_output):
        print(f"Error: Failed to write {fasta_output}", file=sys.stderr)
        sys.exit(1)
    print(f"Success: Wrote {fasta_output} with {len(records)} records")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 extract_sequences.py <csv_file> <compleasm_db>", file=sys.stderr)
        sys.exit(1)
    extract_sequences(sys.argv[1], sys.argv[2])