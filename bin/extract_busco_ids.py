#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_busco_ids.py

Updated on Thu Mar 13 2025

This script extracts BUSCO IDs from assembly tables and outputs a CSV with alignments
and their corresponding BUSCO IDs as lists. It also generates a list of all unique
BUSCO IDs and an initial heatmap DataFrame for the Entheome Genome Extraction Pipeline (EGEP).

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
import pandas as pd

def extract_busco_ids(assemblies_file, compleasm_db, tables_file, total_buscos, threshold):
    print(f"DEBUG: compleasm_db={compleasm_db}, total_buscos={total_buscos}, threshold={threshold}", file=sys.stderr)
    
    # Read assemblies and tables
    with open(assemblies_file) as f:
        assemblies = [line.strip() for line in f]
    with open(tables_file) as f:
        tables = [line.strip() for line in f]
    print(f"DEBUG: Found {len(assemblies)} assemblies, {len(tables)} tables", file=sys.stderr)

    # Keep track of original data
    original_assemblies = assemblies.copy()
    original_tables = tables.copy()
    
    while True:
        # Get BUSCO IDs for each assembly
        busco_sets = []
        assembly_buscos = {}
        for i, table in enumerate(tables):
            df = pd.read_csv(table, sep="\t")
            buscos = set(df[df["Status"].isin(["Complete", "Duplicated"])]["# Busco id"])
            busco_sets.append(buscos)
            assembly_buscos[assemblies[i]] = buscos
            print(f"DEBUG: Assembly {assemblies[i]} has {len(buscos)} BUSCOs", file=sys.stderr)

        # Find shared BUSCO IDs
        shared_buscos = set.intersection(*busco_sets) if busco_sets else set()
        shared_count = len(shared_buscos)
        total = int(total_buscos)
        shared_percent = shared_count / total if total > 0 else 0
        print(f"DEBUG: Shared={shared_count}, Total={total}, Percent={shared_percent*100:.1f}%", file=sys.stderr)
        
        # Check threshold
        if shared_percent >= float(threshold) or len(assemblies) <= 1:
            break
            
        # Remove assembly with smallest overlap
        min_overlap = float('inf')
        assembly_to_remove = None
        for asm, buscos in assembly_buscos.items():
            overlap = len(buscos & shared_buscos)
            if overlap < min_overlap:
                min_overlap = overlap
                assembly_to_remove = asm
                
        remove_index = assemblies.index(assembly_to_remove)
        assemblies.pop(remove_index)
        tables.pop(remove_index)
        print(f"Removed {assembly_to_remove} - shared BUSCOs dropped to {shared_count} ({shared_percent*100:.1f}%)", file=sys.stderr)

    # Write shared BUSCO IDs to CSV
    with open(f"{compleasm_db}_busco_shared.csv", "w") as f:
        f.write("assembly,shared_busco_ids\n")
        for asm in original_assemblies:
            f.write(f"{asm},{','.join(shared_buscos)}\n")

    # Calculate and write overlap string
    overlap_str = f"Shared BUSCO ID overlap for {compleasm_db}: {shared_count}/{total} ({shared_percent*100:.1f}%)"
    with open(f"{compleasm_db}_busco_overlap.txt", "w") as f:
        f.write(overlap_str)

    print(f"DEBUG: Wrote overlap: {overlap_str}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python extract_busco_ids.py <assemblies_file> <compleasm_db> <tables_file> <total_buscos> <threshold>")
        sys.exit(1)
    extract_busco_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
