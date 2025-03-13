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

import os
import sys
import pandas as pd


def extract_busco_ids(assemblies_file, tables_file, compleasm_db,
                      busco_dict_output, all_busco_ids_output, heatmap_df_output):
    """Extract BUSCO IDs from tables and output a CSV with alignments and BUSCO IDs.

    Args:
        assemblies_file (str): Path to the file listing assembly paths.
        tables_file (str): Path to the file listing BUSCO table paths.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        busco_dict_output (str): Path to save the BUSCO ID dictionary CSV.
        all_busco_ids_output (str): Path to save the list of all unique BUSCO IDs.
        heatmap_df_output (str): Path to save the heatmap DataFrame CSV.

    Returns:
        None: Outputs are written to the specified files.
    """
    try:
        with open(assemblies_file, "r") as f:
            assemblies = [line.strip() for line in f]
        with open(tables_file, "r") as f:
            tables = [line.strip() for line in f]
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {e.filename}", file=sys.stderr)
        sys.exit(1)

    if len(assemblies) != len(tables):
        print(f"Error: Mismatch between number of assemblies ({len(assemblies)}) and tables ({len(tables)})", file=sys.stderr)
        sys.exit(1)

    busco_id_dict = {}
    all_busco_ids = set()

    # Extract BUSCO IDs from tables using full assembly paths as keys
    for asm, tab in zip(assemblies, tables):
        print(f"Processing assembly: {asm}...", file=sys.stderr)
        try:
            df = pd.read_csv(tab, sep="\t")
            buscos = df[df["Status"].isin(["Complete", "Duplicated"])]["# Busco id"].tolist()
            busco_id_dict[asm] = buscos
            all_busco_ids.update(buscos)
        except FileNotFoundError:
            print(f"Error: Table file not found: {tab}", file=sys.stderr)
            continue
        except pd.errors.EmptyDataError:
            print(f"Error: BUSCO table is empty or malformed: {tab}", file=sys.stderr)
            continue

    if not busco_id_dict:
        print("Error: No valid BUSCO data extracted from tables", file=sys.stderr)
        sys.exit(1)

    # Save all BUSCO IDs
    with open(all_busco_ids_output, "w") as f:
        f.write("\n".join(sorted(all_busco_ids)))
    print("PASS:\tAll unique BUSCO IDs saved.", file=sys.stderr)

    # Build initial heatmap DataFrame (using assembly basenames for readability)
    heatmap_data = pd.DataFrame(0, index=sorted(all_busco_ids), 
                              columns=[os.path.basename(asm) for asm in busco_id_dict.keys()])
    for asm, buscos in busco_id_dict.items():
        col_name = os.path.basename(asm)
        for busco in buscos:
            heatmap_data.loc[busco, col_name] = 1
    heatmap_data.to_csv(heatmap_df_output)
    print("PASS:\tHeatmap DataFrame saved.", file=sys.stderr)

    # Save busco_id_dict as CSV with alignments and BUSCO lists
    dict_df = pd.DataFrame(
        [(asm, ",".join(buscos)) for asm, buscos in busco_id_dict.items()],
        columns=["alignment", "busco_ids"]
    )
    dict_df.to_csv(busco_dict_output, index=False)
    print("PASS:\tBUSCO ID dictionary saved as CSV with alignments.", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Usage: python extract_busco_ids.py <assemblies_file> <tables_file> "
            "<compleasm_db> <busco_dict_output> <all_busco_ids_output> "
            "<heatmap_df_output>",
            file=sys.stderr
        )
        sys.exit(1)
    extract_busco_ids(
        sys.argv[1], sys.argv[2], sys.argv[3],
        sys.argv[4], sys.argv[5], sys.argv[6]
    )
