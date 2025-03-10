#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_busco_ids.py

Created on Sun Mar 9 18:44:50 2025

This script extracts BUSCO IDs from assembly tables, groups them by species, and
generates an initial heatmap for the Entheome Genome Extraction Pipeline (EGEP).
It outputs a BUSCO ID dictionary, all unique BUSCO IDs, and a heatmap DataFrame.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import sys
import pandas as pd


def extract_busco_ids(assemblies_file, tables_file, compleasm_db, in_genus,
                      busco_dict_output, all_busco_ids_output, heatmap_df_output,
                      heatmap_png_output):
    """Extract BUSCO IDs from tables and generate an initial heatmap.

    Reads assembly and table files, extracts complete and duplicated BUSCO IDs,
    groups them by species with an "in" or "out" label based on genus, and creates
    a presence/absence DataFrame for heatmap generation.

    Args:
        assemblies_file (str): Path to the file listing assembly paths.
        tables_file (str): Path to the file listing BUSCO table paths.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        in_genus (str): Genus prefix for grouping assemblies (e.g., "Ps").
        busco_dict_output (str): Path to save the BUSCO ID dictionary CSV.
        all_busco_ids_output (str): Path to save the list of all unique BUSCO IDs.
        heatmap_df_output (str): Path to save the heatmap DataFrame CSV.
        heatmap_png_output (str): Path to save the initial heatmap PNG.

    Returns:
        None: Outputs are written to the specified files.

    Raises:
        FileNotFoundError: If input files cannot be found.
        pd.errors.EmptyDataError: If a BUSCO table is empty or malformed.
    """
    try:
        with open(assemblies_file, "r") as f:
            assemblies = [line.strip() for line in f]
        with open(tables_file, "r") as f:
            tables = [line.strip() for line in f]
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {e.filename}")
        sys.exit(1)

    if len(assemblies) != len(tables):
        print(f"Error: Mismatch between number of assemblies ({len(assemblies)}) and tables ({len(tables)})")
        sys.exit(1)

    busco_id_dict = {}
    all_busco_ids = set()

    # Extract BUSCO IDs from tables
    for asm, tab in zip(assemblies, tables):
        base_name = os.path.basename(asm)
        if "_" not in base_name:
            print(f"ERROR:\tUnable to parse input fasta for SPECIES_ID: {asm},\n"
                  "Please rename to have the Species and Genus written like Psilocybe cubensis B+ -> "
                  "Ps_cubensis-B+_(rest_of_file_name).fasta")
            continue
        species_id = "_".join(base_name.split("_")[:2])
        print(f"Processing assembly for: {species_id}...")
        try:
            df = pd.read_csv(tab, sep="\t")
            buscos = df[df["Status"].isin(["Complete", "Duplicated"])]["# Busco id"].tolist()
            group = "in" if in_genus in asm else "out"
            busco_id_dict[species_id] = (group, buscos)
            all_busco_ids.update(buscos)
        except FileNotFoundError:
            print(f"Error: Table file not found: {tab}")
            continue
        except pd.errors.EmptyDataError:
            print(f"Error: BUSCO table is empty or malformed: {tab}")
            continue

    if not busco_id_dict:
        print("Error: No valid BUSCO data extracted from tables")
        sys.exit(1)

    # Save all BUSCO IDs
    with open(all_busco_ids_output, "w") as f:
        f.write("\n".join(sorted(all_busco_ids)))
    print("PASS:\tAll unique BUSCO IDs saved.")

    # Build initial heatmap DataFrame
    heatmap_data = pd.DataFrame(0, index=sorted(all_busco_ids), columns=busco_id_dict.keys())
    for species_id, (_, buscos) in busco_id_dict.items():
        for busco in buscos:
            heatmap_data.loc[busco, species_id] = 1
    heatmap_data.to_csv(heatmap_df_output)

    # Save busco_id_dict as CSV
    dict_df = pd.DataFrame(
        [(k, v[0], ",".join(v[1])) for k, v in busco_id_dict.items()],
        columns=["species_id", "group", "busco_ids"]
    )
    dict_df.to_csv(busco_dict_output, index=False)
    print("PASS:\tBUSCO ID dictionary saved.")

    # Generate initial heatmap
    os.system(
        f"plot_heatmap.py {heatmap_df_output} {heatmap_png_output} "
        f"\"BUSCO ID Presence Across {compleasm_db} Assemblies (Before Filtering)\""
    )


if __name__ == "__main__":
    """Command-line entry point for extracting BUSCO IDs and generating a heatmap.

    Expects eight arguments: assemblies file, tables file, Compleasm DB, genus prefix,
    BUSCO dictionary output, all BUSCO IDs output, heatmap DataFrame output, and
    heatmap PNG output.
    """
    if len(sys.argv) != 9:
        print(
            "Usage: python extract_busco_ids.py <assemblies_file> <tables_file> "
            "<compleasm_db> <in_genus> <busco_dict_output> <all_busco_ids_output> "
            "<heatmap_df_output> <heatmap_png_output>"
        )
        sys.exit(1)
    extract_busco_ids(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
        sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]
    )