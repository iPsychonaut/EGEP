#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_busco_overlaps.py

Updated on Thu Mar 13 2025

This script filters assemblies based on BUSCO overlap and outputs a CSV with alignments
and a shared BUSCO ID list. It also generates a filtered heatmap DataFrame and assembly list
for the Entheome Genome Extraction Pipeline (EGEP).

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import sys
import pandas as pd


def filter_busco_overlaps(busco_dict_csv, all_busco_ids_file, assemblies_file,
                          compleasm_db, overlap_threshold, shared_busco_output,
                          filtered_assemblies_output, heatmap_df_output):
    """Filter assemblies by BUSCO overlap and output a CSV with shared BUSCO IDs.

    Args:
        busco_dict_csv (str): Path to the BUSCO ID dictionary CSV (alignment, busco_ids).
        all_busco_ids_file (str): Path to the file listing all unique BUSCO IDs.
        assemblies_file (str): Path to the file listing assembly paths.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        overlap_threshold (float): Quantile threshold for filtering low overlaps.
        shared_busco_output (str): Path to save the filtered shared BUSCO IDs text file.
        filtered_assemblies_output (str): Path to save the filtered assembly list.
        heatmap_df_output (str): Path to save the filtered heatmap DataFrame CSV.

    Returns:
        None: Outputs are written to the specified files.
    """
    try:
        dict_df = pd.read_csv(busco_dict_csv)
        busco_id_dict = {
            row["alignment"]: row["busco_ids"].split(",")
            for _, row in dict_df.iterrows()
        }
        with open(all_busco_ids_file, "r") as f:
            all_busco_ids = [line.strip() for line in f]
        with open(assemblies_file, "r") as f:
            assemblies = [line.strip() for line in f]
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {e.filename}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: BUSCO dictionary CSV is empty or malformed: {busco_dict_csv}", file=sys.stderr)
        sys.exit(1)

    # Initial shared BUSCO IDs
    shared_busco_ids = set.intersection(*[set(ids) for ids in busco_id_dict.values()])
    print(f"PASS:\tInitial number of BUSCO IDs shared by ALL assemblies: {len(shared_busco_ids)}", file=sys.stderr)

    # Build presence/absence matrix (using basenames for heatmap readability)
    heatmap_data = pd.DataFrame(0, index=all_busco_ids, 
                              columns=[os.path.basename(asm) for asm in busco_id_dict.keys()])
    for asm, buscos in busco_id_dict.items():
        col_name = os.path.basename(asm)
        for busco in buscos:
            heatmap_data.loc[busco, col_name] = 1

    # Calculate average BUSCO overlap
    overlap_dict = {}
    for asm in busco_id_dict.keys():
        overlaps = []
        asm_buscos = set(busco_id_dict[asm])
        for other_asm in busco_id_dict.keys():
            if other_asm != asm:
                other_buscos = set(busco_id_dict[other_asm])
                overlap = len(asm_buscos & other_buscos)
                overlaps.append(overlap)
        avg_overlap = sum(overlaps) / len(overlaps) if overlaps else 0
        overlap_dict[asm] = avg_overlap

    # Identify outliers
    overlap_series = pd.Series(overlap_dict)
    threshold = overlap_series.quantile(float(overlap_threshold))
    assembly_outliers = overlap_series[overlap_series < threshold].index.tolist()
    print(f"PASS:\tAverage BUSCO overlap threshold (bottom {overlap_threshold}): {threshold}", file=sys.stderr)
    if assembly_outliers:
        print(f"PASS:\tAssemblies with low average BUSCO overlap: {assembly_outliers}", file=sys.stderr)
    else:
        print("PASS:\tNo assemblies with significantly low BUSCO overlap.", file=sys.stderr)

    # Filter out low-overlap assemblies
    filtered_busco_id_dict = {
        k: v for k, v in busco_id_dict.items() if k not in assembly_outliers
    }
    filtered_assemblies = [asm for asm in assemblies if asm in filtered_busco_id_dict]
    print(f"PASS:\tRemoved {len(busco_id_dict) - len(filtered_busco_id_dict)} assemblies with low BUSCO overlap.", file=sys.stderr)

    # Recalculate shared BUSCO IDs after filtering
    shared_busco_ids_filtered = set.intersection(*[set(ids) for ids in filtered_busco_id_dict.values()])
    shared_busco_str = ",".join(sorted(shared_busco_ids_filtered))

    # Save shared BUSCO IDs as text file
    with open(shared_busco_output, "w") as f:
        f.write("\n".join(sorted(shared_busco_ids_filtered)))
    # Save filtered assemblies
    with open(filtered_assemblies_output, "w") as f:
        f.write("\n".join(filtered_assemblies))
    print(
        f"PASS:\tNumber of BUSCO IDs shared by ALL assemblies after filtering: "
        f"{len(shared_busco_ids_filtered)}", file=sys.stderr
    )

    # Generate filtered heatmap DataFrame
    filtered_heatmap_data = heatmap_data.drop(columns=[os.path.basename(asm) for asm in assembly_outliers], errors="ignore")
    filtered_heatmap_data.to_csv(heatmap_df_output)
    print("PASS:\tFiltered heatmap DataFrame saved.", file=sys.stderr)

    # Save CSV with alignments and shared BUSCO IDs
    filtered_df = pd.DataFrame(
        [(asm, shared_busco_str) for asm in filtered_busco_id_dict.keys()],
        columns=["alignment", "shared_busco_ids"]
    )
    filtered_df.to_csv(f"{compleasm_db}_shared_busco_dict.csv", index=False)
    print("PASS:\tShared BUSCO IDs CSV saved.", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 9:
        print(
            "Usage: python filter_busco_overlaps.py <busco_dict_csv> <all_busco_ids_file> "
            "<assemblies_file> <compleasm_db> <overlap_threshold> <shared_busco_output> "
            "<filtered_assemblies_output> <heatmap_df_output>",
            file=sys.stderr
        )
        sys.exit(1)
    filter_busco_overlaps(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
        sys.argv[6], sys.argv[7], sys.argv[8]
    )
