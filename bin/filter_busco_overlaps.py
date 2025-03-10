#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_busco_overlaps.py

Created on Sun Mar 9 18:44:50 2025

This script filters assemblies based on BUSCO overlap, removing those with low
average overlap, and generates a filtered heatmap for the Entheome Genome Extraction
Pipeline (EGEP). It outputs shared BUSCO IDs and filtered assembly lists.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import sys
import pandas as pd


def filter_busco_overlaps(busco_dict_csv, all_busco_ids_file, assemblies_file,
                          compleasm_db, overlap_threshold, shared_busco_output,
                          filtered_assemblies_output, heatmap_df_output,
                          heatmap_png_output):
    """Filter assemblies by BUSCO overlap and generate a filtered heatmap.

    Loads a BUSCO ID dictionary, calculates average BUSCO overlaps, filters out
    assemblies below a quantile threshold, recalculates shared BUSCO IDs, and
    generates a filtered heatmap.

    Args:
        busco_dict_csv (str): Path to the BUSCO ID dictionary CSV.
        all_busco_ids_file (str): Path to the file listing all unique BUSCO IDs.
        assemblies_file (str): Path to the file listing assembly paths.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        overlap_threshold (float): Quantile threshold for filtering low overlaps
                                   (e.g., 0.3 for bottom 30%).
        shared_busco_output (str): Path to save the filtered shared BUSCO IDs.
        filtered_assemblies_output (str): Path to save the filtered assembly list.
        heatmap_df_output (str): Path to save the filtered heatmap DataFrame CSV.
        heatmap_png_output (str): Path to save the filtered heatmap PNG.

    Returns:
        None: Outputs are written to the specified files.

    Raises:
        FileNotFoundError: If input files cannot be found.
        pd.errors.EmptyDataError: If the BUSCO dictionary CSV is empty or malformed.
    """
    try:
        dict_df = pd.read_csv(busco_dict_csv)
        busco_id_dict = {
            row["species_id"]: (row["group"], row["busco_ids"].split(","))
            for _, row in dict_df.iterrows()
        }
        with open(all_busco_ids_file, "r") as f:
            all_busco_ids = [line.strip() for line in f]
        with open(assemblies_file, "r") as f:
            assemblies = [line.strip() for line in f]
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {e.filename}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: BUSCO dictionary CSV is empty or malformed: {busco_dict_csv}")
        sys.exit(1)

    # Initial shared BUSCO IDs
    shared_busco_ids = set.intersection(*[set(ids) for _, ids in busco_id_dict.values()])
    print(f"PASS:\tInitial number of BUSCO IDs shared by ALL species: {len(shared_busco_ids)}")

    # Build presence/absence matrix
    heatmap_data = pd.DataFrame(0, index=all_busco_ids, columns=busco_id_dict.keys())
    for species_id, (_, buscos) in busco_id_dict.items():
        for busco in buscos:
            heatmap_data.loc[busco, species_id] = 1

    # Calculate average BUSCO overlap
    overlap_dict = {}
    for species_id in busco_id_dict.keys():
        overlaps = []
        species_buscos = set(busco_id_dict[species_id][1])
        for other_id in busco_id_dict.keys():
            if other_id != species_id:
                other_buscos = set(busco_id_dict[other_id][1])
                overlap = len(species_buscos & other_buscos)
                overlaps.append(overlap)
        avg_overlap = sum(overlaps) / len(overlaps)
        overlap_dict[species_id] = avg_overlap

    # Identify outliers
    overlap_series = pd.Series(overlap_dict)
    threshold = overlap_series.quantile(float(overlap_threshold))
    assembly_outliers = overlap_series[overlap_series < threshold].index.tolist()
    print(f"Average BUSCO overlap threshold (bottom {overlap_threshold}): {threshold}")
    if assembly_outliers:
        print(f"Assemblies with low average BUSCO overlap: {assembly_outliers}")
    else:
        print("No assemblies with significantly low BUSCO overlap.")

    # Filter out low-overlap assemblies
    filtered_busco_id_dict = {
        k: v for k, v in busco_id_dict.items() if k not in assembly_outliers
    }
    filtered_assemblies = [
        asm for asm in assemblies
        if "_".join(os.path.basename(asm).split("_")[:2]) in filtered_busco_id_dict
    ]
    print(f"Removed {len(busco_id_dict) - len(filtered_busco_id_dict)} assemblies with low BUSCO overlap.")

    # Recalculate shared BUSCO IDs after filtering
    shared_busco_ids_filtered = set.intersection(*[set(ids) for _, ids in filtered_busco_id_dict.values()])
    with open(shared_busco_output, "w") as f:
        f.write("\n".join(sorted(shared_busco_ids_filtered)))
    with open(filtered_assemblies_output, "w") as f:
        f.write("\n".join(filtered_assemblies))
    print(
        f"PASS:\tNumber of BUSCO IDs shared by ALL species after filtering: "
        f"{len(shared_busco_ids_filtered)}"
    )

    # Generate filtered heatmap
    filtered_heatmap_data = heatmap_data.drop(columns=assembly_outliers, errors="ignore")
    filtered_heatmap_data.to_csv(heatmap_df_output)
    os.system(
        f"plot_heatmap.py {heatmap_df_output} {heatmap_png_output} "
        f"\"BUSCO ID Presence Across {compleasm_db} Assemblies (After Filtering)\""
    )


if __name__ == "__main__":
    """Command-line entry point for filtering BUSCO overlaps and generating a heatmap.

    Expects nine arguments: BUSCO dict CSV, all BUSCO IDs file, assemblies file,
    Compleasm DB, overlap threshold, shared BUSCO output, filtered assemblies output,
    heatmap DataFrame output, and heatmap PNG output.
    """
    if len(sys.argv) != 10:
        print(
            "Usage: python filter_busco_overlaps.py <busco_dict_csv> <all_busco_ids_file> "
            "<assemblies_file> <compleasm_db> <overlap_threshold> <shared_busco_output> "
            "<filtered_assemblies_output> <heatmap_df_output> <heatmap_png_output>"
        )
        sys.exit(1)
    filter_busco_overlaps(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
        sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9]
    )