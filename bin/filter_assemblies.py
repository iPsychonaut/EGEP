#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_assemblies.py

Updated on Wed Mar 13 2025

This script filters genome assemblies based on BUSCO completeness scores, selecting
the best assembly per species and removing outliers using interquartile range (IQR)
analysis. It retains assemblies with completeness > 80% even if they are outliers,
matches summaries and tables by assembly basename in their paths, and saves removal
reasons to a CSV.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import os
import re
import sys
import json
import pandas as pd


def parse_summary(summary_file):
    """Parse a BUSCO summary file to extract single, duplicated, and total BUSCO counts."""
    s_count = None
    d_count = None
    total_buscos = None
    try:
        with open(summary_file, "r") as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("S:"):
                match = re.search(r"S:\S+%,\s*(\d+)", line)
                if match:
                    s_count = int(match.group(1))
            elif line.startswith("D:"):
                match = re.search(r"D:\S+%,\s*(\d+)", line)
                if match:
                    d_count = int(match.group(1))
            elif line.startswith("N:"):
                match = re.search(r"N:(\d+)", line)
                if match:
                    total_buscos = int(match.group(1))
    except FileNotFoundError:
        print(f"Error: Summary file not found: {summary_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading {summary_file}: {e}", file=sys.stderr)
    return s_count, d_count, total_buscos


def filter_assemblies(assembly_file, compiled_busco_file):
    """Filter assemblies, keep those with >80% completeness despite outliers, and save removal reasons.

    Args:
        assembly_file (str): Path to file containing list of assembly paths (from Nextflow Process 2).
        compiled_busco_file (str): Path to compiled BUSCO files dictionary (from Nextflow Process 3).

    Returns:
        None: Outputs filtered assemblies and tables to files, removed assemblies to CSV.
    """
    # Read assembly list
    with open(assembly_file, "r") as f:
        assembly_paths = [line.strip() for line in f if line.strip()]

    if not assembly_paths:
        print("Error: No assemblies provided", file=sys.stderr)
        sys.exit(1)

    # Read compiled BUSCO file
    try:
        with open(compiled_busco_file, "r") as f:
            busco_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: Failed to parse {compiled_busco_file} as JSON: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Compiled BUSCO file not found: {compiled_busco_file}", file=sys.stderr)
        sys.exit(1)

    # Extract all summary and table paths
    all_summary_paths = []
    all_table_paths = []
    compleasm_dbs = list(busco_data.keys())
    for db in compleasm_dbs:
        if "summaries" in busco_data[db] and "tables" in busco_data[db]:
            all_summary_paths.extend(busco_data[db]["summaries"])
            all_table_paths.extend(busco_data[db]["tables"])
        else:
            print(f"Warning: Missing 'summaries' or 'tables' for {db} in {compiled_busco_file}", file=sys.stderr)

    # Create a DataFrame by matching assembly basenames to summary/table paths
    data = []
    for asm in assembly_paths:
        asm_basename = os.path.basename(asm).replace(".fasta", "")
        matching_summaries = [s for s in all_summary_paths if asm_basename in s]
        matching_tables = [t for t in all_table_paths if asm_basename in t]
        if not matching_summaries or not matching_tables:
            print(f"Warning: No matching summary or table found for {asm}", file=sys.stderr)
            continue
        # Take the first match for simplicity
        data.append({
            "assembly": asm,
            "summary": matching_summaries[0],
            "table": matching_tables[0]
        })

    if not data:
        print("Error: No valid assembly-summary-table mappings found", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(data)

    # Track removed assemblies and reasons
    removed_records = []

    # Filter assemblies by species, keeping the best
    processed_species = set()
    assemblies_to_remove = []
    for _, row in df.iterrows():
        assembly_file = row["assembly"]
        base_name = os.path.basename(assembly_file)
        if "_" not in base_name:
            print(f"Warning: Skipping malformed filename: {assembly_file}", file=sys.stderr)
            removed_records.append({"assembly": assembly_file, "reason": "Malformed filename"})
            assemblies_to_remove.append(assembly_file)
            continue
        species_id = "_".join(base_name.split("_")[:2])
        if species_id in processed_species:
            continue
        processed_species.add(species_id)

        species_df = df[df["assembly"].str.contains(species_id)]
        if len(species_df) > 1:
            busco_data = []
            for _, s_row in species_df.iterrows():
                asm, summ = s_row["assembly"], s_row["summary"]
                s_count, d_count, total_buscos = parse_summary(summ)
                if s_count is None or d_count is None or total_buscos is None:
                    print(f"Warning: Could not parse BUSCO counts or total for {summ}", file=sys.stderr)
                    removed_records.append({"assembly": asm, "reason": "Failed BUSCO parsing"})
                    assemblies_to_remove.append(asm)
                    continue
                completeness = s_count + d_count
                completeness_percent = (completeness / total_buscos) * 100 if total_buscos > 0 else 0
                busco_data.append({
                    "assembly": asm,
                    "summary": summ,
                    "table": s_row["table"],
                    "S": s_count,
                    "D": d_count,
                    "completeness": completeness,
                    "completeness_percent": completeness_percent,
                    "total_buscos": total_buscos
                })

            if not busco_data:
                continue

            species_busco_df = pd.DataFrame(busco_data)
            species_busco_df_sorted = species_busco_df.sort_values(by="completeness", ascending=False)
            best_assembly = species_busco_df_sorted.iloc[0]["assembly"]
            for _, s_row in species_busco_df_sorted.iterrows():
                if s_row["assembly"] != best_assembly:
                    assemblies_to_remove.append(s_row["assembly"])
                    removed_records.append({
                        "assembly": s_row["assembly"],
                        "reason": f"Not best for species {species_id} (completeness: {s_row['completeness']}/{s_row['total_buscos']}, {s_row['completeness_percent']:.1f}%)"
                    })

    # Apply initial filtering
    filtered_df = df[~df["assembly"].isin(assemblies_to_remove)]

    # Remove outliers based on completeness IQR, but keep if >80%
    final_busco_data = []
    for _, row in filtered_df.iterrows():
        asm, summ, tab = row["assembly"], row["summary"], row["table"]
        s_count, d_count, total_buscos = parse_summary(summ)
        if s_count is None or d_count is None or total_buscos is None:
            removed_records.append({"assembly": asm, "reason": "Failed BUSCO parsing in final check"})
            continue
        completeness = s_count + d_count
        completeness_percent = (completeness / total_buscos) * 100 if total_buscos > 0 else 0
        final_busco_data.append({
            "assembly": asm,
            "summary": summ,
            "table": tab,
            "completeness": completeness,
            "completeness_percent": completeness_percent,
            "total_buscos": total_buscos
        })

    if not final_busco_data:
        print("Error: No valid BUSCO data available after initial filtering", file=sys.stderr)
        sys.exit(1)

    final_busco_df = pd.DataFrame(final_busco_data)
    q1 = final_busco_df["completeness"].quantile(0.25)
    q3 = final_busco_df["completeness"].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    outliers = final_busco_df[
        (final_busco_df["completeness"] < lower_bound) |
        (final_busco_df["completeness"] > upper_bound)
    ]

    # Filter outliers, but keep those with completeness > 80%
    if not outliers.empty:
        for _, row in outliers.iterrows():
            if row["completeness_percent"] <= 80:
                removed_records.append({
                    "assembly": row["assembly"],
                    "reason": f"Outlier (completeness: {row['completeness']}/{row['total_buscos']}, {row['completeness_percent']:.1f}%)"
                })
            else:
                print(f"Note: Keeping outlier {row['assembly']} with completeness {row['completeness']}/{row['total_buscos']} ({row['completeness_percent']:.1f}% > 80%)", file=sys.stderr)

    # Keep assemblies unless they’re outliers with completeness <= 80%
    final_busco_df_no_outliers = final_busco_df[
        ~((final_busco_df["completeness"] < lower_bound) | 
          (final_busco_df["completeness"] > upper_bound)) |
        (final_busco_df["completeness_percent"] > 80)
    ]

    if final_busco_df_no_outliers.empty:
        print("Error: No assemblies remain after outlier removal", file=sys.stderr)
        sys.exit(1)

    # Write filtered results to files
    with open("filtered_assemblies.txt", "w") as f:
        f.write("\n".join(final_busco_df_no_outliers["assembly"].tolist()))
    with open("filtered_tables.txt", "w") as f:
        f.write("\n".join(final_busco_df_no_outliers["table"].tolist()))

    # Save removed assemblies to CSV
    if removed_records:
        removed_df = pd.DataFrame(removed_records)
        removed_df.to_csv("removed_assemblies.csv", index=False)
    else:
        pd.DataFrame(columns=["assembly", "reason"]).to_csv("removed_assemblies.csv", index=False)


if __name__ == "__main__":
    """Command-line entry point for filtering assemblies."""
    if len(sys.argv) != 3:
        print("Usage: python filter_assemblies.py <assembly_list_file> <compiled_busco_file>", file=sys.stderr)
        sys.exit(1)
    assembly_file = sys.argv[1]
    compiled_busco_file = sys.argv[2]
    filter_assemblies(assembly_file, compiled_busco_file)
