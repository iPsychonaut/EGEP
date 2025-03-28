#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_assemblies.py

Updated on Fri Mar 14 2025

Filters genome assemblies based on BUSCO completeness scores for a specific compleasm_db
and outputs the total number of possible BUSCOs.

@author: ian.bollinger@entheome.org
"""

import os
import re
import sys
import json
import pandas as pd

def parse_summary(summary_file):
    """Parse a BUSCO summary file to extract counts."""
    s_count = d_count = total_buscos = None
    try:
        with open(summary_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("S:"):
                    s_count = int(re.search(r"S:\S+%,\s*(\d+)", line).group(1))
                elif line.startswith("D:"):
                    d_count = int(re.search(r"D:\S+%,\s*(\d+)", line).group(1))
                elif line.startswith("N:"):
                    total_buscos = int(re.search(r"N:(\d+)", line).group(1))
    except Exception as e:
        print(f"Error reading {summary_file}: {e}", file=sys.stderr)
    return s_count, d_count, total_buscos

def filter_assemblies(assembly_file, compiled_busco_file, compleasm_db):
    # Read assembly list
    with open(assembly_file) as f:
        assembly_paths = [line.strip() for line in f if line.strip()]

    # Read BUSCO data
    with open(compiled_busco_file) as f:
        busco_data = json.load(f)

    # Get summaries and tables for this compleasm_db
    summaries = busco_data[compleasm_db]["summaries"]
    tables = busco_data[compleasm_db]["tables"]

    # Match assemblies to summaries and tables
    data = []
    total_buscos = None
    for asm in assembly_paths:
        asm_base = os.path.basename(asm).replace(".fasta", "")
        summary = next((s for s in summaries if asm_base in s), None)
        table = next((t for t in tables if asm_base in t), None)
        if summary and table:
            s_count, d_count, total = parse_summary(summary)
            if total is not None and total_buscos is None:  # Set once, assuming consistent total
                total_buscos = total
            if s_count is not None and d_count is not None and total:
                completeness = (s_count + d_count) / total * 100
                data.append({
                    "assembly": asm,
                    "summary": summary,
                    "table": table,
                    "completeness": completeness
                })

    final_df = pd.DataFrame(data)

    # Write outputs
    with open(f"{compleasm_db}_filtered_assemblies.txt", "w") as f:
        f.write("\n".join(final_df["assembly"].tolist()))
    with open(f"{compleasm_db}_filtered_tables.txt", "w") as f:
        f.write("\n".join(final_df["table"].tolist()))
    with open(f"{compleasm_db}_removed_assemblies.csv", "w") as f:
        f.write("assembly,reason\n")
    with open(f"{compleasm_db}_total_buscos.txt", "w") as f:
        f.write(str(total_buscos or 0))

    print(f"Total BUSCOs for {compleasm_db}: {total_buscos or 0}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_assemblies.py <assembly_list> <busco_file> <compleasm_db>")
        sys.exit(1)
    filter_assemblies(sys.argv[1], sys.argv[2], sys.argv[3])