#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_buscos.py

Created on Sun Mar 9 18:44:50 2025

This script processes a list of FASTA files to check or generate BUSCO files using Compleasm,
producing individual and compiled output lists.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys
import subprocess
import json
from pathlib import Path
import pandas as pd


def find_file(filename, folder=None):
    """
    Search the filesystem for a file with the given name.

    Walks the filesystem from the root (based on the OS) and returns the absolute path
    of the first occurrence of the file, excluding paths containing "$RECYCLE.BIN".

    Parameters:
        filename (str): The file name to search for.
        folder (str, optional): The starting directory for the search. If None, it defaults based on OS.

    Returns:
        str or None: The absolute path if found; otherwise, None.
    """
    if pd.notna(folder):
        root_directory = folder        
    else:
        root_directory = "/"
    for root, dirs, files in os.walk(root_directory):
        # Skip paths containing "$RECYCLE.BIN"
        if "$RECYCLE.BIN" in root:
            continue
        if "ncbi" in dirs:
            continue

        if filename in files:
            return os.path.join(root, filename)
    return None


def generate_buscos(assembly_list_path, compleasm_db, cpu_threads):
    """Generate BUSCO files for assemblies and compile output lists.

    Args:
        assembly_list_path (str): Path to the file containing the list of FASTA files.
        compleasm_db (str): Compleasm database identifier (e.g., "basidiomycota").
        cpu_threads (int): Number of CPU threads to use for Compleasm.

    Returns:
        None: Writes output to "${compleasm_db}_busco_files.txt" and a temp file "temp_${compleasm_db}_compiled_busco_files.json".

    Raises:
        FileNotFoundError: If the assembly list file or required tools are missing.
        RuntimeError: If Compleasm fails or no BUSCO files are generated.
    """
    # Validate inputs
    if not os.path.exists(assembly_list_path):
        print(f"Error: Assembly list file not found: {assembly_list_path}")
        sys.exit(1)

    # Read the assembly list
    with open(assembly_list_path, "r") as f:
        assembly_paths = [line.strip() for line in f if line.strip()]
    
    if not assembly_paths:
        print(f"Error: Assembly list is empty: {assembly_list_path}")
        sys.exit(1)

    # Prepare output files
    output_file = f"{compleasm_db}_busco_files.txt"
    temp_compiled_file = f"temp_{compleasm_db}_compiled_busco_files.json"  # Updated to .json
    
    # Lists to store file paths for this compleasm_db
    summary_files = []
    table_files = []

    # Process each assembly
    for assembly_path in assembly_paths:
        if not os.path.exists(assembly_path):
            print(f"Warning: Assembly file not found, skipping: {assembly_path}")
            continue

        # Derive the base assembly name and potential Compleasm output directory
        base_name = os.path.basename(assembly_path).replace(".fasta", "")
        parent_dir = os.path.dirname(assembly_path)
        
        # Check for existing Compleasm output directory based on the base assembly name
        expected_output_dir = os.path.join(parent_dir, f"{base_name}_{compleasm_db}_odb10_compleasm")
        summary_file = os.path.join(expected_output_dir, "summary.txt")
        table_file = os.path.join(expected_output_dir, f"{compleasm_db}_odb10", "full_table_busco_format.tsv")

        # If both files exist, use them; otherwise, run Compleasm
        if os.path.exists(summary_file) and os.path.exists(table_file):
            print(f"Found existing BUSCO files for {assembly_path} with {compleasm_db}")
            summary_files.append(summary_file)
            table_files.append(table_file)
        else:
            print(f"Generating BUSCO files for {assembly_path} with {compleasm_db}")
            os.makedirs(expected_output_dir, exist_ok=True)
            compleasm_path = find_file("compleasm.py")
            print(compleasm_path)
            cmd = [
                "python3", compleasm_path, "run",
                "-a", assembly_path,
                "-o", expected_output_dir,
                "-t", str(cpu_threads),
                "-l", f"{compleasm_db}_odb10"
            ]
            try:
                result = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
                print(f"Compleasm output: {result.stderr}")
            except subprocess.CalledProcessError as e:
                print(f"Error: Compleasm failed for {assembly_path} with {compleasm_db}: {e.stderr}")
                sys.exit(1)

            # Verify output files were created
            if not (os.path.exists(summary_file) and os.path.exists(table_file)):
                print(f"Error: Compleasm failed to generate expected output for {assembly_path} with {compleasm_db}")
                sys.exit(1)
            
            summary_files.append(summary_file)
            table_files.append(table_file)

    # Write individual output file
    with open(output_file, "w") as f:
        for summary, table in zip(summary_files, table_files):
            f.write(f"{summary}\n{table}\n")
    
    if not summary_files:
        print(f"Error: No BUSCO files generated or found for {compleasm_db}")
        sys.exit(1)

    # Prepare JSON data in the specified format
    compiled_data = {
        "summaries": summary_files,
        "tables": table_files
    }

    # Write to a temporary compiled JSON file
    with open(temp_compiled_file, "w") as f:
        json.dump(compiled_data, f, indent=4)

    print(f"PASS:\tProcessed BUSCO files for {compleasm_db}. Paths written to {output_file} and {temp_compiled_file}")


if __name__ == "__main__":
    """Command-line entry point for generating BUSCO files.

    Expects three arguments: assembly list file, compleasm database, and CPU threads.
    """
    if len(sys.argv) != 4:
        print("Usage: python generate_buscos.py <assembly_list_path> <compleasm_db> <cpu_threads>")
        sys.exit(1)
    assembly_list_path = sys.argv[1]
    compleasm_db = sys.argv[2]
    cpu_threads = int(sys.argv[3])
    generate_buscos(assembly_list_path, compleasm_db, cpu_threads)
