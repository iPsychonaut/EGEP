#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_assemblies.py

Updated on Wed Mar 13 2025

This script recursively searches a directory for .fasta files containing either an
in-group or out-group string in their filenames and writes the full paths to a file.
It skips sub-folders with 'compleasm' or 'quast' in their names.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import os
import sys


def find_fasta_files(base_path, in_group, out_group):
    """Recursively find .fasta files containing in_group or out_group in their names,
    skipping sub-folders with 'compleasm' or 'quast' in their names.

    Args:
        base_path (str): Directory to search recursively.
        in_group (str): String to match in filenames (e.g., "Psilocybe").
        out_group (str): String to match in filenames (e.g., "Agrocybe").

    Returns:
        None: Writes the list of matching file paths to "in-<in_group>_out-<out_group>_assemblies.txt".

    Raises:
        FileNotFoundError: If the base_path does not exist.
        ValueError: If no matching .fasta files are found.
    """
    if not os.path.exists(base_path):
        print(f"Error: Base path does not exist: {base_path}")
        sys.exit(1)

    # List to store matching .fasta file paths
    fasta_files = []

    # Recursively walk through the directory
    for root, dirs, files in os.walk(base_path):
        # Skip directories containing 'compleasm' or 'quast' in their name
        if "compleasm" in root.lower() or "quast" in root.lower():
            continue

        for file in files:
            # Check if the file ends with .fasta and contains in_group or out_group
            if file.endswith(".fasta") and (in_group in file or out_group in file):
                full_path = os.path.join(root, file)
                fasta_files.append(full_path)

    if not fasta_files:
        print(f"Error: No .fasta files found in {base_path} containing '{in_group}' or '{out_group}'")
        sys.exit(1)

    # Construct the output filename with in_group and out_group
    output_filename = f"in-{in_group}_out-{out_group}_assemblies.txt"

    # Write the results to a file
    with open(output_filename, "w") as f:
        f.write("\n".join(fasta_files))


if __name__ == "__main__":
    """Command-line entry point for finding .fasta files.

    Expects three arguments: the base path, in-group string, and out-group string.
    """
    if len(sys.argv) != 4:
        print("Usage: python find_assemblies.py <base_path> <in_group> <out_group>")
        sys.exit(1)
    base_path = sys.argv[1]
    in_group = sys.argv[2]
    out_group = sys.argv[3]
    find_fasta_files(base_path, in_group, out_group)