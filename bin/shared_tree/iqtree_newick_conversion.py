# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 14:36:33 2025

iqtree_newick_conversion.py

@author: theda
"""

from Bio import Phylo
from io import StringIO
import sys

def extract_newick_from_iqtree(iqtree_file):
    output_file = iqtree_file.replace(".iqtree",".nwk")
    # Read the .iqtree file
    with open(iqtree_file, 'r') as f:
        lines = f.readlines()
        print(f"Total lines in file: {len(lines)}")  # Debug: Check if file is read
        for i, line in enumerate(lines):
            # Flexible search: case-insensitive, ignores extra spaces
            if "tree in newick format" in line.lower().strip():
                print(f"Found 'Tree in newick format:' at line {i}")
                # Skip any empty lines after the header
                j = i + 1
                newick_tree = ""
                while j < len(lines):
                    line = lines[j].strip()
                    if line:  # If the line is not empty
                        newick_tree = line
                        break
                    j += 1
                print("Extracted Newick tree:", newick_tree)  # Inspect the string
                # Validate the Newick string
                if not newick_tree:
                    raise ValueError("No Newick tree found after 'Tree in newick format:'.")
                if not newick_tree.endswith(';'):
                    raise ValueError("Extracted Newick tree does not end with ';', might be invalid.")
                break
        else:
            raise ValueError("Newick tree not found in the file.")

    # Parse the Newick string into a tree object using Bio.Phylo
    tree = Phylo.read(StringIO(newick_tree), 'newick')

    # Save the tree to a new file in Newick format
    Phylo.write(tree, output_file, 'newick')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 iqtree_newick_conversion.py <iqtree_file>", file=sys.stderr)
        sys.exit(1)
    extract_newick_from_iqtree(sys.argv[1])