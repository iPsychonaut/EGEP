#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_product_dict.py

Created on Mon Mar 24 2025

This script generates a product dictionary from a FASTA file, mapping gene IDs to their
product names, and saves it as a JSON file.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys, os
from Bio import SeqIO
import json


def generate_product_dict(fasta_file, output_json):
    """Generate a product dictionary from a FASTA file and save it as JSON.

    Parses the FASTA file to extract gene IDs and product names from sequence headers,
    creating a dictionary that is then written to a JSON file.

    Parameters:
        fasta_file (str): Path to the input FASTA file.
        output_json (str): Path to the output JSON file.

    Returns:
        None: Writes the product dictionary to the specified JSON file.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
    """
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        sys.exit(1)

    product_dict = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            parts = record.id.split("|")
            if len(parts) >= 3:
                gene_id = parts[2].split("_")[0]  #e.g., "PSIK"
                desc_parts = record.description.split(" ", 3)
                if len(desc_parts) > 3:
                    product_name = desc_parts[3].split(" OS=")[0]  #e.g., "4-hydroxytryptamine kinase"
                    product_dict[gene_id] = product_name

    with open(output_json, "w") as out:
        json.dump(product_dict, out)

    print(f"PASS:\tGenerated product dictionary saved to {output_json}")


if __name__ == "__main__":
    """Command-line entry point for generating the product dictionary.

    Expects two arguments: input FASTA file and output JSON file.
    """
    if len(sys.argv) != 3:
        print("Usage: python generate_product_dict.py <fasta_file> <output_json>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    output_json = sys.argv[2]
    generate_product_dict(fasta_file, output_json)