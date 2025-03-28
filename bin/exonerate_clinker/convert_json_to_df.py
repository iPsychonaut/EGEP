#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_json_to_df.py

Created on Mon Mar 24 2025

This script loads a JSON file and converts its nested contents into a pandas DataFrame.
It flattens the JSON structure by iterating over contig entries and their gene records,
combining gene-specific and contig-level information. Additionally, it uses Biopython to 
determine the length of the contig (from an assembly FASTA file) referenced as the best_contig,
trim the region based on a specified buffer, and write the trimmed region to a new FASTA file.
The new sequence ID includes the basename of the input FASTA file with the extension removed.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys, os
import pandas as pd
import json
import logging
from Bio import SeqIO  # Use Biopython to read FASTA sequences
from Bio.SeqRecord import SeqRecord

def json_to_dataframe(fasta_file, json_file, trim_buffer):
    """Load a JSON file and convert its contents into a pandas DataFrame.
    
    The JSON file is expected to have a nested structure similar to:
    
    {
      "best_contig": "contig_1136_pilon_1",
      "contig_1136_pilon_1": {
        "genes": {
          "P0DPA9": {"start": 1235, "end": 1483, "bit_score": 84.3},
          "P0DPA8": {"start": 12543, "end": 13385, "bit_score": 392.0},
          ...
        },
        "contig_length": 284602,
        "region": {
          "raw_start": 1235,
          "raw_end": 16172,
          "buffered_start": 6203,
          "buffered_end": 11203,
          "buffer_size": 5000,
          "final_length": 5001
        }
      }
    }
    
    Parameters:
        fasta_file (str): Path to the assembly FASTA file.
        json_file (str): Path to the input JSON file.
        trim_buffer (int): The number of bases to trim from the start and end coordinates.
    
    Returns:
        pd.DataFrame: A DataFrame containing rows for each gene along with its associated contig info.
    
    Raises:
        FileNotFoundError: If the JSON or FASTA file does not exist.
    """
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger(__name__)

    if not os.path.exists(json_file):
        logger.error(f"Error: JSON file not found: {json_file}")
        sys.exit(1)
    
    if not os.path.exists(fasta_file):
        logger.error(f"Error: FASTA file not found: {fasta_file}")
        sys.exit(1)
    
    # Load the JSON data
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Retrieve the best_contig from the JSON
    best_contig = data.get("best_contig", None)
    if not best_contig:
        logger.error("Error: 'best_contig' not specified in JSON file")
        sys.exit(1)
    
    # Parse the assembly FASTA and get the record corresponding to best_contig
    try:
        fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        if best_contig not in fasta_records:
            logger.error(f"Error: Contig '{best_contig}' not found in FASTA file")
            sys.exit(1)
        contig_record = fasta_records[best_contig]
        fasta_length = len(contig_record.seq)
        logger.info(f"Length of best_contig ({best_contig}): {fasta_length}")
    except Exception as e:
        logger.error(f"Error reading FASTA file {fasta_file}: {e}")
        sys.exit(1)
    
    rows_list = []

    # Iterate over keys in the JSON data; skip the 'best_contig' key
    for contig_id, contig_data in data.items():
        if contig_id == "best_contig":
            continue
        if not isinstance(contig_data, dict):
            continue
        
        genes = contig_data.get("genes", {})
        contig_length = contig_data.get("contig_length", None)
        region = contig_data.get("region", {})

        for gene_id, gene_info in genes.items():
            # Create a combined record for each gene
            row = {
                "best_contig": best_contig,
                "contig_id": contig_id,
                "gene_id": gene_id,
                "start": gene_info.get("start"),
                "end": gene_info.get("end"),
                "bit_score": gene_info.get("bit_score"),
                "contig_length": contig_length
            }
            # Flatten region information into the row, if available
            for key, value in region.items():
                row[key] = value
            rows_list.append(row)
    
    df = pd.DataFrame(rows_list)

    # Determine buffered coordinates using trim_buffer
    lowest_start = df['start'].min()  # Finds the lowest value in the "start" column
    print("The lowest value in the 'start' column is:", lowest_start)
    buffered_start = lowest_start - trim_buffer
    if buffered_start < 1:
        buffered_start = 1

    highest_end = df['end'].max()
    print("The highest value in the 'end' column is:", highest_end)
    buffered_end = highest_end + trim_buffer

    if buffered_end > fasta_length:
        buffered_end = fasta_length

    logger.info(f"Buffered coordinates: start = {buffered_start}, end = {buffered_end}")
    
    # Create a FASTA file of the trimmed region from the best contig.
    # Adjust for 0-indexing: subtract 1 from start.
    trimmed_seq = contig_record.seq[buffered_start - 1:buffered_end]
    # Extract the basename of the input FASTA file and remove the extension.
    fasta_basename = os.path.basename(fasta_file)
    fasta_basename_no_ext = fasta_basename.replace(".fasta", "")
    # Create a new SeqRecord with the desired ID.
    trimmed_record = SeqRecord(
        trimmed_seq, 
        id=f"{fasta_basename_no_ext}_trimmed", 
        description=f"Trimmed region from {buffered_start} to {buffered_end}"
    )
    # Build output FASTA file path in the same folder as the JSON file.
    json_folder = os.path.dirname(os.path.abspath(json_file))
    output_fasta = os.path.join(json_folder, f"{fasta_basename_no_ext}_trimmed.fasta")
    # Write the trimmed region to the output FASTA file.
    SeqIO.write(trimmed_record, output_fasta, "fasta")
    logger.info(f"Trimmed FASTA file saved to: {output_fasta}")

    return df


if __name__ == "__main__":
    """Command-line entry point for converting a JSON file into a pandas DataFrame.

    Expects three arguments: the assembly FASTA file, the JSON file, and the trim buffer value.
    """
    if len(sys.argv) != 4:
        print("Usage: python convert_json_to_df.py <fasta_file> <json_file> <trim_buffer>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    json_file = sys.argv[2]
    try:
        trim_buffer = int(sys.argv[3])
    except ValueError:
        print("Error: trim_buffer must be an integer")
        sys.exit(1)

    df = json_to_dataframe(fasta_file, json_file, trim_buffer)
    print(df)
