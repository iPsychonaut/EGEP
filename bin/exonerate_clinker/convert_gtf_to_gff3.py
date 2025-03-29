#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_gtf_to_gff3.py

Created on Mon Mar 24 2025

This script converts an Exonerate GTF file to GFF3 format, adding mRNA, exon, and CDS features,
and uses a dynamically loaded product dictionary for gene annotations.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys, os
import pandas as pd
import json
import logging


def parse_attributes(attributes_str):
    """Parse GTF attributes string into a dictionary.

    Parameters:
        attributes_str (str): Semicolon-separated attributes string from GTF.

    Returns:
        dict: Dictionary of attribute key-value pairs.
    """
    attributes_dict = {}
    try:
        for attribute in attributes_str.split(';'):
            if attribute:
                key, value = attribute.strip().split('=')
                attributes_dict[key.strip()] = value.strip()
    except Exception:
        pass
    return attributes_dict


def calculate_phase(start, previous_end):
    """Calculate phase for CDS features based on start and previous end positions.

    Parameters:
        start (int): Start position of the current CDS feature.
        previous_end (int): End position of the previous CDS feature.

    Returns:
        int: Calculated phase (0, 1, or 2).
    """
    if previous_end < 0:
        return 0
    length = start - previous_end - 1
    return (3 - (length % 3)) % 3


def convert_gtf_to_gff3(gtf_file, gff3_file, product_dict_file):
    """Convert a GTF file to GFF3 format with hierarchical features.

    Reads a GTF file, adds mRNA features, updates attributes with a product dictionary,
    calculates CDS phase, and writes the result as a GFF3 file.

    Parameters:
        gtf_file (str): Path to the input GTF file.
        gff3_file (str): Path to the output GFF3 file.
        product_dict_file (str): Path to the JSON file containing the product dictionary.

    Returns:
        None: Writes the converted data to the specified GFF3 file.

    Raises:
        FileNotFoundError: If the GTF or product dictionary file does not exist.
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger(__name__)

    if not os.path.exists(gtf_file):
        logger.error(f"Error: GTF file not found: {gtf_file}")
        sys.exit(1)
    if not os.path.exists(product_dict_file):
        logger.error(f"Error: Product dictionary file not found: {product_dict_file}")
        sys.exit(1)

    # Load product_dict from JSON
    with open(product_dict_file, 'r') as f:
        product_dict = json.load(f)

    # Read GTF file into DataFrame with explicit column names
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    rows_list = []
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    while len(fields) < 9:
                        fields.insert(-1, '.')
                    row_data = dict(zip(gtf_columns, fields))
                    rows_list.append(row_data)
        if not rows_list:
            logger.warning(f"GTF file {gtf_file} is empty. Creating minimal GFF3.")
            with open(gff3_file, 'w') as f:
                f.write("##gff-version 3\n")
            return
    except Exception as e:
        logger.error(f"Failed to read GTF file {gtf_file}: {e}")
        sys.exit(1)
    
    df = pd.DataFrame(rows_list)

    # Add mRNA features and process the DataFrame
    new_rows = []
    for _, row in df.iterrows():
        new_rows.append(row.to_dict())
        if row['feature'] == 'gene':
            new_row = row.to_dict()
            new_row['feature'] = 'mRNA'
            new_rows.append(new_row)

    df = pd.DataFrame(new_rows)
    df['feature'] = df['feature'].str.replace('cds', 'CDS', case=False)

    # Process attributes and calculate phase
    gene_list = []
    current_gene_name = ''
    current_mrna_id = ''
    previous_cds_end = -1

    for index, row in df.iterrows():
        feature_type = row['feature']
        attributes_dict = parse_attributes(row["attributes"])
        
        if feature_type == 'gene':
            current_gene_name = row["attributes"].split(' ; ')[1].split(' ')[-1].split('.')[0]
            gene_list.append(current_gene_name)
            gene_counter = gene_list.count(current_gene_name)
            current_gene_id = f"{current_gene_name}.{gene_counter}"
            attributes_dict['ID'] = current_gene_id
        elif feature_type == 'mRNA':
            current_mrna_id = f'{current_gene_id}-T1'
            attributes_dict['ID'] = current_mrna_id
            attributes_dict['Parent'] = current_gene_id
            attributes_dict['product'] = product_dict.get(current_gene_name, current_gene_name)
        elif feature_type == 'exon':
            exon_id = f'{current_mrna_id}.exon{gene_counter}'
            attributes_dict['ID'] = exon_id
            attributes_dict['Parent'] = current_mrna_id
        elif feature_type == 'CDS':
            cds_id = f'{current_mrna_id}.CDS'
            attributes_dict['ID'] = cds_id
            attributes_dict['Parent'] = current_mrna_id
            start_position = int(row['start'])
            phase = calculate_phase(start_position, previous_cds_end)
            df.at[index, 'frame'] = str(phase)
            previous_cds_end = int(row['end'])
            attributes_dict['frame'] = row['frame']
            attributes_dict['score'] = row['score']

        attributes_str = ';'.join([f'{k}={v}' for k, v in attributes_dict.items()])
        df.at[index, 'attributes'] = attributes_str

    # Filter desired features
    df = df[df['feature'].isin(['gene', 'mRNA', 'exon', 'CDS'])]

    # Write to GFF3 file
    with open(gff3_file, 'w') as f:
        f.write("##gff-version 3\n")
        df.to_csv(f, sep='\t', header=False, index=False, quoting=None)
    logger.info(f"PASS:\tConverted {gtf_file} to {gff3_file}")


if __name__ == "__main__":
    """Command-line entry point for converting GTF to GFF3.

    Expects three arguments: input GTF file, output GFF3 file, and product dictionary JSON file.
    """
    if len(sys.argv) != 4:
        print("Usage: python convert_gtf_to_gff3.py <gtf_file> <gff3_file> <product_dict_file>")
        sys.exit(1)
    gtf_file = sys.argv[1]
    gff3_file = sys.argv[2]
    product_dict_file = sys.argv[3]
    convert_gtf_to_gff3(gtf_file, gff3_file, product_dict_file)