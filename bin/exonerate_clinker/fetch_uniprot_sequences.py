#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_uniprot_sequences.py

Updated on Tue Mar 18 2025

This script retrieves protein sequences from UniProt based on a list of gene accessions,
filters them by an in-group genus, and generates a FASTA file with labeled entries. If
nucleotide sequences are requested, it back-translates amino acid sequences using IUPAC
ambiguous bases. It is designed for integration into a Nextflow pipeline for genomic analysis.

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
import time
import io

# Simple codon table with IUPAC ambiguous bases for back-translation
AMINO_TO_NUC_AMBIGUOUS = {
    'A': 'GCN', 'C': 'TGY', 'D': 'GAY', 'E': 'GAR', 'F': 'TTY', 'G': 'GGN',
    'H': 'CAY', 'I': 'ATH', 'K': 'AAR', 'L': 'YTN', 'M': 'ATG', 'N': 'AAY',
    'P': 'CCN', 'Q': 'CAR', 'R': 'MGN', 'S': 'WSN', 'T': 'ACN', 'V': 'GTN',
    'W': 'TGG', 'Y': 'TAY', '*': 'TRR'
}

def back_translate(aa_seq):
    """Back-translate an amino acid sequence to a nucleotide sequence using IUPAC ambiguous bases."""
    nt_seq = ""
    for aa in aa_seq:
        nt_seq += AMINO_TO_NUC_AMBIGUOUS.get(aa, 'NNN')  # Default to NNN for unknown
    return nt_seq

def map_uniprot_accession(accession):
    """Map an accession to its current UniProtKB primary accession using the ID mapping service."""
    url = "https://rest.uniprot.org/idmapping/run"
    
    # Try UniProtKB_AC-ID first
    params = {"from": "UniProtKB_AC-ID", "to": "UniProtKB", "ids": accession}
    response = requests.post(url, data=params)
    if response.status_code != 200:
        print(f"ERROR: Failed to start ID mapping for {accession}: {response.status_code}", file=sys.stderr)
        return accession
    
    job_id = response.json().get("jobId")
    if not job_id:
        print(f"ERROR: No job ID returned for {accession}", file=sys.stderr)
        return accession
    
    # Poll for results
    result_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    for _ in range(10):
        result_response = requests.get(result_url)
        if result_response.status_code == 200:
            result = result_response.json()
            print(f"DEBUG: Mapping response for {accession}: {result}", file=sys.stderr)
            if "results" in result and result["results"]:
                mapped = result["results"][0]["to"]["primaryAccession"]
                print(f"DEBUG: {accession} mapped to {mapped}", file=sys.stderr)
                return mapped
            elif "failedIds" in result:
                # Try EMBL-GenBank-DDBJ as a fallback
                print(f"WARNING: {accession} not found in UniProtKB mapping, trying EMBL-GenBank-DDBJ", file=sys.stderr)
                params = {"from": "EMBL-GenBank-DDBJ", "to": "UniProtKB", "ids": accession}
                response = requests.post(url, data=params)
                if response.status_code == 200:
                    job_id = response.json().get("jobId")
                    if job_id:
                        result_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
                        for _ in range(10):
                            result_response = requests.get(result_url)
                            if result_response.status_code == 200:
                                result = result_response.json()
                                print(f"DEBUG: EMBL mapping response for {accession}: {result}", file=sys.stderr)
                                if "results" in result and result["results"]:
                                    mapped = result["results"][0]["to"]["primaryAccession"]
                                    print(f"DEBUG: {accession} mapped to {mapped} via EMBL", file=sys.stderr)
                                    return mapped
                            time.sleep(1)
                print(f"WARNING: {accession} not mapped via EMBL-GenBank-DDBJ", file=sys.stderr)
                return accession
            time.sleep(1)
        else:
            print(f"ERROR: Failed to check status for {accession}: {result_response.status_code}", file=sys.stderr)
            return accession
    print(f"ERROR: Mapping timed out for {accession} after 10 attempts", file=sys.stderr)
    return accession

def fetch_uniprot_sequence(accession, in_group, data_type="aa"):
    """Fetch sequence data from UniProt, mapping accession if needed, and filter by in-group genus."""
    mapped_accession = map_uniprot_accession(accession)
    url = f"https://rest.uniprot.org/uniprotkb/{mapped_accession}.fasta"
    response = requests.get(url)
    
    # If fetch fails and mapped_accession differs from accession, try the original accession
    if response.status_code != 200 and mapped_accession != accession:
        print(f"WARNING: Failed to fetch {mapped_accession} (status: {response.status_code}), retrying with original {accession}", file=sys.stderr)
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        response = requests.get(url)
        if response.status_code != 200:
            print(f"ERROR: Failed to fetch {accession} (mapped to {mapped_accession}) from UniProt: {response.status_code}", file=sys.stderr)
            return None
    elif response.status_code != 200:
        print(f"ERROR: Failed to fetch {accession} (mapped to {mapped_accession}) from UniProt: {response.status_code}", file=sys.stderr)
        return None
    
    fasta_content = response.text
    records = list(SeqIO.parse(io.StringIO(fasta_content), "fasta"))
    
    if not records:
        print(f"ERROR: No FASTA data found for {accession} (mapped to {mapped_accession})", file=sys.stderr)
        return None
    
    record = records[0]
    description = record.description.lower()
    if in_group.lower() not in description:
        print(f"WARNING: {in_group} not found in {accession} (mapped to {mapped_accession}) description, skipping", file=sys.stderr)
        return None
    
    if data_type == "nt":
        aa_seq = str(record.seq)
        nt_seq = back_translate(aa_seq)
        record = SeqRecord(Seq(nt_seq), id=record.id, description=record.description + " [back-translated]")
    
    return record

def generate_fasta(primary_genes, secondary_genes, in_group, data_type, output_label):
    """Generate a FASTA file with sequences from UniProt based on gene lists."""
    sequences = []
    output_path = f"{in_group}_{output_label}_{data_type}.fasta"
    
    for gene in primary_genes:
        record = fetch_uniprot_sequence(gene, in_group, data_type)
        if record:
            record.id = f"{record.id}_Primary"
            record.description = f"{record.description} [Primary]"
            sequences.append(record)
        time.sleep(1)
    
    for gene in secondary_genes:
        record = fetch_uniprot_sequence(gene, in_group, data_type)
        if record:
            record.id = f"{record.id}_Secondary"
            record.description = f"{record.description} [Secondary]"
            sequences.append(record)
        time.sleep(1)
    
    if not sequences:
        print(f"ERROR: No valid sequences retrieved for {output_path}", file=sys.stderr)
        sys.exit(1)
    
    with open(output_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    print(f"DEBUG: Wrote {len(sequences)} sequences to {output_path}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python fetch_uniprot_sequences.py <primary_genes> <secondary_genes> <in_group> <data_type> <output_label>", file=sys.stderr)
        sys.exit(1)
    
    primary_genes = sys.argv[1].split(",")
    secondary_genes = sys.argv[2].split(",")
    in_group = sys.argv[3]
    data_type = sys.argv[4]
    output_label = sys.argv[5]
    
    generate_fasta(primary_genes, secondary_genes, in_group, data_type, output_label)