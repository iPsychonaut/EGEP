#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_clinker.py

Created on Mon Mar 24 2025

DESCRIPTION

@author: ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com
"""

import subprocess
import sys 

def run_subprocess_cmd(cmd_list, shell_check):
    """
    Run a subprocess command and log its execution.

    Executes the command (as a string or list) using subprocess.Popen, streams its output,
    and logs whether it succeeded or failed.

    Parameters:
        cmd_list (str or list): The command to execute.
        shell_check (bool): Whether to execute the command through the shell.

    Returns:
        int: The return code of the subprocess.
    """
    if isinstance(cmd_list, str):
        print(f"CMD:\t{cmd_list}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    else:
        print(f"CMD:\t{' '.join(cmd_list)}")
        process = subprocess.Popen(cmd_list, shell=shell_check, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        print(f"NOTE:\tCommand failed with return code {process.returncode}")
    else:
        print(f"PASS:\tSuccessfully processed command: {' '.join(cmd_list) if isinstance(cmd_list, list) else cmd_list}")
    return process.returncode

def plot_clinker_gene_cluster(images_save, sequence_fasta, gff3_clinker_list, out_label, clinker_pid):
    """
    Generates an HTML plot visualizing gene clusters using the clinker tool.
    
    This function takes a list of gene annotations in GFF3 format along with sequence data and utilizes clinker to generate
    an interactive HTML plot illustrating the relationships and alignments between gene clusters.
    
    Parameters:
        images_save (str): Directory where the generated HTML plot will be saved.
        sequence_fasta (str): Path to the FASTA file used as the reference for gene clusters.
        gff3_clinker_list (list): List of paths to GFF3 files containing gene cluster annotations.
        clinker_pid (float): Percentage identity threshold used by clinker for clustering sequences.
    
    Returns:
        str: Path to the generated clinker matrix file used for creating dendrograms or further analysis.
    
    Examples:
        matrix_file = plot_clinker_gene_cluster("/path/to/save/images", "/path/to/sequence.fasta", ["/path/to/cluster1.gff3", "/path/to/cluster2.gff3"], 0.3)
    """    
    print("Generating clinker cluster plot...")
        
    # Generate paths for requisites and save
    save_path = f"{images_save}/{out_label}_clinker_plot.html"
    gene_functions = sequence_fasta.replace(".fasta","_gene_functions.csv")
    color_map = sequence_fasta.replace(".fasta","_color_map.csv")
    
    # Prepare the clinker command
    gff_list_string = ' '.join(gff3_clinker_list)
    clinker_cmd = ["clinker", gff_list_string,
                    "-i", str(clinker_pid),
                    "-p", save_path,
                    "-gf", gene_functions,
                    "-cm", color_map] 

    # Run the clinker command and check returncode for errors
    run_subprocess_cmd(clinker_cmd, shell_check=False)

    print(f"PASS:\tSuccessfully generated clinker cluster plot: {save_path}")

if __name__ == "__main__":
    """Command-line entry point for generating the clinker plot.

    """
    if len(sys.argv) != 3:
        print("Usage: python plot_clinker.py <images_save> <sequence_fasta> <gff3_clinker_list> <out_label> <clinker_pid>")
        sys.exit(1)
    images_save = sys.argv[1]
    sequence_fasta = sys.argv[2]
    gff3_clinker_list = sys.argv[3]
    out_label = sys.argv[4]
    clinker_pid = sys.argv[5]
    plot_clinker_gene_cluster(images_save, sequence_fasta, gff3_clinker_list, out_label, clinker_pid)