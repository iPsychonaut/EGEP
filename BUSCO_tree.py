# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 10:15:58 2025

@author: ian.bollinger@entheome.org/ian.michael.bollinger@gmail.com
"""
# REQUIRED PYTHON LIBS: biopython, modeltest-ng, mafft, raxml-ng
import os
import pandas as pd

# # inputs
CPU_THREADS = 1
RAM_GB = 1
ASSEMLBY_LIST = ["D:/TESTING_SPACE/EGAP_Test_Data/Ps_semilanceata/Ps_semilanceata_EGAP_assembly.fasta",
                 None,
                 None,
                 None,
                 None]
SPECIES_LIST = ["Ps_semilanceata",
                 None,
                 None,
                 None,
                 None]
GROPUING_LIST = ["in",
                 "in",
                 "in",
                 "out",
                 "out"]
COMPLEASM_DB = "basidiomycota"

faa_list = []
for assembly_file in ASSEMLBY_LIST:
    if ".faa" in assembly_file:
        fasta_type = "aa"
        assembly_file = assembly_file.replace(".faa", ".fasta")
    elif ".fna" in assembly_file:
        fasta_type = "nt"
        assembly_file = assembly_file.replace(".fna", ".fasta")
    if "_" in os.path.basename(assembly_file):
        species_id_list = os.path.basename(assembly_file).split("_")
        SPECIES_ID = "_".join(species_id_list[:2])
    else:
        print("Unable to parse input fasta, please rename to have the Species and Genus written like Panaeolus bisporus -> Pa_bisporus_(rest_of_file_name).fasta")

    # # # Process a single assembly's buscos
    
    # # 1) Determine path to the full_table_busco_format.tsv file. This file contains busco ID's for the gene sequences we have later
    
    # # 2) Determine path to the miniprot_output.gff. This file contains a) the sequences that after the first "\t" on the line starting with "##STA"; and b) the ID is nested between the first "\t" and the first "_" in the line starting with "##PAF"
    
    # # 3) For every sequence we want to create a new combined fasta (.faa) file named after the FINAL_ASSEMBLY and COMPLEASM_DB containing each sequence written as follows: f">{sequence_id}\n{sequence}\n"


faa_df = pd.DataFrame()
# # # Process all .faa's for each organism

# # 4) Once we have multiple combined .faa files we need to create a dataframe containing the following for each .faa file: path_to_faa, SPECIES_ID, GROUPING (either in or out)

# # 5) Based on the grouping we want to create a list of all shared busco genes, we can do this by getting the full list of all Sequence ID's in each .faa file as lists and then determine a Shared BUSCO List of sequences that ARE FOUND IN EACH .faa Sequence ID list.

# # 6) Next, for each .faa file we create a copy, rename it replacing the ".faa" with "_Shared_BUSCO_Sequence.faa", then remove any sequences NOT in the Shared BUSCO List of sequences, followed by a removal of all lines containing ID's (with a ">" on them) and then strip all "\n" from it. The goal is a single concatened sequence of ALL Shared Busco sequences that are in the same order for each FINAL_ASSEMBLY/.faa processed. Finally once the concatenated sequence is generated, add the following line as the first line f"{SPECIES_ID}_{COMPLEASM_DB}_Shared_BUSCO_Sequence\n"


# # # Combining Sequences and Model Testing

# # 7) We now have a list of all the Shared_BUSCO_Sequences.faa files, make them into a single combined.faa file (preserving sequence id's)

# # 8) MAFFT align: mafft combined_fasta > output_msa

# # 9) RAxML Tree: modeltest-ng-static -i output_msa -d aa -o output_raxml -p CPU_THREADS needs to help us set the variable DETERMINED_MODEL


# # # RAxML Tree Generation of Shared BUSCO Genes

# # 10) After model-test-ng has finished running, Build ML Starting tree with RAxML-ng: raxml-ng --msa combined_fasta --model DETERMINED_MODEL --prefix T1 --threads CPU_THREADS -r T1 -P 20

# NEED TO DETERMINE HOW OUTPUT IS BEING GENERATED SO WE CAN RENAME --prefix AND -r APPROPRIATELY

# # 11) Make 1000 bootstrap trees: raxml-ng --bootstrap --msa T1.raxml.rba --model TVM --prefix T2 --threads CPU_THREADS --bs-tree 1000 -r T2 -P 20

# NEED TO DETERMINE HOW OUTPUT IS BEING GENERATED SO WE CAN RENAME --prefix AND -r APPROPRIATELY

# # 12) Check the convergence of the boostrapping test: raxml-ng --bsconverge --bs-trees T2.raxml.bootstraps --prefix Test --threads CPU_THREADS --bs-cutoff 0.03

# NEED TO DETERMINE HOW OUTPUT IS BEING GENERATED SO WE CAN RENAME --prefix APPROPRIATELY

# # 13) Map Support values from Bootstrap test to best scoring ML tree: raxml-ng --support --tree T1.raxml.bestTree --bs-trees T2.raxml.bootstraps --prefix T3 --threads CPU_THREADS

# NEED TO DETERMINE HOW OUTPUT IS BEING GENERATED SO WE CAN RENAME --prefix APPROPRIATELY
