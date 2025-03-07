# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 10:15:58 2025

@author: ian.bollinger@entheome.org/ian.michael.bollinger@gmail.com

EGEP (Entheome Genome Extraction Pipeline) is a versatile bioinformatics pipeline
for annotating gene regions and extracting them for graphical and phylogenetic analysis.
"""
# REQUIRED PYTHON LIBS: pandas, biopython, mafft, IQ-tree
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import subprocess


def run_subprocess_cmd(cmd_list, shell_check):
    """
    Run a subprocess command and log its execution.

    Executes the command (as a string or list) using subprocess.Popen, streams its output,
    and logs whether it succeeded or failed.

    Parameters:
        cmd_list (str or list): The command to execute.
        shell_check (bool): Whether to execute the command through the shell.

    Returns:
        process.returncode (int): The return code of the subprocess.
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


def fasta_ext_fixer(assembly_file):
    """
    DESCRIPTION.

    Parameters:
    assembly_file (TYPE): DESCRIPTION.

    Returns:
    assembly_file (TYPE): DESCRIPTION.
    fasta_type (TYPE): DESCRIPTION.

    """
    if ".faa" in assembly_file:
        fasta_type = "aa"
        assembly_file = assembly_file.replace(".faa", ".fasta")
    elif ".fna" in assembly_file:
        fasta_type = "nt"
        assembly_file = assembly_file.replace(".fna", ".fasta")       
    else:
        fasta_type = None
        assembly_file = assembly_file
    return assembly_file, fasta_type


if __name__ == "__main__":
    print("TEST BUSCO TREE")
    

    CPU_THREADS = 16
    
    ASSEMBLY_LIST = [
                     ]

    GROUP_LIST = ["in",
                  "in",
                  "in",
                  "in",
                  "out"]
    COMPLEASM_LIST = ["basidiomycota", "agaricales"]
    
    
    
    for COMPLEASM_DB in COMPLEASM_LIST:
        
        print("Generating BUSCO ID list shared by ALL species...")
        busco_id_dict = {}
        for index, assembly_file in enumerate(ASSEMBLY_LIST):
            assembly_file, fasta_type = fasta_ext_fixer(assembly_file)
            if "_" in os.path.basename(assembly_file):
                species_id_list = os.path.basename(assembly_file).split("_")
                SPECIES_ID = "_".join(species_id_list[:2])
            else:
                print(f"ERROR:\tUnable to parse input fasta for SPECIES_ID: {assembly_file},\n"
                      "Please rename to have the Species and Genus written like Psilocybe cubensis B+ -> "
                      "Ps_cubensis-B+_(rest_of_file_name).fasta")
                continue
            print(f"Processing Final Assembly & GFF files for: {SPECIES_ID}...")
            compleasm_out_dir = assembly_file.replace(".fasta",f"_{COMPLEASM_DB}_odb10_compleasm")
            full_table_busco_format_tsv = os.path.join(compleasm_out_dir, f"{COMPLEASM_DB}_odb10", "full_table_busco_format.tsv")
            full_table_busco_df = pd.read_csv(full_table_busco_format_tsv, sep='\t')
            filtered_df = full_table_busco_df[full_table_busco_df['Status'].isin(['Complete', 'Duplicated'])]
            busco_id_dict[SPECIES_ID] = (GROUP_LIST[index], filtered_df["# Busco id"].tolist())
        print("PASS:\tBUSCO ID dictionary generated.")
        # return busco_id_dict
        
        
        
        shared_busco_ids = None
        for species_id, (group_type, busco_ids) in busco_id_dict.items():
            if shared_busco_ids is None:
                shared_busco_ids = set(busco_ids)
            else:
                shared_busco_ids.intersection_update(busco_ids)
        shared_busco_ids_list = sorted(list(shared_busco_ids))
        print("PASS:\tNumber of BUSCO IDs shared by ALL species.")
        # return shared_busco_ids_list
        
        
        
        print("Generating ordered BUSCO gene sequences for ALL species...")
        busco_id_fasta_list = []
        for index, assembly_file in enumerate(ASSEMBLY_LIST):
            assembly_file, fasta_type = fasta_ext_fixer(assembly_file)
            base_name = os.path.basename(assembly_file)
            if "_" in base_name:
                species_id_list = base_name.split("_")
                SPECIES_ID = "_".join(species_id_list[:2])
            else:
                print(f"ERROR:\tUnable to parse input fasta for SPECIES_ID: {assembly_file},\n"
                      "NOTE:\tPlease rename it (e.g. 'Ps_cubensis-B+_...')")
                continue
            compleasm_out_dir = assembly_file.replace(".fasta", f"_{COMPLEASM_DB}_odb10_compleasm")
            miniprot_output_gff = os.path.join(compleasm_out_dir, f"{COMPLEASM_DB}_odb10", "miniprot_output.gff")
            out_faa = os.path.join(compleasm_out_dir, f"{SPECIES_ID}_shared_busco.fasta")
            if not os.path.isfile(miniprot_output_gff):
                print(f"NOTE:\tNo GFF file found at: {miniprot_output_gff}")
                continue
            with open(miniprot_output_gff, "r") as f:
                lines = f.read().splitlines()
            records_dict = {}
            i = 0
            n = len(lines)
            while i < n:
                line = lines[i].strip()
                if line.startswith("##PAF"):
                    fields = line.split("\t")
                    if len(fields) > 1:
                        paf_id = fields[1].strip().split("_")[0]
                    else:
                        i += 1
                        continue
                    matched_busco = None
                    for busco in shared_busco_ids_list:
                        if busco in paf_id:
                            matched_busco = busco
                            break
                    i += 1
                    if i < n and lines[i].startswith("##STA"):
                        sta_line = lines[i].strip()
                        prot_seq = sta_line[5:].strip()
                        i += 1
                        block_lines = []
                        while i < n and not lines[i].startswith("##"):
                            block_lines.append(lines[i])
                            i += 1
                        if matched_busco is None:
                            block_string = " ".join(block_lines)
                            for busco in shared_busco_ids_list:
                                if busco in block_string:
                                    matched_busco = busco
                                    break
                    if matched_busco is not None:
                            rec = SeqRecord(Seq(prot_seq),
                                            id=paf_id,
                                            description=f"Shared_{COMPLEASM_DB}_BUSCO:{SPECIES_ID}|{matched_busco}")
                            records_dict[matched_busco] = rec
                    else:
                        pass
                else:
                    i += 1
            ordered_records = []
            for busco in shared_busco_ids_list:
                if busco in records_dict:
                    ordered_records.append(records_dict[busco])
            if ordered_records:
                SeqIO.write(ordered_records, out_faa, "fasta")
                busco_id_fasta_list.append(out_faa)
                print(f"PASS:\t[{SPECIES_ID}] Wrote {len(ordered_records)} BUSCO-matching reads to {out_faa}")
            else:
                print(f"[{SPECIES_ID}] No BUSCO matches found in {miniprot_output_gff}")
        # return busco_id_fasta_list
        
        
        print("")
        busco_seq_list = []
        for busco_id_fasta in busco_id_fasta_list:
            busco_id_fasta, fasta_type = fasta_ext_fixer(busco_id_fasta)
            base_name = os.path.basename(busco_id_fasta)
            if "_" in base_name:
                species_id_list = base_name.split("_")
                SPECIES_ID = "_".join(species_id_list[:2])
            else:
                print(f"ERROR:\tUnable to parse input fasta for SPECIES_ID: {busco_id_fasta},\n"
                      "Please rename it (e.g. 'Ps_cubensis-B+_...')")
                continue
            busco_seq_fasta = busco_id_fasta.replace(".fasta","_SEQ.fasta")
            with open(busco_id_fasta, "r") as infile, open(busco_seq_fasta, "w") as outfile:
                outfile.write("".join(line.strip() for line in infile if not line.startswith(">")))
            
            with open(busco_seq_fasta, "r") as infile:
                content = infile.read()
            
            with open(busco_seq_fasta, "w") as outfile:
                outfile.write(f">{SPECIES_ID}-Shared_BUSCO_SEQ\n" + content)
            busco_seq_list.append(busco_seq_fasta)
        print(f"PASS:\tBUSCO Sequence List generated: {busco_seq_list}")
        # return busco_seq_list
        
        
        
        
        combined_fasta_path = f"Combined_{COMPLEASM_DB}_BUSCO_Sequences.fasta"
        with open(combined_fasta_path, "w") as combined_fasta:
            for busco_seq_fasta in busco_seq_list:
                with open(busco_seq_fasta, "r") as infile:
                    combined_fasta.write(infile.read() + "\n")
        print(f"PASS:\tCombined FASTA file created: {combined_fasta_path}")
        # return combined_fasta_path
        
        
        
        combined_records = []
        for fasta_file in busco_seq_list:
            try:
                records = list(SeqIO.parse(fasta_file, "fasta"))
                combined_records.extend(records)
            except Exception as e:
                print(f"Error parsing {fasta_file}: {e}")
        
        combined_fasta_path = f"Combined_{COMPLEASM_DB}_BUSCO_Sequences.fasta"
        SeqIO.write(combined_records, combined_fasta_path, "fasta")
        combined_fasta_path, fasta_type = fasta_ext_fixer(combined_fasta_path)
        base_name = os.path.basename(combined_fasta_path)
        if "_" in base_name:
            species_id_list = base_name.split("_")
            SPECIES_ID = "_".join(species_id_list[:2])
        else:
            print(f"ERROR:\tUnable to parse input fasta for SPECIES_ID: {combined_fasta_path},\n"
                  "Please rename it (e.g. 'Ps_cubensis-B+_...')")
        with open(combined_fasta_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                print(record.id)
        print(f"PASS:\tCombined FASTA file created: {combined_fasta_path}")
        # return combined_fasta_path
        
        
        
        output_msa = f"{os.path.splitext(combined_fasta_path)[0]}_alignment.msa"
        mafft_cmd = f"mafft --thread {CPU_THREADS} --auto {combined_fasta_path} > {output_msa}"
        _ = run_subprocess_cmd(mafft_cmd, shell_check=True)
        print(f"PASS:\tMAFFT alignment completed. Alignment saved to: {output_msa}")
        # return output_msa
        
        
        
        
        tree_prefix = "/mnt/d/EGAP/Combined_basidiomycota_BUSCO_Sequences_alignment_iqtree"
        treefile = tree_prefix + ".treefile"
        if os.path.exists(treefile):
            print("SKIP:\tIQTree tree modeling and making; already exists: {treefile}")
        else:
            iq_tree = ["iqtree", "-s", output_msa, "-m", "MFP", "-T", int(CPU_THREADS), "--prefix", tree_prefix]
            _ = run_subprocess_cmd(iq_tree, shell_check=False)
        print(f"PASS:\tIQTree tree modeling and making. Treefile saved to: {treefile}")
        # return treefile
        
        support_prefix = "/mnt/d/EGAP/Combined_basidiomycota_BUSCO_Sequences_alignment_bootstrap"
        support_file = support_prefix + ".ufboot"
        if os.path.exists(support_file):
            print("SKIP:\tIQTree tree modeling and supporting; already exists: {support_file}")
        else:
            iq_support = ["iqtree", "-s", output_msa, "-m", "MFP", "-T", int(CPU_THREADS), "--prefix", support_prefix, "-B", int(1000), "--bnni", "--boot-trees"]
            _ = run_subprocess_cmd(iq_support, shell_check=False)
        print(f"PASS:\tIQTree tree modeling and supporting. Support file saved to: {support_file}")
        # return support_file
        
        map_prefix = "/mnt/d/EGAP/Combined_basidiomycota_BUSCO_Sequences_alignment_support"
        supported_tree_file = ""
        if os.path.exists(supported_tree_file):
            print("SKIP:\tIQTree tree support mapping; already exists: {supported_tree_file}")
        else:
            iq_map = ["iqtree", "-t", treefile, "--support", support_file, "-T", int(CPU_THREADS), "--prefix", map_prefix]
            _ = run_subprocess_cmd(iq_map, shell_check=False)
        print(f"PASS:\tIQTree tree support mapping. Supported & Mapped Treefile saved to: {supported_tree_file}")
        # return supported_tree_file
