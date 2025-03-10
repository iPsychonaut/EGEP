#!/usr/bin/env nextflow

// Example usage:
// nextflow run main.nf \
// --base_folder "/mnt/d/TESTING_SPACE/EGAP_Test_Data" \
// --compleasm_list "basidiomycota" \
// --in_genus "Ps" \
// --cpu_threads 12 \
// --output_dir "/mnt/d/TESTING_SPACE/EGAP_Test_Data/EGEP_BUSCO" \
// --overlap_threshold 0.3 \
// -with-singularity ${PWD}/bin/egep.sif

// Input parameters with CLI override capability
params.base_folder = params.base_folder ?: "/mnt/d/TESTING_SPACE/EGAP_Test_Data"
params.compleasm_list = params.compleasm_list ?: "basidiomycota"
params.in_genus = params.in_genus ?: "Ps"
params.cpu_threads = params.cpu_threads ?: 12
params.output_dir = params.output_dir ?: "${launchDir}/output"
params.overlap_threshold = params.overlap_threshold ?: 0.3

// Convert compleasm_list from string to list
def compleasm_list = params.compleasm_list.split(',').toList()

// Channels for input data
compleasm_ch = Channel.fromList(compleasm_list)

// Process 1: Filter assemblies based on BUSCO completeness
process filterAssemblies {
    input:
    val compleasm_db from compleasm_ch
    path base_folder from params.base_folder

    output:
    path "assemblies.txt" into assemblies_ch
    path "tables.txt" into tables_ch

    script:
    """
    filter_assemblies.py ${base_folder} ${compleasm_db}
    """
}

// Process 2: Extract initial BUSCO IDs and generate first heatmap
process extractInitialBuscos {
    input:
    path assemblies from assemblies_ch
    path tables from tables_ch
    val compleasm_db from compleasm_ch
    val in_genus from params.in_genus

    output:
    path "busco_id_dict.csv" into busco_dict_ch
    path "heatmap_before.png" into heatmap_before_ch
    path "all_busco_ids.txt" into all_busco_ids_ch

    script:
    """
    extract_busco_ids.py "${assemblies}" "${tables}" "${compleasm_db}" "${in_genus}" "busco_id_dict.csv" "all_busco_ids.txt" "heatmap_data_before.csv" "heatmap_before.png"
    """
}

// Process 3: Filter assemblies by BUSCO overlap and generate second heatmap
process filterBuscosByOverlap {
    input:
    path busco_dict_csv from busco_dict_ch
    path all_busco_ids from all_busco_ids_ch
    path assemblies from assemblies_ch
    val compleasm_db from compleasm_ch
    val overlap_threshold from params.overlap_threshold

    output:
    path "shared_busco_ids.txt" into shared_busco_ch
    path "filtered_assemblies.txt" into filtered_assemblies_ch
    path "heatmap_after.png" into heatmap_after_ch

    script:
    """
    filter_busco_overlaps.py "${busco_dict_csv}" "${all_busco_ids}" "${assemblies}" "${compleasm_db}" ${overlap_threshold} "shared_busco_ids.txt" "filtered_assemblies.txt" "heatmap_data_after.csv" "heatmap_after.png"
    """
}

// Process 4: Extract sequences from GFF files and combine into a single FASTA
process extractSequences {
    input:
    path assemblies from filtered_assemblies_ch
    path shared_busco_ids from shared_busco_ch
    val compleasm_db from compleasm_ch

    output:
    path "combined_${compleasm_db}_busco_sequences.fasta" into fasta_ch

    script:
    """
    extract_sequences.py "${assemblies}" "${shared_busco_ids}" "${compleasm_db}" "combined_${compleasm_db}_busco_sequences.fasta"
    """
}

// Process 5: Align sequences with MAFFT
process alignSequences {
    input:
    path combined_fasta from fasta_ch

    output:
    path "alignment.msa" into alignment_ch

    script:
    """
    mafft --thread ${params.cpu_threads} --auto ${combined_fasta} > alignment.msa
    """
}

// Process 6: Build phylogenetic tree with IQ-Tree
process buildTree {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path alignment from alignment_ch

    output:
    path "*.treefile"

    script:
    """
    iqtree -s ${alignment} -m MFP -T ${params.cpu_threads} --prefix tree
    iqtree -s ${alignment} -m MFP -T ${params.cpu_threads} --prefix bootstrap -B 1000 --bnni --boot-trees
    iqtree -t tree.treefile --support bootstrap.ufboot -T ${params.cpu_threads} --prefix support
    """
}

// Workflow definition
workflow {
    filterAssemblies
        | extractInitialBuscos
        | filterBuscosByOverlap
        | extractSequences
        | alignSequences
        | buildTree
}