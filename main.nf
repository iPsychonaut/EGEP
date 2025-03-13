#!/usr/bin/env nextflow

// Example usage:
// nextflow run main.nf \
// --base_folder "/mnt/d/TESTING_SPACE/EGEP_Test_Data" \
// --compleasm_list "basidiomycota,agaricales" \
// --in_genus "Psilocybe" \
// --out_genus "Agrocybe" \
// --cpu_threads 12 \
// --output_dir "/mnt/d/TESTING_SPACE/EGEP_Test_Data/EGEP_SharedBUSCOTree" \
// --overlap_threshold 0.3 \
// --calib_mean 67.6 \
// --calib_sigma 6 \
// -with-singularity ${PWD}/bin/egep.sif

// Input parameters with CLI override capability
params.base_folder = params.base_folder ?: "/mnt/d/TESTING_SPACE/EGAP_Test_Data"
params.compleasm_list = params.compleasm_list ?: "basidiomycota,agaricales"
params.in_genus = params.in_genus ?: "Psilocybe"
params.out_genus = params.out_genus ?: "Agrocybe"
params.cpu_threads = params.cpu_threads ?: 12
params.output_dir = params.output_dir ?: "/mnt/d/EGAP"
params.overlap_threshold = params.overlap_threshold ?: 0.3
params.calib_mean = params.calib_mean ?: 67.6
params.calib_sigma = params.calib_sigma ?: 6

// Convert compleasm_list to a list and create channel
def compleasm_list = params.compleasm_list.split(',').toList()
compleasm_ch = Channel.fromList(compleasm_list)

// Calculate number of parallel processes and split CPU threads
def num_processes = compleasm_list.size()
def split_cpu_threads = (params.cpu_threads / num_processes).intValue() ?: 1

// Process 1: Find FASTA files using the Python script
process findAssemblyList {
    tag "Finding assembly list"
    publishDir "${params.output_dir}/assembly_lists", mode: 'copy'

    input:
    val base_folder
    val in_genus
    val out_genus

    output:
    path "in-${in_genus}_out-${out_genus}_assemblies.txt"

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/find_assemblies.py "${base_folder}" "${in_genus}" "${out_genus}"
    """
}

// Process 2: Generate BUSCOs using the Python script
process generateNecessaryBuscos {
    tag "Generating/Compiling BUSCOs for ${compleasm_db}"
    publishDir "${params.output_dir}/compleasm_outputs", mode: 'copy'

    input:
    path assembly_list
    val compleasm_db

    output:
    path "${compleasm_db}_busco_files.txt"  // Individual output per compleasm_db (optional)
    path "temp_${compleasm_db}_compiled_busco_files.json", emit: compiled_list  // JSON file per compleasm_db

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/generate_buscos.py "${assembly_list}" "${compleasm_db}" "${split_cpu_threads}"
    """
}

// Process 3: Initialize the compiled list file
process initializeCompiledList {
    tag "Initializing compiled BUSCO list"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path temp_compiled_files

    output:
    path "final_compiled_busco_files.json"

    script:
    """
    #!/usr/bin/env python3
    import json
    import sys

    compiled_data = {}
    for file in "${temp_compiled_files}".split():
        db_name = file.split('temp_')[1].split('_compiled')[0]
        with open(file, 'r') as f:
            data = json.load(f)
            compiled_data[db_name] = data

    with open('final_compiled_busco_files.json', 'w') as f:
        json.dump(compiled_data, f, indent=2)
    """
}

// Process 4: Filter assemblies based on BUSCO completeness
process filterAssemblies {
    tag "Filtering assemblies based on BUSCO completeness"
    publishDir "${params.output_dir}/filtered_outputs", mode: 'copy'

    input:
    path assembly_list
    path compiled_busco_file

    output:
    path "filtered_assemblies.txt"
    path "filtered_tables.txt"
    path "removed_assemblies.csv"

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/filter_assemblies.py "${assembly_list}" "${compiled_busco_file}"
    """
}

// Process 5: Extract initial BUSCO IDs
process extractInitialBuscos {
    tag "Extracting initial BUSCO IDs for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis/${compleasm_db}", mode: 'copy'

    input:
    path assemblies
    path tables
    val compleasm_db

    output:
    path "${compleasm_db}_busco_id_dict.csv", emit: busco_dict
    path "${compleasm_db}_all_busco_ids.txt", emit: all_busco_ids
    path "${compleasm_db}_heatmap_data_before.csv", emit: heatmap_df_before
    path assemblies, emit: assemblies  // Pass-through assemblies

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/extract_busco_ids.py "${assemblies}" "${tables}" "${compleasm_db}" "${compleasm_db}_busco_id_dict.csv" "${compleasm_db}_all_busco_ids.txt" "${compleasm_db}_heatmap_data_before.csv"
    """
}

// Process 6: Filter assemblies by BUSCO overlap
process filterBuscosByOverlap {
    tag "Filtering BUSCOs by overlap for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis/${compleasm_db}", mode: 'copy'

    input:
    val compleasm_db
    path busco_dict_csv
    path all_busco_ids 
    path assemblies
    val overlap_threshold

    output:
    path "${compleasm_db}_shared_busco_ids.txt", emit: shared_busco
    path "${compleasm_db}_filtered_assemblies.txt", emit: filtered_assemblies
    path "${compleasm_db}_heatmap_data_after.csv", emit: heatmap_df_after
    path "${compleasm_db}_shared_busco_dict.csv", emit: shared_busco_dict

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/filter_busco_overlaps.py "${busco_dict_csv}" "${all_busco_ids}" "${assemblies}" "${compleasm_db}" "${overlap_threshold}" "${compleasm_db}_shared_busco_ids.txt" "${compleasm_db}_filtered_assemblies.txt" "${compleasm_db}_heatmap_data_after.csv"
    """
}

// Graphics Generation Process: Plot BUSCO ID Overlap Heatmap
process plotPreFiltBuscoOverlap {
    tag "Generating BUSCO ID Overlap Heatmap for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis/${compleasm_db}", mode: 'copy'

    input:
    val compleasm_db
    val stage
    path heatmap_df_file

    output:
    path "${compleasm_db}_heatmap_${stage}.png", emit: heatmap_png

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/plot_heatmap.py "${heatmap_df_file}" "${compleasm_db}_heatmap_${stage}.png" "${compleasm_db} BUSCO Overlap (${stage})"
    """
}
process plotPostFiltBuscoOverlap {
    tag "Generating BUSCO ID Overlap Heatmap for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis/${compleasm_db}", mode: 'copy'

    input:
    val compleasm_db
    val stage
    path heatmap_df_file

    output:
    path "${compleasm_db}_heatmap_${stage}.png", emit: heatmap_png

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/plot_heatmap.py "${heatmap_df_file}" "${compleasm_db}_heatmap_${stage}.png" "${compleasm_db} BUSCO Overlap (${stage})"
    """
}

// Process 7: Extract sequences from GFF files and combine into a single FASTA
process extractSequences {
    input:
    val compleasm_db), path(assemblies) from filtered_assemblies_ch
    tuple val(compleasm_db_shared), path(shared_busco_ids) from shared_busco_ch

    output:
    tuple val(compleasm_db), path("Combined_${compleasm_db}_BUSCO_Sequences.fasta") into fasta_ch

    script:
    """
    extract_sequences.py "${assemblies}" "${shared_busco_ids}" "${compleasm_db}" "Combined_${compleasm_db}_BUSCO_Sequences.fasta"
    """
}

workflow {
    // Invoke Process 1 and capture its output
    assembly_list_ch = findAssemblyList(params.base_folder, params.in_genus, params.out_genus)
            
    // Invoke Process 2 with the output from Process 1
    generateNecessaryBuscos(assembly_list_ch, compleasm_ch)

    // Collect all temp_${compleasm_db}_compiled_busco_files.json and initialize the list
    generateNecessaryBuscos.out.compiled_list
        .collect()
        .set { all_temp_compiled_files }

    // Invoke Process 3 to initialize the compiled list
    compiled_busco_ch = initializeCompiledList(all_temp_compiled_files)

    // Invoke Process 4 to filter assemblies and destructure multiple outputs
    (filtered_assemblies_ch, filtered_tables_ch, removed_assemblies_ch) = filterAssemblies(assembly_list_ch, compiled_busco_ch)   
    removed_assemblies_ch
        .flatMap { file -> 
            def content = file.text.trim()
            if (!content) {
                return "INITIAL BUSCO OUTLIER REMOVAL REPORT: No assemblies were removed (file is empty)."
            }
            def lines = content.split('\n')
            if (lines.size() <= 1) {
                return "INITIAL BUSCO OUTLIER REMOVAL REPORT: ${lines[0]} (no additional details)."
            }
            def header = lines[0]
            def body = lines[1..-1].join('\n')
            return "INITIAL BUSCO OUTLIER REMOVAL REPORT:\n${header}\n${body}"
        }
        .view { it }
        
    // Invoke Process 5: Extract BUSCO IDs for each Compleasm DB
    extractInitialBuscos(filtered_assemblies_ch, filtered_tables_ch, compleasm_ch)

    // Generate "before" heatmap for each compleasm_db
    plotPreFiltBuscoOverlap(
        compleasm_ch,
        "Pre-filter",
        extractInitialBuscos.out.heatmap_df_before
    )

    // Invoke Process 6: Filter Assemblies based on shred BUSCO_IDs for each Compleasm DB
    filterBuscosByOverlap(
        compleasm_ch,
        extractInitialBuscos.out.busco_dict,
        extractInitialBuscos.out.all_busco_ids,
        extractInitialBuscos.out.assemblies,
        params.overlap_threshold
    )

    // Generate "after" heatmap for each compleasm_db
    plotPostFiltBuscoOverlap(
        compleasm_ch,
        "Post-filter",
        filterBuscosByOverlap.out.heatmap_df_after
    )
    
    // Invoke Process 7: Extract sequences from GFF files and combine into a single FASTA
}
