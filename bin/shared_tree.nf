#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input parameters with CLI override capability
params.no_file = params.no_file ?: "$projectDir/assets/NO_FILE" 
params.base_folder = params.base_folder ?: "/mnt/d/TESTING_SPACE/EGAP_Test_Data"
params.data_type = params.data_type ?: "aa"
params.compleasm_list = params.compleasm_list ?: "basidiomycota,agaricales"
params.in_genus = params.in_genus ?: "Psilocybe"
params.out_genus = params.out_genus ?: "Agrocybe"
params.cpu_threads = params.cpu_threads ?: 10
params.output_dir = params.output_dir ?: "/mnt/d/EGAP"
params.overlap_threshold = params.overlap_threshold ?: 0.3

// Define compleasm_list early with null safety
def compleasm_list = params.compleasm_list ? params.compleasm_list.split(',').toList() : []

// Create channel from compleasm_list
compleasm_ch = Channel.fromList(compleasm_list)

// Define total BUSCO counts (will be populated later)
def busco_totals = [:]
def threshold = params.overlap_threshold

log.info("""
.---.___________________.---.
|[_]|-------------------|[_]|   ;;;;;;;;;;;  ,;;;;;;;,  ;;;;;;;;;;; ;;;;;;;;;,
`---'~~~~~~~~~~~~~~~~~~~`---'  ;%%%%%%%%%% d%%:"":%%%%;%%%%%%%%%%%;%%%%%%%%%%%%
 |||  .--           ,--  |||  |%%%        |%%%        |%%%        |%%%      %%%
 |||  |-            |  _ |||  |@@@        |@@@        |@@@        |@@@      @@@
 |||  `--           `--' |||  |@@@@@@@@@  |@@@        |@@@@@@@@@@ |@@@@@@@@@@@@
 -.|.'-. .-'.   .'-. .'-.|.-  |#########  |###   #####|########## |##########P'
 ||X||||X||||\\ /||||X||||X||  |###@@@@'   |###   #####|###@@@@@'  |###@@@@@P'
 -'|`~-' `-~' v `~-' `~-'|`~  |###        |###   `@###|###        |###
 |||  .--     v     ,--. |||  |888        |?88     88P|888        |888
 |||  |-      V     |__| |||  ?\$\$\$\$\$\$\$\$\$\$\$ ??\$\$\$\$\$\$\$P ?\$\$\$\$\$\$\$\$\$\$\$|\$\$\$
 |||  `--   \\<">/   |    |||  `?\$\$\$\$\$\$\$\$\$\$ `?\$\$\$\$\$\$\$P `\$\$\$\$\$\$\$\$\$\$\$|\$\$\$
.---._____[:::::::]_____.---.   ╔═══════════════════════════════════════════╗
|[_]|------|:::::|------|[_]|   ║    Entheome Genome Extraction Pipeline    ║
`---'~~~~~~|:::::|~~~~~~`---'   ╚═══════════════════════════════════════════╝

              Curated & Maintained by Ian M Bollinger               
                   (ian.bollinger@entheome.org)                     

                          shared_tree.nf

 Directory -> Assembly-List -> Compleasm-BUSCOs -> 
     -> Shared BUSCO Sequences-> MAFFT Alingment-> Boostrap/Concordance IQTree

===============================================================================

    input from   : ${params.base_folder}
    output to    : ${params.output_dir}
    ----------------------------------------------------------------
    run as       : ${workflow.commandLine}
    started at   : ${workflow.start}
    config files : ${workflow.configFiles}
    container    : ${workflow.containerEngine}:${workflow.container}

===============================================================================
""")

// Process tp Find FASTA files
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

// Process to Generate BUSCOs (for BUSCO mode)
process generateNecessaryBuscos {
    tag "Generating/Compiling BUSCOs for ${compleasm_db}"
    publishDir "${params.output_dir}/compleasm_outputs", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"
    cpus { params.cpu_threads }  // Request all available CPUs initially

    input:
    path assembly_list
    val compleasm_db

    output:
    path "${compleasm_db}_busco_files.txt"
    path "temp_${compleasm_db}_compiled_busco_files.json", emit: compiled_list

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/shared_tree/generate_buscos.py "${assembly_list}" "${compleasm_db}" "${task.cpus}"
    """
}

// Process to Initialize compiled list (for BUSCO mode)
process initializeCompiledList {
    tag "Initializing compiled BUSCO list"
    publishDir "${params.output_dir}", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path temp_compiled_files

    output:
    path "final_compiled_busco_files.json"

    script:
    """
    #!/usr/bin/env python3
    import json
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

// Process to Filter assemblies (for BUSCO mode)
process filterAssemblies {
    tag "Filtering assemblies for ${compleasm_db}"
    publishDir "${params.output_dir}/filtered_outputs", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path assembly_list
    val compleasm_db
    path compiled_busco_file

    output:
    path "${compleasm_db}_filtered_assemblies.txt", emit: filtered_assemblies
    path "${compleasm_db}_filtered_tables.txt", emit: filtered_tables
    path "${compleasm_db}_removed_assemblies.csv", emit: removed_assemblies
    path "${compleasm_db}_total_buscos.txt", emit: total_buscos

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/shared_tree/filter_assemblies.py "${assembly_list}" "${compiled_busco_file}" "${compleasm_db}"
    """
}

// Process to Extract initial BUSCO IDs (for BUSCO mode)
process extractInitialBuscos {
    tag "Extracting BUSCO IDs for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path assemblies
    val compleasm_db
    path tables
    path total_buscos_file

    output:
    path "${compleasm_db}_busco_shared.csv", emit: busco_shared
    path "${compleasm_db}_busco_overlap.txt", emit: overlap_string

    script:
    """
    #!/bin/bash
    total=\$(cat "${total_buscos_file}")
    python3 ${PWD}/bin/shared_tree/extract_busco_ids.py "${assemblies}" "${compleasm_db}" "${tables}" "\${total}" "${threshold}"
    """
}

// Process to Extract and concatenate BUSCO sequences (for BUSCO mode)
process extractSequencesFromBusco {
    tag "Extracting BUSCO sequences for ${compleasm_db}"
    publishDir "${params.output_dir}/sequence_alignments/${compleasm_db}_data", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path busco_shared_csv
    val compleasm_db

    output:
    path "${compleasm_db}_busco_shared_sequences.fasta", emit: busco_sequences

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/shared_tree/extract_sequences_from_busco.py "${busco_shared_csv}" "${compleasm_db}"
    """
}

// Process to Infer gene trees for concordance factors (for both modes)
process inferGeneTrees {
    tag "Inferring gene trees for ${db}"
    publishDir "${params.output_dir}/iqtree_outputs/${db}_data/gene_trees", mode: 'copy'
    publishDir "${params.output_dir}/sequence_alignments/${db}_data/trimmed_msas", mode: 'copy', pattern: "locus_files/*_trimmed.msa"
    container "${workflow.projectDir}/bin/entheome.sif"
    cpus { params.cpu_threads }  // Request all available CPUs initially

    input:
    path sequences_file
    val db

    output:
    path "${db}_busco_gene_trees.treefile", emit: gene_trees
    path "locus_files/*_trimmed.msa", emit: trimmed_msas
    path "${db}_gene_trees_debug.log", emit: debug_log

    script:
    """
    #!/bin/bash
    echo "DEBUG: Starting with ${sequences_file}" > ${db}_gene_trees_debug.log
    seq_count=\$(grep -c ">" "${sequences_file}")
    echo "DEBUG: Total sequences: \$seq_count" >> ${db}_gene_trees_debug.log
    patterns=\$(grep -v ">" "${sequences_file}" | tr -d '\n' | tr -d ' ' | sort -u | wc -l)
    echo "DEBUG: Unique patterns (pre-split): \$patterns" >> ${db}_gene_trees_debug.log
    
    python3 ${PWD}/bin/shared_tree/split_fasta_loci.py "${sequences_file}" "locus_files" 2>> ${db}_gene_trees_debug.log
    echo "DEBUG: Split into locus files" >> ${db}_gene_trees_debug.log

    tree_count=0
    for fasta in locus_files/*.fasta; do
        if [ -f "\$fasta" ]; then
            locus_name=\$(basename "\$fasta" .fasta)
            echo "DEBUG: Processing \$fasta" >> ${db}_gene_trees_debug.log
            seq_count=\$(grep -c ">" "\$fasta")
            echo "DEBUG: \$fasta has \$seq_count sequences" >> ${db}_gene_trees_debug.log
            mafft --thread ${task.cpus} --auto "\$fasta" > "locus_files/\${locus_name}_aligned.msa" 2>> ${db}_gene_trees_debug.log
            msa_count=\$(grep -c ">" "locus_files/\${locus_name}_aligned.msa")
            echo "DEBUG: Aligned MSA has \$msa_count sequences" >> ${db}_gene_trees_debug.log
            trimal -in "locus_files/\${locus_name}_aligned.msa" -out "locus_files/\${locus_name}_trimmed.msa" -automated1 2>> ${db}_gene_trees_debug.log
            trimmed_count=\$(grep -c ">" "locus_files/\${locus_name}_trimmed.msa")
            echo "DEBUG: Trimmed MSA has \$trimmed_count sequences" >> ${db}_gene_trees_debug.log
            patterns=\$(grep -v ">" "locus_files/\${locus_name}_trimmed.msa" | tr -d '\n' | tr -d '-' | tr -d ' ' | sort -u | wc -l)
            echo "DEBUG: Unique patterns (post-trim): \$patterns" >> ${db}_gene_trees_debug.log
            if [ \$trimmed_count -ge 2 ] && [ \$patterns -gt 0 ]; then
                iqtree -s "locus_files/\${locus_name}_trimmed.msa" -m MFP -T ${task.cpus} --prefix "locus_files/\${locus_name}_tree" 2>> ${db}_gene_trees_debug.log
                if [ -f "locus_files/\${locus_name}_tree.treefile" ]; then
                    tree_count=\$((tree_count + 1))
                fi
            else
                echo "DEBUG: Skipping IQ-TREE for \${locus_name} - insufficient sequences or variation" >> ${db}_gene_trees_debug.log
            fi
        fi
    done

    if [ \$tree_count -gt 0 ]; then
        cat locus_files/*_tree.treefile > "${db}_busco_gene_trees.treefile" 2>> ${db}_gene_trees_debug.log
        echo "DEBUG: Combined \$tree_count trees" >> ${db}_gene_trees_debug.log
    else
        echo "ERROR: No valid gene trees generated - insufficient sequences or variation" >> ${db}_gene_trees_debug.log
        exit 1
    fi
    """
}

// Process to Concatenate trimmed alignments
process concatenateTrimmedAlignments {
    tag "Concatenating trimmed alignments for ${db}"
    publishDir "${params.output_dir}/sequence_alignments/${db}_data", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path trimmed_msas
    val db

    output:
    path "${db}_concatenated.msa", emit: concatenated_msa

    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/shared_tree/concatenate_alignments.py "${trimmed_msas}" "${db}_concatenated.msa"
    """
}

// Process to Trim alignment with TrimAl
process trimAlignmentWithTrimal {
    tag "Trimming alignment for ${db}"
    publishDir "${params.output_dir}/sequence_alignments/${db}_data", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"
    cache 'deep'

    input:
    path msa_file
    val db

    output:
    path "${db}_alignment_trimmed.msa", emit: trimmed_msa_file
    path "${db}_trim_debug.log", emit: debug_log

    script:
    def outputFile = "${db}_alignment_trimmed.msa"
    """
    #!/bin/bash
    echo "DEBUG: Trimming ${msa_file}" > ${db}_trim_debug.log
    seq_count=\$(grep -c ">" "${msa_file}")
    echo "DEBUG: Input sequences: \$seq_count" >> ${db}_trim_debug.log
    patterns=\$(grep -v ">" "${msa_file}" | tr -d '\n' | tr -d '-' | tr -d ' ' | sort -u | wc -l)
    echo "DEBUG: Unique patterns (pre-trim): \$patterns" >> ${db}_trim_debug.log
    if [ -f "${params.output_dir}/sequence_alignments/${db}_data/${outputFile}" ]; then
        echo "DEBUG: Output ${outputFile} already exists, skipping trimming" >> ${db}_trim_debug.log
        ln -s "${params.output_dir}/sequence_alignments/${db}_data/${outputFile}" "${outputFile}"
    else
        trimal -in "${msa_file}" -out "${outputFile}" -automated1 || {
            echo "WARNING: trimal failed, copying input as output to continue pipeline" >> ${db}_trim_debug.log
            cp "${msa_file}" "${outputFile}"
        }
    fi
    trimmed_count=\$(grep -c ">" "${outputFile}")
    echo "DEBUG: Trimmed sequences: \$trimmed_count" >> ${db}_trim_debug.log
    trimmed_patterns=\$(grep -v ">" "${outputFile}" | tr -d '\n' | tr -d '-' | tr -d ' ' | sort -u | wc -l)
    echo "DEBUG: Unique patterns (post-trim): \$trimmed_patterns" >> ${db}_trim_debug.log
    """
}

// Process to Build IQ-TREE with concordance - Removed dummy tree
process runIqTree {
    tag "Building tree for ${db}"
    publishDir "${params.output_dir}/iqtree_outputs/${db}_data", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path alignment_file
    path gene_trees_file
    val db

    output:
    path "${db}_iqtree.nwk", emit: tree_file

    script:
    """
    #!/bin/bash
    seq_count=\$(grep -c ">" "${alignment_file}")
    echo "DEBUG: Input sequences: \$seq_count" >&2
    if [ "\$seq_count" -ge 3 ]; then
        iqtree -s "${alignment_file}" -m MFP -T ${task.cpus} -B 1000 --bnni --prefix "${db}_iqtree" --boot-trees
        iqtree -t "${db}_iqtree.treefile" --gcf "${gene_trees_file}" -s "${alignment_file}" --scf 100 --prefix "${db}_iqtree_concord"
        cp "${db}_iqtree_concord.cf.tree" "${db}_iqtree.nwk"
    else
        echo "ERROR: Too few sequences (\$seq_count < 3), cannot build a meaningful tree" >&2
        exit 1
    fi
    """
}

// Workflow
workflow {
    no_file = file(params.no_file)

    // Step 1: Find assemblies
    assembly_list_ch = findAssemblyList(params.base_folder, params.in_genus, params.out_genus)

    // BUSCO Path
    // Step 2B-1: Generate BUSCOs
    busco_results = generateNecessaryBuscos(assembly_list_ch, compleasm_ch)

    // Step 2B-2: Compile BUSCO lists
    compiled_busco_ch = initializeCompiledList(busco_results.compiled_list.collect())

    // Step 2B-3: Filter assemblies
    filter_results = filterAssemblies(assembly_list_ch, compleasm_ch, compiled_busco_ch)

    // Step 2B-4: Extract BUSCO IDs
    busco_ids_ch = extractInitialBuscos(filter_results.filtered_assemblies,
                                        compleasm_ch,
                                        filter_results.filtered_tables,
                                        filter_results.total_buscos)

    // Step 2F-1
    
    // Step 2F-2
    
    // Step 2F-3
    
    // Step 2F-4

    // Step 3: Extract and concatenate sequences
    busco_sequences_ch = extractSequencesFromBusco(busco_ids_ch.busco_shared, compleasm_ch)

    // Step 4: Infer gene trees and get trimmed MSAs
    gene_trees_ch = inferGeneTrees(busco_sequences_ch.busco_sequences, compleasm_ch)
    trimmed_msas_ch = gene_trees_ch.trimmed_msas

    // Step 5: Concatenate trimmed alignments
    concatenated_msa_ch = concatenateTrimmedAlignments(trimmed_msas_ch, compleasm_ch)

    // Step 6: Build IQ-TREE with concatenated alignment and gene trees
    iqtree_ch = runIqTree(concatenated_msa_ch.concatenated_msa, gene_trees_ch.gene_trees, compleasm_ch)
}