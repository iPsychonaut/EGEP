#!/usr/bin/env nextflow

// Input parameters
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

// Define total BUSCO counts (will be populated from filterAssemblies output)
def busco_totals = [:]  // Empty map, filled later from total_buscos_ch
def threshold = params.overlap_threshold

// Calculate CPU split
def num_processes = compleasm_list.size()
def split_cpu_threads = (params.cpu_threads / num_processes).intValue() ?: 1

// Process 1: Find FASTA files
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

// Process 2: Generate BUSCOs
process generateNecessaryBuscos {
    tag "Generating/Compiling BUSCOs for ${compleasm_db}"
    publishDir "${params.output_dir}/compleasm_outputs", mode: 'copy'
    input:
    path assembly_list
    val compleasm_db
    output:
    path "${compleasm_db}_busco_files.txt"
    path "temp_${compleasm_db}_compiled_busco_files.json", emit: compiled_list
    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/generate_buscos.py "${assembly_list}" "${compleasm_db}" "${split_cpu_threads}"
    """
}

// Process 3: Initialize compiled list
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

// Process 4: Filter assemblies
process filterAssemblies {
    tag "Filtering assemblies for ${compleasm_db}"
    publishDir "${params.output_dir}/filtered_outputs", mode: 'copy'
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
    python3 ${PWD}/bin/filter_assemblies.py "${assembly_list}" "${compiled_busco_file}" "${compleasm_db}"
    """
}

// Process 5: Extract initial BUSCO IDs
process extractInitialBuscos {
    tag "Extracting BUSCO IDs for ${compleasm_db}"
    publishDir "${params.output_dir}/busco_analysis", mode: 'copy'
    input:
    path assemblies
    val compleasm_db
    path tables
    path total_buscos_file  // Added to get dynamic total_buscos
    output:
    path "${compleasm_db}_busco_shared.csv", emit: busco_shared
    path "${compleasm_db}_busco_overlap.txt", emit: overlap_string
    script:
    """
    #!/bin/bash
    echo "DEBUG: compleasm_db=${compleasm_db}"
    total=\$(cat "${total_buscos_file}")
    echo "DEBUG: total_buscos=\${total}"
    echo "DEBUG: assemblies=${assemblies}"
    echo "DEBUG: tables=${tables}"
    python3 ${PWD}/bin/extract_busco_ids.py "${assemblies}" "${compleasm_db}" "${tables}" "\${total}" "${threshold}"
    """
}

// Process 6: Extract and concatenate BUSCO sequences
process extractBuscoSequences {
    tag "Extracting BUSCO sequences for ${compleasm_db}"
    publishDir "${params.output_dir}/sequence_alignments/${compleasm_db}_data", mode: 'copy'
    input:
    path busco_shared_csv
    val compleasm_db
    output:
    path "${compleasm_db}_busco_shared_sequences.fasta", emit: busco_sequences
    script:
    """
    #!/bin/bash
    python3 ${PWD}/bin/extract_sequences.py "${busco_shared_csv}" "${compleasm_db}"
    """
}

// Process 7: Align sequences with MAFFT
process alignSequencesWithMafft {
    tag "Aligning sequences for ${compleasm_db} (This can take a while; the larger the dataset, the longer it takes.)"
    publishDir "${params.output_dir}/sequence_alignments/${compleasm_db}_data", mode: 'copy'
    input:
    path busco_sequences_file
    val compleasm_db
    output:
    path "${compleasm_db}_busco_alignment.msa", emit: msa_file
    script:
    """
    mafft --thread ${split_cpu_threads} --auto "${busco_sequences_file}" > "${compleasm_db}_busco_alignment.msa"
    """
}

// Process 8: 
process trimAlignmentWithTrimal {
    tag "Trimming alignment for ${compleasm_db}"
    publishDir "${params.output_dir}/sequence_alignments/${compleasm_db}_data", mode: 'copy'
    input:
    path msa_file
    val compleasm_db
    output:
    path "${compleasm_db}_busco_alignment_trimmed.msa", emit: trimmed_msa_file
    script:
    """
    trimal -in "${msa_file}" \
           -out "${compleasm_db}_busco_alignment_trimmed.msa" \
           -gt 0.8 \
           -st 0.1 \
           -cons 60
    """
}

// Process 9: Build IQTree using Model Finder
process runIqTree {
    tag "Building tree for ${compleasm_db} (This can take a while; the larger the dataset, the longer it takes.)"
    publishDir "${params.output_dir}/iqtree_outputs/${compleasm_db}_data", mode: 'copy'
    input:
    path msa_file
    val compleasm_db
    output:
    path "${compleasm_db}_busco_iqtree.*", emit: iqtree_file
    script:
    """
    iqtree -s "${msa_file}" -m MFP -T ${split_cpu_threads} -B 1000 --bnni --prefix "${compleasm_db}_busco_iqtree" --boot-trees
    """
}

// Process 10: Convert IQTree file into Newick File Format
process convertIqTree {
    tag "Converting IQTree output into Nexus format"
    publishDir "${params.output_dir}/iqtree_outputs/${compleasm_db}_data", mode: 'copy'
    input:
    path iqtree_file
    val compleasm_db
    output:
    path "${compleasm_db}_busco_iqtree.nwk", emit: newick_file
    script:
    """
    python3 ${PWD}/iqtree_newick_conversion.py "${iqtree_file}" "${compleasm_db}_busco_iqtree.nwk"
    """
}


// UNTESTED PROCESSES BUT THE DESIRED FUTURE PATHS, CODE MAY NOT FUNCTION OR PERFORM INTENDED TASK AS IT IS ROUGH PLACEHOLDER

// Process _: Extracting Best Model from IQTree output
process extractBestModel {
    tag "Extracting best model for ${compleasm_db}"
    publishDir "${params.output_dir}/iqtree_outputs/${compleasm_db}_data", mode: 'copy'
    input:
    path iqtree_files
    val compleasm_db
    output:
    path "${compleasm_db}_best_model.txt", emit: best_model
    script:
    """
    # Extract best-fit model from IQ-TREE log
    grep "Best-fit model:" "${compleasm_db}_busco_iqtree.log" | awk '{print \$NF}' > temp_model.txt || echo "JTT" > temp_model.txt
    
    # Map to BEAST2-compatible model (simplified example)
    MODEL=\$(cat temp_model.txt | sed 's/+.*//')  # Strip rate variation (e.g., +G, +I)
    case \$MODEL in
        "JTT"|"LG"|"WAG") echo "\$MODEL" > "${compleasm_db}_best_model.txt" ;;
        *) echo "JTT" > "${compleasm_db}_best_model.txt" ;;  # Default to JTT if unsupported
    esac
    """
}

// Process _:
process convertToNexus {
    tag "Converting alignment to Nexus for ${compleasm_db}"
    publishDir "${params.output_dir}/sequence_alignments/${compleasm_db}_data", mode: 'copy'
    input:
    path msa_file
    val compleasm_db
    output:
    path "${compleasm_db}_busco_alignment.nexus", emit: nexus_file
    script:
    """
    # Use a simple Python script or tool to convert FASTA to Nexus
    python3 - <<EOF
from Bio import AlignIO
alignment = AlignIO.read("${msa_file}", "fasta")
with open("${compleasm_db}_busco_alignment.nexus", "w") as f:
    AlignIO.write(alignment, f, "nexus")
EOF
    """
}

// Process _: Run BEAST2
process runBEAST2 {
    tag "Running BEAST2 for ${compleasm_db}"
    publishDir "${params.output_dir}/beast2_outputs/${compleasm_db}_data", mode: 'copy'
    input:
    path nexus_file
    path newick_file
    path best_model_file
    val compleasm_db
    output:
    path "${compleasm_db}_beast.xml"
    path "${compleasm_db}_beast.log"
    path "${compleasm_db}_beast.trees"
    path "${compleasm_db}_beast_mcc.tree"
    script:
    """
    # Read the best model
    BEST_MODEL=\$(cat "${best_model_file}")

    # Read the Newick tree
    TREE=\$(cat "${newick_file}")

    # Generate BEAST2 XML
    cat > "${compleasm_db}_beast.xml" << 'EOF'
<?xml version="1.0" standalone="yes"?>
<beast version="2.7" namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.speciation:beast.evolution.substitutionmodel:beast.math.distributions">
  <data id="alignment" spec="Alignment" fileName="${nexus_file}"/>
  
  <run id="mcmc" spec="MCMC" chainLength="50000000">
    <state>
      <tree id="Tree.t:alignment" spec="beast.evolution.tree.Tree">
        <taxonset id="taxa" spec="TaxonSet">
          <alignment idref="alignment"/>
        </taxonset>
        <trait id="date" spec="beast.evolution.tree.TraitSet" traitname="date-backward" value="0">
          <taxa idref="taxa"/>
        </trait>
        <newick id="startingTree" spec="beast.util.TreeParser" newick="\${TREE}"/>
      </tree>
      <parameter id="clockRate.c:alignment" spec="RealParameter" value="1.0"/>
    </state>
    
    <distribution id="prior" spec="util.CompoundDistribution">
      <distribution id="YuleModel.t:alignment" spec="beast.evolution.speciation.YuleModel">
        <birthRate id="birthRate.t:alignment" spec="RealParameter" value="0.1" lower="0.0" upper="10.0"/>
      </distribution>
      <distribution id="Agrocybe_Psilocybe_calib" spec="beast.math.distributions.Prior">
        <x id="tmrca.Agrocybe" spec="beast.evolution.tree.MRCAPrior" tree="@Tree.t:alignment">
          <taxonset id="Agrocybe" spec="TaxonSet">
            <taxon idref="Agrocybe"/> <!-- Adjust taxon name based on your data -->
          </taxonset>
        </x>
        <distr spec="beast.math.distributions.Normal">
          <parameter name="mean" value="${params.calib_mean}"/>
          <parameter name="sigma" value="${params.calib_sigma}"/>
        </distr>
      </distribution>
    </distribution>
    
    <operator id="clockScaler" spec="ScaleOperator" parameter="@clockRate.c:alignment" scaleFactor="0.5" weight="1.0"/>
    <operator id="treeScaler" spec="TreeScaler" tree="@Tree.t:alignment" scaleFactor="0.5" weight="1.0"/>
    <operator id="uniformOperator" spec="Uniform" tree="@Tree.t:alignment" weight="10.0"/>
    <logger fileName="${compleasm_db}_beast.log" logEvery="5000">
      <log idref="clockRate.c:alignment"/>
      <log idref="YuleModel.t:alignment"/>
    </logger>
    <logger fileName="${compleasm_db}_beast.trees" logEvery="5000">
      <log idref="Tree.t:alignment"/>
    </logger>
  </run>
  
  <siteModel spec="SiteModel" id="siteModel">
    <substModel spec="beast.evolution.substitutionmodel.\${BEST_MODEL}"/>
    <parameter id="proportionInvariant" estimate="true" lower="0.0" upper="1.0" value="0.0"/>
    <frequencies id="empiricalFreqs" spec="Frequencies" data="@alignment"/>
  </siteModel>
  
  <branchRateModel id="StrictClock" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:alignment"/>
</beast>
EOF

    # Run BEAST2
    beast "${compleasm_db}_beast.xml"

    # Generate MCC tree (10% burn-in)
    treeannotator -burnin 10 -heights ca "${compleasm_db}_beast.trees" "${compleasm_db}_beast_mcc.tree"
    """
}

// END PLACEHOLDER SECTION


workflow {
    // Step 1: Find assemblies
    assembly_list_ch = findAssemblyList(params.base_folder, params.in_genus, params.out_genus)
    
    // Step 2: Generate BUSCOs
    busco_results = generateNecessaryBuscos(assembly_list_ch, compleasm_ch)
    
    // Step 3: Collect compiled files
    compiled_files_ch = busco_results.compiled_list.collect()
    compiled_busco_ch = initializeCompiledList(compiled_files_ch)
    
    // Step 4: Filter assemblies
    filter_results = filterAssemblies(assembly_list_ch, compleasm_ch, compiled_busco_ch)
    filtered_assemblies_ch = filter_results.filtered_assemblies
    filtered_tables_ch = filter_results.filtered_tables
    removed_assemblies_ch = filter_results.removed_assemblies
    total_buscos_ch = filter_results.total_buscos
    
    // Display removed assemblies
    removed_assemblies_ch
        .view { file -> 
            def content = file.text.trim()
            if (!content) "INITIAL BUSCO OUTLIER REMOVAL REPORT: No assemblies removed."
            else "INITIAL BUSCO OUTLIER REMOVAL REPORT:\n${content}"
        }
    
    // Step 5: Extract BUSCO IDs and overlap
    // Pair inputs explicitly without tuples
    assemblies_with_db = compleasm_ch
        .join(filtered_assemblies_ch.map { file -> [file.name.split('_')[0], file] })
    tables_with_db = compleasm_ch
        .join(filtered_tables_ch.map { file -> [file.name.split('_')[0], file] })
    totals_with_db = compleasm_ch
        .join(total_buscos_ch.map { file -> [file.name.split('_')[0], file] })
    
    busco_results = extractInitialBuscos(
        assemblies_with_db.map { db, file -> file },
        assemblies_with_db.map { db, file -> db },
        tables_with_db.map { db, file -> file },
        totals_with_db.map { db, file -> file }
    )
    busco_shared_ch = busco_results.busco_shared
    overlap_string_ch = busco_results.overlap_string
    
    // Print overlap
    overlap_string_ch.view { file -> file.text.trim() }
    
    // Step 6: Extract sequences
    csv_with_db = compleasm_ch
        .join(busco_shared_ch.map { file -> [file.name.split('_')[0], file] })
    busco_sequences_results = extractBuscoSequences(
        csv_with_db.map { db, file -> file },
        csv_with_db.map { db, file -> db }
    )
    compiled_busco_sequences = busco_sequences_results.busco_sequences
       
    // Step 7: MAFFT alignment of compiled shared buscos
    msa_results = alignSequencesWithMafft(compiled_busco_sequences, compleasm_ch)
    
    // Step 8: Clean the alignment with TrimAl
    trimmed_msa_results = trimAlignmentWithTrimal(msa_results.msa_file, compleasm_ch)
    
    // Step 9: Buld IQtree for the compiled shared busco alignments
    iqtree_results = runIqTree(trimmed_msa_results.trimmed_msa_file, compleasm_ch)
    
    // Step 10: Convert to Newick
    newick_results = convertIqTree(iqtree_results.iqtree_files, compleasm_ch)

    // Step 11: Extract best model
    model_results = extractBestModel(iqtree_results.iqtree_files, compleasm_ch)

    // Step 12: Convert alignment to Nexus
    nexus_results = convertToNexus(trimmed_msa_results.msa_file, compleasm_ch)

    // Step 13: Run BEAST2
    beast_results = runBEAST2(nexus_results.nexus_file, newick_results.newick_file, model_results.best_model, compleasm_ch)

    // Step 14: Print summary
    printRunSummary(
        total_buscos_ch,
        overlap_string_ch,
        iqtree_results.iqtree_files,
        beast_results.collect(),
        compleasm_ch
    ).view()
}
