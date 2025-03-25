#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Folder input
params.base_folder = params.base_folder ?: "/mnt/d/TESTING_SPACE/EGEP_Test_Data"

// Processing input
params.primary_genes = params.primary_genes ?: null
params.secondary_genes = params.secondary_genes ?: null
params.data_type = params.data_type ?: "aa"
params.genes_csv = params.genes_csv ?: null
params.in_genus = params.in_genus ?: "Genus_Not_Provided"
params.out_label = params.out_label ?: "Unspecified-Genes"

// Compute resources
params.ram_gb = params.ram_gb ?: 40
params.cpu_threads = params.cpu_threads ?: 12

// Set other inputs
params.output_dir = params.output_dir ?: "${params.base_folder}/EGEP_exonerate_clinker"
params.exonerate_pid = params.exonerate_pid ?: 0.9
params.clinker_pid = params.clinker_pid ?: 0.9
params.organism_kingdom = params.organism_kingdom ?: "Funga"
params.max_intron = params.max_intron ?: null

// Process to fetch UniProt sequences if genes_csv is missing
process fetchUniProtSequences {
    tag "Fetching UniProt sequences"
    publishDir "${params.output_dir}/uniprot_sequences", mode: 'copy'
    
    input:
    val primary_genes
    val secondary_genes
    val in_genus
    val data_type
    val out_label
    
    output:
    path "${in_genus}_${out_label}_${data_type}.fasta"
    
    when:
    params.genes_csv == null && params.primary_genes != null && params.secondary_genes != null
    
    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/fetch_uniprot_sequences.py \
        "${primary_genes}" \
        "${secondary_genes}" \
        "${in_genus}" \
        "${data_type}" \
        "${out_label}"
    """
}

// Process to find FASTA files
process findAssemblyList {
    tag "Finding assembly list"
    publishDir "${params.output_dir}/assembly_lists", mode: 'copy'
    
    input:
    val base_folder
    val in_genus
    val out_label  // Also used to exclude output directories in find_assemblies.py
    
    output:
    stdout
    path "in-${in_genus}_out-${out_label}_assemblies.txt"
    
    script:
    """
    python3 ${workflow.projectDir}/bin/find_assemblies.py "${base_folder}" "${in_genus}" "${out_label}" > in-${in_genus}_out-${out_label}_assemblies.txt
    cat in-${in_genus}_out-${out_label}_assemblies.txt  # Just the paths
    """
}

// Process to BLAST input genes against assemblies to find the best contig
process blastAndFindBestContig {
    tag "BLASTing assembly: ${assembly_file.simpleName}"
    publishDir "${params.output_dir}/blast_results", mode: 'copy'

    input:
    path assembly_file
    path genes_fasta
    val primary_genes
    val data_type

    output:
    tuple path(assembly_file), path("${assembly_file.simpleName}_blast_hits.json")
    path "${assembly_file.simpleName}_blast_debug.log"

    script:
    def blast_cmd = (data_type == "aa") ? "tblastn" : "blastn"
    def genes_list = primary_genes.split(',').join(' ')
    """
    makeblastdb -in ${assembly_file} -dbtype nucl -out ${assembly_file.baseName}
    ${blast_cmd} -query ${genes_fasta} -db ${assembly_file.baseName} -out ${assembly_file.simpleName}_blast.out -outfmt 6
    python3 ${workflow.projectDir}/bin/exonerate_clinker/find_best_contig.py \
        ${assembly_file.simpleName}_blast.out \
        "${genes_list}" \
        ${assembly_file.simpleName}_blast_hits.json \
        > ${assembly_file.simpleName}_blast_debug.log 2>&1
    """
}

// Process to trim the contig to a consensus region
process extractTrimmedRegion {
    tag "Extracting trimmed region from: ${assembly_file.simpleName}"
    publishDir "${params.output_dir}/trimmed_fastas", mode: 'copy'

    input:
    tuple path(assembly_file), path(blast_hits_json)
    val out_label

    output:
    path "${assembly_file.simpleName}_${out_label}_trimmed.fasta", optional: true

    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/extract_trimmed_region.py \
        ${assembly_file} \
        ${blast_hits_json} \
        ${assembly_file.simpleName}_${out_label}_trimmed.fasta || \
        echo "Warning: No BLAST hits aligned for ${assembly_file.simpleName}" >&2
    """
}

// Process to extract the best contig FASTA from the assembly
process extractBestContigFasta {
    tag "Extracting contig from: ${assembly_file.simpleName}"
    publishDir "${params.output_dir}/contig_fastas", mode: 'copy'

    input:
    tuple path(assembly_file), path(best_contig_dict)
    val out_label

    output:
    path "${assembly_file.simpleName}_${out_label}_best_contig.fasta", optional: true

    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/extract_best_contig.py \
        ${assembly_file} \
        ${best_contig_dict} \
        ${assembly_file.simpleName}_${out_label}_best_contig.fasta || \
        echo "Skipping ${assembly_file.simpleName}: No valid contig ID found in ${best_contig_dict}"
    """
}

// Process to run Exonerate analysis on the best contig FASTA (Step 4)
process exonerateContigAnalysis {
    tag "Exonerating contig: ${contig_fasta.simpleName}"
    publishDir "${params.output_dir}/exonerate_results/contig", mode: 'copy'

    input:
    path contig_fasta
    path genes_fasta
    val exonerate_pid
    val organism_kingdom
    val max_intron
    val data_type

    output:
    path "${contig_fasta.simpleName}.gtf"

    script:
    def model = (data_type == "nt") ? "cdna2genome" : "protein2genome:bestfit"
    def maxintron = max_intron ?: (organism_kingdom == "Funga" ? 200 : (organism_kingdom in ["Flora", "Fauna"] ? 2000 : 200))
    """
    if [ -f ${contig_fasta.simpleName}.gtf ] && [ \$(stat -c%s ${contig_fasta.simpleName}.gtf) -gt 1000 ]; then
        echo "Skipping Exonerate: ${contig_fasta.simpleName}.gtf already exists and is >1KB"
        exit 0
    fi

    exonerate --model ${model} \
              --percent ${exonerate_pid} \
              --showalignment yes \
              --showtargetgff yes \
              --showvulgar no \
              -E \
              --query ${genes_fasta} \
              --maxintron ${maxintron} \
              --target ${contig_fasta} \
              > ${contig_fasta.simpleName}.gtf
    """
}

// Process to generate product_dict from the genes FASTA file
process generateProductDict {
    tag "Generating product_dict from: ${genes_fasta.simpleName}"
    publishDir "${params.output_dir}/product_dict", mode: 'copy'

    input:
    path genes_fasta

    output:
    path "${genes_fasta.simpleName}_product_dict.json"

    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/generate_product_dict.py \
        ${genes_fasta} \
        ${genes_fasta.simpleName}_product_dict.json
    """
}

// Process to convert Exonerate GTF to GFF3 format for contig (Step 5)
process convertContigGtfToGff3 {
    tag "Converting contig GTF to GFF3: ${gtf_file.simpleName}"
    publishDir "${params.output_dir}/gff3_results/contig", mode: 'copy'

    input:
    path gtf_file
    path product_dict_json

    output:
    path "${gtf_file.simpleName}.gff3"

    script:
    """
    if [ -f ${gtf_file.simpleName}.gff3 ]; then
        echo "Skipping conversion: ${gtf_file.simpleName}.gff3 already exists"
        cp ${gtf_file.simpleName}.gff3 ${gtf_file.simpleName}.gff3
    else
        python3 ${workflow.projectDir}/bin/exonerate_clinker/convert_gtf_to_gff3.py \
            ${gtf_file} \
            ${gtf_file.simpleName}.gff3 \
            ${product_dict_json}
    fi
    """
}

// Process to extract cluster subsequence from contig FASTA
process extractCluster {
    tag "Extracting cluster from: ${contig_fasta.simpleName}"
    publishDir "${params.output_dir}/cluster_fastas", mode: 'copy'

    input:
    tuple path(gff3_file), path(contig_fasta)

    output:
    path "${contig_fasta.simpleName}_cluster.fasta"

    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/extract_cluster.py \
        ${gff3_file} \
        ${contig_fasta} \
        ${contig_fasta.simpleName}_cluster.fasta
    """
}

// Process to run Exonerate analysis on the cluster FASTA (Step 7)
process exonerateClusterAnalysis {
    tag "Exonerating cluster: ${contig_fasta.simpleName}"
    publishDir "${params.output_dir}/exonerate_results/cluster", mode: 'copy'

    input:
    path contig_fasta
    path genes_fasta
    val exonerate_pid
    val organism_kingdom
    val max_intron
    val data_type

    output:
    path "${contig_fasta.simpleName}.gtf"

    script:
    def model = (data_type == "nt") ? "cdna2genome" : "protein2genome:bestfit"
    def maxintron = max_intron ?: (organism_kingdom == "Funga" ? 200 : (organism_kingdom in ["Flora", "Fauna"] ? 2000 : 200))
    """
    if [ -f ${contig_fasta.simpleName}.gtf ] && [ \$(stat -c%s ${contig_fasta.simpleName}.gtf) -gt 1000 ]; then
        echo "Skipping Exonerate: ${contig_fasta.simpleName}.gtf already exists and is >1KB"
        exit 0
    fi

    exonerate --model ${model} \
              --percent ${exonerate_pid} \
              --showalignment yes \
              --showtargetgff yes \
              --showvulgar no \
              -E \
              --query ${genes_fasta} \
              --maxintron ${maxintron} \
              --target ${contig_fasta} \
              > ${contig_fasta.simpleName}.gtf
    """
}

// Process to convert Exonerate GTF to GFF3 format for cluster (Step 8)
process convertClusterGtfToGff3 {
    tag "Converting cluster GTF to GFF3: ${gtf_file.simpleName}"
    publishDir "${params.output_dir}/gff3_results/cluster", mode: 'copy'

    input:
    path gtf_file
    path product_dict_json

    output:
    path "${gtf_file.simpleName}.gff3"

    script:
    """
    if [ -f ${gtf_file.simpleName}.gff3 ]; then
        echo "Skipping conversion: ${gtf_file.simpleName}.gff3 already exists"
        cp ${gtf_file.simpleName}.gff3 ${gtf_file.simpleName}.gff3
    else
        python3 ${workflow.projectDir}/bin/exonerate_clinker/convert_gtf_to_gff3.py \
            ${gtf_file} \
            ${gtf_file.simpleName}.gff3 \
            ${product_dict_json}
    fi
    """
}

// Process to extract amino acid sequences from GTF
process extractAASequences {
    tag "Extracting AA sequences from: ${gtf_file.simpleName}"
    publishDir "${params.output_dir}/aa_sequences", mode: 'copy'

    input:
    path gtf_file

    output:
    path "${gtf_file.simpleName}_aa.fasta"

    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/extract_aa_sequences.py \
        ${gtf_file} \
        ${gtf_file.simpleName}_aa.fasta
    """
}

// Process to generateClinkerPlot process
process generateClinkerPlot {
    tag "Generating clinker plot"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path genes_fasta
    path fasta_gff3_pairs  // List of [fasta, gff3] pairs flattened into a single list
    val clinker_pid
    val out_label

    output:
    path "${out_label}_clinker_plot.html", optional: true  // Mark as optional to avoid pipeline failure
    path "${out_label}_clinker_matrix.csv", optional: true

    script:
    """
    # Log inputs for debugging
    echo "DEBUG: Genes FASTA: ${genes_fasta}"
    echo "DEBUG: FASTA-GFF3 pairs: ${fasta_gff3_pairs}"

    # Filter GFF3 files explicitly
    GFF3_FILES=()
    for file in ${fasta_gff3_pairs}; do
        if [[ "\$file" == *.gff3 ]]; then
            GFF3_FILES+=("\$file")
        fi
    done

    # Check if GFF3 files exist
    if [ \${#GFF3_FILES[@]} -eq 0 ]; then
        echo "ERROR: No GFF3 files found in input. Available files: ${fasta_gff3_pairs}"
        exit 1
    else
        echo "DEBUG: Found GFF3 files: \${GFF3_FILES[@]}"
    fi

    # Run clinker with verbose output
    clinker \${GFF3_FILES[@]} \
        -i ${clinker_pid} \
        -p ${out_label}_clinker_plot.html \
        > clinker_log.txt 2>&1

    # Check if HTML file was created
    if [ -f "${out_label}_clinker_plot.html" ]; then
        echo "SUCCESS: Generated clinker plot: ${out_label}_clinker_plot.html"
    else
        echo "ERROR: Clinker failed to generate ${out_label}_clinker_plot.html. See clinker_log.txt"
        cat clinker_log.txt
        exit 1
    fi

    # Move matrix file if it exists
    if [ -f "${out_label}_clinker_plot_matrix.csv" ]; then
        mv ${out_label}_clinker_plot_matrix.csv ${out_label}_clinker_matrix.csv
        echo "SUCCESS: Matrix file moved to ${out_label}_clinker_matrix.csv"
    else
        echo "WARNING: No matrix file generated by clinker"
    fi
    """
}

// Workflow
workflow {
    // Step 1: Set up genes_csv_ch and data_type
    if (params.genes_csv) {
        genes_csv_ch = Channel.fromPath(params.genes_csv, checkIfExists: true)
        data_type = params.genes_csv =~ /_(nt|aa)\.fasta$/ ? (params.genes_csv =~ /_(nt|aa)\.fasta$/)[0][1] : 'aa'
    } else if (params.primary_genes && params.secondary_genes && params.data_type) {
        genes_csv_ch = fetchUniProtSequences(params.primary_genes, params.secondary_genes, params.in_genus, params.data_type, params.out_label)
        data_type = params.data_type
    } else {
        error "Missing required inputs: provide genes_csv or primary_genes, secondary_genes, and data_type"
    }

    // Step 2: Generate product dict
    product_dict_ch = generateProductDict(genes_csv_ch)

    // Step 3: Find assembly files
    assembly_list_ch = findAssemblyList(params.base_folder, params.in_genus, params.out_label)
    assembly_files_ch = assembly_list_ch[0]
        .splitText()
        .map { it.trim() }
        .filter { it.endsWith('.fasta') && file(it).exists() && !it.contains(params.out_label) }

    // Step 4: BLAST each assembly file and find hit locations
    blast_results_ch = blastAndFindBestContig(
        assembly_files_ch,
        genes_csv_ch.first(),
        params.primary_genes,
        data_type
    )
    blast_hits_ch = blast_results_ch[0]  // Tuple of [assembly_file, blast_hits_json]
    debug_logs_ch = blast_results_ch[1]  // Debug log files
    blast_hits_ch.view { "Blast Hits JSON: ${it[1]}" }

    // Step 5: Extract trimmed region FASTA
    trimmed_fasta_ch = extractTrimmedRegion(blast_hits_ch, params.out_label)
    trimmed_fasta_ch.view { "Trimmed FASTA: $it" }

    // Step 6: Run Exonerate on trimmed region - collect all outputs
    gtf_ch = exonerateClusterAnalysis(
        trimmed_fasta_ch,
        genes_csv_ch.first(),
        params.exonerate_pid,
        params.organism_kingdom,
        params.max_intron,
        data_type
    ).collect()  // Collect all GTF files into a single list
    gtf_ch.view { "All GTFs collected: $it" }

    // Step 7: Convert GTF to GFF3 - process all GTFs
    gff3_ch = convertClusterGtfToGff3(gtf_ch.flatten(), product_dict_ch.first())
    gff3_ch.view { "GFF3: $it" }

    // Step 8: Extract amino acid sequences - process all GTFs
    aa_sequences_ch = extractAASequences(gtf_ch.flatten())
    aa_sequences_ch.view { "AA Sequences: $it" }

    // Step 9: Pair trimmed FASTA with GFF3 files
    fasta_gff3_ch = trimmed_fasta_ch.join(gff3_ch)
    fasta_gff3_ch.view { "FASTA-GFF3 Pair: $it" }

    // Step 10: Generate clinker plot
    clinker_plot_ch = generateClinkerPlot(genes_csv_ch, fasta_gff3_ch.collect(), params.clinker_pid, params.out_label)
}

// Step 11: Dendrogram: Uses a MAFFT alignment of the trimmed region fasta to then build a ML IQTree with Bootstrap output