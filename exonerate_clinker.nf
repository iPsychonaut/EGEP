#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Folder input
params.base_folder = params.base_folder ?: "/mnt/d/TESTING_SPACE/EGEP_Test_Data"

// Processing input
params.primary_genes = params.primary_genes ?: null
params.secondary_genes = params.secondary_genes ?: null
params.data_type = params.data_type ?: "aa"
params.genes_fasta = params.genes_fasta ?: null
params.in_genus = params.in_genus ?: "Genus_Not_Provided"
params.out_dir = params.out_dir ?: "Unspecified-Genes"

// Compute resources
params.ram_gb = params.ram_gb ?: 40
params.cpu_threads = params.cpu_threads ?: 12

// Other inputs
params.output_dir = params.output_dir ?: "${params.base_folder}/EGEP_exonerate_clinker"
params.exonerate_pid = params.exonerate_pid ?: 0.9
params.clinker_pid = params.clinker_pid ?: 0.9
params.organism_kingdom = params.organism_kingdom ?: "Funga"
params.max_intron = params.max_intron ?: nulla
params.max_trimmed_length = params.max_trimmed_length ?: 5000  // Maximum allowed trimmed sequence length (bp)

// Print ASCII art before pipeline runs
log.info """\
      EEEEEEEEEEEEEEE      GGGGGGGGG     EEEEEEEEEEEEEEE  PPPPPPPPPPPPPP
      EEEEEEEEEEEEEEE    GGGGGGGGGGGGG   EEEEEEEEEEEEEEE  PPPPPP    PPPPP
      EEEE              GGGGG      GGGG  EEEE             PPPP        PPPP
      EEEEEEEEE        GGGG              EEEEEEEEE        PPPPPP    PPPPP
      EEEEEEEEE        GGGG     GGGGGG   EEEEEEEEE        PPPPPPPPPPPPPP
      EEEE              GGGGG      GGGG  EEEE             PPPP  
      EEEEEEEEEEEEEEE    GGGGGGGGGGGGG   EEEEEEEEEEEEEEE  PPPP
      EEEEEEEEEEEEEEE      GGGGGGGGG     EEEEEEEEEEEEEEE  PPPP0
                                                                    
               ╔═══════════════════════════════════╗                 
               ║Entheome Genome Extraction Pipeline║                 
               ╚═══════════════════════════════════╝                 
    
              Curated & Maintained by Ian M Bollinger               
                   (ian.bollinger@entheome.org)                     

                        exonerate_clinker.nf
                                  
        Directory -> FASTA-List -> Exonerate-GFFs -> Clinker Plot

         ==========================
         input from   : ${params.base_folder}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """

// Process to fetch UniProt sequences if genes_fasta is missing
process fetchUniProtSequences {
    tag "Fetching UniProt sequences"
    publishDir "${params.out_dir}/uniprot_sequences", mode: 'copy'
    
    input:
    val primary_genes
    val secondary_genes
    val in_genus
    val data_type
    val out_dir
    
    output:
    path "output.fasta", emit: fasta_file
    
    script:
    """
    python3 ${workflow.projectDir}/bin/exonerate_clinker/fetch_uniprot_sequences.py "${primary_genes}" "${secondary_genes}" "${in_genus}" "${data_type}" "${out_dir}"
    mv ${in_genus}_${out_dir}_${data_type}.fasta output.fasta
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

// Process to find FASTA files
process findAssemblyList {
    tag "Finding assembly list"
    publishDir "${params.output_dir}/assembly_lists", mode: 'copy'
    
    input:
    val base_folder
    val in_genus
    val out_dir  // Also used to exclude output directories in find_assemblies.py
    
    output:
    stdout emit: assembly_paths          // Named output for stdout
    path "in-${in_genus}_out-${out_dir}_assemblies.txt", emit: assembly_file  // Named output for file

    script:
    """
    python3 ${workflow.projectDir}/bin/find_assemblies.py "${base_folder}" "${in_genus}" "${out_dir}" > in-${in_genus}_out-${out_dir}_assemblies.txt
    cat in-${in_genus}_out-${out_dir}_assemblies.txt  # Just the paths
    """
}

// Process to BLAST input genes and extract trimmed region in one step
process blastAndFindBestContig {
    tag "BLASTing and trimming: ${assembly_file.simpleName}"
    publishDir "${params.output_dir}/blast_results", mode: 'copy', pattern: "*.json"
    publishDir "${params.output_dir}/clinker_inputs", mode: 'copy', pattern: "*_trimmed.fasta"

    input:
    path assembly_file
    path genes_fasta
    val primary_genes
    val secondary_genes
    val data_type
    val max_trimmed_length
    val out_dir

    output:
    path "${assembly_file.simpleName}_blast_hits.json", emit: blast_hits
    path "${assembly_file.simpleName}_${out_dir}_trimmed.fasta", optional: true, emit: trimmed_fasta
    path "${assembly_file.simpleName}_blast_debug.log", emit: debug_log

    script:
    def blast_cmd = (data_type == "aa") ? "tblastn" : "blastn"
    def genes_list = primary_genes.split(',').join(' ')
    def secondary_list = secondary_genes ? secondary_genes.split(',').join(' ') : ''
    """
    # Create BLAST database
    makeblastdb -in ${assembly_file} -dbtype nucl -out ${assembly_file.baseName} || {
        echo "ERROR: makeblastdb failed" >&2
        exit 1
    }

    # Run BLAST with explicit -outfmt 6
    ${blast_cmd} -query ${genes_fasta} -db ${assembly_file.baseName} -out ${assembly_file.simpleName}_blast.out -outfmt 6 || {
        echo "ERROR: ${blast_cmd} failed" >&2
        exit 1
    }

    # Check BLAST output exists and has content
    if [ ! -s ${assembly_file.simpleName}_blast.out ]; then
        echo "WARNING: BLAST output is empty for ${assembly_file.simpleName}" >&2
        touch ${assembly_file.simpleName}_blast_hits.json  # Empty JSON for downstream compatibility
        echo '{"best_contig": "None"}' > ${assembly_file.simpleName}_blast_hits.json
        exit 0  # Graceful exit, trimmed_fasta is optional
    fi

    # Analyze BLAST and extract trimmed region
    python3 ${workflow.projectDir}/bin/exonerate_clinker/find_and_extract_best_contig.py \
        ${assembly_file} \
        ${assembly_file.simpleName}_blast.out \
        "${genes_list}" \
        "${secondary_list}" \
        ${assembly_file.simpleName}_blast_hits.json \
        ${assembly_file.simpleName}_${out_dir}_trimmed.fasta \
        --buffer 5000 \
        --max_length ${max_trimmed_length} \
        > ${assembly_file.simpleName}_blast_debug.log 2>&1 || {
            echo "ERROR: Python script failed, see debug log" >&2
            cat ${assembly_file.simpleName}_blast_debug.log >&2
            exit 2
        }
    """
}


// Process to run Exonerate analysis on a list of contig FASTAs using exonerate_fasta_list.py
process exonerateContigAnalysis {
    tag "Exonerating contigs with ${cpu_threads} threads"
    publishDir "${params.output_dir}/exonerate_results/contig", mode: 'copy'

    input:
    path contig_fastas  // Collection of FASTA files
    path genes_fasta
    val exonerate_pid
    val organism_kingdom
    val max_intron
    val data_type
    val cpu_threads

    output:
    path "*.gtf"  // Collects all generated GTF files

    script:
    def maxintron_arg = max_intron ?: 'nan'  // Pass 'nan' if not specified
    def fasta_list = contig_fastas.collect { it.name }.join(',')  // Use filenames only
    """
    echo "Running exonerate_fasta_list.py with inputs: ${fasta_list}"
    python3 ${workflow.projectDir}/bin/exonerate_clinker/exonerate_fasta_list.py \
        "${fasta_list}" \
        "${genes_fasta}" \
        "${exonerate_pid}" \
        "${organism_kingdom}" \
        "${maxintron_arg}" \
        "${cpu_threads}"
    ls -l *.gtf || echo "No GTF files found after running script"
    """
}

// Process to convert Exonerate GTF to GFF3 format for contig (Step 5)
process convertContigGtfToGff3 {
    tag "Converting contig GTF to GFF3: ${gtf_file.simpleName}"
    publishDir "${params.output_dir}/clinker_inputs", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

    input:
    path gtf_file
    path product_dict_json

    output:
    path "${gtf_file.simpleName}.gff3", optional: true
    path "${gtf_file.simpleName}_convert.log", emit: debug_log

    script:
    """
    if [ -f ${gtf_file.simpleName}.gff3 ]; then
        echo "Skipping conversion: ${gtf_file.simpleName}.gff3 already exists" > ${gtf_file.simpleName}_convert.log
        cp ${gtf_file.simpleName}.gff3 ${gtf_file.simpleName}.gff3
    else
        python3 ${workflow.projectDir}/bin/exonerate_clinker/convert_gtf_to_gff3.py \
            ${gtf_file} \
            ${gtf_file.simpleName}.gff3 \
            ${product_dict_json} \
            > ${gtf_file.simpleName}_convert.log 2>&1 || echo "Conversion failed; see log" >&2
    fi
    """
}

// Process to extract amino acid sequences from GTF
process extractAASequences {
    tag "Extracting AA sequences from: ${gtf_file.simpleName}"
    publishDir "${params.output_dir}/aa_sequences", mode: 'copy'
    container "${workflow.projectDir}/bin/entheome.sif"

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
    val clinker_pid
    val out_dir
    val gff3_ch

    output:
    path "${out_dir}_clinker_plot.html"  // Only the plot, not optional

    script:
    """
    # Run clinker with staged GFF3 file names only, redirect output for debugging
    clinker -i ${clinker_pid} \
            -p ${out_dir}_clinker_plot.html \
            ${params.output_dir}/clinker_inputs/*.gff3
    """
}


// Workflow
workflow {
    // Step 1: Set up genes_fasta_ch and data_type
    if (params.genes_fasta) {
        if (!file(params.genes_fasta).exists()) {
            error "Provided genes_fasta does not exist: ${params.genes_fasta}"
        }
        genes_fasta_ch = Channel.fromPath(params.genes_fasta)
        data_type = params.genes_fasta =~ /_(nt|aa)\.fasta$/ ? (params.genes_fasta =~ /_(nt|aa)\.fasta$/)[0][1] : 'aa'
        println "Using genes_fasta: ${params.genes_fasta}"
    } else if (params.primary_genes && params.secondary_genes && params.data_type) {
        println "Fetching UniProt sequences"
        genes_fasta_ch = fetchUniProtSequences(
            params.primary_genes,
            params.secondary_genes,
            params.in_genus,
            params.data_type,
            params.out_dir
        ).fasta_file
        data_type = params.data_type
    } else {
        error "Must provide either --genes_fasta or --primary_genes, --secondary_genes, and --data_type"
    }

    // Step 2: Generate product dict
    product_dict_ch = generateProductDict(genes_fasta_ch)

    // Step 3: Find assembly files
    assembly_list_ch = findAssemblyList(params.base_folder, params.in_genus, params.out_dir)
    assembly_files_ch = assembly_list_ch.assembly_paths
        .splitText()
        .map { it.trim() }
        .filter { it.endsWith('.fasta') && file(it).exists() && !it.contains(params.out_dir) }

    // Step 4: BLAST and trim in one go
    blast_results_ch = blastAndFindBestContig(
        assembly_files_ch,
        genes_fasta_ch.first(),
        params.primary_genes,
        params.secondary_genes,
        data_type,
        params.max_trimmed_length,
        params.out_dir
    )
    blast_hits_ch = blast_results_ch.blast_hits
    trimmed_fasta_ch = blast_results_ch.trimmed_fasta
    debug_logs_ch = blast_results_ch.debug_log

    // Step 5: Collect trimmed FASTAs and proceed
    trimmed_fasta_list_ch = trimmed_fasta_ch
        .filter { it && file(it).exists() }
        .collect()

    gtf_ch = exonerateContigAnalysis(
        trimmed_fasta_list_ch,
        genes_fasta_ch.first(),
        params.exonerate_pid,
        params.organism_kingdom,
        params.max_intron,
        data_type,
        params.cpu_threads
    )

    // Step 6: Convert GTF to GFF3 and separate outputs
    convert_results_ch = convertContigGtfToGff3(gtf_ch.flatten(), product_dict_ch.first())
    gff3_ch = convert_results_ch[0]
        .filter { file(it).exists() }
        .filter { file(it).readLines().size() > 1 }
        .collect()
        .ifEmpty { error "No non-empty GFF3 files generated from convertContigGtfToGff3" }

    // Step 10: Generate clinker plot with paired GFF3 and FASTA
    clinker_plot_ch = generateClinkerPlot(
        params.clinker_pid,
        params.out_dir,
        gff3_ch
    )
}
