#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

def helpMessage() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_block = params.monochrome_logs ? '' : "\033[3m";
  c_ul = params.monochrome_logs ? '' : "\033[4m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_bul = c_bold + c_ul;
  is_viruses = (params.taxids == 10239) ? " (Viruses)" : ""
  log.info"""
  =${c_dim}=================================================================${c_reset}
  ${c_blue+c_bold}${workflow.manifest.name}${c_reset}   ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}==================================================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bul}Usage:${c_reset}
  Given some barcoded and demultiplexed reads, the typical command for running the pipeline is as follows:
  
    nextflow run ${workflow.manifest.name} \\
      ${c_red}--reads "${params.reads}"${c_reset} \\
      ${c_green}--outdir ${params.outdir}${c_reset} \\
      --refs "refs/*.fasta" \\
      -profile singularity # recommended to run with Singularity

  ${c_yellow+c_bold+c_block}NOTE:${c_yellow} For best results, please ensure you have ${c_bul}Singularity${c_yellow} installed prior to running this workflow.${c_dim}(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)${c_reset}

  Note: 
  The argument supplied to "--reads" must be quoted if using "*" and other 
  characters and symbols that could be shell expanded!

  ${c_bul}Mandatory Options:${c_reset}
    ${c_red}--reads${c_reset}   Input reads directory and pattern (default: ${c_red}"${params.reads}"${c_reset})
    --refs      Reference genomes multiFASTA files (one reference genome per file!) (default: "${params.refs}")
  ${c_bul}Amplicon Sequencing Options:${c_reset}
    --bedfile        BED format file with amplicon sequencing primers info (optional). 
                     Produced as output from PrimalScheme.

  ${c_bul}Consensus Generation Options:${c_reset}
    --low_coverage   Low coverage threshold (default=${params.low_coverage}).
                     Replace consensus sequence positions below this depth
                     threshold with a low coverage character 
                     (see ${c_dim}--low_cov_char${c_reset})
    --no_coverage    No coverage threshold (default=${params.no_coverage}).
                     Replace consensus sequence positions with less than or 
                     equal this depth with a no coverage character 
                     (see ${c_dim}--no_cov_char${c_reset})
    --low_cov_char   Low coverage character (default="${params.low_cov_char}")
    --no_cov_char    No coverage character (default="${params.no_cov_char}")

  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}

  ${c_bul}Other Options:${c_reset}
    ${c_green}--outdir${c_reset}          The output directory where the results will be saved
                      (default: ${c_green}${params.outdir}${c_reset})
    -w/--work-dir     The temporary directory where intermediate data will be 
                      saved (default: ${workflow.workDir})
    -profile          Configuration profile to use. [standard, singularity, 
                      conda, slurm] (default '${workflow.profile}')
    --tracedir        Pipeline run info output directory (default: 
                      ${params.tracedir})

  Note: 
  It is recommended that this workflow be executed with Singularity using the 
  Singularity profile (`-profile singularity`) for maximum reproducibility and
  ease of execution on different platforms.
  """.stripIndent()
}

//=============================================================================
// Help info
//=============================================================================
// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

//=============================================================================
// Check user input params
//=============================================================================
if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

//=============================================================================
// LOG EXECUTION START PARAMS
//=============================================================================
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Ref FASTAs'] = params.refs
if(params.bedfile) {
  summary['Primer Scheme'] = params.bedfile
}
summary['Consensus No Coverage'] = "<=${params.no_coverage}X positions replaced with '${params.no_cov_char}'"
summary['Consensus Low Coverage'] = "<${params.low_coverage}X positions replaced with '${params.low_cov_char}'"
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
summary['Command-Line']   = workflow.commandLine
summary['Nextflow version'] = workflow.nextflow.version
log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
log.info "========================================="

//=============================================================================
// PROCESSES
//=============================================================================

process SNIPPY {
  tag "$sample VS $ref_name"
  publishDir "${params.outdir}/snippy/", pattern: "${sample}-VS-${ref_name}", mode: 'copy'

  input:
    tuple sample,
          path(reads),
          path(ref)
  output:
    tuple sample,
          path(ref),
          path("${sample}-VS-${ref_name}/")

  script:
  ref_name = file(ref).getBaseName()
  """
  snippy \\
    --cpus ${task.cpus} \\
    --outdir ${sample}-VS-${ref_name} \\
    --prefix $sample \\
    --ref $ref \\
    --R1 ${reads[0]} \\
    --R2 ${reads[1]} 
  """
}

process MAP_STATS {
  tag "$sample VS $ref_name"
  publishDir "${params.outdir}/mapping_stats", 
      mode: 'copy', 
      pattern: "*.{tsv,flagstat,idxstats,stats}"

  input:
    tuple val(sample),
          path(ref_fasta),
          path(snippy_outdir)

  output:
  tuple val(sample),
        path(ref_fasta),
        path(snippy_outdir),
        path(depths), emit: depths
  path '*.{flagstat,idxstats,stats}', emit: mqc
  script:
  ref_name = ref_fasta.getBaseName()
  prefix = "${sample}-VS-${ref_name}"
  depths = "${prefix}-depths.tsv"
  flagstat = "${prefix}.flagstat"
  idxstats = "${prefix}.idxstats"
  stats = "${prefix}.stats"
  bam = "${snippy_outdir}/${sample}.bam"
  """
  samtools flagstat $bam > $flagstat
  samtools depth -a -d 0 $bam | perl -ne 'chomp \$_; print "${sample}\t\$_\n"' > $depths
  samtools idxstats $bam > $idxstats
  samtools stats $bam > $stats
  """
}

// Filter for ALT allele variants that have greater depth than REF and
// that have greater than 2X coverage 
process BCF_FILTER {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/variants",
    pattern: "*.filtered.vcf",
    mode: 'copy'
  publishDir "${params.outdir}/variants",
    pattern: "*.raw.vcf",
    saveAs: { "${sample}-VS-${ref_name}.raw.vcf" },
    mode: 'copy'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(snippy_outdir),
          path(depths)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(snippy_outdir),
          path(depths),
          path(filt_vcf)
  script:
  ref_name = ref_fasta.getBaseName()
  vcf = "${snippy_outdir}/${sample}.raw.vcf"
  filt_vcf = "${sample}-VS-${ref_name}.filtered.vcf"
  """
  bcftools filter \\
    -e 'FORMAT/AO < FORMAT/RO' \\
    $vcf \\
    -Ov \\
    -o $filt_vcf
  """
}

process CONSENSUS {
  tag "$sample - $ref_name"
  publishDir "${params.outdir}/consensus", 
    pattern: "*.consensus.fasta",
    mode: 'copy'

  input:
    tuple val(sample),
          path(ref_fasta),
          path(snippy_outdir),
          path(depths),
          path(filt_vcf)
  output:
    tuple val(sample),
          path(ref_fasta),
          path(snippy_outdir),
          path(depths),
          path(consensus)

  script:
  ref_name = ref_fasta.getBaseName()
  consensus = "${sample}-${ref_name}.consensus.fasta"
  """
  vcf_consensus_builder \\
    -v $filt_vcf \\
    -d $depths \\
    -r $ref_fasta \\
    -o $consensus \\
    --low-coverage $params.low_coverage \\
    --no-coverage $params.no_coverage \\
    --low-cov-char $params.low_cov_char \\
    --no-cov-char $params.no_cov_char \\
    --sample-name $sample
  """
}

process COVERAGE_PLOT {
  publishDir "${params.outdir}/plots", 
    pattern: '*.pdf',
    mode: 'copy'

  input:
  tuple val(sample),
        path(ref_fasta),
        path(snippy_outdir),
        path(depths),
        path(filt_vcf)
  
  output:
  path("*.pdf")

  script:
  ref_name = ref_fasta.getBaseName()
  plot_base_filename = "coverage_plot-${sample}-VS-${ref_name}"
  """
  plot_coverage.py -d $depths -o ${plot_base_filename}.pdf
  plot_coverage.py -d $depths -o ${plot_base_filename}-log_scale.pdf --log-scale-y
  plot_coverage.py -d $depths -v $filt_vcf -o ${plot_base_filename}-with_variants.pdf
  plot_coverage.py -d $depths -v $filt_vcf -o ${plot_base_filename}-with_variants-log_scale.pdf --log-scale-y
  plot_coverage.py --no-highlight -d $depths -o ${plot_base_filename}-no_low_cov_highlighting.pdf
  plot_coverage.py --no-highlight --log-scale-y -d $depths -o ${plot_base_filename}-log_scale-no_low_cov_highlighting.pdf
  plot_coverage.py --no-highlight -d $depths -v $filt_vcf -o ${plot_base_filename}-with_variants-no_low_cov_highlighting.pdf
  plot_coverage.py --no-highlight --log-scale-y -d $depths -v $filt_vcf -o ${plot_base_filename}-with_variants-log_scale-no_low_cov_highlighting.pdf
  """
}

process MULTIQC {
  publishDir "${params.outdir}", mode: params.publish_dir_mode,
      saveAs: { filename ->
                "multiqc/$filename"
              }
  input:
  path('samtools/*')

  output:
  path "*multiqc_report.html"
  path "*_data"

  script:
  """
  multiqc .
  """
}

//=============================================================================
// WORKFLOW
//=============================================================================
workflow {

  // Reference genomes input channel
  Channel.fromPath( params.refs )
    .set { ch_refs }

  // Map each sample's reads against each reference genome sequence
  // - Combine each ref seq with each sample's reads
  // - map reads against ref
  if (params.single_end) {
    Channel.fromPath(params.reads)
      .map { [file(it).getBaseName(), it] }
      .set { ch_reads }
  } else {
    ch_reads = Channel.fromFilePairs(
        params.reads,
        checkIfExists: true)
      .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}*.fastq.gz\"' (quotes around reads path required if using `*` and other characters expanded by the shell!)"}
      .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1] ] }
  }
  ch_reads | combine(ch_refs) | SNIPPY | MAP_STATS

  MAP_STATS.out.depths \
    | filter { 
      // Filter for alignments that did have some reads mapping to the ref genome
      depth_linecount = file(it[3]).readLines().size()
      if (depth_linecount == 1) {
        println "No reads from \"${it[0]}\" mapped to reference ${it[1]}"
      }
      depth_linecount > 2
    } \
    | BCF_FILTER \
    | CONSENSUS

  COVERAGE_PLOT(BCF_FILTER.out)
  MULTIQC(MAP_STATS.out.mqc.collect().ifEmpty([]))
}

//=============================================================================
// INTROSPECTION
//=============================================================================
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
