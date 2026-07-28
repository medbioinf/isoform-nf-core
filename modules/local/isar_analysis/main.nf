process ISAR_ANALYSIS {
    tag "isar"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path analysis_script
    path samplesheet
    path contrasts
    path quant_dirs
    path gtf
    path transcript_fasta

    output:
    path "isar_analysis", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def contrast_arg = contrasts ? "--contrasts ${contrasts}" : ""
    def extra_args = args ? "${args}" : ""
    """
    mkdir -p salmon_quant
    for quant_dir in ${quant_dirs}; do
        ln -s "\$(readlink -f "\$quant_dir")" "salmon_quant/\$(basename "\$quant_dir")"
    done

    Rscript ${analysis_script} \\
        --samplesheet ${samplesheet} \\
        --quant-dir salmon_quant \\
        --gtf ${gtf} \\
        --transcript-fasta ${transcript_fasta} \\
        --outdir isar_analysis \\
        ${contrast_arg} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --top-n ${params.isar_top_n} \\
        ${extra_args}
    """

    stub:
    """
    mkdir -p isar_analysis
    touch isar_analysis/design_matrix.csv
    touch isar_analysis/comparisons.csv
    touch isar_analysis/analysis_notes.txt
    touch isar_analysis/sessionInfo.txt
    """
}
