process ISAR_CONTRAST_SUMMARY {
    tag "isar_contrast_summary"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path summary_script
    path isar_results

    output:
    path "isar_contrast_summary", emit: results
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), emit: versions_ggplot2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? " \\\n        ${task.ext.args}" : ''
    """
    Rscript ${summary_script} \\
        --isar-dir ${isar_results} \\
        --outdir isar_contrast_summary \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --top-n ${params.isar_contrast_summary_top_n}${args}
    """

    stub:
    """
    mkdir -p isar_contrast_summary
    touch isar_contrast_summary/isar_contrast_summary_notes.txt
    touch isar_contrast_summary/significant_isoform_switches_per_comparison.csv
    touch isar_contrast_summary/isoform_switch_intersections.csv
    touch isar_contrast_summary/isoform_switch_intersection_members.csv
    """
}
