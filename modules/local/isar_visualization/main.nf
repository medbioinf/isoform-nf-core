process ISAR_VISUALIZATION {
    tag "isar_visualization"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path visualization_script
    path isar_results

    output:
    path "isar_visualization", emit: results
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), emit: versions_ggplot2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${visualization_script} \\
        --isar-dir ${isar_results} \\
        --outdir isar_visualization \\
        --top-n ${params.isar_visualization_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        ${args}
    """

    stub:
    """
    mkdir -p isar_visualization
    touch isar_visualization/visualization_notes.txt
    touch isar_visualization/top_isoform_candidates.csv
    touch isar_visualization/top_gene_summary.csv
    """
}
