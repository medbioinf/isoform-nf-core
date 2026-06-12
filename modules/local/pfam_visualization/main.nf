process PFAM_VISUALIZATION {
    tag "pfam_visualization"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path visualization_script
    path pfam_import

    output:
    path "pfam_visualization", emit: results
    tuple val("${task.process}"), val('r-ggplot2'), eval("Rscript -e \"cat(as.character(packageVersion('ggplot2')))\""), emit: versions_ggplot2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${visualization_script} \\
        --pfam-dir ${pfam_import} \\
        --outdir pfam_visualization \\
        --top-n ${params.pfam_visualization_top_n} \\
        ${args}
    """

    stub:
    """
    mkdir -p pfam_visualization
    touch pfam_visualization/pfam_visualization_notes.txt
    touch pfam_visualization/pfam_top_domain_change_candidates.csv
    """
}
