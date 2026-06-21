process IUPRED2A_IMPORT {
    tag "iupred2a_import"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path import_script
    path isar_results
    path iupred2a_results

    output:
    path "iupred2a_import", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${import_script} \\
        --isar-dir ${isar_results} \\
        --iupred2a-results ${iupred2a_results} \\
        --outdir iupred2a_import \\
        ${args}
    """

    stub:
    """
    mkdir -p iupred2a_import
    touch iupred2a_import/iupred2a_import_notes.txt
    touch iupred2a_import/iupred2a_idr_summary.csv
    """
}
