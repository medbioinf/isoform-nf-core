process DEEPTMHMM_IMPORT {
    tag "deeptmhmm_import"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path import_script
    path isar_results
    path deeptmhmm_results

    output:
    path "deeptmhmm_import", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${import_script} \\
        --isar-dir ${isar_results} \\
        --deeptmhmm-results ${deeptmhmm_results} \\
        --outdir deeptmhmm_import \\
        ${args}
    """

    stub:
    """
    mkdir -p deeptmhmm_import
    touch deeptmhmm_import/deeptmhmm_import_notes.txt
    touch deeptmhmm_import/deeptmhmm_summary.csv
    """
}
