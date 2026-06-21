process SIGNALP_IMPORT {
    tag "signalp_import"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path import_script
    path isar_results
    path signalp_results

    output:
    path "signalp_import", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${import_script} \\
        --isar-dir ${isar_results} \\
        --signalp-results ${signalp_results} \\
        --outdir signalp_import \\
        ${args}
    """

    stub:
    """
    mkdir -p signalp_import
    touch signalp_import/signalp_import_notes.txt
    touch signalp_import/signalp_summary.csv
    """
}
