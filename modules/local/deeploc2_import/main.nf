process DEEPLOC2_IMPORT {
    tag "deeploc2_import"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path import_script
    path isar_results
    path deeploc2_results

    output:
    path "deeploc2_import", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${import_script} \\
        --isar-dir ${isar_results} \\
        --deeploc2-results ${deeploc2_results} \\
        --outdir deeploc2_import \\
        ${args}
    """

    stub:
    """
    mkdir -p deeploc2_import
    touch deeploc2_import/deeploc2_import_notes.txt
    touch deeploc2_import/deeploc2_summary.csv
    """
}
