process PFAM_IMPORT {
    tag "pfam_import"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path import_script
    path isar_results
    path pfam_results

    output:
    path "pfam_import", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${import_script} \\
        --isar-dir ${isar_results} \\
        --pfam-results ${pfam_results} \\
        --outdir pfam_import \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        ${args}
    """

    stub:
    """
    mkdir -p pfam_import
    touch pfam_import/pfam_import_notes.txt
    touch pfam_import/pfam_domain_summary.csv
    """
}
