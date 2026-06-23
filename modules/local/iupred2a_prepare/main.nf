process IUPRED2A_PREPARE {
    tag "iupred2a_prepare"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path prepare_script
    path isar_results

    output:
    path "iupred2a_prepare", emit: results
    path "iupred2a_prepare/*_AA.fasta", optional: true, emit: aa_fasta
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? " \\\n        ${task.ext.args}" : ''
    """
    Rscript ${prepare_script} \\
        --isar-dir ${isar_results} \\
        --outdir iupred2a_prepare \\
        --top-n ${params.iupred2a_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --target-genes "${params.annotated_switch_genes ?: ''}"${args}
    """

    stub:
    """
    mkdir -p iupred2a_prepare
    touch iupred2a_prepare/iupred2a_prepare_notes.txt
    touch iupred2a_prepare/iupred2a_top_candidates.csv
    printf ">stub_isoform\\nMA\\n" > iupred2a_prepare/isoform_iupred2a_candidates_AA.fasta
    """
}
