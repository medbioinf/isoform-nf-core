process DEEPLOC2_PREPARE {
    tag "deeploc2_prepare"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path prepare_script
    path isar_results

    output:
    path "deeploc2_prepare", emit: results
    path "deeploc2_prepare/*_AA.fasta", optional: true, emit: aa_fasta
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? " \\\n        ${task.ext.args}" : ''
    """
    Rscript ${prepare_script} \\
        --isar-dir ${isar_results} \\
        --outdir deeploc2_prepare \\
        --top-n ${params.deeploc2_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --target-genes "${params.annotated_switch_genes ?: ''}"${args}
    """

    stub:
    """
    mkdir -p deeploc2_prepare
    touch deeploc2_prepare/deeploc2_prepare_notes.txt
    touch deeploc2_prepare/deeploc2_top_candidates.csv
    printf ">stub_isoform\\nMKTLLAILAVATLALA\\n" > deeploc2_prepare/isoform_deeploc2_candidates_AA.fasta
    """
}
