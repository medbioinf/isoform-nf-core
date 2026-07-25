process ANNOTATED_SWITCH_PLOTS {
    tag "annotated_switch_plots"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path plot_script
    path annotated_rds
    path comparisons

    output:
    path "annotated_switch_plots", emit: results
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genes = params.annotated_switch_genes ?: '-'
    def condition1 = params.annotated_switch_condition1 ?: ''
    def condition2 = params.annotated_switch_condition2 ?: ''
    def plotTopology = (params.annotated_switch_plot_topology instanceof String ? params.annotated_switch_plot_topology.toBoolean() : params.annotated_switch_plot_topology) ? 'true' : 'false'
    def args = task.ext.args ?: ''
    """
    Rscript ${plot_script} \\
        ${annotated_rds} \\
        annotated_switch_plots \\
        ${params.annotated_switch_top_n} \\
        "${genes}" \\
        "${condition1}" \\
        "${condition2}" \\
        ${plotTopology} \\
        ${params.isar_qvalue_cutoff} \\
        ${params.isar_dif_cutoff} \\
        ${comparisons} \\
        ${args}
    """

    stub:
    """
    mkdir -p annotated_switch_plots
    touch annotated_switch_plots/annotated_switch_plot_comparisons.csv
    touch annotated_switch_plots/annotated_switch_plot_notes.txt
    """
}
