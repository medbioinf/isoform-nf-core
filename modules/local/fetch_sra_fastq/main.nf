process FETCH_SRA_FASTQ {
    tag "${meta.sra_run}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37aacd127aa32161d8b38a83efb18df01a8ab1d769a93e88f80342d27801b548/data' :
        'community.wave.seqera.io/library/sra-tools_pigz:4a694d823f6f7fcf' }"

    input:
    tuple val(meta), val(run_accession)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('sratools'), eval("fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+' | head -n 1"), topic: versions, emit: versions_sratools
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"), topic: versions, emit: versions_pigz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${run_accession}"
    """
    prefetch \\
        ${run_accession} \\
        --output-directory sra_cache \\
        --max-size ${params.sra_prefetch_max_size}

    SRA_PATH=\$(find sra_cache -type f \\( -name "${run_accession}.sra" -o -name "${run_accession}.sralite" \\) -print -quit)
    if [ -z "\${SRA_PATH}" ]; then
        echo "Unable to locate prefetched archive for ${run_accession}" >&2
        exit 1
    fi

    fasterq-dump \\
        --threads ${task.cpus} \\
        --skip-technical \\
        --split-files \\
        --outdir fastq \\
        ${args} \\
        "\${SRA_PATH}"

    if [ -f "fastq/${run_accession}_1.fastq" ] && [ -f "fastq/${run_accession}_2.fastq" ]; then
        mv "fastq/${run_accession}_1.fastq" "${prefix}_1.fastq"
        mv "fastq/${run_accession}_2.fastq" "${prefix}_2.fastq"
    elif [ -f "fastq/${run_accession}_1.fastq" ]; then
        mv "fastq/${run_accession}_1.fastq" "${prefix}.fastq"
    elif [ -f "fastq/${run_accession}.fastq" ]; then
        mv "fastq/${run_accession}.fastq" "${prefix}.fastq"
    else
        echo "No FASTQ files produced for ${run_accession}" >&2
        exit 1
    fi

    pigz --no-name --processes ${task.cpus} *.fastq
    """

    stub:
    def prefix = "${meta.id}_${run_accession}"
    """
    echo '' | gzip > ${prefix}.fastq.gz
    """
}
