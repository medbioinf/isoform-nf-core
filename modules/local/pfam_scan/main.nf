process PFAM_SCAN {
    tag "pfam_scan"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pfam_scan:1.6--hdfd78af_5' :
        'quay.io/biocontainers/pfam_scan:1.6--hdfd78af_5' }"

    input:
    path aa_fasta
    path pfam_db

    output:
    path "pfam_scan.out", emit: results
    path "pfam_scan.log", emit: log
    tuple val("${task.process}"), val('pfam_scan'), val('1.6'), emit: versions_pfam_scan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    db_dir="\$(readlink -f ${pfam_db})"
    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for Pfam scanning." > pfam_scan.log
        touch pfam_scan.out
        exit 0
    fi

    if [ ! -f "\$db_dir/Pfam-A.hmm" ]; then
        echo "Missing Pfam-A.hmm in: \$db_dir" >&2
        exit 1
    fi

    if [ ! -f "\$db_dir/Pfam-A.hmm.h3f" ] || [ ! -f "\$db_dir/Pfam-A.hmm.h3i" ] || [ ! -f "\$db_dir/Pfam-A.hmm.h3m" ] || [ ! -f "\$db_dir/Pfam-A.hmm.h3p" ]; then
        echo "Indexing Pfam-A.hmm with hmmpress" | tee pfam_scan.log
        hmmpress "\$db_dir/Pfam-A.hmm" 2>&1 | tee -a pfam_scan.log
    else
        echo "Using pre-indexed Pfam database: \$db_dir" > pfam_scan.log
    fi

    pfam_scan.pl \\
        -fasta "\$fasta" \\
        -dir "\$db_dir" \\
        -cpu ${task.cpus} \\
        -outfile pfam_scan.out \\
        ${args} \\
        2>&1 | tee -a pfam_scan.log
    """

    stub:
    """
    touch pfam_scan.out
    echo "stub pfam_scan" > pfam_scan.log
    """
}
