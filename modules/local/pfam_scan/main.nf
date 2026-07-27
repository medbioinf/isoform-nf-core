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

    source_db_dir="\$(readlink -f ${pfam_db})"
    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for Pfam scanning." > pfam_scan.log
        touch pfam_scan.out
        exit 0
    fi

    if [ ! -f "\$source_db_dir/Pfam-A.hmm" ]; then
        echo "Missing Pfam-A.hmm in: \$source_db_dir" >&2
        exit 1
    fi

    if [ -f "\$source_db_dir/Pfam-A.hmm.h3f" ] && [ -f "\$source_db_dir/Pfam-A.hmm.h3i" ] && [ -f "\$source_db_dir/Pfam-A.hmm.h3m" ] && [ -f "\$source_db_dir/Pfam-A.hmm.h3p" ]; then
        db_dir="\$source_db_dir"
        echo "Using pre-indexed Pfam database: \$db_dir" > pfam_scan.log
    else
        # Never write hmmpress indexes into the supplied reference directory. Only
        # the large HMM file is copied; ancillary database files remain symlinked.
        db_dir="\$PWD/pfam_db_indexed"
        mkdir -p "\$db_dir"
        cp "\$source_db_dir/Pfam-A.hmm" "\$db_dir/Pfam-A.hmm"
        find "\$source_db_dir" -mindepth 1 -maxdepth 1 ! -name 'Pfam-A.hmm' ! -name 'Pfam-A.hmm.h3f' ! -name 'Pfam-A.hmm.h3i' ! -name 'Pfam-A.hmm.h3m' ! -name 'Pfam-A.hmm.h3p' -exec ln -s {} "\$db_dir/" \\;
        echo "Indexing a task-local copy of Pfam-A.hmm with hmmpress" | tee pfam_scan.log
        hmmpress "\$db_dir/Pfam-A.hmm" 2>&1 | tee -a pfam_scan.log
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
