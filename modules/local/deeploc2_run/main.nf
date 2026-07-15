process DEEPLOC2_RUN {
    tag "deeploc2"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/hannharris/deeploc2.1@sha256:4ed94de91083c332cc7405b99f2de5c0a94381cd8bce9e2fe81f7d5b32f97759'
    containerOptions '--entrypoint='

    input:
    path aa_fasta
    path convert_script

    output:
    path "deeploc2/deeploc2_isar.csv", emit: results
    path "deeploc2/results_*.csv", emit: raw
    path "deeploc2/deeploc2.log", emit: log
    tuple val("${task.process}"), val('deeploc2'), val('2.1@sha256:4ed94de91083'), emit: versions_deeploc2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p deeploc2
    export HOME="\$PWD/deeploc2_home"
    export MPLCONFIGDIR="\$PWD/deeploc2_mpl"
    mkdir -p "\$HOME" "\$MPLCONFIGDIR"

    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for DeepLoc2." > deeploc2/deeploc2.log
        printf 'Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome\\n' > deeploc2/deeploc2_isar.csv
        touch deeploc2/results_empty.csv
        exit 0
    fi

    deeploc2 \\
        -f "\$fasta" \\
        -o deeploc2 \\
        -m Fast \\
        -d cpu \\
        ${args} \\
        > deeploc2/deeploc2.log \\
        2>&1

    raw_csv="\$(find deeploc2 -maxdepth 1 -name 'results_*.csv' | sort | head -n 1)"
    if [ -z "\$raw_csv" ]; then
        echo "DeepLoc2 did not write a results_*.csv file." >&2
        exit 1
    fi

    python3 ${convert_script} "\$raw_csv" deeploc2/deeploc2_isar.csv
    """

    stub:
    """
    mkdir -p deeploc2
    export HOME="\$PWD/deeploc2_home"
    export MPLCONFIGDIR="\$PWD/deeploc2_mpl"
    mkdir -p "\$HOME" "\$MPLCONFIGDIR"

    cat > deeploc2/results_stub.csv <<'EOF'
Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome,Peripheral,Transmembrane,Lipid anchor
stub_isoform,Cytoplasm,,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
EOF
    cat > deeploc2/deeploc2_isar.csv <<'EOF'
Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome
stub_isoform,Cytoplasm,,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
EOF
    echo "stub deeploc2" > deeploc2/deeploc2.log
    """
}
