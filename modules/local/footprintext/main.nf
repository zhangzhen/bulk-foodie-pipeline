process FOOTPRINTEXT {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(cut_counts_file), path(mappable_file), path(DHS_bed_file)
    val scripts_dir

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    footprinting_run_all_4321.sh \\
        $cut_counts_file \\
        $mappable_file \\
        $DHS_bed_file \\
        ${prefix}.bed \\
        ${meta.id} \\
        $scripts_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octave: \$(octave -version | head -n1 | cut -d' ' -f2)
        python: \$(python --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octave: \$(octave -version | head -n1 | cut -d' ' -f2)
        python: \$(python --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
} 