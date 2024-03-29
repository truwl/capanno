// nf-core DSL2 modules are provided by Truwl to teach internal Nextflow process mechanics. See https://github.com/nf-core/modules for instructions on how to use these in your recipes
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BEDTOOLS_MASKFASTA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    tuple val(meta), path(bed)
    path  fasta

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bedtools \\
        maskfasta \\
        $options.args \\
        -fi $fasta \\
        -bed $bed \\
        -fo ${prefix}.fa
    bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
    """
}
