// nf-core DSL2 modules are provided by Truwl to teach internal Nextflow process mechanics. See https://github.com/nf-core/modules for instructions on how to use these in your recipes
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bcftools=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
    }

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path  "*.version.txt"        , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    bcftools merge -Oz \\
        --output ${prefix}.vcf.gz \\
        $options.args \\
        *.vcf.gz
    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
