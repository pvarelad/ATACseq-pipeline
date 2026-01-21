#!/usr/bin/env nextflow

process TAGDIR {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output: 
    tuple val(sample_id), path("tagdir_${sample_id}"), emit: tagdir

    shell: 
    """
    makeTagDirectory tagdir_${sample_id} ${bam_file}
    """

    stub:
    """
    mkdir ${sample_id}_tags
    """
}