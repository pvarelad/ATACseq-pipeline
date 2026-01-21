#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    label 'process_low'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), emit: flagstat

    shell:
    """
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt
    """

    stub:
    """
    touch ${sample_id}_flagstat.txt
    """
}