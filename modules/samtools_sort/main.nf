#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
    label 'process_medium'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted

    shell:
    """ 
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.stub.sorted.bam
    """
}