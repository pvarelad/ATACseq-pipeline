#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    label 'process_high'
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    tuple val(genome_name), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    bowtie2 -p ${task.cpus} -x ${index}/${genome_name} -U ${reads} | \
        samtools view -bS - > ${sample_id}.bam
    """

    stub:
    """
    touch ${sample_id}.bam
    """
}