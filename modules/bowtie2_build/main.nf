#!/usr/bin/env nextflow

process BOWTIE2_BUILD {
    label 'process_high'
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir params.outdir, mode:'copy'

    input:
    path genome

    output:
    tuple val(genome.baseName), path("bowtie2_index"), emit: index

    script:
    """
    mkdir bowtie2_index
    bowtie2-build --threads ${task.cpus} ${genome} bowtie2_index/${genome.baseName}
    """

    stub:
    """
    mkdir bowtie2_index
    """
}