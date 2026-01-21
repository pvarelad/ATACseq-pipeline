#!/usr/bin/env nextflow

process FASTQC {
    label 'process_low'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip

    shell:
    """
    fastqc -t $task.cpus $reads
    """

    stub:
    """
    touch stub_${sample_id}_fastqc.zip
    touch stub_${sample_id}_fastqc.html
    """
}