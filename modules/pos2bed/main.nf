#!/usr/bin/env nextflow

process POS2BED {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'

    input: 
    tuple val(sample_id), path(homer_txt)

    output: 
    tuple val(sample_id), path("${sample_id}_peaks.bed")


    shell:
    """
    pos2bed.pl ${homer_txt} > ${sample_id}_peaks.bed
    """

    stub:
    """
    touch ${homer_txt.baseName}.bed
    """
}