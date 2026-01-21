#!/usr/bin/env nextflow

process FINDPEAKS {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'
 input:
    tuple val(sample_id), path(tagdir)
    
    output:
    tuple val(sample_id), path("${sample_id}_peaks.txt"), emit: peaks
    path "${sample_id}_peaks.txt"
    
    script:
    """
    findPeaks ${tagdir} \
        -style factor \
        -o ${sample_id}_peaks.txt
    """
}