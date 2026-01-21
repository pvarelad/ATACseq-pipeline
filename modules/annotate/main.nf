#!/usr/bin/env nextflow

process ANNOTATE {
    label 'process_low'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    path(filtered_peaks)
    path(gtf)
    path(genome)

    output:
    path("annotated_peaks.txt")

    script:
    """
    annotatePeaks.pl \
        ${filtered_peaks} \
        ${genome} \
        -gtf ${gtf} \
        > annotated_peaks.txt
    """

    stub:
    """
    touch annotated_peaks.txt
    """
}