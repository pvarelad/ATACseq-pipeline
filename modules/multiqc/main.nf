#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'
    container 'ghcr.io/bf528/multiqc:latest'
    publishDir 'results/', mode: 'copy'

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc ${qc_files.join(' ')} --outdir . --filename multiqc_report.html -f
    """
    
    stub:
    """
    touch multiqc_report.html
    """
}