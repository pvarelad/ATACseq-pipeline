#!/usr/bin/env nextflow

process REMOVE_MT {
    label 'process_high'
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.noMT.bam"), path("${sample_id}.noMT.bam.bai")
    
    script:
    """
    # Remove chrM (mitochondrial DNA)
    samtools idxstats ${bam} | cut -f 1 | grep -v chrM | grep -v MT | \
    xargs samtools view -@ 4 -b ${bam} > ${sample_id}.noMT.bam
    
    samtools index ${sample_id}.noMT.bam
    """
    }