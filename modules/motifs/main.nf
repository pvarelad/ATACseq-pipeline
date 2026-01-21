#!/usr/bin/env nextflow

process MOTIFS {
    label 'process_high'
    container 'ghcr.io/bf528/homer_samtools:latest'
    publishDir params.outdir, mode:'copy'
    
    cpus 16
    memory '64.GB'
    time '12.h'

    input: 
    path(bed)
    path(fasta)

    output:
    path("motifs/*")

    script:
    """
    samtools faidx ${fasta}
    
    # Fix duplicate IDs and take top 5000 peaks
    # Sort by score (column 5), take top peaks, assign unique IDs
    sort -k5,5nr ${bed} | head -5000 | \
    awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, "peak_"NR, \$5, \$6}' > top_peaks.bed
    
    echo "Created top 5000 peaks with unique IDs"
    head top_peaks.bed
    
    findMotifsGenome.pl \
        top_peaks.bed \
        ${fasta} \
        motifs \
        -size 200 \
        -mask \
        -p ${task.cpus}
    
    # Verify completion
    ls -lh motifs/*.html
    """
}