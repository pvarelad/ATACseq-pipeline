#!/usr/bin/env nextflow 
nextflow.enable.dsl=2

include {FASTQC} from './modules/fastqc'
include {BOWTIE2_BUILD} from './modules/bowtie2_build'
include {BOWTIE2_ALIGN} from './modules/bowtie2_align'
include {SAMTOOLS_SORT} from './modules/samtools_sort'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {REMOVE_MT} from './modules/remove_mt'
include {BAMCOVERAGE} from './modules/bamcoverage'
include {MULTIQC} from './modules/multiqc'
include {TAGDIR} from './modules/tagdir'
include {FINDPEAKS} from './modules/homer_findpeaks'
include {POS2BED} from './modules/pos2bed'
include {ANNOTATE} from './modules/annotate'
include {MOTIFS} from './modules/motifs'

    workflow {
    
    read_ch = Channel.fromFilePairs('fastq/*.fastq.gz', size: 1)
    

    FASTQC(read_ch)

    BOWTIE2_BUILD(params.genome) 
    BOWTIE2_ALIGN(read_ch, BOWTIE2_BUILD.out.index)

    SAMTOOLS_SORT(BOWTIE2_ALIGN.out)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.sorted)
    
    REMOVE_MT(SAMTOOLS_SORT.out)
    
    multiqc_ch = FASTQC.out.zip.map{ it[1] }.flatten() \
    .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map{ it[1] })
    .collect()
    
    MULTIQC(multiqc_ch)

    TAGDIR(REMOVE_MT.out)

    FINDPEAKS(TAGDIR.out)

    BAMCOVERAGE(REMOVE_MT.out)

    POS2BED(FINDPEAKS.out.peaks)

    peaks_ch = Channel.fromPath('significant_peaks.bed')

    ANNOTATE(peaks_ch, params.gtf,params.genome)
   
    MOTIFS(peaks_ch,params.genome)

    }