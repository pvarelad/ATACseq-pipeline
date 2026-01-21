# ATAC-seq Nextflow Pipeline (Single-End Reads)

This repository contains a **Nextflow** pipeline for processing **single-end ATAC-seq** data, from raw FASTQ files to peak calling, annotation, and motif analysis. The workflow performs quality control, alignment, filtering, peak calling with HOMER, and downstream analyses, with results summarized using MultiQC.

---

## Workflow Overview

The pipeline performs the following steps:

1. **Quality control** of raw FASTQ files (FastQC)
2. **Genome index building** (Bowtie2)
3. **Read alignment** to the reference genome (Bowtie2)
4. **Sorting BAM files** (SAMtools)
5. **Alignment statistics** (SAMtools flagstat)
6. **Mitochondrial read removal**
7. **Tag directory creation** (HOMER)
8. **Peak calling** (HOMER findPeaks)
9. **Conversion of peaks to BED format**
10. **Peak annotation** (HOMER annotatePeaks.pl)
11. **Motif discovery** (HOMER findMotifsGenome.pl)
12. **Coverage track generation** (bamCoverage)
13. **Summary reporting** (MultiQC)

---

## Directory Structure

```text
.
├── main.nf
├── nextflow.config
├── modules/
│   ├── fastqc.nf
│   ├── bowtie2_build.nf
│   ├── bowtie2_align.nf
│   ├── samtools_sort.nf
│   ├── samtools_flagstat.nf
│   ├── remove_mt.nf
│   ├── bamcoverage.nf
│   ├── tagdir.nf
│   ├── homer_findpeaks.nf
│   ├── pos2bed.nf
│   ├── annotate.nf
│   └── motifs.nf
├── fastq/
│   └── *.fastq.gz
└── results/
```

---

## Input Requirements

### FASTQ files

* Single-end ATAC-seq reads
* Compressed FASTQ format: `.fastq.gz`
* Located in the `fastq/` directory

Example:

```text
fastq/sample1.fastq.gz
fastq/sample2.fastq.gz
```

The pipeline reads inputs using:

```groovy
Channel.fromFilePairs('fastq/*.fastq.gz', size: 1)
```

---

## Required Parameters

The following parameters must be provided, either via `nextflow.config` or on the command line:

| Parameter  | Description                             |
| ---------- | --------------------------------------- |
| `--genome` | Path to the reference genome FASTA file |
| `--gtf`    | Gene annotation file (GTF format)       |
| `--outdir` | Output directory for results (optional) |

Example `nextflow.config` snippet:

```groovy
params {
    genome = '/path/to/genome.fa'
    gtf    = '/path/to/annotation.gtf'
    outdir = 'results'
}
```

---

## Software Requirements

This pipeline relies on the following tools:

* Nextflow (>= 22.x)
* FastQC
* Bowtie2
* SAMtools
* HOMER
* deepTools (bamCoverage)
* MultiQC

It is recommended to run the pipeline using **containers** (Docker or Singularity) or **Conda environments**, depending on how the individual modules are configured.

---

## Running the Pipeline

From the root directory of the project, run:

```bash
nextflow run main.nf \
  --genome /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  -resume
```

The `-resume` flag allows Nextflow to reuse completed steps if the pipeline is interrupted.

---

## Output Description

Results are organized by analysis step and sample. Typical outputs include:

* **FastQC reports** (`*_fastqc.html`, `*_fastqc.zip`)
* **Aligned BAM files** (sorted, mitochondrial reads removed)
* **Alignment statistics** (`flagstat` files)
* **Coverage tracks** (`.bw` files)
* **Peak files** (HOMER peaks and BED format)
* **Annotated peaks** (gene and genomic feature annotations)
* **Motif analysis results**
* **MultiQC report** summarizing QC and alignment metrics

---

## Notes

* The pipeline is designed specifically for **single-end ATAC-seq** data.
* Mitochondrial reads are removed prior to peak calling.
* Peak calling, annotation, and motif analysis are performed using **HOMER**.
* The pipeline assumes one FASTQ file per sample.

---

## Troubleshooting

* Ensure FASTQ filenames are unique and correctly placed in the `fastq/` directory.
* Verify that the reference genome and GTF files are compatible (same genome build).
* If errors occur, check `.nextflow.log` for detailed diagnostics.

---

## Author

Pipeline developed for ATAC-seq analysis using Nextflow.

