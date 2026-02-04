# ATAC-seq Nextflow Pipeline (Single-End Reads)

This pipeline performs comprehensive ATAC-seq analysis from raw FASTQ files to differential chromatin accessibility and motif discovery. Built with Nextflow and containerized tools via Singularity, it ensures complete reproducibility across different computing environments. The pipeline is optimized for identifying open chromatin regions and analyzing differential accessibility between experimental conditions.

## Pipeline Workflow

The pipeline consists of the following major steps:

1. **Quality Control** - Assessment of raw sequencing data
2. **Read Alignment** - Alignment to reference genome
3. **Post-alignment Processing** - Sorting, indexing, filtering, and mitochondrial read removal
4. **Peak Calling** - Identification of open chromatin regions
5. **Differential Accessibility Analysis** - Identification of condition-specific accessible regions
6. **Peak Annotation** - Genomic feature assignment
7. **Motif Enrichment** - Discovery of enriched transcription factor binding motifs
8. **Signal Visualization** - Coverage track generation

## Requirements

### Software Dependencies

All tools are provided through containerized environments (`ghcr.io/bf528`) executed via Singularity:

- [Nextflow](https://www.nextflow.io/) (≥21.0)
- [Singularity](https://sylabs.io/singularity/) (≥3.7)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.12.1)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.5.4)
- [SAMtools](http://www.htslib.org/) (v1.21)
- [MultiQC](https://multiqc.info/) (v1.25)
- [HOMER](http://homer.ucsd.edu/homer/) (v5.1)
- [deepTools](https://deeptools.readthedocs.io/) (v3.5.5)
- R (≥4.0)
  - DiffBind

### Reference Data Requirements

- Mouse reference genome: GRCm39
- Gene annotations for GRCm39

### Input Requirements

- Single-end FASTQ files (gzipped or uncompressed)
- Sample metadata with experimental conditions
- Multiple biological replicates per condition (recommended)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/atacseq-pipeline.git
cd atacseq-pipeline

# Install Nextflow (if not already installed)
curl -s https://get.nextflow.io | bash

# Install Singularity (if not already installed)
# See: https://sylabs.io/guides/3.0/user-guide/installation.html
```

## Usage

### Basic Usage

```bash
nextflow run main.nf \
  --input samples.csv \
  --genome GRCm39 \
  --outdir results \
  -profile singularity
```

### Input Format

Create a CSV file (`samples.csv`) with the following format:

```csv
sample,fastq,condition,replicate
sample1,/path/to/sample1.fastq.gz,control,1
sample2,/path/to/sample2.fastq.gz,control,2
sample3,/path/to/sample3.fastq.gz,treatment,1
sample4,/path/to/sample4.fastq.gz,treatment,2
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to input CSV file | Required |
| `--genome` | Reference genome (GRCm39, GRCm38) | `GRCm39` |
| `--outdir` | Output directory | `./results` |
| `--remove_mito` | Remove mitochondrial reads | `true` |
| `--mito_name` | Mitochondrial chromosome name | `chrM` |

### Advanced Options

```bash
nextflow run main.nf \
  --input samples.csv \
  --genome GRCm39 \
  --remove_mito true \
  --outdir results \
  -profile singularity \
  -resume
```

## Pipeline Details

### 1. Quality Control (FastQC)

Raw single-end FASTQ files are assessed for:
- Per-base sequence quality
- Adapter content
- GC content distribution
- Sequence duplication levels
- Overrepresented sequences

All analyses use default parameters.

### 2. Alignment (Bowtie2)

Filtered reads are aligned to the GRCm39 mouse reference genome using Bowtie2 with default parameters optimized for ATAC-seq data.

```bash
bowtie2 -x genome_index -U input.fastq.gz -S output.sam
```

### 3. Post-alignment Processing (SAMtools)

Aligned reads are processed through multiple quality control steps:

#### 3a. Sorting and Indexing
- **Sorting**: Coordinate-based sorting of alignments
- **Indexing**: BAM file indexing for rapid access

#### 3b. Filtering
- Removal of unmapped reads
- Removal of low-quality alignments
- Filtering of non-primary alignments

#### 3c. Mitochondrial Read Removal
Mitochondrial reads are removed to reduce technical noise, as they can represent:
- Non-nuclear DNA fragments
- Potential sources of bias in accessibility analysis
- High-abundance artifacts that obscure true signal

```bash
samtools view -h input.bam | \
  grep -v chrM | \
  samtools view -b > output_filtered.bam
```

#### 3d. Alignment Statistics
Quality metrics are generated using the `flagstat` utility:
- Total reads aligned
- Properly paired reads
- Mapping quality distribution
- Duplicate rates

All SAMtools operations use default parameters.

### 4. Quality Report Aggregation (MultiQC)

FastQC, Bowtie2, and SAMtools metrics are aggregated into a single comprehensive HTML report for easy visualization across all samples.

### 5. Peak Calling (HOMER)

Open chromatin regions are identified using a multi-step process:

#### 5a. Tag Directory Creation
Filtered BAM files are converted into HOMER tag directories:

```bash
makeTagDirectory tag_dir/ input.bam
```

#### 5b. Peak Calling
Accessible chromatin regions are identified using HOMER with default parameters:

```bash
findPeaks tag_dir/ -style factor -o peaks.txt
```

#### 5c. Format Conversion
Peak files are converted to BED format for downstream analyses and visualization.

### 6. Differential Accessibility Analysis (DiffBind)

DiffBind in R is used to identify differentially accessible regions between conditions:

- **Count matrix generation**: Read counts within peaks across all samples
- **Normalization**: Library size and composition normalization
- **Statistical testing**: Identification of significantly differential peaks
- **FDR correction**: Multiple testing correction
- **Visualization**: MA plots, PCA, correlation heatmaps

### 7. Peak Annotation (HOMER)

Differentially accessible peaks are annotated with genomic features:
- Promoter regions (TSS ± defined distance)
- Introns
- Exons
- Intergenic regions
- Distance to nearest transcription start site (TSS)
- Associated gene symbols

```bash
annotatePeaks.pl peaks.bed GRCm39 > annotations.txt
```

### 8. Motif Enrichment Analysis (HOMER)

Motif discovery and enrichment analysis is performed to identify transcription factor binding sites:

```bash
findMotifsGenome.pl peaks.bed GRCm39 output_dir/ -size 200
```

Analysis includes:
- *De novo* motif discovery
- Known motif enrichment
- Comparison to motif databases
- P-value and fold-enrichment calculations

### 9. Signal Visualization (deepTools)

Genome-wide chromatin accessibility signal tracks are generated as BigWig files using bamCoverage with default parameters:

```bash
bamCoverage -b input.bam -o output.bw
```

BigWig files enable:
- Visualization in genome browsers (IGV, UCSC)
- Signal comparison across samples
- Metagene profiling
- Heatmap generation

## Output Structure

```
results/
├── fastqc/                      # Raw read quality reports
│   ├── sample1_fastqc.html
│   └── sample1_fastqc.zip
├── aligned/                     # Alignment files
│   ├── sample1.bam
│   ├── sample1.bam.bai
│   └── sample1.flagstat
├── filtered/                    # Filtered BAM files (mito removed)
│   ├── sample1_filtered.bam
│   └── sample1_filtered.bam.bai
├── multiqc/                     # Aggregated QC report
│   └── multiqc_report.html
├── tag_directories/             # HOMER tag directories
│   └── sample1/
├── peaks/                       # Peak calling results
│   ├── sample1_peaks.txt
│   └── sample1_peaks.bed
├── diffbind/                    # Differential accessibility
│   ├── differential_peaks.csv
│   ├── normalized_counts.csv
│   ├── pca_plot.pdf
│   ├── ma_plot.pdf
│   └── correlation_heatmap.pdf
├── annotations/                 # Peak annotations
│   └── annotated_peaks.txt
├── motifs/                      # Motif enrichment results
│   ├── homerResults.html
│   ├── knownResults.txt
│   └── motif_files/
└── bigwig/                      # Coverage tracks
    └── sample1.bw
```

## Example Dataset

This pipeline was developed and validated using:
- **Data source**: NCBI GEO (GSE266583)
- **Organism**: Mouse (*Mus musculus*)
- **Genome**: GRCm39
- **Design**: Single-end sequencing with biological replicates

## Reproducibility

All analyses are performed using containerized tools from `ghcr.io/bf528` executed through Singularity, ensuring:
- Consistent software versions across environments
- Portable execution on HPC clusters and cloud platforms
- Complete reproducibility of results
- Simplified dependency management

## Best Practices

### Quality Control Thresholds
- Minimum read quality: Q20
- Adapter contamination: <5%
- Library complexity: Check fragment size distribution
- Alignment rate: >70% for ATAC-seq
- Mitochondrial contamination: <20% before filtering

### Experimental Design Recommendations
- Minimum 2-3 biological replicates per condition
- 20-50 million reads per sample recommended
- Include appropriate controls
- Balance sample processing to minimize batch effects

### Peak Calling Considerations
- ATAC-seq typically generates 40,000-150,000 peaks
- Lower peak counts may indicate poor sample quality
- Higher peak counts may indicate background noise

### Mitochondrial Read Removal Rationale
Removing mitochondrial reads is critical for ATAC-seq because:
- Mitochondria are highly abundant and transcriptionally active
- Mitochondrial fragments can represent 20-80% of reads
- These reads can mask true chromatin accessibility signal
- Removal improves signal-to-noise ratio

## Biological Interpretation

### Understanding ATAC-seq Results

**Open Chromatin Regions**: Peaks indicate regions of accessible chromatin where transcription factors can bind.

**Differential Accessibility**: Changes in peak intensity between conditions suggest:
- Altered transcription factor binding
- Chromatin remodeling
- Changes in gene regulatory activity

**Motif Enrichment**: Enriched motifs indicate:
- Transcription factors likely driving observed changes
- Regulatory networks active in each condition
- Potential mechanistic insights

### Tool Citations

- **FastQC**: Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.
- **Bowtie2**: Langmead & Salzberg (2012) Fast gapped-read alignment with Bowtie 2. Nature Methods.
- **SAMtools**: Danecek et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience.
- **MultiQC**: Ewels et al. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics.
- **HOMER**: Heinz et al. (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Molecular Cell.
- **DiffBind**: Ross-Innes et al. (2012) Differential oestrogen receptor binding is associated with clinical outcome in breast cancer. Nature.
- **deepTools**: Ramírez et al. (2016) deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research.

### Data Source

Data analyzed with this pipeline:
- **GEO Accession**: GSE266583
- **NCBI GEO**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266583

### ATAC-seq Method Citation

- **ATAC-seq**: Buenrostro et al. (2013) Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. Nature Methods.

## Troubleshooting

### Common Issues

**Issue**: High mitochondrial read percentage (>50%)
```bash
# Solution: Check sample quality and ATAC-seq protocol
# May indicate poor nuclei isolation or cell lysis issues
```

**Issue**: Low peak count (<20,000 peaks)
- Check alignment rates
- Verify sequencing depth is sufficient
- Assess fragment size distribution
- Review ATAC-seq library preparation

**Issue**: Singularity container fails to pull
```bash
# Solution: Manually pull the container
singularity pull docker://ghcr.io/bf528/atacseq-tools:latest
```

**Issue**: Memory errors during peak calling
```bash
# Solution: Increase memory in nextflow.config
process {
    withName: 'HOMER_FINDPEAKS' {
        memory = '32 GB'
    }
}
```

**Issue**: Poor replicate concordance
- Check QC metrics for individual samples
- Verify consistent sample processing
- Consider biological variability
- Review batch effects

### Quality Control Red Flags

- Alignment rate <60%
- Mitochondrial reads >80% before filtering
- Extreme GC bias
- High duplication rates (>50%)
- Unusual fragment size distribution


