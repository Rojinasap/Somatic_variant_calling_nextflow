# Somatic Variant Calling Pipeline - Nextflow

A comprehensive Nextflow pipeline for detecting somatic mutations from tumor-normal paired sequencing data.

## Overview

This pipeline performs end-to-end somatic variant calling analysis from paired tumor-normal FASTQ files. It includes quality control, alignment, post-processing, variant calling, and filtering steps.

## Pipeline Steps

1. **Quality Control** - FastQC analysis of raw sequencing reads
2. **Alignment** - BWA-MEM alignment to reference genome
3. **BAM Processing** - Sorting, marking duplicates, and base quality score recalibration (BQSR)
4. **Variant Calling** - Somatic variant detection using GATK Mutect2
5. **Variant Filtering** - Apply filters to remove false positives
6. **Report Generation** - MultiQC aggregation of quality metrics

## Requirements

### Software Dependencies

- Nextflow (>= 21.04.0)
- FastQC
- BWA
- SAMtools
- GATK (>= 4.0)
- MultiQC

Alternatively, use Docker or Singularity containers (see Execution Profiles below).

### Input Data

- Paired-end FASTQ files for tumor and matched normal samples
- Reference genome (FASTA format with `.fai` and `.dict` indexes)
- Optional: Known sites VCF for BQSR (e.g., dbSNP)
- Optional: Germline resource VCF for Mutect2
- Optional: Panel of normals VCF for Mutect2

## Quick Start

### 1. Prepare Samplesheet

Create a CSV file with the following format:

```csv
patient,sample,type,fastq_1,fastq_2
patient1,tumor1,tumor,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
patient1,normal1,normal,/path/to/normal_R1.fastq.gz,/path/to/normal_R2.fastq.gz
```

**Columns:**
- `patient`: Patient identifier (samples with same patient ID will be paired)
- `sample`: Unique sample identifier
- `type`: Either "tumor" or "normal"
- `fastq_1`: Path to Read 1 FASTQ file
- `fastq_2`: Path to Read 2 FASTQ file

An example samplesheet is provided: `samplesheet.csv`

### 2. Run the Pipeline

**Basic command:**

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --reference /path/to/reference.fa \
  --outdir results
```

**With optional resources:**

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --reference /path/to/reference.fa \
  --known_sites /path/to/dbsnp.vcf.gz \
  --germline_resource /path/to/af-only-gnomad.vcf.gz \
  --panel_of_normals /path/to/pon.vcf.gz \
  --outdir results
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input samplesheet CSV file |
| `--reference` | Path to reference genome FASTA file |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `results` |
| `--reference_fai` | Path to reference index (.fai) | `${reference}.fai` |
| `--reference_dict` | Path to reference dictionary (.dict) | Auto-detected |
| `--known_sites` | Path to known sites VCF for BQSR | None |
| `--known_sites_idx` | Path to known sites VCF index | None |
| `--germline_resource` | Path to germline resource VCF | None |
| `--germline_resource_idx` | Path to germline resource index | None |
| `--panel_of_normals` | Path to panel of normals VCF | None |
| `--panel_of_normals_idx` | Path to panel of normals index | None |

## Execution Profiles

The pipeline supports multiple execution profiles:

### Standard (Local Execution)

```bash
nextflow run main.nf -profile standard --input samplesheet.csv --reference genome.fa
```

### Docker

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --reference genome.fa
```

### Singularity

```bash
nextflow run main.nf -profile singularity --input samplesheet.csv --reference genome.fa
```

### SLURM Cluster

```bash
nextflow run main.nf -profile slurm --input samplesheet.csv --reference genome.fa
```

### SGE Cluster

```bash
nextflow run main.nf -profile sge --input samplesheet.csv --reference genome.fa
```

## Output Structure

```
results/
├── fastqc/                    # FastQC quality control reports
├── aligned/                   # Aligned BAM files
├── markduplicates/           # Duplicate marking metrics
├── bqsr/                     # Base recalibration tables
├── variants/                 # Raw variant calls (VCF)
├── variants_filtered/        # Filtered variant calls
├── multiqc/                  # Aggregated QC report
└── pipeline_info/            # Pipeline execution reports
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

## Resource Configuration

The pipeline automatically configures resources for each process. You can customize these in `nextflow.config` or override them at runtime:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --reference genome.fa \
  --max_cpus 16 \
  --max_memory 64.GB
```

## Resume Capability

Nextflow supports resuming interrupted runs:

```bash
nextflow run main.nf -resume --input samplesheet.csv --reference genome.fa
```

## Troubleshooting

### Common Issues

1. **Missing reference indexes**: Ensure `.fai` and `.dict` files exist for your reference genome
   ```bash
   samtools faidx reference.fa
   gatk CreateSequenceDictionary -R reference.fa
   ```

2. **Memory errors**: Increase memory allocation in `nextflow.config` for specific processes

3. **Container issues**: Ensure Docker/Singularity is properly installed and configured

### Getting Help

```bash
nextflow run main.nf --help
```

## Citation

If you use this pipeline in your research, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
- **GATK**: McKenna, A., et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297-1303.
- **BWA**: Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25(14), 1754-1760.

## License

This pipeline is distributed under the MIT License.

## Contact

For questions or issues, please open an issue on the GitHub repository.
