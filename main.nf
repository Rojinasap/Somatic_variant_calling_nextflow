#!/usr/bin/env nextflow

/*
 * Somatic Variant Calling Pipeline
 * 
 * This pipeline performs somatic variant calling from tumor-normal paired samples
 * 
 * Pipeline steps:
 * 1. Quality control (FastQC)
 * 2. Alignment (BWA-MEM)
 * 3. BAM post-processing (sorting, deduplication, base quality score recalibration)
 * 4. Somatic variant calling (Mutect2)
 * 5. Variant filtering
 */

nextflow.enable.dsl=2

// Parameters
params.input = null
params.outdir = 'results'
params.reference = null
params.reference_dict = null
params.reference_fai = null
params.known_sites = null
params.known_sites_idx = null
params.germline_resource = null
params.germline_resource_idx = null
params.panel_of_normals = null
params.panel_of_normals_idx = null
params.help = false

// Help message
def helpMessage() {
    log.info """
    ================================================================
    Somatic Variant Calling Pipeline
    ================================================================
    
    Usage:
      nextflow run main.nf --input samplesheet.csv --reference genome.fa [options]
    
    Required arguments:
      --input               Path to input samplesheet CSV file
      --reference           Path to reference genome FASTA file
      --outdir              Output directory (default: results)
    
    Optional arguments:
      --reference_dict      Path to reference dictionary file
      --reference_fai       Path to reference index file
      --known_sites         Path to known sites VCF for BQSR
      --known_sites_idx     Path to known sites VCF index
      --germline_resource   Path to germline resource VCF for Mutect2
      --germline_resource_idx Path to germline resource VCF index
      --panel_of_normals    Path to panel of normals VCF for Mutect2
      --panel_of_normals_idx Path to panel of normals VCF index
    
    Input samplesheet format (CSV):
      patient,sample,type,fastq_1,fastq_2
      patient1,tumor1,tumor,/path/to/tumor_R1.fastq.gz,/path/to/tumor_R2.fastq.gz
      patient1,normal1,normal,/path/to/normal_R1.fastq.gz,/path/to/normal_R2.fastq.gz
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input) {
    log.error "Error: --input parameter is required"
    helpMessage()
    exit 1
}

if (!params.reference) {
    log.error "Error: --reference parameter is required"
    helpMessage()
    exit 1
}

/*
 * PROCESS: Quality Control with FastQC
 */
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(patient), val(sample_id), val(type), path(fastq_1), path(fastq_2)
    
    output:
    tuple val(patient), val(sample_id), val(type), path("*_fastqc.{html,zip}"), emit: fastqc_results
    
    script:
    """
    fastqc -t 2 ${fastq_1} ${fastq_2}
    """
}

/*
 * PROCESS: Alignment with BWA-MEM
 */
process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(patient), val(sample_id), val(type), path(fastq_1), path(fastq_2)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(patient), val(sample_id), val(type), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}' \
        ${reference} ${fastq_1} ${fastq_2} | \
        samtools view -@ ${task.cpus} -Sb - > ${sample_id}.bam
    """
}

/*
 * PROCESS: Sort BAM files
 */
process SORT_BAM {
    tag "$sample_id"
    
    input:
    tuple val(patient), val(sample_id), val(type), path(bam)
    
    output:
    tuple val(patient), val(sample_id), val(type), path("${sample_id}.sorted.bam"), emit: sorted_bam
    
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

/*
 * PROCESS: Mark Duplicates
 */
process MARK_DUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/markduplicates", mode: 'copy', pattern: "*.metrics.txt"
    
    input:
    tuple val(patient), val(sample_id), val(type), path(bam)
    
    output:
    tuple val(patient), val(sample_id), val(type), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: dedup_bam
    path "${sample_id}.metrics.txt", emit: metrics
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.dedup.bam \
        -M ${sample_id}.metrics.txt \
        --CREATE_INDEX true
    """
}

/*
 * PROCESS: Base Quality Score Recalibration (BQSR)
 */
process BASE_RECALIBRATION {
    tag "$sample_id"
    publishDir "${params.outdir}/bqsr", mode: 'copy', pattern: "*.recal_data.table"
    
    input:
    tuple val(patient), val(sample_id), val(type), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path known_sites
    path known_sites_idx
    
    output:
    tuple val(patient), val(sample_id), val(type), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bam.bai"), emit: recal_bam
    path "${sample_id}.recal_data.table", emit: recal_table
    
    script:
    known_sites_arg = known_sites.name != 'NO_FILE' ? "--known-sites ${known_sites}" : ""
    """
    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${reference} \
        ${known_sites_arg} \
        -O ${sample_id}.recal_data.table
    
    gatk ApplyBQSR \
        -I ${bam} \
        -R ${reference} \
        --bqsr-recal-file ${sample_id}.recal_data.table \
        -O ${sample_id}.recal.bam
    
    samtools index ${sample_id}.recal.bam
    """
}

/*
 * PROCESS: Somatic Variant Calling with Mutect2
 */
process MUTECT2 {
    tag "$patient"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(patient), val(tumor_id), path(tumor_bam), path(tumor_bai), val(normal_id), path(normal_bam), path(normal_bai)
    path reference
    path reference_fai
    path reference_dict
    path germline_resource
    path germline_resource_idx
    path panel_of_normals
    path panel_of_normals_idx
    
    output:
    tuple val(patient), path("${patient}.vcf.gz"), path("${patient}.vcf.gz.tbi"), emit: vcf
    path "${patient}.vcf.gz.stats", emit: stats
    
    script:
    germline_arg = germline_resource.name != 'NO_FILE' ? "--germline-resource ${germline_resource}" : ""
    pon_arg = panel_of_normals.name != 'NO_FILE' ? "--panel-of-normals ${panel_of_normals}" : ""
    """
    gatk Mutect2 \
        -R ${reference} \
        -I ${tumor_bam} \
        -I ${normal_bam} \
        -tumor ${tumor_id} \
        -normal ${normal_id} \
        ${germline_arg} \
        ${pon_arg} \
        -O ${patient}.vcf.gz
    """
}

/*
 * PROCESS: Filter Mutect2 Calls
 */
process FILTER_MUTECT_CALLS {
    tag "$patient"
    publishDir "${params.outdir}/variants_filtered", mode: 'copy'
    
    input:
    tuple val(patient), path(vcf), path(vcf_idx)
    path reference
    path reference_fai
    path reference_dict
    
    output:
    tuple val(patient), path("${patient}.filtered.vcf.gz"), path("${patient}.filtered.vcf.gz.tbi"), emit: filtered_vcf
    
    script:
    """
    gatk FilterMutectCalls \
        -R ${reference} \
        -V ${vcf} \
        -O ${patient}.filtered.vcf.gz
    """
}

/*
 * PROCESS: MultiQC - Aggregate QC reports
 */
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    
    script:
    """
    multiqc .
    """
}

/*
 * WORKFLOW
 */
workflow {
    // Parse input samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            tuple(row.patient, row.sample, row.type, file(row.fastq_1), file(row.fastq_2))
        }
        .set { input_samples }
    
    // Prepare reference files
    reference_ch = Channel.fromPath(params.reference)
    reference_fai_ch = params.reference_fai ? Channel.fromPath(params.reference_fai) : Channel.fromPath("${params.reference}.fai")
    reference_dict_ch = params.reference_dict ? Channel.fromPath(params.reference_dict) : Channel.fromPath(params.reference.replaceAll(/\.fa(sta)?$/, '.dict'))
    
    // Prepare known sites (optional)
    known_sites_ch = params.known_sites ? Channel.fromPath(params.known_sites) : Channel.value(file('NO_FILE'))
    known_sites_idx_ch = params.known_sites_idx ? Channel.fromPath(params.known_sites_idx) : Channel.value(file('NO_FILE'))
    
    // Prepare germline resource (optional)
    germline_resource_ch = params.germline_resource ? Channel.fromPath(params.germline_resource) : Channel.value(file('NO_FILE'))
    germline_resource_idx_ch = params.germline_resource_idx ? Channel.fromPath(params.germline_resource_idx) : Channel.value(file('NO_FILE'))
    
    // Prepare panel of normals (optional)
    panel_of_normals_ch = params.panel_of_normals ? Channel.fromPath(params.panel_of_normals) : Channel.value(file('NO_FILE'))
    panel_of_normals_idx_ch = params.panel_of_normals_idx ? Channel.fromPath(params.panel_of_normals_idx) : Channel.value(file('NO_FILE'))
    
    // Quality control
    FASTQC(input_samples)
    
    // Alignment
    BWA_MEM(
        input_samples,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch
    )
    
    // Sort BAM
    SORT_BAM(BWA_MEM.out.bam)
    
    // Mark duplicates
    MARK_DUPLICATES(SORT_BAM.out.sorted_bam)
    
    // Base quality recalibration
    BASE_RECALIBRATION(
        MARK_DUPLICATES.out.dedup_bam,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch,
        known_sites_ch,
        known_sites_idx_ch
    )
    
    // Group tumor-normal pairs by patient
    BASE_RECALIBRATION.out.recal_bam
        .branch {
            tumor: it[2] == 'tumor'
            normal: it[2] == 'normal'
        }
        .set { sample_types }
    
    // Join tumor and normal samples for the same patient
    tumor_normal_pairs = sample_types.tumor
        .map { patient, sample_id, type, bam, bai -> tuple(patient, sample_id, bam, bai) }
        .combine(
            sample_types.normal.map { patient, sample_id, type, bam, bai -> tuple(patient, sample_id, bam, bai) },
            by: 0
        )
    
    // Somatic variant calling
    MUTECT2(
        tumor_normal_pairs,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch,
        germline_resource_ch,
        germline_resource_idx_ch,
        panel_of_normals_ch,
        panel_of_normals_idx_ch
    )
    
    // Filter variants
    FILTER_MUTECT_CALLS(
        MUTECT2.out.vcf,
        reference_ch,
        reference_fai_ch,
        reference_dict_ch
    )
    
    // Aggregate QC reports
    MULTIQC(
        FASTQC.out.fastqc_results.map { it[3] }.collect().ifEmpty([])
            .mix(MARK_DUPLICATES.out.metrics.collect().ifEmpty([]))
    )
}

workflow.onComplete {
    log.info """
    ================================================================
    Pipeline completed!
    ================================================================
    Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:  ${workflow.duration}
    Output:    ${params.outdir}
    ================================================================
    """.stripIndent()
}
