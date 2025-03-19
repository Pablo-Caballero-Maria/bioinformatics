#!/bin/bash

# filepath: variant_calling_pipeline.sh
# Variant calling pipeline using GATK and VEP for SRR files

set -e  # Exit on error
set -o pipefail  # Pipe will exit with last non-zero status

# Default parameters - modify as needed
THREADS=4
REFERENCE="/path/to/reference.fasta"
KNOWN_VARIANTS="/path/to/known_variants.vcf.gz"
OUTPUT_DIR="./output"

# Create output directories
mkdir -p ${OUTPUT_DIR}/{fastq,aligned,processed,variants,logs,metrics}

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a ${OUTPUT_DIR}/logs/pipeline.log
}

log "Starting variant calling pipeline"

# Process command line arguments
SRR_LIST=$1
if [ -z "$SRR_LIST" ]; then
    echo "Usage: $0 <srr_accessions_file>"
    echo "Example: $0 srr_list.txt"
    exit 1
fi

# Download and prepare reference if needed
if [ ! -f ${REFERENCE}.fai ]; then
    log "Indexing reference genome"
    samtools faidx ${REFERENCE}
    bwa index ${REFERENCE}
    gatk CreateSequenceDictionary -R ${REFERENCE}
fi

# Process each SRR accession
cat ${SRR_LIST} | while read SRR; do
    log "Processing sample: ${SRR}"
    
    # Download SRR files if they don't exist
    if [ ! -f ${OUTPUT_DIR}/fastq/${SRR}_1.fastq.gz ]; then
        log "Downloading ${SRR}"
        fasterq-dump ${SRR} -O ${OUTPUT_DIR}/fastq -e ${THREADS}
        pigz -p ${THREADS} ${OUTPUT_DIR}/fastq/${SRR}*.fastq
    fi
    
    # Quality control with FastQC
    log "Running FastQC"
    fastqc -o ${OUTPUT_DIR}/logs -t ${THREADS} ${OUTPUT_DIR}/fastq/${SRR}*.fastq.gz
    
    # Trim reads with Trimmomatic
    log "Trimming reads"
    trimmomatic PE -threads ${THREADS} \
        ${OUTPUT_DIR}/fastq/${SRR}_1.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_2.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_1.trimmed.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_1.unpaired.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_2.trimmed.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_2.unpaired.fastq.gz \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # Align with BWA-MEM
    log "Aligning reads"
    bwa mem -t ${THREADS} -R "@RG\tID:${SRR}\tSM:${SRR}\tPL:ILLUMINA" \
        ${REFERENCE} \
        ${OUTPUT_DIR}/fastq/${SRR}_1.trimmed.fastq.gz \
        ${OUTPUT_DIR}/fastq/${SRR}_2.trimmed.fastq.gz | \
        samtools view -Sb - | \
        samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/aligned/${SRR}.bam -
    
    # Mark duplicates
    log "Marking duplicates"
    gatk MarkDuplicatesSpark \
        -I ${OUTPUT_DIR}/aligned/${SRR}.bam \
        -O ${OUTPUT_DIR}/processed/${SRR}.markdup.bam \
        -M ${OUTPUT_DIR}/metrics/${SRR}.markdup_metrics.txt \
        --spark-master local[${THREADS}]
    
    # Index BAM
    log "Indexing BAM file"
    samtools index ${OUTPUT_DIR}/processed/${SRR}.markdup.bam
    
    # Base quality score recalibration
    log "Base quality score recalibration"
    gatk BaseRecalibrator \
        -R ${REFERENCE} \
        -I ${OUTPUT_DIR}/processed/${SRR}.markdup.bam \
        --known-sites ${KNOWN_VARIANTS} \
        -O ${OUTPUT_DIR}/processed/${SRR}.recal_data.table
    
    gatk ApplyBQSR \
        -R ${REFERENCE} \
        -I ${OUTPUT_DIR}/processed/${SRR}.markdup.bam \
        --bqsr-recal-file ${OUTPUT_DIR}/processed/${SRR}.recal_data.table \
        -O ${OUTPUT_DIR}/processed/${SRR}.recalibrated.bam
    
    # Collect alignment metrics
    log "Collecting alignment metrics"
    gatk CollectAlignmentSummaryMetrics \
        -R ${REFERENCE} \
        -I ${OUTPUT_DIR}/processed/${SRR}.recalibrated.bam \
        -O ${OUTPUT_DIR}/metrics/${SRR}.alignment_metrics.txt
    
    gatk CollectInsertSizeMetrics \
        -I ${OUTPUT_DIR}/processed/${SRR}.recalibrated.bam \
        -O ${OUTPUT_DIR}/metrics/${SRR}.insert_size_metrics.txt \
        -H ${OUTPUT_DIR}/metrics/${SRR}.insert_size_histogram.pdf
    
    # Variant calling
    log "Calling variants with HaplotypeCaller"
    gatk HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${OUTPUT_DIR}/processed/${SRR}.recalibrated.bam \
        -O ${OUTPUT_DIR}/variants/${SRR}.raw.vcf.gz \
        -ERC GVCF
    
    # Select SNPs and Indels
    log "Selecting variants"
    gatk SelectVariants \
        -R ${REFERENCE} \
        -V ${OUTPUT_DIR}/variants/${SRR}.raw.vcf.gz \
        -select-type SNP \
        -O ${OUTPUT_DIR}/variants/${SRR}.snps.vcf.gz
    
    gatk SelectVariants \
        -R ${REFERENCE} \
        -V ${OUTPUT_DIR}/variants/${SRR}.raw.vcf.gz \
        -select-type INDEL \
        -O ${OUTPUT_DIR}/variants/${SRR}.indels.vcf.gz
    
    # Annotate variants with VEP
    log "Annotating variants with VEP"
    vep --cache --offline \
        --species homo_sapiens \
        --assembly GRCh38 \
        --input_file ${OUTPUT_DIR}/variants/${SRR}.snps.vcf.gz \
        --output_file ${OUTPUT_DIR}/variants/${SRR}.snps.annotated.vcf \
        --vcf --symbol --terms SO --sift b --polyphen b \
        --hgvs --regulatory --force_overwrite
    
    vep --cache --offline \
        --species homo_sapiens \
        --assembly GRCh38 \
        --input_file ${OUTPUT_DIR}/variants/${SRR}.indels.vcf.gz \
        --output_file ${OUTPUT_DIR}/variants/${SRR}.indels.annotated.vcf \
        --vcf --symbol --terms SO \
        --hgvs --regulatory --force_overwrite
    
    log "Completed processing sample: ${SRR}"
done

log "Pipeline completed successfully"
