#!/bin/bash
# filepath: rna_seq_pipeline.sh

# RNA-seq analysis pipeline
# Usage: ./rna_seq_pipeline.sh <reference_genome> <gtf_file> <output_dir>

set -e  # Exit on error

# Check input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference_genome> <gtf_file> <output_dir>"
    exit 1
fi

REFERENCE_GENOME=$1  # Reference genome in FASTA format
GTF_FILE=$2          # Gene annotation file
OUTPUT_DIR=$3        # Output directory

# Create output directory structure
mkdir -p ${OUTPUT_DIR}/{fastq,qc,trim,aligned,sorted,counts,results}

# Log file setup
LOG_FILE="${OUTPUT_DIR}/pipeline_log.txt"
exec > >(tee -a ${LOG_FILE}) 2>&1

echo "Starting RNA-seq pipeline at $(date)"

# Convert SRR files to FASTQ format
echo "Converting SRA to FASTQ..."
for SRR_FILE in *.sra; do
    SAMPLE=$(basename ${SRR_FILE} .sra)
    fastq-dump --split-files --gzip ${SRR_FILE} -O ${OUTPUT_DIR}/fastq/
    echo "Converted ${SRR_FILE} to FASTQ"
done

# Quality control with FastQC
echo "Running FastQC on raw reads..."
fastqc -o ${OUTPUT_DIR}/qc/ ${OUTPUT_DIR}/fastq/*.fastq.gz
multiqc ${OUTPUT_DIR}/qc/ -o ${OUTPUT_DIR}/qc/multiqc_report

# Trim reads for quality and adapters
echo "Trimming reads..."
for FQ1 in ${OUTPUT_DIR}/fastq/*_1.fastq.gz; do
    SAMPLE=$(basename ${FQ1} _1.fastq.gz)
    FQ2=${OUTPUT_DIR}/fastq/${SAMPLE}_2.fastq.gz
    
    # Check if paired-end or single-end
    if [ -f "$FQ2" ]; then
        # Paired-end trimming
        trimmomatic PE -phred33 ${FQ1} ${FQ2} \
            ${OUTPUT_DIR}/trim/${SAMPLE}_1_trimmed.fastq.gz ${OUTPUT_DIR}/trim/${SAMPLE}_1_unpaired.fastq.gz \
            ${OUTPUT_DIR}/trim/${SAMPLE}_2_trimmed.fastq.gz ${OUTPUT_DIR}/trim/${SAMPLE}_2_unpaired.fastq.gz \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        # Single-end trimming
        trimmomatic SE -phred33 ${FQ1} ${OUTPUT_DIR}/trim/${SAMPLE}_trimmed.fastq.gz \
            ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    fi
    echo "Trimmed ${SAMPLE}"
done

# Index the reference genome for BWA
echo "Indexing reference genome for BWA..."
bwa index ${REFERENCE_GENOME}

# Align reads to reference genome using BWA
echo "Aligning reads to reference genome..."
for SAMPLE in $(ls ${OUTPUT_DIR}/trim/*_1_trimmed.fastq.gz 2>/dev/null | sed 's/_1_trimmed.fastq.gz//g' | xargs -n 1 basename); do
    # Check if paired-end
    if [ -f "${OUTPUT_DIR}/trim/${SAMPLE}_2_trimmed.fastq.gz" ]; then
        bwa mem -t 8 ${REFERENCE_GENOME} \
            ${OUTPUT_DIR}/trim/${SAMPLE}_1_trimmed.fastq.gz \
            ${OUTPUT_DIR}/trim/${SAMPLE}_2_trimmed.fastq.gz \
            > ${OUTPUT_DIR}/aligned/${SAMPLE}.sam
    else
        bwa mem -t 8 ${REFERENCE_GENOME} \
            ${OUTPUT_DIR}/trim/${SAMPLE}_1_trimmed.fastq.gz \
            > ${OUTPUT_DIR}/aligned/${SAMPLE}.sam
    fi
    echo "Aligned ${SAMPLE}"
done

# Handle any single-end samples
for FQ in ${OUTPUT_DIR}/trim/*_trimmed.fastq.gz; do
    # Skip files that are part of paired-end data
    if [[ $FQ != *_1_trimmed.fastq.gz && $FQ != *_2_trimmed.fastq.gz ]]; then
        SAMPLE=$(basename ${FQ} _trimmed.fastq.gz)
        bwa mem -t 8 ${REFERENCE_GENOME} ${FQ} > ${OUTPUT_DIR}/aligned/${SAMPLE}.sam
        echo "Aligned ${SAMPLE} (single-end)"
    fi
done

# Convert SAM to sorted BAM and index
echo "Converting SAM to sorted BAM..."
for SAM in ${OUTPUT_DIR}/aligned/*.sam; do
    SAMPLE=$(basename ${SAM} .sam)
    samtools view -bS ${SAM} | samtools sort -o ${OUTPUT_DIR}/sorted/${SAMPLE}.sorted.bam
    samtools index ${OUTPUT_DIR}/sorted/${SAMPLE}.sorted.bam
    # Remove SAM file to save space
    rm ${SAM}
    echo "Processed ${SAMPLE} to sorted BAM"
done

# Count reads per gene using featureCounts
echo "Counting reads per gene..."
featureCounts -a ${GTF_FILE} -o ${OUTPUT_DIR}/counts/gene_counts.txt \
    -T 8 -p ${OUTPUT_DIR}/sorted/*.sorted.bam

# Generate count matrix for differential expression analysis
echo "Generating count matrix..."
Rscript - <<EOF
library(edgeR)
counts <- read.delim("${OUTPUT_DIR}/counts/gene_counts.txt", skip=1)
# Extract gene IDs and count columns
gene_ids <- counts[,1]
count_matrix <- counts[,7:ncol(counts)]
row.names(count_matrix) <- gene_ids
# Save the count matrix
write.table(count_matrix, file="${OUTPUT_DIR}/results/count_matrix.txt", 
            sep="\t", quote=FALSE, row.names=TRUE)
cat("Count matrix saved to ${OUTPUT_DIR}/results/count_matrix.txt\n")
EOF

echo "RNA-seq analysis pipeline completed at $(date)"
echo "Results can be found in ${OUTPUT_DIR}/results"
echo "Log file: ${LOG_FILE}"
