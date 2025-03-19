#!/bin/bash
# filepath: 16S_processing_pipeline.sh

# Exit on error
set -e

# Script for processing 16S rRNA sequencing data from SRR files
# Usage: ./16S_processing_pipeline.sh input_directory output_directory

# Check if required arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=4  # Adjust based on your system

# Create necessary directories
mkdir -p ${OUTPUT_DIR}/{fastq,qc,trimmed,merged,asv,taxonomy,diversity}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check required tools
echo "Checking required tools..."
required_tools=(fastq-dump fastqc trimmomatic qiime)
for tool in "${required_tools[@]}"; do
    if ! command_exists $tool; then
        echo "Error: $tool is not installed. Please install it before running this script."
        exit 1
    fi
done

# Convert SRR files to FASTQ
echo "Converting SRR files to FASTQ..."
for srr_file in ${INPUT_DIR}/*.sra; do
    base_name=$(basename "$srr_file" .sra)
    if [ ! -f "${OUTPUT_DIR}/fastq/${base_name}_1.fastq.gz" ]; then
        fastq-dump --split-files --gzip "$srr_file" -O "${OUTPUT_DIR}/fastq/"
    fi
done

# Quality control with FastQC
echo "Running FastQC..."
mkdir -p ${OUTPUT_DIR}/qc/fastqc
fastqc ${OUTPUT_DIR}/fastq/*.fastq.gz -o ${OUTPUT_DIR}/qc/fastqc -t $THREADS

# Trim low quality reads and adapters
echo "Trimming reads..."
for forward in ${OUTPUT_DIR}/fastq/*_1.fastq.gz; do
    reverse=${forward/_1./_2.}
    base_name=$(basename "$forward" _1.fastq.gz)
    
    if [ -f "$reverse" ]; then
        # Paired-end trimming
        trimmomatic PE -threads $THREADS \
            "$forward" "$reverse" \
            ${OUTPUT_DIR}/trimmed/${base_name}_1_trimmed.fastq.gz ${OUTPUT_DIR}/trimmed/${base_name}_1_unpaired.fastq.gz \
            ${OUTPUT_DIR}/trimmed/${base_name}_2_trimmed.fastq.gz ${OUTPUT_DIR}/trimmed/${base_name}_2_unpaired.fastq.gz \
            ILLUMINACLIP:/path/to/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        # Single-end trimming
        trimmomatic SE -threads $THREADS \
            "$forward" \
            ${OUTPUT_DIR}/trimmed/${base_name}_trimmed.fastq.gz \
            ILLUMINACLIP:/path/to/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    fi
done

# Initialize QIIME 2 environment
echo "Initializing QIIME 2 analysis..."
cd ${OUTPUT_DIR}

# Import data into QIIME2 format
echo "Importing data..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${OUTPUT_DIR}/trimmed \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ${OUTPUT_DIR}/qiime2/demux.qza

# Generate visualization of sequence quality
qiime demux summarize \
  --i-data ${OUTPUT_DIR}/qiime2/demux.qza \
  --o-visualization ${OUTPUT_DIR}/qiime2/demux.qzv

# Denoise with DADA2 (generates ASVs)
echo "Denoising sequences and generating ASVs..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${OUTPUT_DIR}/qiime2/demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --o-table ${OUTPUT_DIR}/asv/table.qza \
  --o-representative-sequences ${OUTPUT_DIR}/asv/rep-seqs.qza \
  --o-denoising-stats ${OUTPUT_DIR}/asv/denoising-stats.qza \
  --p-n-threads $THREADS

# Taxonomic classification using Silva database
echo "Assigning taxonomy..."
qiime feature-classifier classify-sklearn \
  --i-classifier /path/to/silva-138-99-515-806-nb-classifier.qza \
  --i-reads ${OUTPUT_DIR}/asv/rep-seqs.qza \
  --o-classification ${OUTPUT_DIR}/taxonomy/taxonomy.qza

# Generate taxonomy visualization
qiime metadata tabulate \
  --m-input-file ${OUTPUT_DIR}/taxonomy/taxonomy.qza \
  --o-visualization ${OUTPUT_DIR}/taxonomy/taxonomy.qzv

# Basic diversity analysis
echo "Generating diversity metrics..."

# Create phylogenetic tree for diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${OUTPUT_DIR}/asv/rep-seqs.qza \
  --o-alignment ${OUTPUT_DIR}/diversity/aligned-rep-seqs.qza \
  --o-masked-alignment ${OUTPUT_DIR}/diversity/masked-aligned-rep-seqs.qza \
  --o-tree ${OUTPUT_DIR}/diversity/unrooted-tree.qza \
  --o-rooted-tree ${OUTPUT_DIR}/diversity/rooted-tree.qza

# Calculate core diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${OUTPUT_DIR}/diversity/rooted-tree.qza \
  --i-table ${OUTPUT_DIR}/asv/table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file /path/to/sample-metadata.tsv \
  --output-dir ${OUTPUT_DIR}/diversity/core-metrics-results

# Generate alpha diversity visualization
qiime diversity alpha-group-significance \
  --i-alpha-diversity ${OUTPUT_DIR}/diversity/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /path/to/sample-metadata.tsv \
  --o-visualization ${OUTPUT_DIR}/diversity/faith-pd-group-significance.qzv

# Generate beta diversity visualization
qiime diversity beta-group-significance \
  --i-distance-matrix ${OUTPUT_DIR}/diversity/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /path/to/sample-metadata.tsv \
  --m-metadata-column your_condition_column \
  --o-visualization ${OUTPUT_DIR}/diversity/unweighted-unifrac-condition-significance.qzv \
  --p-pairwise

# Generate taxa bar plots
qiime taxa barplot \
  --i-table ${OUTPUT_DIR}/asv/table.qza \
  --i-taxonomy ${OUTPUT_DIR}/taxonomy/taxonomy.qza \
  --m-metadata-file /path/to/sample-metadata.tsv \
  --o-visualization ${OUTPUT_DIR}/taxonomy/taxa-bar-plots.qzv

echo "16S analysis complete. Results are in ${OUTPUT_DIR}"
