#!/bin/bash

BAM_FILE="$1"

# samtools view "$BAM_FILE" | awk '{print $3}' | grep -E 'NC_0000(0[1-9]|1[0-9]|2[0-4])\.[0-9]+' | sort -u > unique_chromosomes_bam.txt

samtools view "$BAM_FILE" | awk '{print $3}' | sort -u > unique_chromosomes_bam.txt
