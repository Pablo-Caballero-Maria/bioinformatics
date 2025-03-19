#!/bin/bash

GTF_FILE="$1"

awk '!/^#/ {print $1}' "$GTF_FILE" | sort -Vu >unique_chromosomes_gtf.txt
