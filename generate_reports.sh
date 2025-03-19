#!/usr/bin/bash
for file in *.fastq; do
  echo "Generating report of: $file"
  fastqc $file
done
