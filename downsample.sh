#!/usr/bin/bash
for file in *.fastq; do
  echo "Downsampling $file"
  seqtk sample -s100 $file 3000000 >downsampled/$file
done
