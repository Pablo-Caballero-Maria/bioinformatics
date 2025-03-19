for b in *_unmapped.bam; do
  sample=$(basename $b _unmapped.bam)

  echo "Transforming bam to fastq for $sample"
  samtools fastq -1 ${sample}_R1_001.fastq -2 ${sample}_R2_001.fastq -n $b
  echo "Done :)"
done
