for fq in *_R1_001.fastq; do
  sample=$(basename $fq _R1_001.fastq)

  echo "Mapping $sample"

  bwa mem GCF_000001405.26_GRCh38_genomic.fna \
    ${sample}_R1_001.fastq ${sample}_R2_001.fastq |
    samtools view -bS -o ${sample}.bam

  echo "Done :)"
done
