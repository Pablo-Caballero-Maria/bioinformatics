for b in *.bam; do
  sample=$(basename $b .bam)

  echo "Getting unmapped reads from $sample"
  samtools view -b -f 12 $b > ${sample}_unmapped.bam
  echo "Done :)"
done
