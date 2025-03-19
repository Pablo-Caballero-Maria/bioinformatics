list="SRR_list.txt"
while read line; do
  if [[ $line != "#"* ]]; then
    echo "Downloading $line"
    fastq-dump --split-files $line
  fi
done <$list
