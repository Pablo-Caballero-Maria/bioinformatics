#!/bin/bash

GTF_FILE=$1
BAM_IDS=$2
OUT_FILE=$3

# Process the GTF file line by line
while IFS= read -r line; do
  # Skip comments
  if [[ $line =~ ^# ]]; then
    echo "$line" >>"$OUT_FILE"
    continue
  fi

  # Get the first field (chromosome identifier)
  chromosome_id_gtf=$(echo "$line" | awk '{print $1}')
  if [[ $chromosome_id_gtf == "MT" ]]; then
    nc_id="NC_012920.1"
    echo "$line" | awk -v nc="$nc_id" -v OFS="\t" '{ $1 = nc; print }' >>"$OUT_FILE"
    continue
  fi

  # Convert X and Y to numbers
  if [[ $chromosome_id_gtf == "X" ]]; then
    chromosome_number="23"
  elif [[ $chromosome_id_gtf == "Y" ]]; then
    chromosome_number="24"
  else
    chromosome_number=$chromosome_id_gtf
  fi

  # Find the corresponding NC identifier (add as many zeros as needed)
  printf -v chromosome_padded "%02d" "$chromosome_number"
  nc_id=$(grep -oP "NC_0000${chromosome_padded}\.\d+" "$BAM_IDS")
  # Replace the chromosome identifier with the NC identifier
  echo "$line" | awk -v nc="$nc_id" -v OFS="\t" '{ $1 = nc; print }' >>"$OUT_FILE"

done <"$GTF_FILE"
