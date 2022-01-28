#!/usr/bin/env bash

## This script identifies the sequence context for each position using bedtools
## Sequence contexts are given in the + strand.

module load bio/BEDTools/2.27.1-foss-2018b

function get_seqs {
  bedtools getfasta \
    -fi "${fi}" \
    -bed ../outputs/"${x}".bed \
    -fo ../outputs/"${x}"_contexts.tsv \
    -tab
}

fi="/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
regions=( "near_splice" "coding" )

for x in "${regions[@]}"
do
  echo "Retrieving ${x} sequences"
  get_seqs
done
