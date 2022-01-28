#!/usr/bin/env bash

## This script indexes the "all_snv" VCFs with tabix

# Load the necessary modules
module load bio/BCFtools/1.9-foss-2019b

# Run the bgzip and tabix commands in a function
function index {
  bgzip < "../outputs/${vcf}.vcf" > "../outputs/${vcf}.vcf.gz"
  tabix -p vcf -f "../outputs/${vcf}.vcf.gz"
}

# VCFs of interest:
ns_chr="near_splice_all_snvs_chr"
ns_no_chr="near_splice_all_snvs_no_chr_prefix"
c_chr="coding_all_snvs_chr"
c_no_chr="coding_all_snvs_no_chr_prefix"

vcfs=( ${c_chr} ${c_no_chr} ${ns_chr} ${ns_no_chr} )

# Index each VCF with tabix
for vcf in ${vcfs[@]}
do
  echo "Compressing and indexing ${vcf}.vcf"
  index
done
