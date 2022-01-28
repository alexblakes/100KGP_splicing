#!/usr/bin/env bash

## This script indexes the "all_snv" VCFs with tabix

# Load the necessary modules
module load bio/BCFtools/1.9-foss-2019b

function index {
  bgzip < "../outputs/${vcf}.vcf" > "../outputs/${vcf}.vcf.gz"
  tabix -p vcf -f "../outputs/${vcf}.vcf.gz"
}

# VCFs of interest:
chr="branch_all_snvs_chr"
no_chr="branch_all_snvs_no_chr_prefix"

vcfs=( ${chr} ${no_chr} )

# Index each VCF with tabix
for vcf in ${vcfs[@]}
do
  echo "Compressing and indexing ${vcf}.vcf"
  index
done
