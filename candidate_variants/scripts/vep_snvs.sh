#!/usr/bin/env bash

#BSUB -q inter
#BSUB -P re_gecip_machine_learning
#BSUB -J "vep_snvs"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000
#BSUB -e ../outputs/%J_err.txt
#BSUB -o ../outputs/%J_out.txt

# This script annotates GEL DNMs with VEP, in order to find the HGNC ID of the
# affected gene. This will allow merging with the G2P data.

module load bio/VEP/99.1-foss-2019a-Perl-5.28.1

vep \
  --input_file ../outputs/near_splice_and_branch_dnms_chr.vcf \
  --output_file ../outputs/near_splice_and_branch_dnms_vep.tsv \
  --tab \
  --species homo_sapiens \
  --assembly GRCh38 \
  --offline \
  --cache \
  --dir_cache ${CACHEDIR} \
  --cache_version 99 \
  --fork 4 \
  --buffer_size 100000 \
  --per_gene \
  --no_stats \
  --symbol \
  --fields "Uploaded_variation,Gene,Feature,HGNC_ID,Consequence"
