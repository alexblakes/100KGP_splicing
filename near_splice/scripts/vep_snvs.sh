#!/usr/bin/env bash

#BSUB -q short
#BSUB -P re_gecip_machine_learning
#BSUB -J "vep_snvs[1-22]"
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000
#BSUB -e ../outputs/vep/hpc_err/%J_%I_err.txt
#BSUB -o ../outputs/vep/hpc_out/%J_%I_out.txt

# This script annotates all the coding SNVs in the unaffected parents with VEP.

module load bio/VEP/99.1-foss-2019a-Perl-5.28.1

vep \
  --input_file ../outputs/vep/in/coding_vep_input_${LSB_JOBINDEX}.vcf \
  --output_file ../outputs/vep/out/coding_vep_output_${LSB_JOBINDEX}.tsv \
  --tab \
  --species homo_sapiens \
  --assembly GRCh38 \
  --offline \
  --cache \
  --dir_cache ${CACHEDIR} \
  --cache_version 99 \
  --fork 4 \
  --buffer_size 100000 \
  --pick \
  --no_stats
