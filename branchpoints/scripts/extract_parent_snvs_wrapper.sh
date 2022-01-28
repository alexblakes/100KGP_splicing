#!/usr/bin/env bash

## This is a wrapper script for the main extract_parent_snvs.sh script
## The array job takes ~3 days to run.

dir="../outputs/MAPS_snvs"
dir_hpc_out="${dir}/hpc_out"
dir_hpc_err="${dir}/hpc_err"
dir_out="${dir}/out"

# Empty the output directories
rm -r ${dir}; \
mkdir ${dir} && \
mkdir ${dir_hpc_out} && \
mkdir ${dir_hpc_err} && \
mkdir ${dir_out}

# Submit the array job to the HPC
bsub < extract_parent_snvs.sh

# Wait for the job to complete before advancing
bwait -w 'ended(MAPS_snvs_)'

# Concatenate the outputs
cat ../outputs/MAPS_snvs/out/* > ../outputs/unaff_parents_branch_snvs.tsv
