#!/usr/bin/env bash

## This script sends the VEP annotation script to an array job on the HPC
bsub < vep_snvs.sh
bwait -w 'ended(vep_snvs)'

# Concatenate the VEP output:
cat ../outputs/vep/out/* > ../outputs/unaff_parents_coding_snvs_vep_out.tsv
