#!/usr/bin/env bash

## Create and empty the output directories for running VEP on the coding SNVs
rm -r ../outputs/vep ; \
mkdir ../outputs/vep && \
mkdir ../outputs/vep/in && \
mkdir ../outputs/vep/out && \
mkdir ../outputs/vep/hpc_out && \
mkdir ../outputs/vep/hpc_err
