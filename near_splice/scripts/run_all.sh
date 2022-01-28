#!/usr/bin/env bash

## This bash script runs all of the scripts in the near_splice directory.

## Reformat and annotate the GENCODE .gtf
python3 exon_filter.py # Filters exons from GENCODE
python3 near_splice_sites.py # Annotates near-splice positions
python3 coding_sites.py # Annotates coding positions
python3 positions_to_bed.py # Converts the above outputs to .bed format

## phyloP
python3 phylop.py # Get phylop scores for every near-splice position

## SpliceAI
# NB sequences are annotated against the + strand, because downstream all
# variants in VCFs are + stranded.
bash get_fasta.sh # Get sequence contexts with bedtools
# The below script is memory intensive, using ~36GB memory & takes ~10 mins to
# run
python3 all_possible_variants_to_vcf.py # Annotate all possible SNVs
bash index_all_snv_vcfs.sh # Index VCFs with tabix
# The next script takes ~30 mins to run.
bash extract_spliceai_scores.sh # Extract pre-computed SpliceAI scores.
python3 spliceai_stats.py # Get summary stats for the SpliceAI scores

## MAPS
python3 unaffected_parents.py # Identify unaffected parents
python3 positions_for_bcftools.py # Get tab-separated positions file
# The command below took ~3 days to run (large array job to HPC)
# The wrapper script waits for the array job to complete before proceeding
bash extract_parent_snvs_wrapper.sh # Run the "wrapper" script
bash vep_directories.sh # Empty the VEP input/output directories
python3 vep_input_format.py # Reformat the coding SNVs to VCF format
bash vep_wrapper.sh # Send the VEP script to an array job
python3 tidy_near_splice_and_coding_snvs.py # Collate the SNVs
# NB variants are given on the forward strand.
# The mutability of reverse-complemented contexts (+ and - strand) is the same.
python3 MAPS.py # Calculate MAPS
