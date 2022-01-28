#!/usr/bin/env bash

## This bash script runs all of the scripts relevant to the branchpoint data.

## Preliminary
python3 branchpoints_bed.py # Convert branchpoint positions to .bed

## phyloP scores
python3 phylop_branch.py # Get phylop scores for every near-splice position

## SpliceAI
# NB sequences are annotated against the + strand, because downstream all
# variants in VCFs are + stranded.
bash get_fasta_branch.sh # Get sequence contexts with bedtools
python3 all_possible_variants_to_vcf.py # Annotate all possible SNVs
bash index_all_snv_vcfs.sh # Index VCFs with tabix
# The next script takes ~30 mins to run.
bash extract_spliceai_scores.sh # Extract pre-computed SpliceAI scores.
python3 spliceai_stats.py # Get summary stats for the SpliceAI scores

## MAPS
python3 unaffected_parents.py # Identify unaffected parents
python3 positions_for_bcftools.py # Get tab-separated positions file
bash extract_parent_snvs_wrapper.sh # Run the "wrapper" script
python3 tidy_branch_and_coding_snvs.py # Collate branchpoint and coding SNVs
python3 MAPS.py # Calculate MAPS for branchpoint positions

## Position-weight matrices
python3 pwm.py
