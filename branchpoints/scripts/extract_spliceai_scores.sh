#!/usr/bin/env bash

## This script extracts SpliceAI scores for each SNV.

# Load BCFtools
module load bio/BCFtools/1.9-foss-2019b

exome="/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/OLD/genome_scores-173156985/genome_scores-ds.4fed5ab77f51400ca5c74b9dab27a059/spliceai_scores.exome.hg38.vcf.gz"
genome="/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.masked.snv.hg38.vcf.gz"

chr="../outputs/branch_all_snvs_chr.vcf.gz"
no_chr="../outputs/branch_all_snvs_no_chr_prefix.vcf.gz"

# Intersect our positions of interest with the pre-computed SpliceAI file.
bcftools isec \
	-n=2 \
	-w1 \
	${genome} \
	${no_chr} \
	> ../outputs/branch_spliceai_scores.vcf

###########################################################################################################################
