#!/usr/bin/env bash

## This script extracts SpliceAI scores for each SNV.

#BSUB -q short
#BSUB -P re_gecip_machine_learning
#BSUB -J "extract_spliceai_scores"
#BSUB -R "rusage[mem=4000]"
#BSUB -M 4000
#BSUB -e ../outputs/%J_err.txt
#BSUB -o ../outputs/%J_out.txt

# Load BCFtools
module load bio/BCFtools/1.9-foss-2019b

genome="/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.masked.snv.hg38.vcf.gz"
no_chr="../outputs/near_splice_and_branch_dnms_no_chr_prefix.vcf"

bgzip < ${no_chr} > "${no_chr}.gz"
tabix -p vcf -f "${no_chr}.gz"

bcftools isec \
	-w1 \
	-n=2 \
	${genome} \
	"${no_chr}.gz" \
	-o ../outputs/near_splice_and_branch_dnms_spliceai_scores.vcf

################################################################################
