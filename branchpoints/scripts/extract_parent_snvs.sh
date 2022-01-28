#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_machine_learning
#BSUB -J "MAPS_snvs_[1-1371]"
#BSUB -R "rusage[mem=500]"
#BSUB -M 1000
#BSUB -e ../outputs/MAPS_snvs/hpc_err/%J_%I_err.txt
#BSUB -o ../outputs/MAPS_snvs/hpc_out/%J_%I_out.txt

# This scripts identifies SNVs in 26,660 unaffected parents which fall within
# 5bp of a putative splicing branchpoint.

module load bio/BCFtools/1.9-foss-2019b

function extract_snvs {
  bcftools view \
    -v snps \
    -S ${unaffected_parent_samples} \
    -O u \
    -R "../outputs/branch_regions_for_bcftools.bed" \
    ${agg_vcf_chunk} | \
  bcftools query \
    -i "FILTER='PASS' & GT='alt'" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AN\t%INFO/AC\t[%SAMPLE,%GT,%DP,%GQ,%AD;]\n' # | \
}

# Get a VCF chunk, depending on the job index
agg_vcf_chunks="../data/aggV2_vcf_file_paths.tsv"
agg_vcf_chunk=`sed -n "${LSB_JOBINDEX}p" ${agg_vcf_chunks} | awk '{print $(NF-1)}'`

unaffected_parent_samples="../outputs/unaffected_RD_parents.tsv"

extract_snvs > "../outputs/MAPS_snvs/out/branch_snvs_${LSB_JOBINDEX}"
wait
echo "All done."
