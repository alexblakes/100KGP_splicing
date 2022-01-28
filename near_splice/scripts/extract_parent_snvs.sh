#!/usr/bin/env bash

#BSUB -q medium
#BSUB -P re_gecip_machine_learning
#BSUB -J "MAPS_snvs_[1-1371]"
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000
#BSUB -e ../outputs/MAPS_snvs/hpc_err/%J_%I_err.txt
#BSUB -o ../outputs/MAPS_snvs/hpc_out/%J_%I_out.txt

# This script extracts all SNVs in near-splice and coding regions for 26,660
# unaffected parents in GEL.

module load bio/BCFtools/1.9-foss-2019b

function extract_snvs {
  bcftools view \
    -v snps \
    -S ${unaffected_parent_samples} \
    -O u \
    -T "../outputs/${region}_positions_for_bcftools.tsv" \
    ${agg_vcf_chunk} | \
  bcftools query \
    -i "FILTER='PASS' & GT='alt'" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AN\t%INFO/AC\t[%SAMPLE,%GT,%DP,%GQ,%AD;]\n' # | \
}

# Get a VCF chunk, depending on the job index
agg_vcf_chunks="../data/aggV2_vcf_file_paths.tsv"
agg_vcf_chunk=`sed -n "${LSB_JOBINDEX}p" ${agg_vcf_chunks} | awk '{print $(NF-1)}'`

unaffected_parent_samples="../outputs/unaffected_RD_parents.tsv"

regions=( "near_splice" "coding" )

for region in ${regions[@]}
do
  extract_snvs > "../outputs/MAPS_snvs/out/${region}/${region}_snvs_${LSB_JOBINDEX}" &
done
wait
echo "All done."
