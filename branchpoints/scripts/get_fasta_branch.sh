# This script annotates the sequence context around each putative splicing
# branchpoint.

module load bio/BEDTools/2.27.1-foss-2018b

bedtools getfasta \
-fo ../outputs/branch_contexts.tsv \
-fi "/public_data_resources/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" \
-bed ../outputs/branchpoints.bed \
-tab
