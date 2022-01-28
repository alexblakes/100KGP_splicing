#!/usr/bin/env bash

### This script runs all of the scripts in the diagnostic_SpliceAI_vars
### directory

python3 extract_dnms.py
bsub < vep_snvs.sh
bash extract_spliceai_scores.sh
bwait -w 'ended(extract_spliceai_scores)'
python3 spliceai_stats.py
python3 dnms_in_g2p_genes.py
python3 participant_data.py
