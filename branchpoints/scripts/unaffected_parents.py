"""
This script gets sample names (Plate Key IDs) for unaffected parents in the rare
diseases cohort.

Only participants who have consented to use of their data, and who are
represented in the GEL aggV2 VCF are included.
"""

# Import relevant modules
import numpy as np
import pandas as pd

# Collate and merge the relevant data
platekeys = "../data/rare_disease_analysis_2021-10-22_12-52-21.tsv"
affection_status = "../data/rare_diseases_pedigree_member_2021-10-22_12-48-13.tsv"
consented_samples="/gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/docs/main_programme_v12_samples.txt"

samples = pd.read_csv(platekeys, sep="\t")
unaff_participants = pd.read_csv(affection_status, sep="\t")
consented_samples = pd.read_csv(consented_samples, sep="\t", header=None).squeeze()

df = samples.merge(unaff_participants)

# 26660 unaffected parents are represented in the AggV2 VCF:
df = df[df["Affection Status"] == "Unaffected"]\
    .drop_duplicates(subset="Plate Key")\
    .loc[df["Plate Key"].isin(consented_samples)]

# Write to output
df["Plate Key"].to_csv(
    "../outputs/unaffected_RD_parents.tsv",
    sep="\t",
    index=False,
    header=None)
