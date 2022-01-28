"""
This script reformats the unaffected parental coding SNVs to VCF in order to
run VEP
"""

import numpy as np
import pandas as pd

df = pd.read_csv(
    "../outputs/unaff_parents_coding_snvs.tsv",
    sep="\t",
    header=None,
    names=["chrom","pos","ref","alt","filter","an","ac","sample"],
    usecols=["chrom","pos","ref","alt"]
    )

# convert to .vcf format
df["id"] = "."
df["qual"] = "."
df["filter"] = "."
df["info"] = "."

df = df[["chrom","pos","id","ref","alt","qual","filter","info"]]\
    .sort_values(by=["chrom","pos"])

for i, frame in enumerate(np.array_split(df, 22)):
    frame.to_csv(
        f"../outputs/vep/in/coding_vep_input_{i+1}.vcf",
        sep="\t",
        index=False,
        header=False
        )
