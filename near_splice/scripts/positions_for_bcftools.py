"""
This script gets a tab-delimited positions file for near-splice and coding
positions, for use by bcftools as part of the MAPS SNV filtering exercise.
"""

# Import relevant modules
import numpy as np
import pandas as pd

def positions_for_bcftools(region):
    """ Reformat the positions file
    """
    df = pd.read_csv(
        f"../outputs/{region}_positions.tsv",
        sep="\t"
        )

    df = df[["chrom","pos"]].sort_values(by=["chrom","pos"])

    df.to_csv(
        f"../outputs/{region}_positions_for_bcftools.tsv",
        sep="\t",
        index=False,
        header=False
        )

ns, cd = "near_splice", "coding"
for region in [ns, cd]: positions_for_bcftools(region)
