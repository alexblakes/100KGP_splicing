"""
This script creates a tab-delimited .bed file of positions within 5bp of a
branchpoint, for use by bcftools as part of the MAPS SNV filtering process.
"""

# Import relevant modules
import numpy as np
import pandas as pd

def positions_for_bcftools():
    """ Reformat the positions file
    """
    df = pd.read_csv(
        f"../outputs/branchpoints.bed",
        sep="\t",
        header=None,
        names=["chrom","start","end","site","score","strand"]
        )

    df0 = df[df.site==0].copy()
    df0["start"] = df0.start - 4
    df0["end"] = df0.end + 4
    df0 = (df0[["chrom","start","end"]]
        .sort_values(by=["chrom","start"])
    )

    df0.to_csv(
        f"../outputs/branch_regions_for_bcftools.bed",
        sep="\t",
        index=False,
        header=False
        )

if __name__ == '__main__':
    positions_for_bcftools()
