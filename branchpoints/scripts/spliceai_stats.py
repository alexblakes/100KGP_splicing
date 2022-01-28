"""
This script generates summary statistics for the SpliceAI scores of
near-splice variants.
"""

# Load the relevant modules
import numpy as np
import pandas as pd
import os

def load_spliceai_data(path):
    """ Read a SpliceAI vcf to memory.
    """
    df = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        header=None,
        names=["chrom","pos","id","ref","alt","qual","filter","info"],
        usecols=["chrom","pos","ref","alt","info"]
        )

    return df

def tidy_spliceai(df):
    """ Tidy the SpliceAI vcf, and get the near-splice annotations for each
    position.
    """
    # Extract delta score (ds) and delta position (dp)
    scores = [x.split(";")[0] for x in df.iloc[:,-1].copy()]
    ds =  pd.DataFrame([x.split("|")[2:6] for x in scores]).astype(float)
    dp =  pd.DataFrame([x.split("|")[6:] for x in scores]).astype(int)

    # Calculate the probability of any splicing impact
    ds_any = (1-(1-ds).product(axis=1)).rename("DS_any")

    # Tidy data
    df = pd.concat([df, ds, dp, ds_any], axis=1)\
        .drop("info", axis=1)

    df.columns = ["chrom", "pos", "ref", "alt", "DS_AG", "DS_AL", "DS_DG",\
        "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL", "DS_any"]

    df["chrom"] = "chr" + df["chrom"].astype(str) # For merging, next.

    return df

def merge_near_splice_positions(df):
    """ Retreive branchpoint annotations for each variant
    """
    branch = pd.read_csv(
        "../outputs/branchpoints.bed",
        sep="\t",
        header=None,
        names=["chrom","start","end","site","branch_score","strand"]
        )
    branch["pos"] = branch.end - 1
    branch = branch[["chrom","pos","site","branch_score"]]

    df = branch.merge(df)

    return df

def spliceai_stats(df, type):
    """ Calculate summary statistics of SpliceAI scores at near-splice positions
    """
    # Get the greatest splicing impact of any SNV at a given position.
    pos_group = df.groupby(["site","chrom","pos","ref"])
    ds_any_max = pos_group.DS_any.max()

    # Calculate summary stats by near-splice position
    site_group = ds_any_max.groupby(level=["site"])
    stats = site_group.agg(["count","mean","std","sem"]).reset_index()
    stats["ci_upper"] = stats["mean"] + (stats["sem"] * 1.96)
    stats["ci_lower"] = stats["mean"] - (stats["sem"] * 1.96)

    stats["type"] = type

    return stats

if __name__ == "__main__":
    in_path = "../outputs/branch_spliceai_scores.vcf"
    output = "../outputs/branch_spliceai_scores.tsv"
    stats_out = "../stats/branch_spliceai.tsv"

    if os.path.exists(output) & os.path.exists(stats_out):
        df = pd.read_csv(output, sep="\t")

    else:
        df = load_spliceai_data(in_path)\
            .pipe(tidy_spliceai)\
            .pipe(merge_near_splice_positions)

        df.to_csv(output, sep="\t", index=False)

    df85 = df[df.branch_score >=0.85].copy()

    stats = pd.concat([spliceai_stats(df, "All"), spliceai_stats(df85, ">=0.85")])

    stats.to_csv(stats_out, sep="\t", index=False)
