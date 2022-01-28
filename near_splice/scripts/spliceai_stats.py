"""
This script generates summary statistics for the SpliceAI scores of
near-splice variants.
"""

# Load the relevant modules
import numpy as np
import pandas as pd
import os

def main():
    """ Runs all the functions in this script
    """
    in_path = "../outputs/near_splice_spliceai_scores/0000.vcf"
    output = "../outputs/near_splice_spliceai_scores.tsv"
    stats_out = "../stats/near_splice_spliceai.tsv"

    if os.path.exists(output) & os.path.exists(stats_out):
        df = pd.read_csv(output, sep="\t")

    else:
        df = load_spliceai_data(in_path)\
            .pipe(tidy_spliceai)\
            .pipe(merge_near_splice_positions, path=output)\
            .pipe(near_splice_spliceai_stats, path=stats_out)

    return df

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

def merge_near_splice_positions(df, path):
    """ Retreive near-splice annotations for each variant
    """
    ns = pd.read_csv("../outputs/near_splice_positions.tsv", sep="\t")
    df = ns.merge(df)
    df.to_csv(path, sep="\t", index=False) # Save the entire dataframe

    return df

def near_splice_spliceai_stats(df, path):
    """ Calculate summary statistics of SpliceAI scores at near-splice positions
    """
    # Get the greatest splicing impact of any SNV at a given position.
    pos_group = df.groupby(["region","site","chrom","pos","ref"])
    ds_any_max = pos_group.DS_any.max()

    # Calculate summary stats by near-splice position
    site_group = ds_any_max.groupby(level=["region", "site"])
    stats = site_group.agg(["count","mean","std","sem"])
    stats["ci_upper"] = stats["mean"] + (stats["sem"] * 1.96)
    stats["ci_lower"] = stats["mean"] - (stats["sem"] * 1.96)

    # Write to output
    stats = stats.reset_index()
    stats.to_csv(path, sep="\t", index=False)

    return df

if __name__ == "__main__":
    df = main()
