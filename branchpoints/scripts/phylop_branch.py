"""
This script annotates all near-splice positions with phyloP scores
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import pyBigWig

def load_data(data):
    """ Read the branchpoint positions into memory
    """
    df = pd.read_csv(
        data,
        sep="\t",
        header=None,
        names=["chrom","start","end","site","score","strand"]
        )

    return df

def phylop(row):
    """ Get the phyloP score for a given position.
    """
    return pbw.values(row["chrom"], row["start"]+1, row["end"]-1)[0]

def get_phylop_scores(df):
    """ Get phyloP scores for every position of interest
    """
    print("Extracting phyloP scores")

    df["phylop"] = df.apply(phylop, axis=1)

    return df

def stats(df, type):
    """ Generate summary statistics and save to .tsv
    """
    grouped = df.groupby(["site"])
    stats = grouped.phylop.agg(["count", "mean", "std", "sem"]).reset_index()

    stats["ci_upper"] = stats["mean"] + (stats["sem"] * 1.96)
    stats["ci_lower"] = stats["mean"] - (stats["sem"] * 1.96)

    stats["type"] = type

    return stats

if __name__ == "__main__":
    bw = "/public_data_resources/phylop100way/hg38.phyloP100way.bw"
    pbw = pyBigWig.open(bw)
    data = "../outputs/branchpoints.bed"

    df = load_data(data)\
        .pipe(get_phylop_scores)

    df85 = df[df.score >= 0.85].copy()

    phylop_stats = pd.concat([stats(df, "All"), stats(df85, ">=0.85")])

    phylop_stats.to_csv(
        "../stats/branchpoint_phylop.tsv",
        sep="\t",
        index=False
        )
