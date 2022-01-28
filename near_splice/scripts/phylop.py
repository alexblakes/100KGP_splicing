"""
This script annotates all near-splice positions with phyloP scores
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import pyBigWig

def main():
    """ Run all functions in this script
    """
    data = "../outputs/near_splice_positions.tsv"

    df = load_data(data)\
        .pipe(get_phylop_scores)\
        .pipe(stats)

def load_data(data):
    """ Read the near-splice positions into memory
    """
    df = pd.read_csv(
        data,
        sep="\t"
        )

    return df

def phylop(row):
    """ Get the phyloP score for a given variant.
    """
    return pbw.values(row["chrom"], row["pos"]-1, row["pos"])[0]

def get_phylop_scores(df):
    """ Get phyloP scores for every position of interest
    """
    print("Extracting phyloP scores")

    df["phylop"] = df.apply(phylop, axis=1)

    return df

def stats(df):
    """ Generate summary statistics and save to .tsv
    """
    grouped = df.groupby(["region","site"])
    stats = grouped.phylop.agg(["count", "mean", "std", "sem"])

    stats["ci_upper"] = stats["mean"] + (stats["sem"] * 1.96)
    stats["ci_lower"] = stats["mean"] - (stats["sem"] * 1.96)

    stats.reset_index()\
        .to_csv("../stats/near_splice_phylop.tsv", sep="\t", index=False)

bw = "/public_data_resources/phylop100way/hg38.phyloP100way.bw"
pbw = pyBigWig.open(bw)

if __name__ == "__main__":
    main()
