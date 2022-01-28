"""
This script converts earlier outputs to .bed format.
This is so that we can determine the reference allele and sequence context at
each position using bedtools.
"""

# Import the relevant modules
import numpy as np
import pandas as pd

def main(in_path, out_path):
    """ Runs all the functions in this script
    """
    df = load_data(in_path)\
        .pipe(convert_to_bed)

    df.to_csv(out_path, sep="\t", index=False, header=False)

    return df

def load_data(data):
    """ Read the coding / near-splice positions into memory
    """

    print(f"Loading {data}")

    df = pd.read_csv(
        data,
        sep="\t",
        usecols=["chrom", "pos", "strand"],
        dtype = {"chrom":"category", "strand":"category"}
        )

    return df

def convert_to_bed(df):
    """ Convert to .bed format. We are interested in sequence context too, so
    will include the immediately adjacent positions.
    """

    print("Converting to bed.")

    df["start"] = df["pos"] - 2
    df["end"] = df["pos"] + 1
    df["name"] = "."
    df["score"] = "."
    df = df[["chrom", "start", "end", "name", "score", "strand"]]\
        .sort_values(by=["chrom","start"])

    return df

splicing = ["../outputs/near_splice_positions.tsv", "../outputs/near_splice.bed"]
coding = ["../outputs/coding_positions.tsv", "../outputs/coding.bed"]

if __name__ == "__main__":
    for x in [splicing, coding]:
        main(*x)
