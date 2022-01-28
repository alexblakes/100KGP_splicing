"""
This script annotates all of the coding positions in GENCODE.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import near_splice_sites as nss

def main():
    """ Run all functions in this script.
    Write the output to .tsv and return the dataframe.
    These functions are imported from near_splice_sites.py:
        load_data
        autosomes_only
    """
    data = "../outputs/coding_exons.tsv"
    df = nss.load_data(data)\
        .pipe(nss.autosomes_only)\
        .pipe(coding_sites)

    output = "../outputs/coding_positions.tsv"

    df.to_csv(output, sep="\t", index=False)

    return df

def coding_sites(df):
    """ Annotate all of the coding positions in GENCODE.
    """
    # Identify the coding positions
    coding_positions = lambda x: list(range(x.start, x.end + 1))

    df["pos"] = df.apply(coding_positions, axis=1)

    # Subset to the relevant columns
    # Explode so that each position is represented in one row
    # Keep only one representation for each position
    df = df[["chrom","pos","strand"]]\
        .explode("pos", ignore_index=True)\
        .drop_duplicates(subset=["chrom","pos"])

    print(f"There are {len(df)} distinct coding positions in GENCODE.")

    return df

if __name__ == "__main__":
    main()
