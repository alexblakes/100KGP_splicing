""" This script calculates a position-weight matrix for the branchpoint and
near-splice positions. Only forward stranded positions are used.
"""

import numpy as np
import pandas as pd

def get_branch():
    """ Get reference alleles at branchpoint positions.
    """
    df = pd.read_csv("../outputs/branch_contexts.tsv", sep="\t", header=None)

    # create dummy dataframes for simplicity
    a = df[0].str.split(":", expand=True)
    b = a[1].str.split("-", expand=True)

    # extract the relevant information
    chrom = a[0].rename("chrom")
    start = b[0].astype(int).rename("start")
    end = b[1].astype(int).rename("end")
    pos = (b[1].astype(int) - 1).rename("pos")
    ref = df[1].str.slice(1,2).rename("ref")

    df = (pd.concat([chrom, pos, ref, start, end], axis=1)
        .sort_values(by=["chrom","pos"])
        .assign(region="Branchpoint")
    )
    bp = pd.read_csv(
        "../outputs/branchpoints.bed",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "site", "score", "strand"],
        usecols = ["chrom", "start", "end", "site", "strand"]
        )
    df = (df.merge(bp)
        .query("strand == '+'")
        .loc[:, ["region","site","ref"]]
    )
    return df

def get_near_splice():
    """ Get reference alleles at near-splice positions.
    """
    df = (pd.read_csv(
        "/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/near_splice/outputs/near_splice_spliceai_scores.tsv",
        sep="\t",
        usecols = ["chrom","pos","region","site","ref","strand"]
        )
        .drop_duplicates()
        .query("strand == '+'")
        .drop(["chrom","pos","strand"], axis=1)
        .replace({"donor":"Donor", "acceptor":"Acceptor"})
    )
    return df

df = pd.concat([get_branch(), get_near_splice()])

# Calculate position-weight matrices for these positions.
pwm = (df.groupby(["region","site"])
    .ref.value_counts()
    .unstack(level="ref")
    .fillna(0)
    .reset_index()
    .astype({"A":int, "T":int, "C":int, "G":int})
)
pwm.to_csv("../stats/pwm.tsv", sep="\t", index=False)
