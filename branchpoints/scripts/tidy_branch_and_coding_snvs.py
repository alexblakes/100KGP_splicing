"""
This script collates the branchpoint and coding SNVs.
It gives one definitive consequence to each SNV.
Variants are labelled as "nonsense", "missense", "synonymous", or a given a
branchpoint annotation.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import os

def get_allele_counts(path):
    """ Read the allele counts for SNVs in the unaffected parents.
    """

    names = ["chrom", "pos", "ref", "alt", "filter", "an", "ac", "sample"]

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=names,
        usecols=["chrom", "pos", "ref", "alt", "ac"]
        )\
        .drop_duplicates(subset=["chrom", "pos", "ref", "alt"], keep=False)

    return df

def get_branch_csqs(path):
    """ Get the branchpoint annotations for each branchpoint variant.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom","start","end","csq","score","strand"],
        usecols=["chrom","end","csq","score"],
        )

    df["region"] = "Branchpoint"
    df["pos"] = df.end - 1
    df = df.drop("end", axis=1)

    return df

def get_branch_contexts(path):
    """ Get the sequence context for each genomic position.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["span","context"]
        )

    # List comprehensions are faster than string methods here:
    a = pd.DataFrame([x.split(":") for x in df.span])
    b = pd.DataFrame([x.split("-") for x in a[1]]).astype(int)

    chrom = a[0].rename("chrom")
    pos = (b[0] + 2).rename("pos")

    df = pd.concat([chrom, pos, df.context], axis=1)\
        .drop_duplicates()\
        .drop_duplicates(subset=["chrom","pos"], keep=False)

    return df

def get_coding_csqs(path):
    df = pd.read_csv(
        path,
        sep="\t",
        )
    coding_csqs = ['missense_variant', 'synonymous_variant', 'stop_gained']
    df = df[df.csq.isin(coding_csqs)]
    df["region"] = "Coding"

    return df

def combine_annotations():
    """ Run the functions above to retrieve the annotations.
    Merge these annotations together.
    """
    branch_snvs = get_allele_counts("../outputs/unaff_parents_branch_snvs.tsv")
    branch_csqs = get_branch_csqs("../outputs/branchpoints.bed")
    branch_ctxt = get_branch_contexts("../outputs/branch_contexts.tsv")

    # Merge the data, to retrieve consequences and context annotations for each
    # allele
    branch = branch_snvs.merge(branch_csqs).merge(branch_ctxt)

    # Concatenate the coding and near-splice variants
    coding = get_coding_csqs("/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/near_splice/outputs/unaff_parents_allele_counts.tsv")

    df = (pd.concat([coding, branch], ignore_index=True)
        .drop_duplicates(["chrom","pos","ref","alt"], keep=False)
    )
    return df

if __name__ == '__main__':
    output = "../outputs/unaff_parents_allele_counts.tsv"
    if os.path.exists(output):
        df = pd.read_csv(
            output,
            sep="\t",
            dtype={"csq":"object"}
            )
    else:
        df = combine_annotations()
        df.to_csv(output, sep="\t", index=False)
