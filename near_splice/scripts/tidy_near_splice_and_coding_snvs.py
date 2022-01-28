"""
This script collates the near-splice and coding SNVs.
It gives one definitive consequence to each SNV.
Variants are labelled as "nonsense", "missense", "synonymous", or a given a
near-splice annotation.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import os

def read_allele_counts(region):
    """ Read the allele counts for SNVs in the unaffected parents.
    """

    names = ["chrom", "pos", "ref", "alt", "filter", "an", "ac", "sample"]

    df = pd.read_csv(
        f"../outputs/unaff_parents_{region}_snvs.tsv",
        sep="\t",
        header=None,
        names=names,
        usecols=["chrom", "pos", "ref", "alt", "ac"]
        )\
        .drop_duplicates(subset=["chrom", "pos", "ref", "alt"], keep=False)

    return df

def read_coding_csqs():
    """ Get the VEP consequence annotations for each coding variant
    """

    names = (["variant", "a", "b", "ensg", "enst", "d", "csq", "e", "f", "g", "h", "i",
        "j", "k", "l", "m", "n"])

    df = pd.read_csv(
        "../outputs/unaff_parents_coding_snvs_vep_out.tsv",
        sep="\t",
        comment="#",
        header=None,
        names=names,
        usecols=["variant", "csq"]
        )

    # Reformat the data
    var = df.variant.str.split("_", expand=True)

    chrom = var[0].rename("chrom")
    pos = var[1].astype(int).rename("pos")
    ref = var[2].str.split("/").str[0].rename("ref")
    alt = var[2].str.split("/").str[1].rename("alt")

    # Keep only SNVs with the following selected consequences
    consequences = ['synonymous_variant', 'missense_variant', 'stop_gained']
    csq = df.csq.str.split(",").explode().rename("csq")
    csq = csq[csq.isin(consequences)]

    df = pd.concat([chrom, pos, ref, alt, csq], axis=1, join="inner")\
        .drop_duplicates() # 3,039,395 variants (78 duplicated variants)

    return df

def read_near_splice_csqs():
    """ Get the consequence annotations for each near-splice variant.
    """

    df = pd.read_csv("../outputs/near_splice_positions.tsv", sep="\t")
    df["csq"] = df.region + "_" + df.site.astype(str)
    df = df[["chrom", "pos", "csq"]]

    return df

def read_contexts(region):
    """ Get the sequence context for each genomic position.
    """

    df = pd.read_csv(
        f"../outputs/{region}_contexts.tsv",
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

def combine_annotations():
    """ Run the functions above to retrieve the annotations.
    Merge these annotations together.
    """
    cd_snvs = read_allele_counts("coding")
    ns_snvs = read_allele_counts("near_splice")

    cd_csqs = read_coding_csqs()
    ns_csqs = read_near_splice_csqs()

    cd_ctxt = read_contexts("coding")
    ns_ctxt = read_contexts("near_splice")

    # Merge the data, to retrieve consequences and context annotations for each
    # allele
    cd = cd_snvs.merge(cd_csqs).merge(cd_ctxt)
    ns = ns_snvs.merge(ns_csqs).merge(ns_ctxt)

    # Concatenate the coding and near-splice variants
    df = pd.concat([cd, ns])

    return df

def reduce_annotations(df):
    """ Reduce the data so that each allele has one definitive annotation.
    Where a synonymous variant is found in a near-splice position, label it as
    near-splice.
    Where a missense/nonsense variant is found in a near-splice position, label
    it as missense/nonsense.
    """
    # Keep those variants which are not (duplicated and synonymous) and those
    # variants which are not (duplicated and (not missense and not stop_gained)
    mask1 = df.duplicated(subset=["chrom","pos","ref","alt"], keep=False)
    mask2 = df.csq == "synonymous_variant"
    df = df[~(mask1 & mask2)]

    mask3 = df.duplicated(subset=["chrom","pos","ref","alt"], keep=False)
    mask4 = df.csq != "missense_variant"
    mask5 = df.csq != "stop_gained"
    df = df[~(mask3 & (mask4 & mask5))]

    return df

if __name__ == '__main__':
    output = "../outputs/unaff_parents_allele_counts.tsv"

    if os.path.exists(output):
        df = pd.read_csv(output, sep="\t")
    else:
        df = combine_annotations()\
            .pipe(reduce_annotations)
        df.to_csv(output, sep="\t", index=False)
