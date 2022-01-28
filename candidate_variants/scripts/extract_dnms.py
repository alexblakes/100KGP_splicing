""" This script identifies "high-confidence" DNMs overlapping the near-splice
and branchpoint positions of interest. """

import numpy as np
import pandas as pd

def get_dnms():
    """ Read the "stringent" GEL DNMs into memory.
    """
    dnms = (pd.read_csv(
        "../data/denovo_flagged_variants_2021-11-18_16-55-59.tsv",
        sep="\t",
        usecols=["Trio Id", "Assembly","Chrom","Position","Reference","Alternate"],
        )
        .rename(columns={
            "Chrom":"chrom",
            "Position":"pos",
            "Reference":"ref",
            "Alternate":"alt",
            "Trio Id":"trio_id" ###
            })
    )
    dnms = (dnms[
        (dnms.Assembly == "GRCh38") &
        (dnms.ref.str.len() == 1) &
        (dnms.alt.str.len() == 1) &
        (~dnms.chrom.isin(["X", "Y"]))
        ]
        .copy()
        .drop("Assembly", axis=1)
    )
    dnms["chrom"] = "chr" + dnms.chrom.astype(str)

    return dnms

def get_dnm_cohort():
    """ Get participant information for the GEL DNMs.
    """
    chrt = pd.read_csv(
        "../data/denovo_cohort_information_2021-11-18_16-57-09.tsv",
        sep="\t",
        usecols=["Trio Id", "Plate Key", "Participant Id", "Member",
            "Affection Status", "Assembly"],
        )
    chrt.columns=(["trio_id", "platekey", "participant_id", "member",
        "affection_status", "assembly"])
    chrt = (chrt[
        (chrt.member=="Offspring") &
        (chrt.assembly=="GRCh38")
        ]
        .drop(["member","assembly"], axis=1)
    )
    return chrt

def get_near_splice_positions():
    """ Get the near-splice positions
    """
    ns = pd.read_csv(
    "/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/near_splice/outputs/near_splice_positions.tsv",
    sep="\t",
    )
    return ns

def get_branch_positions():
    """ Get the branchpoint positions
    """
    branch = pd.read_csv(
        "/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/branchpoints/outputs/branchpoints.bed",
        sep="\t",
        header=None,
        names=["chrom","start","end","site","score","strand"]
        )
    branch["pos"] = branch.end - 1
    branch["region"] = "branch"
    branch = branch[["chrom","pos","strand","region","site"]]

    return branch

def to_vcf(df):
    df["id"] = "."
    df["qual"] = "."
    df["filter"] = "."
    df["info"] = "."

    df = (df[["chrom","pos","id","ref","alt","qual","filter","info"]]
        .drop_duplicates(["chrom","pos","ref","alt"])
        .sort_values(by=["chrom","pos"])
    )
    df_no_chr = df.copy()
    df_no_chr["chrom"] = df_no_chr.chrom.str.slice(3) # Remove "chr" prefix

    header = open("../data/vcf_38_header.txt").read()

    with open("../outputs/near_splice_and_branch_dnms_chr.vcf", "w") as output:
        output.write(header)
        df.to_csv(output, index=False, sep="\t", header=False)

    with open("../outputs/near_splice_and_branch_dnms_no_chr_prefix.vcf", "w") as output:
        output.write(header)
        df_no_chr.to_csv(output, index=False, sep="\t", header=False)

if __name__ == "__main__":
    # Find DNMs overlapping branchpoint / near-splice positions.
    dnms = get_dnms().merge(get_dnm_cohort())
    branch = get_branch_positions().merge(dnms)
    ns = get_near_splice_positions().merge(dnms)

    df = (pd.concat([ns, branch])
        .drop_duplicates(["chrom","pos","participant_id"])
    )   # 251 positions seen in both branchpoints and near-splice regions.
        # Keep only near-splice (splice acceptor) annotation for these.

    df.to_csv(
        "../outputs/near_splice_and_branch_dnms.tsv",
        sep="\t",
        index=False,
    )
    to_vcf(df) # Save as VCF
