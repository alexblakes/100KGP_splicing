"""
This script annotates DNMs with exit questionnaire, phenotype, and tiering data.
"""

import numpy as np
import pandas as pd

def get_dnms():
    """ Get the DNMs of interest.
    """
    dnms = pd.read_csv(
        "../outputs/dnms_in_g2p_genes.tsv",
        sep="\t",
        dtype={"pos":"Int64"}
        )
    return dnms

def get_eq(dnms):
    """Get exit questionnaire data for each participant"""
    eq = (pd.read_csv(
        "../data/gmc_exit_questionnaire_2021-11-18_17-01-01.tsv",
        sep="\t",
        header=0,
        usecols=["Participant Id", "Case Solved Family", "Acmg Classification",
            "Assembly", "Chromosome", "Position", "Reference", "Alternate"],
        dtype={"Position":"Int64"}
        )
        .set_axis(["participant_id", "case_solved", "acmg", "assembly",
                    "chrom", "pos", "ref", "alt"], axis=1)
        .drop(23952) # An old instance of a now-solved case
    )
    eq = eq[eq.assembly!="GRCh37"]
    eq_na = eq[eq.chrom.isna()].copy()
    eq_chr = eq[~eq.chrom.isna()].copy()
    eq_chr.loc[~eq_chr.chrom.str.startswith("chr"), "chrom"] = "chr" + eq.chrom.astype(str)
    eq = pd.concat([eq_na, eq_chr])

    exact_match = dnms.merge(eq).drop_duplicates()

    df = (dnms[~dnms.participant_id.isin(exact_match.participant_id)]
        .merge(
            eq.drop("acmg", axis=1),
            on="participant_id",
            how="left",
            suffixes=[None,"_y"]
            )
        .loc[:,:"assembly"]
        .append(exact_match)
        .sort_values(by="acmg")
        .drop_duplicates()
        .drop_duplicates(["chrom","pos","ref","alt","participant_id"])
    )
    #
    return df

def get_hpo():
    """Get HPO terms."""
    hpo = (pd.read_csv(
        "../data/rare_diseases_participant_phen_2021-11-18_17-00-10.tsv",
        sep = "\t",
        usecols = ["Participant Id", "Hpo Term"]
        )
        .set_axis(["participant_id", "hpo"], axis=1)
        .groupby("participant_id")
            ["hpo"]
            .apply(",".join)
            .reset_index()
    )
    return hpo

def get_tiers(df):
    """Get tiering data."""
    tiers = (pd.read_csv(
        "../data/tiering_data_2021-11-18_17-21-18.tsv",
        sep="\t",
        usecols=["Participant Id","Assembly","Chromosome","Position","Reference","Alternate","Tier"]
        )
        .set_axis(["participant_id","assembly","chrom","pos","ref","alt","tier"], axis=1)
        .dropna()
    )
    tiers["chrom"] = "chr" + tiers.chrom.astype(str)
    tiers["tier"] = tiers.tier.str.slice(4).astype(int)
    tiers["tier"] = (tiers
        .groupby(["participant_id","chrom","pos","ref","alt"])
        .tier
        .transform(min) # If conflicting, give the "highest" tier.
    )
    tiers = (tiers[tiers.assembly=="GRCh38"]
        .drop("assembly", axis=1)
        .drop_duplicates()
    )
    max_tier = (tiers
        .groupby("participant_id")
        .tier.min()
        .rename("max_tier")
        .reset_index()
    )
    return df.merge(tiers, how="left").merge(max_tier, how="left")

if __name__ == "__main__":
    # Annotate DNMs with exit questionnaire, phenotype, and tiering data.
    df = (get_dnms()
        .pipe(get_eq)
        .merge(get_hpo())
        .pipe(get_tiers)
    )
    df["region"] = pd.Categorical(df.region, ["branch","acceptor","donor"], ordered=True)
    df = df.sort_values(by=["region","site","DS_any"])

    # Decide whether our prioritised variant was causative.
    # The order of these lines matters!
    df["eq_outcome"] = np.nan

    m1 = (df.case_solved=="yes")
    m2 = (df.case_solved=="yes") & (df.max_tier.isna())
    m3 = ((df.case_solved=="yes") &
          (df.acmg.isin(["pathogenic_variant", "likely_pathogenic_variant"])))
    m4 = (df.case_solved.isin(["no","partially"]))
    m5 = (df.case_solved.isna()) | (df.case_solved == "unknown")

    df.loc[m1, "eq_outcome"] = "Case solved, different variant"
    df.loc[m2, "eq_outcome"] = "Case solved, variant unknown"
    df.loc[m3, "eq_outcome"] = "Case solved, same variant"
    df.loc[m4, "eq_outcome"] = "Case not solved"
    df.loc[m5, "eq_outcome"] = "No data"

    df.to_csv("../outputs/candidate_dnms.tsv", sep="\t", index=False)
