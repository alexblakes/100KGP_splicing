"""
This script generates summary statistics for the SpliceAI scores of
the GEL DNMs.
"""

# Load the relevant modules
import numpy as np
import pandas as pd

def get_spliceai_data():
    """ Read the SpliceAI output vcf to memory.
    """
    df = pd.read_csv(
        "../outputs/near_splice_and_branch_dnms_spliceai_scores.vcf",
        comment="#",
        sep="\t",
        header=None,
        names=["chrom","pos","id","ref","alt","qual","filter","info"],
        usecols=["chrom","pos","ref","alt","info"],
        dtype={"chrom":str}
        ) \
        .drop_duplicates(["chrom","pos","ref","alt"])

    return df

def tidy_spliceai(df):
    """ Tidy the SpliceAI VCF, splitting the delta scores and getting DS_any and
    DS_max annotations scores.
    """
    # Extract delta score (ds) and delta position (dp)
    scores = [x.split(";")[0] for x in df.iloc[:,-1].copy()]
    ds =  pd.DataFrame([x.split("|")[2:6] for x in scores]).astype(float)
    ds.columns = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]

    ds_any = np.round((1-(1-ds).product(axis=1)), 3) # Probability of any splicing impact
    ds_max = ds.max(axis=1) # Maximum predicated splicing impact
    ds_max_type = ds.idxmax(axis=1) # Type of maximal splicing impact
    ds_diff = np.round((ds_any - ds_max), 3)

    # Tidy data
    df = df.drop("info", axis=1)
    df = pd.concat([df, ds, ds_any, ds_max, ds_max_type, ds_diff], axis=1)
    df.columns = ["chrom", "pos", "ref", "alt", "DS_AG", "DS_AL", "DS_DG",\
        "DS_DL", "DS_any", "DS_max", "DS_max_type", "DS_diff"]
    df = df[["chrom","pos","ref","alt","DS_any","DS_max","DS_max_type"]]

    df["chrom"] = "chr" + df["chrom"].astype(str) # Re-introduce the "chr-" prefix

    return df

if __name__ == "__main__":
    sai = get_spliceai_data()\
        .pipe(tidy_spliceai)\

    dnms = pd.read_csv("../outputs/near_splice_and_branch_dnms.tsv", sep="\t")

    df = dnms.merge(sai, how="left")

    df.to_csv(
        "../outputs/near_splice_and_branch_dnms_spliceai.tsv",
        sep="\t",
        index=False,
        )
