"""
This script finds DNMs which fall within known dominant LoF disease genes.
"""

# Import relevant modules
import numpy as np
import pandas as pd

# Load variant data to memory
def get_dnms():
    """Get the GEL DNMs of interest."""
    dnms = pd.read_csv(
        "../outputs/near_splice_and_branch_dnms_spliceai.tsv",
        sep="\t",
        )
    dnms = (dnms[["chrom", "pos", "ref", "alt", "region", "site", "platekey",
        "participant_id", "affection_status", "DS_any", "DS_max",
        "DS_max_type"]])

    return dnms

def get_vep_output():
    """Get the VEP output."""
    vep = pd.read_csv(
        "../outputs/near_splice_and_branch_dnms_vep.tsv",
        comment="#",
        sep="\t",
        header=None,
        names=["variant","ensg","enst","hgnc_id","csq"],
        na_values="-",
        )

    variant = vep.variant.str.split("_", expand=True)
    vep["chrom"] = variant[0]
    vep["pos"] = variant[1].astype(int)
    vep["ref"] = variant[2].str.slice(0,1)
    vep["alt"] = variant[2].str.slice(2)
    vep["hgnc_id"] = vep["hgnc_id"].str.slice(5).astype(float)

    vep = vep[["chrom","pos","ref","alt","ensg","enst","hgnc_id"]]

    return vep

def read_g2p(path):
    """Get G2P data.
    """
    df = (pd.read_csv(
        path,
        dtype={"hgnc id": float},
        usecols=["gene symbol","DDD category","allelic requirement",
            "mutation consequence","panel","hgnc id"]
        )
        .set_axis(
            labels=["symbol","category","allele_requirement",
                "g2p_mutation_effect","g2p_panel", "hgnc_id"],
            axis=1,
            )
        )

    return df

def get_g2p():
    """Combine and filter the individual G2P annotation files.
    """
    paths = ([
        "../data/DDG2P_27_10_2021.csv",
        "../data/EyeG2P_27_10_2021.csv",
        "../data/SkinG2P_27_10_2021.csv",
        ])
    g2p = pd.concat([read_g2p(p) for p in paths])
    g2p = g2p[g2p.category.isin(["confirmed","probable"])]
    g2p = g2p[g2p.allele_requirement=="monoallelic"]
    g2p = g2p[g2p.g2p_mutation_effect=="loss of function"]
    g2p = g2p.drop_duplicates("hgnc_id")
    return g2p

# Merge the data to find DNMs within known RD genes
if __name__ == "__main__":
    dnms = get_dnms()
    vep = get_vep_output()
    g2p = get_g2p()

    df = dnms.merge(vep)\
        .merge(g2p)

    df.to_csv("../outputs/dnms_in_g2p_genes.tsv", sep="\t", index=False)
