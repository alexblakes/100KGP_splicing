"""
This script calculates chi-squared p-values for the near-splice and
branchpoint MAPS statistics.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportions_chisquare

def get_branch():
    """ Get the MAPS scores for branchpoint variants.
    """
    br = (pd.read_csv("../stats/branch_maps_output.tsv", sep="\t")
        .rename(columns={"csq":"site"})
        .query("~site.isin(['Synonymous', 'Missense', 'Nonsense'])")
    )
    return br

def get_near_splice():
    """ Get MAPS scores for near-splice variants.
    """
    ns = (pd.read_csv("/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/near_splice/stats/maps_output.tsv",
        sep="\t"
        )
        .drop("xlab", axis=1)
        .query("region != 'Coding'")
        .assign(set = "all")
    )
    return ns

def chi_sq_p(row):
    """ Calculate chi-squared p-values for the MAPS score of each position
    against the MAPS score for synonymous variants.
    """
    syn_singletons = 481163
    syn_alleles = 914195
    count = pd.Series([row.ns_norm, syn_singletons])
    nobs = pd.Series([row.n_alleles, syn_alleles])

    p = proportions_chisquare(count = count, nobs = nobs)[1]

    return p

if __name__ == "__main__":
    df = pd.concat([get_branch(), get_near_splice()]).reset_index(drop=True)
    df["maps_norm"] = df.maps + 0.526344 #Synonymous ps_raw
    df["ns_norm"] = (df.maps_norm * df.n_alleles).apply(lambda x: int(round(x)))
    df["chi_sq_p"] = df.apply(chi_sq_p, axis=1)
    df.to_csv("../stats/maps_chi_square_test.tsv", sep="\t", index=False)
