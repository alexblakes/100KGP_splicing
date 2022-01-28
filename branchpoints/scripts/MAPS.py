"""
This script calculates the mutability-adjusted proportion of singletons, MAPS.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv(
    "../outputs/unaff_parents_allele_counts.tsv",
    sep="\t",
    dtype={"csq":"object"}
    )

a = [x[0] for x in df.context]
b = [x[2] for x in df.context]
df["alt_context"] = a + df.alt + b

mu_snp = pd.read_csv("../data/forSanger_1KG_mutation_rate_table.txt", sep=" ")

df = df.merge(
    mu_snp,
    left_on=["context","alt_context"],
    right_on=["from","to"]
    )\
    .drop(["from","to"], axis=1)

df = df[["csq","ac","context","alt_context","mu_snp","region","score"]]
df["set"] = "all"
df85 = df[df.score >= 0.85].copy()
df85["set"] = "85"

# Construct a linear model describing how the proportion of singletons varies
# with mutability.
# Build the model on synonymous variants:
syn = df[df.csq=="synonymous_variant"]
syn_g = syn.groupby(["context", "alt_context"])

# Get the relevant summary statistics
n_singletons = syn_g.ac.apply(lambda x: (x==1).sum()).rename("n_singleton")
n_alleles = syn_g.ac.count().rename("n_alleles")
ps_raw = (n_singletons/n_alleles).rename("ps_raw")
mu_snp = syn_g.mu_snp.max().rename("mu_snp")
syn_ps = pd.concat([n_singletons, n_alleles, ps_raw, mu_snp], axis=1)

# Construct the linear model. Weighted least-squares.
maps_model = smf.wls(
    "ps_raw ~ mu_snp",
    data=syn_ps,
    weights=syn_ps.n_alleles
    )

res = maps_model.fit()

def maps_adjust(df):
    df_g = df.groupby("csq")

    # Mean mutability across each group, weighted by number of alleles
    n_singletons = df_g.apply(lambda x: (x.ac==1).sum()).rename("n_singletons")
    n_alleles = df_g.ac.count().rename("n_alleles")
    ps_raw = (n_singletons/n_alleles).rename("ps_raw") # Unadjusted ps
    se = np.sqrt((ps_raw * (1 - ps_raw))/n_alleles).rename("se")
    mu = df_g.mu_snp.mean().rename("mu_snp") # Weighted mutability
    ps_pred = res.predict(mu).rename("ps_pred") # Expected proportion of singletons
    maps = (ps_raw - ps_pred).round(6).rename("maps") # Calculate MAPS
    ci_upper = (maps + 1.96 * se).rename("ci_upper")
    ci_lower = (maps - 1.96 * se).rename("ci_lower")
    region = df_g.region.first()
    _set = df_g.set.first()

    # Reformat the data index for easier plotting later
    maps_out = (pd.concat([
        region, _set, n_singletons, n_alleles, ps_raw, se, mu, ps_pred, maps, ci_upper,
        ci_lower],
        axis=1,
        )
        .reset_index()
    )

    return maps_out

maps_out = (pd.concat([maps_adjust(df), maps_adjust(df85)], ignore_index=True)
    .replace({
        "missense_variant":"Missense",
        "stop_gained":"Nonsense",
        "synonymous_variant":"Synonymous"
        })
)

maps_out["csq"] = pd.Categorical(maps_out.csq, ["Synonymous","Missense","Nonsense","-5","-4","-3","-2","-1","0","1","2","3","4","5"])
maps_out["region"] = pd.Categorical(maps_out.region, ["Coding","Branchpoint"])
maps_out["set"] = pd.Categorical(maps_out.set, ["all","85"])

maps_out = maps_out.sort_values(by=["region","set","csq"])

maps_out.to_csv("../stats/branch_maps_output.tsv", sep="\t", index=False)
