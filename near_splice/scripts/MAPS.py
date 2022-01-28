"""
This script calculates the mutability-adjusted proportion of singletons, MAPS.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

df = pd.read_csv("../outputs/unaff_parents_allele_counts.tsv", sep="\t")

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

df = df[["csq","ac","context","alt_context","mu_snp"]]
df_g = df.groupby("csq")

# Construct a linear model describing how the proportion of singletons varies
# with mutability.
# Build the model on synonymous variants:
syn = df_g.get_group("synonymous_variant")
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

# Mean mutability across each group, weighted by number of alleles
n_singletons = df_g.apply(lambda x: (x.ac==1).sum()).rename("n_singletons")
n_alleles = df_g.ac.count().rename("n_alleles")
ps_raw = (n_singletons/n_alleles).rename("ps_raw") # Unadjusted ps
se = np.sqrt((ps_raw * (1 - ps_raw))/n_alleles).rename("se")
mu = df_g.mu_snp.mean().rename("mu_snp") # Weighted mutability
ps_pred = res.predict(exog=mu).rename("ps_pred") # Expected proportion of singletons
maps = (ps_raw - ps_pred).round(6).rename("maps") # Calculate MAPS
ci_upper = (maps + 1.96 * se).rename("ci_upper")
ci_lower = (maps - 1.96 * se).rename("ci_lower")

# Reformat the data index for easier plotting later
labels = maps.index.rename("csq").to_frame()
labels.loc[:, "region"] = "Coding"
labels.loc[labels.csq.str.contains("acceptor"), "region"] = "Acceptor"
labels.loc[labels.csq.str.contains("donor"), "region"] = "Donor"
labels.loc[labels.region.isin(["Acceptor", "Donor"]), "site"] = labels.csq.str.split("_").str[1]

labels["xlab"] = labels.site
labels.loc[labels.region.isin(["Acceptor", "Donor"]), "xlab"] = labels.region.str.slice(0,1) + " " + labels.site
labels.loc[labels.csq.str.contains("s"), "xlab"] = labels.csq.str.split("_").str[0].str.capitalize()
labels["xlab"] = labels.xlab.str.replace("Stop", "Nonsense")

maps_out = pd.concat([labels.loc[:, "region":], n_singletons, n_alleles, ps_raw,
    se, mu, ps_pred, maps, ci_upper, ci_lower], axis=1).reset_index(drop=True)

maps_out.to_csv("../stats/maps_output.tsv", sep="\t", index=False)
