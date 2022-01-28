""" This script combines the phyloP, SpliceAI, and MAPS data for near-splice
and branchpoint positions, so they can be plotted together.
"""

import numpy as np
import pandas as pd

dir = "/re_gecip/machine_learning/AlexBlakes/near_splice/write_up/paper/code/"

ns_phylop = (pd.read_csv(
    dir + "near_splice/stats/near_splice_phylop.tsv",
    sep="\t",
    usecols = ["region","site","mean","sem"]
    )
    .assign(score = "phyloP")
    .assign(subset = "All")
    .replace({"acceptor":"Acceptor", "donor":"Donor"})
)
ns_sai = (pd.read_csv(
    dir + "near_splice/stats/near_splice_spliceai.tsv",
    sep="\t",
    usecols = ["region","site","mean","sem"]
    )
    .assign(score = "SpliceAI")
    .assign(subset = "All")
    .replace({"acceptor":"Acceptor", "donor":"Donor"})
)
ns_maps = (pd.read_csv(
    dir + "near_splice/stats/maps_output.tsv",
    sep="\t",
    usecols = ["region","site","maps","se"]
    )
    .assign(score = "MAPS")
    .assign(subset = "All")
    .rename(columns={"maps":"mean", "se":"sem"})
    .dropna()
    .astype({"site":"int"})
)
b_phylop = (pd.read_csv(
    dir + "branchpoints/stats/branchpoint_phylop.tsv",
    sep="\t",
    usecols = ["site","mean","sem", "type"]
    )
    .assign(region = "Branchpoint")
    .assign(score = "phyloP")
    .rename(columns={"type":"subset"})
    .replace({">=0.85":"85"})
)
b_sai = (pd.read_csv(
    dir + "branchpoints/stats/branch_spliceai.tsv",
    sep="\t",
    usecols = ["site","mean","sem", "type"]
    )
    .assign(region = "Branchpoint")
    .assign(score = "SpliceAI")
    .rename(columns={"type":"subset"})
    .replace({">=0.85":"85"})
)
b_maps = (pd.read_csv(
    dir + "branchpoints/stats/branch_maps_output.tsv",
    sep="\t",
    usecols = ["region","csq","maps","se", "set"]
    )
    .assign(score = "MAPS")
    .rename(columns={"maps":"mean", "se":"sem", "csq":"site", "set":"subset"})
    .replace({"all":"All"})
)
b_maps = (b_maps[b_maps.region == "Branchpoint"]
    .astype({"site":"int"})
)
df = (pd.concat([ns_phylop, ns_sai, ns_maps, b_phylop, b_sai, b_maps])
    [["score","region","subset","site","mean","sem"]]
)
df.to_csv(
    "../stats/near_splice_branch_stats_combined.tsv",
    sep="\t",
    index=False
    )
