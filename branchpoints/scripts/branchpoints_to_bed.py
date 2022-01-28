""" This script reformats branchpoints and adjacent positions into .bed format.
"""

# Import relevant modules
import numpy as np
import pandas as pd

# Open the relevant files
df = (pd.read_csv(
    "../data/branchpoints_38_liftover.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "pos", "name", "score", "strand"]
    )
    .sort_values(by=["chrom","start","score"])
    .drop_duplicates(["chrom","start"], keep="last") # Keep highest score
)

# Remove positions not on autosomes
df = df[~df.chrom.str.contains("_alt")]
df = df[~df.chrom.isin(["chrX","chrY"])]

# Annotate positions 5 bases either side of the branchpoint
df["site"] = [range(-5,6)]*len(df)
df = df.explode("site")
df["site"] = df.site.astype(int)

# Get strand-specific annotation of positions near branchpoints
df_fwd = df[df["strand"] == "+"].copy()
df_rev = df[df["strand"] == "-"].copy()

df_fwd["pos"] = df_fwd["pos"] + df_fwd["site"]
df_rev["pos"] = df_rev["pos"] - df_rev["site"]

df = (pd.concat([df_fwd,df_rev])
    .drop_duplicates()
    .drop_duplicates(subset=["chrom","pos"], keep=False) # No ambiguous sites
    .sort_values(by=["chrom","pos"])
)

# Convert to .be format and write to output
df["start"] = df["pos"]-2
df["end"] = df["pos"]+1

df = df[["chrom","start","end","site","score","strand"]]

df.to_csv("../outputs/branchpoints.bed", sep="\t", index=False, header=False)
