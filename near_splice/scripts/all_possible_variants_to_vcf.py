""" This script annotates all possible SNVs at near-splice and coding positions
and writes them to a VCF output.
NB this script is rather memory intensive, requiring up to 36GB memory.
"""

# Import the relevant modules
import numpy as np
import pandas as pd
import os

def main(region):
    """ Run all functions in this script.
    """
    print(f"Processing {region} data.")

    data = f"../outputs/{region}_contexts.tsv"

    df = read_data(data)\
        .pipe(vcf_format)\
        .pipe(chrom_format)

def read_data(data):
    """ Read the data to memory.
    """
    df = pd.read_csv(data, sep="\t", header=None)
    return df

def vcf_format(df):
    """ Annotate every possible SNV at each position.
    Reformat the data to .vcf format.
    Starting format line 0 (bedtools output): chr1:926001-926004(+)  TCC
    """
    print("Reformatting to vcf.")
    # create dummy dataframes for simplicity
    a = df[0].str.split(":", expand=True)
    b = a[1].str.split("-", expand=True)

    # extract the relevant information
    chrom = a[0].rename("chrom")
    pos = (b[0].astype(int) + 2).rename("pos")
    ref = df[1].str.slice(1,2).rename("ref")

    df = pd.concat([chrom, pos, ref], axis=1)\
        .sort_values(by=["chrom","pos"])

    # identify all possible SNVs
    df["alt"] = pd.Series([["A","T","C","G"]] * len(df))
    df = df.explode("alt", ignore_index=True)
    df = df[df.ref != df.alt]

    # convert to .vcf format
    df["id"] = "."
    df["qual"] = "."
    df["filter"] = "."
    df["info"] = "."

    df = df[["chrom","pos","id","ref","alt","qual","filter","info"]]

    return df

def chrom_format(df):
    """ The GRCh38 format usually has a "chr" prefix before the chromosome name.
    However, the pre-computed GRCh38 SpliceAI scores do not.
    This function creates a second dataframe, without the "chr" prefix, and
    passes both dfs to write_vcf.
    """
    df_no_chr = df.copy()
    df_no_chr["chrom"] = df_no_chr["chrom"].str.slice(3) # Remove "chr" prefix

    # Dummy variables
    chr = [df, "chr"]
    no_chr = [df_no_chr, "no_chr_prefix"]

    for x in [chr, no_chr]: write_vcfs(*x)

def write_vcfs(df, chr_format="chr"):
    """ Write to .vcf
    """
    print(f"Writing {chr_format} to output.")

    # Define the .vcf header
    head = open("../data/vcf_38_header.txt")
    header = head.read()

    # Write split dataframe to output
    with open(f"../outputs/{region}_all_snvs_{chr_format}.vcf", "w") as output:
        output.write(header)
        df.to_csv(output, index=False, sep="\t", header=False)

if __name__ == "__main__":
    for region in ["near_splice", "coding"]: main(region)
