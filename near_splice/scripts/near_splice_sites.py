"""
This script annotates the near splice positions for every exon of interest from
GENCODE v29. Every genomic position of interest will be annotated for its
near-splice position (A-25 to A+10 and D-10 to D+10)

NB identical exons for multiple transcripts are present.

Features on chrX and chrY are not included.
"""

# Import the relevant modules
import numpy as np
import pandas as pd

c = "category" # We often convert dtypes to categories below

def main():
    """ Run all functions in the script.
    """

    df = load_data()\
    .pipe(autosomes_only)\
    .pipe(extract_attributes)\
    .pipe(splice_junctions)\
    .pipe(splice_positions)\

    output = "../outputs/near_splice_positions.tsv"

    df.to_csv(output, sep="\t", index=False)

    return df

def load_data(data="../outputs/coding_exons.tsv"):
    """ Read the filtered GENCODE exons into memory
    """
    df = pd.read_csv(
        data,
        sep="\t",
        usecols=["chrom", "start", "end", "strand", "attr"],
        dtype = {"chrom":c, "strand":c}
        )

    return df

def autosomes_only(df):
	""" Keeps only autosomal features
	"""
	print("Dropping features on sex chromosomes")
	sex_chroms = ["chrX", "chrY"]
	df = df[~df.chrom.isin(sex_chroms)]
	return df

def extract_attributes(df):
    """ Extract relevant identifiers from the attributes column
    """
    attr = df.attr.str.split("; ", expand=True).iloc[:,:11]

    ensg = attr[0].str.slice(9, -1).rename("ensg") # gene
    enst = attr[1].str.slice(15, -1).rename("enst") # transcript
    ense = attr[7].str.slice(9, -1).rename("ense") # exon
    exon_number = attr[6].str.slice(12).astype(int).rename("exon_number")

    # Merge the attributes back into the original dataframe
    df = pd.concat([df, ensg, enst, ense, exon_number], axis=1)\
        .drop("attr", axis=1)

    print(f"There are {df.ense.nunique()} unique CDS exons in {df.enst.nunique()} unique transcripts in the GENCODE annotation")

    return df

def splice_junctions(df):
    """ Identifies the splice junctions between every feature in a transcript.
    Ensure to subset the MANE annotation to the feature of interest,
    e.g. "CDS", before applying this function.
	Annotate exon boundaries with a "region" (acceptor or donor).
    """
    # Identify the extreme ends of the exons, which are not
    # involved in splicing.
    grouped = df.groupby("enst")
    df.loc[:,"start_min"] = grouped["start"].transform(min)
    df.loc[:,"end_max"] = grouped["end"].transform(max)

    # The function evaluates start and end positions separately.
    # Extreme starts and ends are excluded.
    mask_start = df["start"] == df["start_min"]
    mask_end = df["end"] == df["end_max"]

    df_starts = df[~mask_start].copy()
    df_starts.loc[:,"pos"] = df_starts["start"]

    df_ends = df[~mask_end].copy()
    df_ends.loc[:,"pos"] = df_ends["end"]

    # The strand determines whether a start or end is a splice acceptor or
	# donor.
    df_starts.loc[:,"region"] = np.where(df_starts["strand"]=="+","acceptor","donor")
    df_ends.loc[:,"region"] = np.where(df_ends["strand"]=="+","donor","acceptor")

    # Start and end annotations are combined.
    df = pd.concat([df_starts, df_ends])\
        .drop(["end_max","start_min"], axis=1)\
        .sort_values(["chrom","pos"])\
        .drop_duplicates()

    return df

def splice_positions(df):
    """ Annotate the near-splice positions genome-wide.
    """
    # Find the "site" for each position (its distance from the splice junction)
    a_sites = pd.Series([list(range(-25,11))]*len(df)) # acceptors
    d_sites = pd.Series([list(range(-10,11))]*len(df)) # donors

    df.loc[df.region=="acceptor", "site"] = a_sites
    df.loc[df.region=="donor", "site"] = d_sites

    df = df.explode("site", ignore_index=True)

    df = df.astype({"ensg":c, "enst":c, "ense":c, "region":c, "site":int})

    # The splicing "site" is dictated by the strand of the transcript
    df.loc[df.strand == "+", "pos"] = df.pos + df.site
    df.loc[df.strand == "-", "pos"] = df.pos - df.site

    # Drop unnecessary columns
    df = df[["chrom", "pos", "strand", "region", "site"]]

    # Drop any splice positions with ambiguous annotation, and drop duplicates
    print("Dropping ambiguous positions.")
    df = df.drop_duplicates(subset=["chrom", "pos", "region", "site"])
    df = df.drop_duplicates(subset=["chrom","pos"], keep=False)

    print(f"There are {len(df)} valid near-splice positions")

    return df

if __name__ == "__main__":
    main()
