"""
This script filters for exons meeting certain criteria within the GENCODE
v29 comprehensive gene annotation (GRCh38)

Filtering criteria:
- feature = "CDS"
- gene_type = "protein_coding"
- transcript_type = "protein_coding"
- annotation != "level 3" (automated annotation)
- tag = "CCDS", "appris_principal", "appris_candidate_longest", "appris_candidate", or "exp_conf"
"""

# Import the relevant modules
import numpy as np
import pandas as pd

def main():
	""" Run every function in this script.
	"""
	gencode_path = "/public_data_resources/GENCODE/v29/GRCh38/gencode.v29.annotation.gtf"
	output_path = "../outputs/coding_exons.tsv"

	df = read_gtf(gencode_path)\
		.pipe(cds_filter, feature="CDS")\
		.pipe(write_output, path=output_path)

	return df

def read_gtf(path):
	""" Read a .gtf file to memory.
	"""
	names = (["chrom", "source", "feature", "start", "end", "score", "strand",
		"phase","attr"])

	df = pd.read_csv(
		path,
		sep="\t",
		comment="#",
		header=None,
		names = names
		)

	return df

def cds_filter(df, feature):
	""" Filter for CDS exons meeting certain criteria
	"""
	df = df[df.feature==feature]

	tags = (['tag "CCDS"',
		'tag "appris_principal_1"',
	 	'tag "appris_candidate_longest"',
		'tag "appris_candidate"',
		'tag "exp_conf"'])

	tags = "|".join(tags) # For OR filtering in mask4, below.

	mask1 = df.attr.str.contains('gene_type "protein_coding"')
	mask2 = df.attr.str.contains('transcript_type "protein_coding"')
	mask3 = df.attr.str.contains('; level 3;')
	mask4 = df.attr.str.contains(tags)

	df = df[mask1 & mask2 & ~mask3 & mask4]

	print(f"{len(df)} {feature} features meet the filtering criteria.") # n = 401314

	return df

def write_output(df, path):
	""" Write the dataframe to .tsv
	"""
	df.to_csv(path, sep="\t", index=False)

if __name__ == "__main__":
	main()
