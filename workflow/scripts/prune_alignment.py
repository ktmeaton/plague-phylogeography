#!/usr/bin/env python3

# Usage
# ./prune_alignment.py \
#   --metadata ../../results/metadata/all/metadata.tsv \
#   --matrix ../../results/snippy_multi/all/chromosome/filter5/snippy-multi.snps.dist \
#   --aln ../../results/snippy_multi/all/chromosome/filter5/snippy-multi.snps.aln \
#   --outdir ../../results/snippy_multi/all/chromosome/filter5/prune/

import os
import pandas as pd
import copy
import argparse
from Bio import AlignIO, Align, SeqIO

# Command-line argument capture
parser = argparse.ArgumentParser(
    description="Prune alignment for phylogeography.", add_help=True
)

parser.add_argument(
    "--metadata",
    help="Path to the metadata tsv.",
    action="store",
    dest="metadataPath",
    required=True,
)

parser.add_argument(
    "--matrix",
    help="SNP distance matrix.",
    action="store",
    dest="snpMatrixPath",
    required=True,
)

parser.add_argument(
    "--aln", help="Input alignment.", action="store", dest="alnPath", required=True,
)

parser.add_argument(
    "--outdir", help="Output directory.", action="store", dest="outDir", required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())
metadata_path = args["metadataPath"]
snp_matrix_path = args["snpMatrixPath"]
aln_path = args["alnPath"]
out_dir = args["outDir"]

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Load Metadata
metadata_basename = os.path.basename(metadata_path)
metadata_df = pd.read_csv(metadata_path, sep="\t")
metadata_df.fillna("NA", inplace=True)
metadata_df.set_index(metadata_df.columns[0], inplace=True)

# Load SNP Matrix
snp_mat_df = pd.read_csv(snp_matrix_path, sep="\t")
snp_mat_df.set_index(snp_mat_df.columns[0], inplace=True)

# -----------------------------------------------------
# Processing
branch_dict = {}
TIME_WINDOW = 25  # exclude if dates are within this range
ANCIENT_DATE_THRESHOLD = 1900  # samples older than this are considered ancient
GEO = "province"
GEO_ALT = "country"


for rec in metadata_df.iterrows():
    sample = rec[0]
    branch = rec[1]["branch_minor"]
    if branch == "NA":
        print("WARNING: sample {} does not have branch metadata.", format(sample))
        continue
    # Remove the subclade letter
    while branch[-1].isalpha():
        branch = branch[:-1]

    if branch not in branch_dict:
        branch_dict[branch] = {
            GEO: {},
        }

    geo_val = metadata_df[GEO][sample]
    if geo_val == "NA":
        geo_val = metadata_df[GEO_ALT][sample]

    date = metadata_df["date"][sample].lstrip("[").rstrip("]")
    # If it's a range, take the mean
    date_split = [int(d) for d in date.split(":")]
    if len(date_split) > 1:
        date = sum(date_split) / len(date_split)
    date = int(date)

    # Add the country if it hasn't been observed
    if geo_val not in branch_dict[branch][GEO]:
        branch_dict[branch][GEO][geo_val] = {"dates": {date: sample}}
        continue

    # Keep all ancient samples
    if date < ANCIENT_DATE_THRESHOLD:
        branch_dict[branch][GEO][geo_val]["dates"][date] = sample

    else:
        # Compare dates
        i_date = 0
        # By default, assume we're adding the sample
        add_sample = True
        for c_date in branch_dict[branch][GEO][geo_val]["dates"]:
            date_diff = abs(date - c_date)
            # If the date difference is too small, exclude
            if date_diff < TIME_WINDOW:
                # How to resolve ties? Want to minimize terminal branch length
                add_sample = False
            i_date += 1

        if add_sample:
            branch_dict[branch][GEO][geo_val]["dates"][date] = sample

sample_list = []
for branch in branch_dict:
    # print()
    # print(branch)
    for geo_val in branch_dict[branch][GEO]:
        # print("\t", geo_val)
        for date in branch_dict[branch][GEO][geo_val]["dates"]:
            sample = branch_dict[branch][GEO][geo_val]["dates"][date]
            strain = metadata_df["strain"][sample]
            sample_list.append(sample)
            # print("\t\t", date, strain, sample)

# -----------------------------------------------------
# Filter Metadata
filter_df = copy.deepcopy(metadata_df)
for sample in metadata_df.index:
    if sample not in sample_list:
        filter_df.drop(sample, inplace=True)
out_path_filter_df = os.path.join(out_dir, metadata_basename)

filter_df.to_csv(out_path_filter_df, sep="\t")

# -----------------------------------------------------
# Filter Alignment
aln_basename = os.path.basename(aln_path)
aln = AlignIO.read(aln_path, "fasta")
filter_seq = [rec for rec in aln if rec.id in sample_list or rec.id == "Reference"]
filter_aln = Align.MultipleSeqAlignment(filter_seq)

out_path_filter_aln = os.path.join(out_dir, aln_basename)
with open(out_path_filter_aln, "w") as outfile:
    count = SeqIO.write(filter_aln, outfile, "fasta")
    # print("\n",count, "alignments written.")
