#!/usr/bin/env python3

# Usage
# prune_alignment.py \
#   --metadata metadata.tsv \
#   --matrix snippy-multi.snps.dist \
#   --aln snippy-multi.snps.aln \
#   --outdir prune/

import os
import pandas as pd
import copy
import argparse
from Bio import AlignIO, Align, SeqIO
import subprocess

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
# ANCIENT_DATE_THRESHOLD = 1900  # samples older than this are considered ancient
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
    if date != "NA":
        # If it's a range, take the mean
        date_split = [int(d) for d in date.split(":")]
        if len(date_split) > 1:
            date = sum(date_split) / len(date_split)
        date = int(date)

    # Get the shortest pairwise distance (not to iteself)
    snp_diffs = snp_mat_df.loc[sample]
    min_diff = max(snp_diffs)

    for compare_sample in snp_mat_df.columns:
        # skip itself
        if compare_sample == sample:
            continue
        diff = snp_diffs[compare_sample]
        if diff < min_diff:
            min_diff = diff

    # Add the country if it hasn't been observed
    if geo_val not in branch_dict[branch][GEO]:
        branch_dict[branch][GEO][geo_val] = {"dates": {date: {sample: min_diff}}}
        continue

    # Keep all ancient samples
    # if date < ANCIENT_DATE_THRESHOLD:
    #    branch_dict[branch][GEO][geo_val]["dates"][date] = {sample: min_diff}

    else:
        # By default, assume we're adding the sample
        update_status = "add"
        for c_date in branch_dict[branch][GEO][geo_val]["dates"]:
            # if the date is NA, automatically include
            if date == "NA":
                update_status = "add"
                break

            date_diff = abs(date - c_date)
            # If the date difference is too small, exclude
            if date_diff < TIME_WINDOW:
                champion = branch_dict[branch][GEO][geo_val]["dates"][c_date]
                champion_sample = list(champion.keys())[0]
                champion_diff = list(champion.values())[0]
                contender = {sample: min_diff}

                if min_diff < champion_diff:
                    update_status = "replace"
                    break
                update_status = "none"

        # Replacing sample
        if update_status == "replace":
            print("Replacing:", champion, "with", contender)
            # Remove old champion
            branch_dict[branch][GEO][geo_val]["dates"].pop(c_date)
            # Add contender
            branch_dict[branch][GEO][geo_val]["dates"][date] = contender

        # Add new sample
        if update_status == "add":
            branch_dict[branch][GEO][geo_val]["dates"][date] = {sample: min_diff}

sample_list = []
for branch in branch_dict:
    print()
    print(branch)
    for geo_val in branch_dict[branch][GEO]:
        print("\t", geo_val)
        for date in branch_dict[branch][GEO][geo_val]["dates"]:
            sample = list(branch_dict[branch][GEO][geo_val]["dates"][date].keys())[0]
            strain = metadata_df["strain"][sample]
            sample_list.append(sample)
            print("\t\t", date, strain, sample)

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

out_path_aln = os.path.join(out_dir, aln_basename)

# Write temporary file
with open(out_path_aln + ".tmp", "w") as outfile:
    count = SeqIO.write(filter_aln, outfile, "fasta")
    # print("\n",count, "alignments written.")

# Remove constant sites
result = subprocess.run(["snp-sites", "-m", "-o", out_path_aln, out_path_aln + ".tmp"])

# Remove temporary file
os.remove(out_path_aln + ".tmp")
