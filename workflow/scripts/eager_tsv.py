# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import os  # file path searching and creation

# ------------------------------------------------------------------------------#
# Argument Parsing                                                             #
# ------------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=("Create the tsv input file for the nf-core/eager pipeline."),
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--files",
    help="Full path of input fastq files to run with nf-core/eager.",
    action="store",
    dest="filesPath",
    required=True,
)

parser.add_argument(
    "--organism",
    help="Organism name for nf-core/eager.",
    action="store",
    dest="orgName",
    required=True,
)

parser.add_argument(
    "--tsv",
    help="Output tsv for nf-core/eager.",
    action="store",
    dest="tsvPath",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())
files_path = args["filesPath"]
org_name = args["orgName"]
tsv_path = args["tsvPath"]

# EAGER tsv input column names
EAGER_HEADER = (
    "Sample_Name\t"
    "Library_ID\t"
    "Lane\t"
    "Colour_Chemistry\t"
    "SeqType\t"
    "Organism\t"
    "Strandedness\t"
    "UDG_Treatment\t"
    "R1\t"
    "R2\t"
    "BAM"
)

# Default values
sample_name = "NA"
library_id = "NA"
LANE = "1"
# For now, default to HiSeq color chemistry
COLOR_CHEMISTRY = "4"
seq_type = "NA"
ORGANISM = org_name
STRANDEDNESS = "double"
UDG = "none"
R1 = "NA"
R2 = "NA"
BAM = "NA"

sample_file_dict = {}
tsv_file = open(tsv_path, "w")

# Iterate through the different subreads/sra accessions
for file in files_path.split(" "):
    sample_dir = os.path.dirname(file)
    sample_name = os.path.basename(sample_dir)
    sample_files = os.listdir(sample_dir)
    # Iterate through the single or paired fastq files
    for file in sample_files:
        library_id = file.split("_")[0]
        if library_id not in sample_file_dict:
            sample_file_dict[library_id] = []
        file_path = os.path.join(sample_dir, file)
        sample_file_dict[library_id].append(file_path)

tsv_file.write(EAGER_HEADER + "\n")

for library_id in sample_file_dict:
    # Reset file specific values
    seq_type = "NA"
    R1 = "NA"
    R2 = "NA"
    # Single end
    if len(sample_file_dict[library_id]) == 1:
        R1 = sample_file_dict[library_id][0]
        seq_type = "SE"
    elif len(sample_file_dict[library_id]) == 2:
        for file_path in sample_file_dict[library_id]:
            if "_1.fastq.gz" in file_path:
                R1 = file_path
            elif "_2.fastq.gz" in file_path:
                R2 = file_path
            seq_type = "PE"
    tsv_file.write(
        sample_name
        + "\t"
        + library_id
        + "\t"
        + LANE
        + "\t"
        + COLOR_CHEMISTRY
        + "\t"
        + seq_type
        + "\t"
        + ORGANISM
        + "\t"
        + STRANDEDNESS
        + "\t"
        + UDG
        + "\t"
        + R1
        + "\t"
        + R2
        + "\t"
        + BAM
        + "\n"
    )

tsv_file.close()
