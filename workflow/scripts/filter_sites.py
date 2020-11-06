# -----------------------------------------------------------------------------#
# Modules and Packages                                                         #
# -----------------------------------------------------------------------------#
import argparse  # Command-line argument parsing
from Bio import AlignIO  # Multiple Alignment parsing
from Bio import SeqIO  # Nucleotide parsing
from Bio.Seq import Seq  # Nucleotide parsing
from Bio.SeqRecord import SeqRecord  # Sequencing create for writing
import numpy as np  # Storing multiple alignment as array
import logging  # log rather than print

# import sys # Redirect logging info to sys stdout

# ------------------------------------------------------------------------------#
# Argument Parsing                                                             #
# ------------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=(
        "Filter sites in a multi-fasta alignment from snippy for missing data."
    ),
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--fasta",
    help="Path of input multi-fasta alignment.",
    action="store",
    dest="fastaPath",
    required=True,
)

parser.add_argument(
    "--missing",
    help="Percentage of samples that are allowed to be missing data at a site (AGCT).",
    action="store",
    dest="percMissing",
    required=True,
)

parser.add_argument(
    "--output",
    help="Path to the output filtered multi-fasta file.",
    action="store",
    dest="outputPath",
    required=True,
)

parser.add_argument(
    "--log",
    help="Write output messages to a log file.",
    action="store",
    dest="logPath",
    required=True,
)

parser.add_argument(
    "--keep-invariant",
    help="Retain invariant positions",
    action="store_true",
    dest="keepInvariant",
    required=False,
)


# Retrieve user parameters
args = vars(parser.parse_args())
fasta_path = args["fastaPath"]
prop_missing = float(int(args["percMissing"]) / 100)
prop_data = 1 - prop_missing
output_path = args["outputPath"]
log_path = args["logPath"]
keep_invariant = args["keepInvariant"]

# ------------------------------------------------------------------------------#
# Setup                                                                        #
# ------------------------------------------------------------------------------#
output_file = open(output_path, "w")
alignment_out_list = []
num_output_col = 0

# setup logging if a log path was supplied
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    handlers=[
        logging.FileHandler(log_path),
        # logging.StreamHandler(sys.stdout)
    ],
)

if keep_invariant:
    logging.info("Invariant positions will be included.")
else:
    logging.info("Invariant positions will be excluded.")
# ------------------------------------------------------------------------------#
# Processing                                                                   #
# ------------------------------------------------------------------------------#

# Analyze the input file
logging.info("Reading the input fasta file: " + str(fasta_path))
alignment_in = AlignIO.read(open(fasta_path), "fasta")

# Analyze the length of the alignment
logging.info("Reading alignment length...")
alignment_in_len = alignment_in.get_alignment_length()
num_samples = len(alignment_in[:, 1])
logging.info("Alignment length: " + str(alignment_in_len))

# Store the name of the samples for the output file
logging.info("Reading alignment sample names...")
alignment_in_samples = [rec.id for rec in alignment_in]
num_samples = len(alignment_in[:, 1])
logging.info("Number of samples: " + str(num_samples))

# Convert the input alignment into a numpy array to operate on columns
alignment_out_array = np.array([list("") for rec in alignment_in])

# Progress log setup
progress_log_increment = [
    5,
    10,
    15,
    20,
    25,
    30,
    35,
    40,
    45,
    50,
    55,
    60,
    65,
    70,
    75,
    80,
    85,
    90,
    95,
    100,
]
increment_ind = 0

# Iterate over records in the input alignment
logging.info(
    "Analyzing alignment for sites with >=" + str(prop_data) + " samples with data..."
)
for column in range(0, alignment_in_len):
    # check progress log
    prop_process = (column / alignment_in_len) * 100
    if prop_process >= progress_log_increment[increment_ind]:
        logging.info(str(int(prop_process)) + "% Complete")
        increment_ind += 1
    # Check the bases present in the current column/site
    column_seq = Seq(alignment_in[:, column])
    num_data = (
        column_seq.lower().count("g")
        + column_seq.lower().count("t")
        + column_seq.lower().count("a")
        + column_seq.lower().count("c")
    )
    site_prop = num_data / num_samples
    # Check if the amount of missing data passes user parameter
    if site_prop >= prop_data:
        column_seq_array = np.array([list(char) for char in column_seq])
        # If we're keeping invariant sites, write regardless
        if keep_invariant:
            alignment_out_array = np.append(
                alignment_out_array, column_seq_array, axis=1
            )
        else:
            join_seq = "".join(column_seq_array[:, 0])
            # Check if site is variable
            if join_seq.count(join_seq[0]) != len(join_seq):
                alignment_out_array = np.append(
                    alignment_out_array, column_seq_array, axis=1
                )

logging.info("100% Complete")

# Create a list of filtered fasta sequences
logging.info("Creating the filtered multiple alignment..")
for id, seq in zip(alignment_in_samples, alignment_out_array):
    join_seq = "".join(seq)
    format_seq = Seq(join_seq)
    alignment_out_list.append(SeqRecord(format_seq, id=id, description=""))

logging.info("Writing to the output fasta file: " + str(output_path))
SeqIO.write(alignment_out_list, output_file, "fasta")

alignment_out_len = len(alignment_out_array[0])
logging.info("Wrote a multi-fasta alignment of length: " + str(alignment_out_len))

# Clean up
output_file.close()
