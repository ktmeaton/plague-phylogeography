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
    "--keep-singleton",
    help="Keep singleton sites.",
    action="store_true",
    dest="keepSingleton",
    required=False,
)

# Retrieve user parameters
args = vars(parser.parse_args())
fasta_path = args["fastaPath"]
prop_missing = float(int(args["percMissing"]) / 100)
prop_data = 1 - prop_missing
output_path = args["outputPath"]
log_path = args["logPath"]
keep_singleton = args["keepSingleton"]

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
    filename=log_path,
    filemode="w",
    # handlers=[
    #    logging.FileHandler(log_path),
    #    # logging.StreamHandler(sys.stdout)
    # ],
)

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
# Counter to store number of singleton sites
biallelic_singleton_sites = 0
multiallelic_singleton_sites = 0
pseudo_singleton_sites = 0
parsimony_informative_sites = 0
passing_filter_sites = 0
failing_filter_sites = 0
passing_filter_parsimony_sites = 0
passing_filter_singleton_sites = 0

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
    count_g = column_seq.lower().count("g")
    count_t = column_seq.lower().count("t")
    count_a = column_seq.lower().count("a")
    count_c = column_seq.lower().count("c")
    count_list = [count_g, count_t, count_a, count_c]
    # store the variant type (default singleton)
    site_type = "singleton"

    # Remove singleton columns sites
    # Biallelic format [0,0,1,x]
    if count_list.count(0) == 2 and count_list.count(1) == 1:
        biallelic_singleton_sites += 1
        if not keep_singleton:
            continue
    # Multi allelic format [0,1,1,x]
    elif count_list.count(0) == 1 and count_list.count(1) == 2:
        multiallelic_singleton_sites += 1
        if not keep_singleton:
            continue
    # Pseudo allelic format [0,1,x,x] or [1,1,x,x]
    elif count_list.count(1) > 0:
        pseudo_singleton_sites += 1
        if not keep_singleton:
            continue
    else:
        parsimony_informative_sites += 1
        site_type = "parsimony"

    # Decide if this
    num_data = count_g + count_t + count_a + count_c
    site_prop = num_data / num_samples
    # Check if the amount of missing data passes user parameter
    if site_prop >= prop_data:
        passing_filter_sites += 1
        if site_type == "singleton":
            passing_filter_singleton_sites += 1
        elif site_type == "parsimony":
            passing_filter_parsimony_sites += 1

        column_seq_array = np.array([list(char) for char in column_seq])
        join_seq = "".join(column_seq_array[:, 0])
        # Check if site is variable
        if join_seq.count(join_seq[0]) != len(join_seq):
            alignment_out_array = np.append(
                alignment_out_array, column_seq_array, axis=1
            )
    else:
        failing_filter_sites += 1

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
total_singleton_sites = (
    biallelic_singleton_sites + multiallelic_singleton_sites + pseudo_singleton_sites
)
logging.info(
    "Parsimony informative sites: "
    + str(parsimony_informative_sites)
    + " ("
    + str(int(parsimony_informative_sites / alignment_in_len * 100))
    + "%)"
)
logging.info(
    "Biallelic singleton sites: "
    + str(biallelic_singleton_sites)
    + " ("
    + str(int(biallelic_singleton_sites / alignment_in_len * 100))
    + "%)"
)
logging.info(
    "Multiallelic singleton sites: "
    + str(multiallelic_singleton_sites)
    + " ("
    + str(int(multiallelic_singleton_sites / alignment_in_len * 100))
    + "%)"
)
logging.info(
    "Pseudo singleton sites: "
    + str(pseudo_singleton_sites)
    + " ("
    + str(int(pseudo_singleton_sites / alignment_in_len * 100))
    + "%)"
)
logging.info(
    "Total singleton sites: "
    + str(total_singleton_sites)
    + " ("
    + str(int(total_singleton_sites / alignment_in_len * 100))
    + "%)"
)
logging.info(
    "Singleton sites passing missing data filter: {}".format(
        passing_filter_singleton_sites
    )
)
logging.info(
    "Parsimony informative sites passing missing data filter: {}".format(
        passing_filter_parsimony_sites
    )
)

logging.info("Total sites passing missing data filter: " + str(passing_filter_sites))
logging.info("Total sites failing missing data filter: " + str(failing_filter_sites))
logging.info("Wrote a multi-fasta alignment of length: " + str(alignment_out_len))

# Clean up
output_file.close()
