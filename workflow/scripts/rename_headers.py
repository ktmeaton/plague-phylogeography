import argparse
import os

# ------------------------------------------------------------------------------#
# Argument Parsing                                                             #
# ------------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=(
        "Rename headers of a fasta or gff file by removing the version number."
    ),
    add_help=True,
)

# Argument groups for the program

parser.add_argument(
    "--file",
    help="Path to input file.",
    action="store",
    dest="filePath",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())
file_path = args["filePath"]
output_path = file_path + ".tmp"
ext = os.path.splitext(file_path)[1]
with open(file_path) as temp_in_file:
    with open(output_path, "w") as temp_out_file:
        # Header parsing for fasta files
        if ext == ".fna":
            for line in temp_in_file:
                strip_line = line.strip()
                if ">" in line:
                    header = os.path.splitext(strip_line.split(" ")[0])[0]
                    temp_out_file.write(header + "\n")
                else:
                    temp_out_file.write(strip_line + "\n")
        if ext == ".gff":
            for line in temp_in_file:
                if "#" in line:
                    temp_out_file.write(line)
                if "#" not in line:
                    split_line = line.strip().split("\t")
                    header = os.path.splitext(split_line[0])[0]
                    other = "\t".join(split_line[1:])
                    temp_out_file.write(header + "\t" + other + "\n")
os.rename(output_path, file_path)
