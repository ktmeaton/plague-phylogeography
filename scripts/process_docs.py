#!/usr/bin/env python3
"""
@author: Katherine Eaton

Generate process docs for ReadTheDocs from nextflow file.

Usage:
./process_docs.py --nf ../main.nf --rst ../docs/process/process_all.rst

"""

# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#
import argparse  # Command-line argument parsing
import os  # File path checking
import sys  # Generate system exit error

# This program should only be called from the command-line
if __name__ != "__main__":
    quit()

# -----------------------------------------------------------------------------#
#                            Argument Parsing                                  #
# -----------------------------------------------------------------------------#
parser = argparse.ArgumentParser(
    description="Generate process docs for ReadTheDocs from nextflow file.",
    add_help=True,
)

parser.add_argument(
    "--nf",
    help="Path to the nextflow pipeline file.",
    action="store",
    dest="nfPath",
    required=True,
)

parser.add_argument(
    "--rst",
    help="Path to the output rst docs file.",
    action="store",
    dest="rstPath",
    required=True,
)

# Retrieve user parameters
args = vars(parser.parse_args())

nf_path = args["nfPath"]
rst_path = args["rstPath"]

# -----------------------------------------------------------------------------#
#                                Error Catching                                #
# -----------------------------------------------------------------------------#
# Check if NF file exists
if not os.path.exists(nf_path):
    print("An error occurred while trying to open", nf_path)
    sys.exit(1)

# -----------------------------------------------------------------------------#
#                            Constants and Variables                           #
# -----------------------------------------------------------------------------#
rst_file = open(rst_path, "w")
H1_CHAR = "-"
H2_CHAR = "*"
H3_CHAR = "-"
TABLE_CHAR = "="
TABLE_COL_WIDTH = 40
# -----------------------------------------------------------------------------#
#                               Processing                                     #
# -----------------------------------------------------------------------------#
with open(nf_path, "r") as nf_file:
    for line in nf_file:
        line = line.strip()
        # Skip everything that is not process
        if line.startswith("process"):
            # Parse the process lines
            split_process = line.split(" ")
            process_name = split_process[1].strip("{")
            format_process_name = process_name.replace("_", " ").title()
            rst_file.write(
                "\n"
                + format_process_name
                + "\n"
                + H3_CHAR * len(format_process_name)
                + "\n\n"
            )
            line = nf_file.readline().strip()
            # Begin the process docstring
            if line == ("/*"):
                # Write the process description.
                process_description = nf_file.readline().strip()
                rst_file.write(process_description + "\n\n")
                # Process the IO docs
                line = nf_file.readline().strip()

                io_doc_exists = False
                while not line == ("*/"):
                    # Blank lines signal table line
                    # If line is blank and no io docstring found yet
                    if not line and not io_doc_exists:
                        line = nf_file.readline().strip()
                        continue
                    # If line is blank and io docstring has been found
                    if not line and io_doc_exists:
                        rst_file.write(
                            TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + "\n\n"
                        )
                        line = nf_file.readline().strip()
                        continue
                    # Process a subsection of the docstring
                    if line in ["Input:", "Output:", "Publish:"]:
                        io_doc_exists = True
                        io_section = line
                        # Write top border
                        rst_file.write(
                            TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + "\n"
                        )
                        # Write input column headers
                        rst_file.write(
                            io_section
                            + " " * (TABLE_COL_WIDTH - len(io_section) + 1)
                            + "Type"
                            + " " * (TABLE_COL_WIDTH - len("Type") + 1)
                            + "Description"
                            + " " * (TABLE_COL_WIDTH - len("Description") + 1)
                            + "\n"
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + " "
                            + TABLE_CHAR * TABLE_COL_WIDTH
                            + "\n"
                        )
                        line = nf_file.readline().strip()
                        continue
                    io_split = line.split("(")
                    # Retrieve the input value
                    io_name = io_split[0]
                    # Retrieve the input type
                    io_split = io_split[1].split(")")
                    io_type = io_split[0]
                    # Retrieve the input description
                    io_split = io_split[1].split(": ")
                    io_desc = io_split[1]
                    # Figure out where the put the process links
                    if "process" in io_desc:
                        # process :ref:`ncbimeta_db_update<NCBImeta DB Update>`
                        io_desc_split = io_desc.split("process ")
                        io_process_name = io_desc_split[1].strip(".")
                        io_process_link = (
                            "process :ref:`"
                            + io_process_name
                            + "<"
                            + io_process_name.title()
                            + ">`"
                        )
                        io_desc = io_desc_split[0] + io_process_link
                    rst_file.write(
                        io_name
                        + " " * (TABLE_COL_WIDTH - len(io_name) + 1)
                        + io_type
                        + " " * (TABLE_COL_WIDTH - len(io_type) + 1)
                        + io_desc
                        + " " * (TABLE_COL_WIDTH - len(io_desc) + 1)
                        + "\n"
                    )
                    # Read in the next IO line
                    line = nf_file.readline().strip()
                # Read past the ending docstring
                if line == "*/":
                    line = nf_file.readline().strip()

        # Process the script code
        if line == "script:" or line == "shell:":
            script_type = line.strip(":")
            # Skip the current line which is """ or '''
            line = nf_file.readline().strip()
            # Check if the immediate next line is """ or '''
            line = nf_file.readline().strip()
            if line != "'''" and line != '"""':
                rst_file.write("**" + script_type + "**::" + "\n\n")
                while line != "'''" and line != '"""':
                    rst_file.write("\t" + line + "\n")
                    line = nf_file.readline().strip()
