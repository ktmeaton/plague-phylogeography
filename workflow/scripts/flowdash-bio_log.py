import os
import platform  # Get node name
from dotenv import load_dotenv  # Reading env files

# import logging
# logging.basicConfig(level=logging.DEBUG)

# --------------------------------------------------------------------------#
# General Setup
# --------------------------------------------------------------------------#
load_dotenv()
token = os.getenv("FLOWDASH_BIO_TOKEN")

if not token:
    print(
        "Please generate an API token for flowdash-bio using your username as: ",
        "\ncurl -X GET -u username https://flowdash-bio-stage.herokuapp.com/api/tokens",
        "\nand then export FLOWDASH_BIO_TOKEN='yourtoken' to a .env file.",
    )
    exit(-1)

# FLOWDASH_BIO_URL=
node_name = platform.node()


# If this script was imported by snakemake, use the default log_handler function
if __name__ != "__main__":

    def log_handler(msg):
        """A log handler for slack integration., adapted from snakemake"""
        # print(msg)
        # --------------------------------------------------------------------------#
        # Job Start
        # --------------------------------------------------------------------------#
        if msg["level"] == "run_info":
            print(msg["msg"])
        # --------------------------------------------------------------------------#
        # Error message
        # --------------------------------------------------------------------------#
        if msg["level"] == "error":
            print(msg["msg"])

        # --------------------------------------------------------------------------#
        # Finished message
        # --------------------------------------------------------------------------#
        if msg["level"] == "progress" and msg["done"] == msg["total"]:
            print(msg["msg"])

        # --------------------------------------------------------------------------#
        # Report Upload
        # --------------------------------------------------------------------------#
        if msg["level"] == "info" and "Report created" in msg["msg"]:
            report_file = msg["msg"].replace("Report created: ", "").rstrip(".")
            print(report_file)


# test_msg = {"level": "progress", "done": 1, "total" : 1, "msg": "Test message."}
# log_handler(test_msg)
