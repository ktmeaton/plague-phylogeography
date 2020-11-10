import os
import platform  # Get node name
from dotenv import load_dotenv  # Reading env files
import requests
from requests.auth import HTTPBasicAuth

# import logging
# logging.basicConfig(level=logging.DEBUG)

# --------------------------------------------------------------------------#
# General Setup
# --------------------------------------------------------------------------#
load_dotenv()


# Static Token
flowdash_bio_token = os.getenv("FLOWDASH_BIO_TOKEN")

# Dynamic Token
if not flowdash_bio_token:
    flowdash_bio_username = os.getenv("FLOWDASH_BIO_USERNAME")
    flowdash_bio_password = os.getenv("FLOWDASH_BIO_PASSWORD")
    if not flowdash_bio_username or not flowdash_bio_password:
        print(
            "Please export one of the following in your .env file: ",
            "\nFLOWDASH_BIO_TOKEN",
        )
        exit(-1)

    response = requests.get(
        url="https://flowdash-bio-stage.herokuapp.com/api/tokens",
        auth=HTTPBasicAuth(flowdash_bio_username, flowdash_bio_password),
    )

    flowdash_bio_token = response.json()["token"]

FLOWDASH_BIO_URL = (
    "https://flowdash-bio-stage.herokuapp.com/api/workflows/attr?"
    "node={}&total_jobs={}&completed_jobs={}&running_jobs={}&failed_jobs={}"
)
flowdash_bio_headers = {"Authorization": "Bearer %s" % flowdash_bio_token}

# Default values
data = {
    "node": platform.node(),
    "total_jobs": 0,
    "completed_jobs": 0,
    "running_jobs": 0,
    "failed_jobs": 0,
}

# If this script was imported by snakemake, use the default log_handler function
if __name__ != "__main__":

    def log_handler(msg):
        """A log handler for slack integration., adapted from snakemake"""
        # print(msg)
        api_method = ""
        # --------------------------------------------------------------------------#
        # Job Start
        # --------------------------------------------------------------------------#
        if msg["level"] == "run_info":
            struct_msg = [line.strip().split("\t") for line in msg["msg"].split("\n")]
            # Total jobs will be the last element
            total_jobs = int(struct_msg[-1][0])
            data["total_jobs"] = total_jobs
            api_method = "POST"

        # --------------------------------------------------------------------------#
        # Error message
        # --------------------------------------------------------------------------#
        if msg["level"] == "error":
            print(msg)

        # --------------------------------------------------------------------------#
        # Finished message
        # --------------------------------------------------------------------------#
        elif msg["level"] == "progress" and msg["done"] == msg["total"]:
            data["total_jobs"] = msg["total"]
            data["completed_jobs"] = msg["done"]
            api_method = "PUT"

        # --------------------------------------------------------------------------#
        # Progress message
        # --------------------------------------------------------------------------#
        elif msg["level"] == "progress":
            data["total_jobs"] = msg["total"]
            data["completed_jobs"] = msg["done"]
            api_method = "PUT"

        # --------------------------------------------------------------------------#
        # Report Upload
        # --------------------------------------------------------------------------#
        elif msg["level"] == "info" and "Report created" in msg["msg"]:
            report_file = msg["msg"].replace("Report created: ", "").rstrip(".")
            print(report_file)

        query_url = FLOWDASH_BIO_URL.format(
            data["node"],
            data["total_jobs"],
            data["completed_jobs"],
            data["running_jobs"],
            data["failed_jobs"],
        )
        if api_method == "POST":
            # Post to the flowdash-bio API
            result = requests.post(url=query_url, headers=flowdash_bio_headers)
            print(api_method)
            print(query_url)
            print(result)
        elif api_method == "PUT":
            # Post to the flowdash-bio API
            result = requests.put(url=query_url, headers=flowdash_bio_headers)
            print(api_method)
            print(query_url)
            print(result)


# test_msg = {"level": "progress", "done": 1, "total" : 1, "msg": "Test message."}
# log_handler(test_msg)
