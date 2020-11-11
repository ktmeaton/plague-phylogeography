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

# --------------------------------------------------------------------------#
# API URLS
FLOWDASH_BIO_URL_BASE = "https://flowdash-bio.herokuapp.com"
FLOWDASH_BIO_URL = (
    FLOWDASH_BIO_URL_BASE
    + "/api/workflows"
    + "?node={}&total_jobs={}&completed_jobs={}&running_jobs={}&failed_jobs={}"
)

# For debugging, use ngrok tunnels
try:
    result = requests.get("http://localhost:4040/api/tunnels")
    if result.status_code == 200:
        data = result.json()
        public_url = data["tunnels"][0]["public_url"]
        FLOWDASH_BIO_URL_BASE = public_url
        FLOWDASH_BIO_URL = (
            FLOWDASH_BIO_URL_BASE
            + "/api/workflows"
            + "?node={}&total_jobs={}&completed_jobs={}&running_jobs={}&failed_jobs={}"
        )
except requests.exceptions.ConnectionError:
    pass

# --------------------------------------------------------------------------#
# API Tokens
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

    token_url = FLOWDASH_BIO_URL_BASE + "/api/tokens"
    response = requests.get(
        url=token_url, auth=HTTPBasicAuth(flowdash_bio_username, flowdash_bio_password),
    )

    flowdash_bio_token = response.json()["token"]

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
        # Job Restart
        # --------------------------------------------------------------------------#
        if msg["level"] == "info" and "Trying to restart" in msg["msg"]:
            # search for failed workflows
            query_url = (
                FLOWDASH_BIO_URL_BASE
                + "/api/workflows?node={}&status=Failed".format(data["node"])
            )
            result = requests.get(url=query_url, headers=flowdash_bio_headers)
            workflows = result.json()["workflows"]
            # Grab latest workflow (-1)
            workflow_id = list(workflows.keys())[-1]
            # Substract a "failed_job" if retrying
            data["failed_jobs"] = workflows[workflow_id]["failed_jobs"] - 1
            data["total_jobs"] = workflows[workflow_id]["total_jobs"]
            data["completed_jobs"] = workflows[workflow_id]["completed_jobs"]
            data["running_jobs"] = workflows[workflow_id]["running_jobs"]
            api_method = "PUT"

        # --------------------------------------------------------------------------#
        # Error message
        # --------------------------------------------------------------------------#
        if msg["level"] == "error":
            # first search for running workflows
            query_url = (
                FLOWDASH_BIO_URL_BASE
                + "/api/workflows?node={}&status=Running".format(data["node"])
            )
            result = requests.get(url=query_url, headers=flowdash_bio_headers)
            workflows = result.json()["workflows"]
            try:
                workflow_id = list(workflows.keys())[-1]
                data["failed_jobs"] = 1
            except IndexError:
                # instead try failed workflows
                query_url = (
                    FLOWDASH_BIO_URL_BASE
                    + "/api/workflows?node={}&status=Failed".format(data["node"])
                )
                result = requests.get(url=query_url, headers=flowdash_bio_headers)
                workflows = result.json()["workflows"]
                workflow_id = list(workflows.keys())[-1]
                data["failed_jobs"] = workflows[workflow_id]["failed_jobs"] + 1

            data["total_jobs"] = workflows[workflow_id]["total_jobs"]
            data["completed_jobs"] = workflows[workflow_id]["completed_jobs"]
            data["running_jobs"] = workflows[workflow_id]["running_jobs"]
            api_method = "PUT"

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
            # report_file = msg["msg"].replace("Report created: ", "").rstrip(".")
            # print(report_file)
            pass

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
        elif api_method == "PUT":
            # Post to the flowdash-bio API
            result = requests.put(url=query_url, headers=flowdash_bio_headers)


# test_msg = {"level": "progress", "done": 1, "total" : 1, "msg": "Test message."}
# log_handler(test_msg)
