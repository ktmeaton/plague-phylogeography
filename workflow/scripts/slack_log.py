import os
from slack import WebClient
import platform  # Get system name
from dotenv import load_dotenv  # Reading env files
from datetime import datetime  # Logging time
import argparse  # Parse arguments if this script is main

# import logging
# logging.basicConfig(level=logging.DEBUG)

# --------------------------------------------------------------------------#
# General Setup
# --------------------------------------------------------------------------#
load_dotenv()
token = os.getenv("SLACK_TOKEN")

if not token:
    print(
        "Please create an API token for your Slack App and",
        "export SLACK_TOKEN='yourtoken' to a .env file.",
    )
    exit(-1)

client = WebClient(token=token)
node_name = platform.node()
now = datetime.now()
current_time = now.strftime("%H:%M")

# If this script was not imported by snakemake, perform a custom slack operation
if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True,)
    parser.add_argument(
        "--rules",
        help="Rule(s) successfully completed by snakemake in csv format.",
        action="store",
        dest="rulesList",
        required=True,
    )
    parser.add_argument(
        "--complete",
        help="Include flag is snakemake rules are complete.",
        action="store_true",
        dest="rulesComplete",
        required=False,
    )

    args = vars(parser.parse_args())
    rules_list = args["rulesList"].split(",")
    rules_complete = args["rulesComplete"]

    # If we're only beginning
    if not rules_complete:
        client.chat_postMessage(
            channel="snakemake",
            blocks=[
                {"type": "divider"},
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": ":heavy_check_mark: Snakemake workflow on *"
                        + node_name
                        + "* is starting the following rules:\n"
                        + "\n".join(["\t\t - " + rule for rule in rules_list]),
                    },
                },
            ],
        )
    if rules_complete:
        client.chat_postMessage(
            channel="snakemake",
            blocks=[
                {"type": "divider"},
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": ":heavy_check_mark: Snakemake workflow on *"
                        + node_name
                        + "* completed the following rules:\n"
                        + "\n".join(["\t\t - " + rule for rule in rules_list]),
                    },
                },
            ],
        )


# If this script was imported by snakemake, use the default log_handler function
if __name__ != "__main__":

    def log_handler(msg):
        """A log handler for slack integration., adapted from snakemake"""
        print(msg)
        # --------------------------------------------------------------------------#
        # Error message
        # --------------------------------------------------------------------------#
        if msg["level"] == "error":
            client.chat_postMessage(
                channel="snakemake",
                blocks=[
                    {"type": "divider"},
                    {
                        "type": "section",
                        "text": {
                            "type": "mrkdwn",
                            "text": ":x: Snakemake workflow on *"
                            + node_name
                            + "* exited with an error at "
                            + "*"
                            + current_time
                            + "*.",
                        },
                    },
                    {"type": "section", "text": {"type": "mrkdwn", "text": msg["msg"]}},
                ],
            )

        # --------------------------------------------------------------------------#
        # Finished message
        # --------------------------------------------------------------------------#
        if msg["level"] == "progress" and msg["done"] == msg["total"]:
            client.chat_postMessage(
                channel="snakemake",
                blocks=[
                    {"type": "divider"},
                    {
                        "type": "section",
                        "text": {
                            "type": "mrkdwn",
                            "text": ":heavy_check_mark: Snakemake workflow on *"
                            + node_name
                            + "* completed sucessfully at "
                            + "*"
                            + current_time
                            + "*.",
                        },
                    },
                ],
            )


# test_msg = {"level": "progress", "done": 1, "total" : 1, "msg": "Test message."}
# log_handler(test_msg)
