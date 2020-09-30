import os
from slack import WebClient
import platform  # Get system name
from dotenv import load_dotenv  # Reading env files
from datetime import datetime  # Logging time

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


# If this script was imported by snakemake, use the default log_handler function
if __name__ != "__main__":

    def log_handler(msg):
        """A log handler for slack integration., adapted from snakemake"""
        # print(msg)
        # --------------------------------------------------------------------------#
        # Job Start
        # --------------------------------------------------------------------------#
        if msg["level"] == "run_info":
            # print(msg['msg'])
            client.chat_postMessage(
                channel="snakemake",
                blocks=[
                    {"type": "divider"},
                    {
                        "type": "section",
                        "text": {
                            "type": "mrkdwn",
                            "text": ":cyclone: Snakemake workflow on *"
                            + node_name
                            + "* started running the following rules at "
                            + "*"
                            + current_time
                            + "*.",
                        },
                    },
                    {"type": "section", "text": {"type": "mrkdwn", "text": msg["msg"]}},
                ],
            )
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
