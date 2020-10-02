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


if __name__ == "__main__":
    client.chat_postMessage(
        channel="snakemake",
        blocks=[
            {"type": "divider"},
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": ":exclamation: *"
                    + node_name
                    + "* is testing slack integration at "
                    + "*"
                    + current_time
                    + "*.",
                },
            },
        ],
    )
