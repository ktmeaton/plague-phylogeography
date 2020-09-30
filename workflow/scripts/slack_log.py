import os
from slack import WebClient
from slack.errors import SlackApiError
import platform  # Get system name
from dotenv import load_dotenv  # Reading env files
from datetime import datetime  # Logging time

# import logging
# logging.basicConfig(level=logging.DEBUG)


def log_handler(msg):
    """A log handler for slack integration., adapted from snakemake"""
    # --------------------------------------------------------------------------#
    # Token Setup
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
    # --------------------------------------------------------------------------#
    # Compose Message
    # --------------------------------------------------------------------------#
    node_name = platform.node()
    now = datetime.now()
    current_time = now.strftime("%H:%M")
    # --------------------------------------------------------------------------#
    # Error message
    # --------------------------------------------------------------------------#
    if msg["level"] == "error":
        try:
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
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    # --------------------------------------------------------------------------#
    # Finished message
    # --------------------------------------------------------------------------#
    print(msg)
    if msg["level"] == "progress" and msg["done"] == msg["total"]:
        try:
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
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'


# test_msg = {"level": "progress", "done": 1, "total" : 1, "msg": "Test message."}
# log_handler(test_msg)
