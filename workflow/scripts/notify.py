import requests
import os
import json

ntfy_url = os.getenv("NTFY_TOPIC")

def log_handler(msg):

    if msg["level"] == "error":
        requests.post(
            ntfy_url,
            data = "Workflow failed :(".encode(encoding = "utf-8")
        )

    if msg["level"] == "progress" and msg["done"] == msg["total"]:
        requests.post(
            ntfy_url,
            data = "Workflow finished :)".encode(encoding = "utf-8")
        )
