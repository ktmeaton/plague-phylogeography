#!/usr/bin/env python3

import asyncio
from pyppeteer import launch
import time

# ------------------------------------------
# URL Construction

clade = "1.ORI"
dataset = "main/beast/geo/" + clade
geo_res = "province"
color = "continent"
panel = "map"
transmissions = "show"
repo = "plague-phylogeography-projects"
base_url = "https://nextstrain.org/community/ktmeaton"
branch = "main"

project = base_url + "/" + repo + "@" + branch + "/" + dataset
options = [
    "r=" + geo_res,
    "c=" + color,
    "d=" + panel,
    "sidebar=closed",
    "transmissions=" + transmissions,
    "legend=closed",
    "p=full",
    "onlyPanels",
]
query = "?" + "&".join(options)
url = project + query

output = "{}_{}_{}.pdf".format(clade, panel, color)

print(url)
print()

wait = 10
width = 720
height = 720
options = {
    "path": output,
    "format": "A0",
    "printBackground": True,
    "landscape": True,
}


async def main():
    browser = await launch()
    page = await browser.newPage()

    await page.setViewport({"width": width, "height": height, "deviceScaleFactor": 4})
    print("Waiting {} seconds for Nextstrain to load page...".format(wait))

    # Navigate to page
    await page.goto(url)
    # Change to landscape
    time.sleep(wait)
    print("Taking a pdf screenshot...")
    await page.emulateMedia("screen")
    await page.pdf(options)
    print("Saved screenshot to " + output)
    await browser.close()


asyncio.get_event_loop().run_until_complete(main())
