#!/usr/bin/env python3

import asyncio
from pyppeteer import launch
import time

# ------------------------------------------
# URL Construction

clade = "all"
# dataset = "main/beast/geo/" + clade
dataset = "main/ml/" + clade
geo_res = "province"
color = "population_inferred"
panel = "map"
transmissions = "hide"
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

url = "https://nextstrain.org/community/ktmeaton/plague-phylogeography-projects@main/main/ml/all?branchLabel=Branch%20Support%20Conf%20Char&branches=hide&d=tree&l=scatter&legend=closed&onlyPanels&p=full&scatterX=date_mean&sidebar=closed&transmissions=hide"
output = "{}_{}_{}.pdf".format("all", "scatter", "time")


print(url)
print()

wait = 10
width = 1920
height = 1080
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
