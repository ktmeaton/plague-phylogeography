#!/usr/bin/env python3

import asyncio
from pyppeteer import launch
import time

# dataset = (
#  "https://nextstrain.org/community/" +
#  "ktmeaton/plague-phylogeography-projects@main/main/full/all"
# )
# query = "?d=tree&legend=open&m=div&onlyPanels&p=full&sidebar=closed"

# dataset = "http://localhost:4000/1.PRE"
# query = "?c=date_mean&d=map&m=div&onlyPanels&p=full&sidebar=closed&tl=blank"
# url = dataset + query
# url = "https://nextstrain.org/community/ktmeaton/plague-phylogeography-projects@3a11d4ecd0dc8cdf1354c3599503b71553b75c12/main/full/1.PRE?c=date_mean&label=Biosample%20Accession:NA&m=div&tl=country"

# Russia Geocoding Compare
url = "https://nextstrain.org/community/ktmeaton/plague-phylogeography-projects@main/main/full/all?c=countryn&d=map&f_country=Russia&onlyPanels&p=full&sidebar=closed&transmissions=hide"
wait = 10
width = 450
height = 450


async def main():
    browser = await launch()
    page = await browser.newPage()
    await page.setViewport({"width": width, "height": height, "deviceScaleFactor": 4})
    print("Waiting {} seconds for Nextstrain to load page...".format(wait))
    await page.goto(url)
    time.sleep(wait)
    print("Taking a pdf screenshot...")
    await page.emulateMedia("screen")
    await page.pdf({"path": "screenshot.pdf"})
    await browser.close()


asyncio.get_event_loop().run_until_complete(main())
