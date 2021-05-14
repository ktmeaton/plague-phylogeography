import asyncio
from pyppeteer import launch
import time

# dataset = (
#  "https://nextstrain.org/community/" +
#  "ktmeaton/plague-phylogeography-projects@main/main/full/all"
# )
# query = "?d=tree&legend=open&m=div&onlyPanels&p=full&sidebar=closed"

dataset = "http://localhost:4000/1.PRE"
query = "?c=date_mean&d=map&m=div&onlyPanels&p=full&sidebar=closed&tl=blank"


url = dataset + query

wait = 5


async def main():
    browser = await launch()
    page = await browser.newPage()
    await page.setViewport({"width": 900, "height": 600, "deviceScaleFactor": 4})
    print("Waiting {} seconds for Nextstrain to load page...".format(wait))
    await page.goto(url)
    time.sleep(wait)
    print("Taking a pdf screenshot...")
    await page.emulateMedia("screen")
    await page.pdf({"path": "screenshot.pdf"})
    await browser.close()


asyncio.get_event_loop().run_until_complete(main())
