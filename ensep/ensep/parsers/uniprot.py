import json
from multiprocessing.pool import ThreadPool
import time
from typing import Iterable

from more_itertools import chunked
import requests
from tqdm import tqdm
from typing_extensions import Final

URL: Final = "https://www.ebi.ac.uk/proteins/api/proteins/{}"

get_protein_accs = """
MATCH (p:Protein)-[:HAS_IDENTIFIER]->(i:Identifier {database:"UniProt"})
WHERE p.sequence IS NULL
RETURN i.accession
"""

update_protein_properties = """
UNWIND $items AS item
MATCH (p:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database:"UniProt", accession:item.accession})
SET p.sequence = item.sequence
SET p.taxonomy = item.taxonomy
SET p.names = item.names
"""

def collect_url(
    url, backoff=1, try_count=0, max_tries=3, response_type="application/json"
):
    if try_count > max_tries:
        raise Exception("Maximum number of tries at obtaining URL exceeded")

    headers = {"accept": response_type}

    try:
        response = requests.get(url, headers=headers)
        return response.text

    except KeyboardInterrupt:
        raise KeyboardInterrupt()
    except Exception as E:
        print(E)
        time.sleep(backoff)
        return collect_url(
            url,
            backoff=backoff * 2,
            try_count=try_count + 1,
            max_tries=max_tries,
            response_type=response_type,
        )


class UniProt:
    def __init__(self, threads: int = 20):
        self.threads: int = threads

    def collect(self, driver):
        with driver.session() as session:
            accs = [i["i.accession"] for i in session.run(get_protein_accs)]
        
        total = sum(1 for _ in accs)
        urls = (URL.format(i) for i in accs)

        items = []

        with ThreadPool(self.threads) as pool:
            for response in tqdm( pool.imap_unordered(collect_url, urls), total=total ):
                if "errorMessage" in response:
                    continue
                r = json.loads(response)

                if "recommendedName" in r["protein"]:
                    names = [r["protein"]["recommendedName"]["fullName"]["value"]]
                    
                elif "submittedName" in r["protein"]:
                    names =  [i["fullName"]["value"] for i in r["protein"]["submittedName"]]
                else:
                    raise Exception()

                org = r["organism"]["taxonomy"]
                seq = r["sequence"]["sequence"]
                acc = r["accession"]

                item = {
                    "taxonomy" : org,
                    "sequence" : seq,
                    "names" : names,
                    "accession": acc
                }

                items.append(item)
        
        with driver.session() as session:
            for chunk in chunked(items, 1_000):
                session.run(update_protein_properties, items=chunk)
