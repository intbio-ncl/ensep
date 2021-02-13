from itertools import product
from pathlib import Path
import xml.etree.cElementTree as et

from more_itertools import chunked

_ie = lambda i: f"{{http://www.ebi.ac.uk/intenz}}{i}"
_cml = lambda i: f"{{http://www.xml-cml.org/schema/cml2/react}}{i}"

get_chebi_ids_query = """
MATCH (i:Identifier {database:"ChEBI"})<-[:HAS_IDENTIFIER]-(c:Compound)
RETURN i.accession
"""

create_reaction_left = """
UNWIND $items as item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "Rhea", accession: item.rhea_acc})
WITH rx, item.left as left
UNWIND left as l
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"ChEBI", accession: l})
MERGE (rx)-[:LEFT]->(c)
"""

create_reaction_right = """
UNWIND $items as item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "Rhea", accession: item.rhea_acc})
WITH rx, item.right as right
UNWIND right as r
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"ChEBI", accession: r})
MERGE (rx)-[:RIGHT]->(c) 
"""

create_protein = """
UNWIND $items as item
MERGE (en:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "UniProt", accession: item.protein})
"""

create_enzyme_reaction = """
UNWIND $items as item
MATCH (rx)-[:HAS_IDENTIFIER]->(:Identifier {database: "Rhea", accession: item.reaction})
MATCH (en:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "UniProt", accession: item.protein})
MERGE (rx)-[:CATALYZED_BY]->(en)
"""

class Intenz:
    def __init__(self, path):
        self.path = Path(path)
        if not self.path.exists():
            raise Exception()

    def parse_enzymes(self, driver):
        infile = self.path / "intenz.xml"

        enzymes = []
        for _, elem in et.iterparse(infile.open(), events=("end",)):
            if elem.tag == _ie("enzyme"):

                try:
                    reactions = [r.attrib["id"] for r in elem.find(_ie("reactions")).findall(_cml("reaction"))]
                except AttributeError:
                    continue # Log that it indicates there are no reactions. Example is 1.1.1.5, which is deprecated

                try:
                    proteins = [
                        link.attrib["accession_number"] 
                        for link in elem.find(_ie("links")).findall(_ie("link"))
                        if link.attrib.get("db") == "UniProt"
                        ]
                except AttributeError:
                    continue # Log that it indicates there are no proteins. Example is 1.1.5.7.

                for r, p in product(reactions, proteins):
                    a = {"reaction": r, "protein": p}
                    enzymes.append(a)

        with driver.session() as session:
            for chunk in chunked(enzymes, 1_000):
                session.run(create_protein, items=chunk)
                session.run(create_enzyme_reaction, items=chunk)



    def parse_reactions(self, driver):
        infile = self.path / "intenz.xml"

        with driver.session() as session:
            chebi_accs = {i["i.accession"] for i in session.run(get_chebi_ids_query)}

        reactions = []
        for _, elem in et.iterparse(infile.open(), events=("end",)):
            if elem.tag == _cml("reaction"):
                reaction_id = elem.attrib["id"]

                reactants = [reac.find(_cml("molecule")) for reac in elem.find(_cml("reactantList"))]
                reactants = [i.find(_cml("identifier")).attrib["value"] for i in reactants]
                reactants = [i.split(":")[1] for i in reactants]
                reactants = [i for i in reactants if i in chebi_accs]

                products = [reac.find(_cml("molecule")) for reac in elem.find(_cml("productList"))]
                products = [i.find(_cml("identifier")).attrib["value"] for i in products]
                products = [i.split(":")[1] for i in products]
                products = [i for i in products if i in chebi_accs]

                a = {"rhea_acc": reaction_id, "left": reactants, "right": products}
                reactions.append(a)
        
        with driver.session() as session:
            for chunk in chunked(reactions, 1_000):
                session.run(create_reaction_left, items=chunk)
                session.run(create_reaction_right, items=chunk)