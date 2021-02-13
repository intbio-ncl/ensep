from collections import defaultdict
import gzip
from pathlib import Path

from more_itertools import chunked
from tqdm import tqdm

add_new_compounds_query = """
UNWIND $items AS item
MATCH (c)-[:HAS_IDENTIFIER]->(i:Identifier {database:"ChEBI", accession: item.chebi_acc})
MERGE (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"KEGG", accession: item.kegg_acc})
"""

create_reaction_left = """
UNWIND $items AS item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "KEGG", accession: item.kegg_acc})
WITH rx, item.left AS left
UNWIND left AS l
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"KEGG", accession: l})
MERGE (rx)-[:LEFT]->(c)
"""

create_reaction_right = """
UNWIND $items AS item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "KEGG", accession: item.kegg_acc})
WITH rx, item.right as right
UNWIND right as r
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"KEGG", accession: r})
MERGE (rx)-[:RIGHT]->(c) 
"""

add_proteins = """
UNWIND $items AS item
MERGE (p:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "UniProt", accession: item.uniprot_acc})
"""

add_enzymes = """
UNWIND $items AS item
MATCH (p:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "UniProt", accession: item.uniprot_acc})
UNWIND item.reactions as r
MATCH (rxn:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "KEGG", accession: r})
MERGE (rxn)-[:CATALYZED_BY]->(p)
"""



class KEGG:
    def __init__(self, path):
        self.path = Path(path)
        if not self.path.exists():
            raise Exception(f"File {self.path} does not exist.")

    def parse_flatfile(self, infile):
        d = defaultdict(list)
        for line in infile.open():
            if line.startswith("///"):
                yield d
                d = defaultdict(list)

            elif line.startswith(" "):
                d[k].append(line.strip())
            
            else:
                stripped_line = line.strip()

                k = stripped_line.split()[0]
                k_len = len(k)
                v = f'{k_len * " "}{stripped_line[k_len:]}'
                d[k].append(v)

    def parse_compounds(self, driver):
        cpd_file = self.path / "kegg_compound_details.txt"
        if not cpd_file.exists():
            raise Exception()

        compounds = []

        for i in self.parse_flatfile(cpd_file):
            kegg_id = i["ENTRY"][0].strip().split()[0]
            links = [ln.strip().partition(": ") for ln in i["DBLINKS"] if "ChEBI:" in ln]
            links = [ln[2] for ln in links if ln]

            if len(links) != 1:
                continue

            a = {"kegg_acc": kegg_id, "chebi_acc": links[0]}
            compounds.append(a)
        
        with driver.session() as session:
            session.run(add_new_compounds_query, items=compounds)

    def parse_reactions(self, driver):
        rxn_file = self.path / "kegg_reaction_details.txt"
        if not rxn_file.exists():
            raise Exception()

        reactions = []

        for i in self.parse_flatfile(rxn_file):
            kegg_id = i["ENTRY"][0].strip().split()[0]

            eq = i["EQUATION"][0].partition(" <=> ")
            if not eq[1]:
                raise Exception()

            left, _ , right = eq
            left = [l.strip() for l in left.split("+")]
            right = [r.strip() for r in right.split("+")]
            a = {"kegg_acc": kegg_id, "left":left, "right":right}  
            reactions.append(a)

        with driver.session() as session:
            session.run(create_reaction_left, items=reactions)

        with driver.session() as session:
            session.run(create_reaction_right, items=reactions)

    def parse_orthology_hits(self, driver, annotations):
        ko_file = self.path / "kegg_ko_details.txt"
        if not ko_file.exists():
            raise Exception()

        ko_rxn = {}

        for i in self.parse_flatfile(ko_file):
            kegg_id = i["ENTRY"][0].strip().split()[0]
            reactions = [line.strip() for line in i["DBLINKS"] if line.strip().startswith("RN: ")]

            if len(reactions) == 0:
                continue
            assert len(reactions) == 1

            reactions = reactions[0].split()[1:]
            ko_rxn[kegg_id] = reactions
        
        items = []

        with gzip.open(annotations, "rt") as f:
            for line in f:
                # Comments line
                if line.startswith("#"):
                    continue
                
                # Insignificant match
                if not line.startswith("*"):
                    continue


                data = line.strip().split()
                data = data[:6] + [" ".join(data[6:])]
                protein = data[1].split("|")[1]
                ko = data[2]
                reactions = ko_rxn.get(ko)
                if not reactions:
                    continue

                item = {
                    "uniprot_acc" : protein,
                    "reactions": reactions
                }

                items.append(item)
        
        with driver.session() as session:

            for chunk in tqdm(chunked(items, 1_000)):
                session.run( add_proteins, items=chunk )
                session.run( add_enzymes, items=chunk )