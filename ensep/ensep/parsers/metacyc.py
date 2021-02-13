from collections import defaultdict
from itertools import product
from pathlib import Path
from pprint import PrettyPrinter

from rdkit import Chem, RDLogger
from typing_extensions import Final

from ensep.common.logger import logger

RDLogger.DisableLog("rdApp.*")

PRINTER: Final = PrettyPrinter()

get_chebi_ids_query = """
MATCH (i:Identifier {database:"ChEBI"})<-[:HAS_IDENTIFIER]-(c:Compound)
RETURN i.accession
"""

get_metacyc_ids_query = """
MATCH (i:Identifier {database:"MetaCyc"})<-[:HAS_IDENTIFIER]-(c:Compound)
RETURN i.accession
"""

add_accession_to_entity_query = """
UNWIND $items as item
MATCH (i:Identifier {database: "ChEBI", accession: item.chebi_acc})<-[:HAS_IDENTIFIER]-(c:Compound)
MERGE (c)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.metacyc_acc})
"""

create_new_compound_query = """
UNWIND $items as item
MERGE (c:Compound {smiles: item.smiles})
MERGE (i:Identifier {database: "MetaCyc", accession: item.metacyc_acc})<-[:HAS_IDENTIFIER]-(c)
"""

create_new_protein_query = """
UNWIND $items as item
MERGE (p:Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "UniProt", accession: item.uniprot_acc})
MERGE (p)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.metacyc_acc})
"""

create_new_complex_query = """
UNWIND $items as item
MERGE (c:Complex)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.metacyc_acc})
WITH c, item.components AS components
UNWIND components as component
MATCH (p: Protein)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: component})
MERGE (c)-[:HAS_COMPONENT]->(p)
"""

create_reaction_left = """
UNWIND $items as item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.metacyc_acc})
WITH rx, item.left as left
UNWIND left as l
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"MetaCyc", accession: l})
MERGE (rx)-[:LEFT]->(c)
"""

create_reaction_right = """
UNWIND $items as item
MERGE (rx:Reaction)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.metacyc_acc})
WITH rx, item.right as right
UNWIND right as r
MATCH (c)-[:HAS_IDENTIFIER]->(:Identifier {database:"MetaCyc", accession: r})
MERGE (rx)-[:RIGHT]->(c) 
"""

create_enzyme_reaction = """
UNWIND $items as item
MATCH (rx)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.reaction})
MATCH (en)-[:HAS_IDENTIFIER]->(:Identifier {database: "MetaCyc", accession: item.enzyme})
MERGE (rx)-[:CATALYZED_BY]->(en)
"""


class MetaCyc:
    def __init__(self, path):
        self.path = Path(path)
        if not self.path.exists():
            raise Exception(f"File {self.path} does not exist.")

    def parse_flatfile(self, infile):
        d = defaultdict(list)
        for line in infile.open(encoding="latin1"):
            if line.startswith("#"):
                continue

            elif line.startswith("//"):
                yield d
                d = defaultdict(list)

            elif line.startswith("/"):
                d[k][-1] = f"{d[k][-1].strip()} {line.strip()[1:]}"

            else:
                k, _, v = line.strip().partition(" - ")
                d[k].append(v)

    def parse_compounds(self, driver):
        logger.info("Parsing compounds from MetaCyc")
        cpd_path = self.path / "data/compounds.dat"
        if not cpd_path.exists():
            raise Exception(f"File {cpd_path} could not be found.")

        logger.info("Obtaining ChEBI accessions of compounds in EnSeP")
        with driver.session() as session:
            chebi_accs = {
                i["i.accession"] for i in session.run(get_chebi_ids_query)
            }

        # Two lists -- one for items with 1 ChEBI xref, and one for items with none.
        chebi_xref = []
        no_chebi_xref = []

        for i in self.parse_flatfile(cpd_path):
            cpd_id = i['UNIQUE-ID'][0]
            # Get accessions of cross-reference to ChEBI
            links = [
                ln[1:-1].split()[1][1:-1]
                for ln in i["DBLINKS"]
                if "CHEBI" in ln
            ]
            # If there's more than one ChEBI cross-reference, skip due to 
            # ambiguity.
            if len(links) > 1:
                logger.warning(
                    f"Multiple ChEBI xrefs for {cpd_id}, skipping"
                    )
                continue

            # Confirm the accession is in our database.
            links = [ln for ln in links if ln in chebi_accs]

            if links:
                a = {"chebi_acc": links[0], "metacyc_acc": cpd_id}
                chebi_xref.append(a)
            else:
                # Get, and parse, SMILES string (if available)
                smiles = i.get("SMILES")
                if not smiles:
                    logger.warning(
                        f"No ChEBI xref / SMILES string for {cpd_id}, skipping"
                    )
                    continue

                m = Chem.MolFromSmiles(smiles[0])
                if not m:
                    logger.warning(
                        f"No ChEBI xref / parsable SMILES string for {cpd_id}, skipping"
                    )
                    continue

                canonical_smiles = Chem.MolToSmiles(m)
                a = {
                    "smiles": canonical_smiles,
                    "metacyc_acc": cpd_id,
                }
                no_chebi_xref.append(a)

        logger.info(
            f"{len(chebi_xref)} compounds in MetaCyc with a single ChEBI cross-reference in EnSeP"
            )
        logger.info(
            f"{len(no_chebi_xref)} compounds in MetaCyc without a ChEBI cross-reference in EnSeP")

        with driver.session() as session:
            session.run(add_accession_to_entity_query, items=chebi_xref)
        with driver.session() as session:
            session.run(create_new_compound_query, items=no_chebi_xref)

        logger.info("Compounds successfully added to EnSeP")

    def parse_proteins(self, driver):
        logger.info("Parsing proteins and complexes from MetaCyc")
        pro_path = self.path / "data/proteins.dat"
        if not pro_path.exists():
            raise Exception(f"File {pro_path} could not be found.")

        complexes = []
        proteins = []

        for i in self.parse_flatfile(pro_path):
            # Determine the type based on the presence of "COMPONENTS" in keys.
            t = "Complex" if "COMPONENTS" in i.keys() else "Protein"

            # If complex, store accession and component accessions.
            if t == "Complex":
                a = {
                    "metacyc_acc": i["UNIQUE-ID"][0],
                    "components": i["COMPONENTS"],
                }
                complexes.append(a)

            # If protein, get UniProt links (because this we where we'll obtain
            # the sequences from). Multiple links / zero link records are dropped.
            else:
                protein_id = i["UNIQUE-ID"][0]
                # [1:-1] here gets rid of bracekets.
                links = [ln[1:-1] for ln in i["DBLINKS"] if "UNIPROT" in ln]
                # [1:-1] here gets rid of quotes.
                links = [ln.split()[1][1:-1] for ln in links]

                if not links:
                    logger.warning(
                        f"No UniProt IDs for {protein_id}, skipping"
                    )
                    continue

                if len(links) > 1:
                    logger.warning(
                        f"Multiple UniProt IDs for {protein_id}, skipping"
                    )
                    continue

                a = {"metacyc_acc": i["UNIQUE-ID"][0], "uniprot_acc": links[0]}
                proteins.append(a)

        logger.info(f"First pass completed -- {len(proteins)} proteins and {len(complexes)} complexes identified")
        
        # Due to our dropping of Proteins, check to make sure the protein
        # components of complexes are going to be stored in EnSeP.
        logger.info(f"Filtering complexes to ensure at least one component present in EnSeP")
        proteins_accs = {pro["metacyc_acc"] for pro in proteins}
        updated_complexes = []
        for cplx in complexes:
            cplx["components"] = [i for i in cplx["components"] if i in proteins_accs]
            if cplx["components"]:
                updated_complexes.append(cplx)
        complexes = updated_complexes
        del proteins_accs
        logger.info(f"Second pass completed -- {len(complexes)} complexes identified")

        logger.info("Adding proteins to EnSeP")
        with driver.session() as session:
            session.run(create_new_protein_query, items=proteins)
        logger.info("Proteins successfully added to EnSeP")
        logger.info("Adding complexes to EnSeP")
        with driver.session() as session:
            session.run(create_new_complex_query, items=complexes)
        logger.info("Complexes successfully added to EnSeP")

    def parse_reactions(self, driver):
        logger.info("Parsing reactions from MetaCyc")

        reac_path = self.path / "data/reactions.dat"
        if not reac_path.exists():
            raise Exception(f"File {reac_path} could not be found.")

        logger.info("Obtaining compounds in EnSeP with MetaCyc IDs")
        with driver.session() as session:
            metacyc_accs = {
                i["i.accession"] for i in session.run(get_metacyc_ids_query)
            }

        reactions = []

        for i in self.parse_flatfile(reac_path):
            rid = i["UNIQUE-ID"][0]
            a = {"metacyc_acc": rid}
            a["left"] = [cpd for cpd in i["LEFT"] if cpd in metacyc_accs]
            a["right"] = [cpd for cpd in i["RIGHT"] if cpd in metacyc_accs]

            if (not a["left"]) or (not a["right"]) :
                logger.warning(
                    f"Reaction {rid} lacking compounds (all left or all right) in EnSeP, skipping"
                    )
                continue

            reactions.append(a)

        logger.info(f"Adding {len(reactions)} reactions to EnSeP")    
        with driver.session() as session:
            session.run(create_reaction_left, items=reactions)
        with driver.session() as session:
            session.run(create_reaction_right, items=reactions)
        logger.info("Reactions successfully added to EnSeP")

    def parse_enzymes(self, driver):
        logger.info("Parsing enzymes from MetaCyc")
        enrx_path = self.path / "data/enzrxns.dat"
        if not enrx_path.exists():
            raise Exception(f"File {enrx_path} could not be found.")

        enzrxns = []

        for i in self.parse_flatfile(enrx_path):
            if not i.get("ENZYME") or not i.get("REACTION"):
                continue

            for e, r in product(i["ENZYME"], i["REACTION"]):
                a = {"enzyme": e, "reaction": r}
                enzrxns.append(a)

        logger.info("Adding enzymes to EnSeP")
        with driver.session() as session:
            session.run(create_enzyme_reaction, items=enzrxns)
        logger.info("Successfully added enzymes to EnSeP")

    def parse(self, driver):
        logger.info("Parsing MetaCyc data files")
        self.parse_compounds(driver)
        self.parse_proteins(driver)
        self.parse_reactions(driver)
        self.parse_enzymes(driver)
