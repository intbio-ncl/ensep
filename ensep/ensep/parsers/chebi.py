import csv
import gzip
from pathlib import Path

from more_itertools import chunked
from rdkit import Chem, RDLogger

from ensep.common.logger import logger

RDLogger.DisableLog("rdApp.*")


def process_structures(reader, driver):
    logger.info("Processing ChEBI structures")
    compounds = []
    r = filter(lambda i: i["TYPE"] == "SMILES", reader)

    for item in r:
        smiles_string = item["STRUCTURE"]
        chebi_id = item["COMPOUND_ID"]

        m = Chem.MolFromSmiles(smiles_string)
        if not m:
            logger.debug(
                f"Unable to parse structure for ChEBI {chebi_id}, skipping"
                )
            continue

        canonical_smiles = Chem.MolToSmiles(m)
        cpd_object = {
            "smiles": canonical_smiles,
            "accession": chebi_id
            }
        compounds.append(cpd_object)

    query = """
    UNWIND $cpds AS d
    MERGE (c:Compound {smiles: d.smiles})
    MERGE (i:Identifier {accession: d.accession, database: "ChEBI"})
    MERGE (i)<-[:HAS_IDENTIFIER]-(c)
    """

    logger.info("Adding compounds to EnSeP")

    with driver.session() as session:
        for chunk in chunked(compounds, 1_000):
            session.run(query, cpds=chunk)

    logger.info(f"Added {len(compounds):,} compounds to EnSeP")


class ChebiParser:
    def __init__(self, infile):
        self.infile = Path(infile)
        if not self.infile.exists():
            raise Exception(f"{infile} given does not exist")

    def parse(self, driver):
        logger.info("Parsing ChEBI data files")
        structure_file = self.infile / "structures.csv.gz"
        with gzip.open(structure_file, "rt") as f:
            reader = csv.DictReader(f)
            process_structures(reader, driver)
