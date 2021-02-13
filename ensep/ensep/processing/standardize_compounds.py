from multiprocessing import Pool

import molvs
from more_itertools import chunked
from rdkit import Chem, RDLogger
import tqdm

from ensep.common.logger import logger

STANDARDIZER = molvs.Standardizer(tautomer_transforms=[])

get_all_smiles_query = """
MATCH (c:Compound)
RETURN c.smiles
"""

add_super_parent_relationsips = """
UNWIND $items AS item
MATCH (base:Compound {smiles: item.base_smiles})
MERGE (parent:Compound {smiles: item.parent_smiles})
MERGE (base)-[:HAS_SUPER_PARENT]->(parent)
"""

add_compound_inchi = """
UNWIND $items as item
MATCH (c:Compound {smiles: item.smiles})
SET c.inchi = item.inchi
"""

def generate_compound_super_parent(smiles):
    m = Chem.MolFromSmiles(smiles)
    if not m:
        logger.warning(f"SMILES string cannot be parsed: {smiles}")
        return None
    try:
        p = STANDARDIZER.super_parent(m)
    except Chem.rdchem.AtomValenceException:
        logger.warning(f"Cannot generate super parent: {smiles}")
        return None
    except Chem.rdchem.KekulizeException:
        logger.warning(f"Cannot generate super parent: {smiles}")
        return None

    p_smiles = Chem.MolToSmiles(p)

    details = {"base_smiles": smiles, "parent_smiles": p_smiles}
    return details

def generate_compound_inchi(smiles):
    m = Chem.MolFromSmiles(smiles)
    if not m:
        logger.warning(f"SMILES string cannot be parsed: {smiles}")
        return None
    inchi = Chem.MolToInchi(m)

    details = {"smiles": smiles, "inchi": inchi}
    return details

def get_compounds(driver):
    logger.debug("Collecting all compounds from EnSeP KB")

    # Set compounds
    compounds = []

    logger.debug("Generating super parents of compounds in EnSeP KB")
    with driver.session() as session, Pool(4) as P:
        cpd_count = sum(1 for _ in session.run(get_all_smiles_query))
        logger.debug(f"Collected {cpd_count:,} compounds from EnSeP KB")
        cpd_smiles = (i["c.smiles"] for i in session.run(get_all_smiles_query))
        
        i = P.imap_unordered(generate_compound_super_parent, cpd_smiles)
        for details in tqdm.tqdm(i, position=0, total=cpd_count):
            if details:
                compounds.append(details)  

    logger.debug(f"Adding {len(compounds):,} super parent relationships")
    with driver.session() as session:
        for chunk in chunked(compounds, 1_000):
            session.run(add_super_parent_relationsips, items=chunk)


    # Reset compounds
    compounds = []

    logger.debug("Generating InChI strings for compounds in EnSeP KB")
    with driver.session() as session, Pool(4) as P:
        cpd_cound = sum(1 for _ in session.run(get_all_smiles_query))
        logger.debug(f"Collected {cpd_count:,} compounds from EnSeP KB")
        cpd_smiles = (i["c.smiles"] for i in session.run(get_all_smiles_query))

        i = P.imap_unordered(generate_compound_inchi, cpd_smiles)
        for details in tqdm.tqdm(i, position=0, total=cpd_count):
            if details:
                compounds.append(details)  

    logger.debug(f"Adding InChI strings to {len(compounds):,} compounds")
    with driver.session() as session:
        session.run(add_compound_inchi, items=compounds)