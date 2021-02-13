from neo4j import GraphDatabase
import neobolt

from . import parsers, processing, common

driver = None
_logger = common.logger.logger

def connect(uri, username, password):
    _logger.info("Creating Neo4j DB driver")
    global driver
    driver = GraphDatabase.driver(uri=uri, auth=(username, password), encrypted=False)

def create_indexes(driver):
    _logger.info("Creating indexes for EnSeP")
    idxs = [":Compound(smiles)", ":Identifier(accession)", ":Identifier(database)"]
    with driver.session() as session:
        for idx in idxs:
            query = f"CREATE INDEX ON {idx}"
            try:
                session.run(query)
            except neobolt.exceptions.ClientError:
                continue
            except neobolt.exceptions.TransientError:
                continue

def delete_db(driver):
    query = """MATCH (n)
    DETACH DELETE n
    """
    with driver.session() as session:
        try:
            session.run(query)
        except neobolt.exceptions.TransientError:
            pass

