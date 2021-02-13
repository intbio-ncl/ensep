from collections import defaultdict

import ensep
import rdkit.Chem as Chem

import msgpack
import tqdm


ensep.connect("bolt://repotrial.bioswarm.net:1687", "neo4j", "ensep")

query_1 = """
MATCH (rxn:Reaction)-[:CATALYZED_BY]->(p:Protein)
MATCH (rxn)-[:LEFT]->(l:Compound)
MATCH (l)-[:HAS_SUPER_PARENT]->(sl:Compound)
MATCH (rxn)-[:RIGHT]->(r:Compound)
MATCH (r)-[:HAS_SUPER_PARENT]->(sr:Compound)
RETURN DISTINCT p.sequence, sl.smiles, sr.smiles
"""

query_2 = """
MATCH (rxn:Reaction)-[:CATALYZED_BY]->(cplx:Complex)
MATCH (cplx)-[:HAS_COMPONENT]->(p:Protein)
MATCH (rxn)-[:LEFT]->(l:Compound)
MATCH (l)-[:HAS_SUPER_PARENT]->(sl:Compound)
MATCH (rxn)-[:RIGHT]->(r:Compound)
MATCH (r)-[:HAS_SUPER_PARENT]->(sr:Compound)
RETURN DISTINCT p.sequence, sl.smiles, sr.smiles
"""


d = defaultdict(list)

c = 0
with ensep.driver.session() as session:
    for i in tqdm.tqdm( session.run(query_1) ):
        a, b = sorted( [i["sl.smiles"], i["sr.smiles"]] )
        d[(a,b)].append( i["p.sequence"] )        


    for i in tqdm.tqdm( session.run(query_2) ):
        a, b = sorted( [i["sl.smiles"], i["sr.smiles"]] )
        d[(a,b)].append( i["p.sequence"] )


d = dict(d)

with open("ensep.msgpack", "wb") as f:
    f.write(msgpack.packb(d))
