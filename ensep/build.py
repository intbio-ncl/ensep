import ensep

import rdkit.Chem as Chem

print("Connect")
ensep.connect("bolt://localhost:7687", "neo4j", "ensep")

ensep.create_indexes(ensep.driver)

# cp = ensep.parsers.chebi.ChebiParser("../datafiles/chebi")
# cp.parse(ensep.driver)

# mp = ensep.parsers.metacyc.MetaCyc("../datafiles/metacyc")
# mp.parse(ensep.driver)

# kegg = ensep.parsers.kegg.KEGG("../datafiles/kegg")
# kegg.parse_compounds(ensep.driver)
# kegg.parse_reactions(ensep.driver)
# kegg.parse_orthology_hits(ensep.driver, "../datafiles/kofamscan_analyses/uniprot_sprot.fasta.out.gz")

# intenz = ensep.parsers.intenz.Intenz("../datafiles/intenz")
# intenz.parse_reactions(ensep.driver)
# intenz.parse_enzymes(ensep.driver)

# Processing
# ensep.processing.standardize_compounds.get_compounds(ensep.driver)

uniprot = ensep.parsers.uniprot.UniProt()
uniprot.collect(ensep.driver)
