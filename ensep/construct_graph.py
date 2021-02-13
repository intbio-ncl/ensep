import subprocess
from itertools import chain, combinations
from multiprocessing import Pool, cpu_count
from pathlib import Path
from pprint import PrettyPrinter
from tempfile import NamedTemporaryFile, TemporaryDirectory

import msgpack
import networkx as nx
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

printer = PrettyPrinter()


def list_to_fasta(lst, outfile):
    seq_records = [
        SeqRecord(id=f"SEQ{idx}", seq=Seq(val), description="")
        for idx, val in enumerate(lst)
    ]
    SeqIO.write(seq_records, outfile, "fasta")


def makeblastdb(infile, db):
    command = [
        "makeblastdb",
        "-in",
        f"{infile.absolute()}",
        "-dbtype",
        "prot",
        "-out",
        f"{db.absolute()}",
        "-hash_index",
    ]
    subprocess.call(command)


def run_blastp(infile, db, outfile):
    command = [
        "blastp",
        "-query",
        f"{infile.absolute()}",
        "-db",
        f"{db.absolute()}",
        "-evalue",
        "1e-3",
        "-out",
        f"{outfile.absolute()}",
        "-outfmt",
        "6",
        "-num_threads",
        f"{cpu_count()}",
    ]
    subprocess.call(command)


def run_needleall(tup):
    aseqs, bseqs = tup

    # Write to file.
    a = NamedTemporaryFile("wt+", suffix=".fasta")
    b = NamedTemporaryFile("wt+", suffix=".fasta")
    SeqIO.write(aseqs, a, "fasta")
    SeqIO.write(bseqs, b, "fasta")

    result_file = NamedTemporaryFile("wt+")

    a.seek(0)
    b.seek(0)
    # Run search.
    command = [
        "needleall",
        "-asequence",
        a.name,
        "-bsequence",
        b.name,
        "-gapopen",
        "10",
        "-gapextend",
        "0.5",
        "-outfile",
        result_file.name,
        "-aformat3",
        "pair",
    ]
    subprocess.call(command)

    # Parse output.
    identities = []
    with open(result_file.name, "r") as f:
        id1, id2, identity = None, None, None
        for line in f:
            if line.strip().startswith("# 1:"):
                id1 = line.strip().split()[2]
            elif line.strip().startswith("# 2:"):
                id2 = line.strip().split()[2]
            elif line.strip().startswith("# Identity"):
                num, den = line.strip().split()[2].split("/")
                identity = int(num) / int(den)

                if identity > 0.2:
                    tup = (id1, id2, identity)
                    identities.append(tup)
    # Tidy up
    a.close()
    b.close()
    result_file.close()

    # Return identities
    return identities


def genreate_needleall_jobs(network: nx.DiGraph, fasta_file):
    for node in network:
        neighbours = list(network.neighbors(node))
        if not neighbours:
            continue
        a_records = [
            i for i in SeqIO.parse(fasta_file, "fasta") if i.id == node
        ]
        b_records = [
            i for i in SeqIO.parse(fasta_file, "fasta") if i.id in neighbours
        ]
        yield a_records, b_records


def run(input):
    with open(input, "rb") as f:
        d = msgpack.unpackb(f.read(), raw=False, use_list=False)

    d = {k: v for k, v in d.items() if any(v)}

    # Proteins
    proteins = set(chain(*d.values()))
    if None in proteins:
        proteins.remove(None)
    proteins = sorted(proteins)

    # Write Proteins to a temproary FASTA file.
    directory = TemporaryDirectory()
    fasta_file = Path(directory.name) / "input.fasta"
    list_to_fasta(proteins, fasta_file)

    blastdb = Path(directory.name) / "db"
    makeblastdb(fasta_file, blastdb)

    blasttab = Path(directory.name) / "results.tab"
    run_blastp(fasta_file, blastdb, blasttab)

    # Parse the results of the BLAST search, to build a graph indicating the
    # searches that need doing with needleall.
    G = nx.DiGraph()

    for record in SeqIO.parse(fasta_file, "fasta"):
        G.add_node(record.id)

    for query in SearchIO.parse(blasttab, "blast-tab"):
        for hit in query:
            if hit.id == query.id:
                continue
            if not G.has_edge(hit.id, query.id):
                G.add_edge(query.id, hit.id)

    # Run needleall searches to construct a new graph.
    H = nx.Graph()
    for idx, val in enumerate(proteins):
        H.add_node(f"SEQ{idx}", sequence=val, graphics={"fill": "#50C878"})

    P = Pool(cpu_count())
    for result in P.imap_unordered(
        run_needleall, genreate_needleall_jobs(G, fasta_file)
    ):
        for id1, id2, identity in result:
            H.add_edge(id1, id2, identity=identity)
    del G

    # Compound processing.
    compounds = sorted({k[0] for k in d.keys()})
    for idx, val in enumerate(compounds):
        H.add_node(f"CPD{idx}", smiles=val, graphics={"fill": "#0F52BA"})

    fps = {
        cpd: AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(cpd), 2, nBits=16384
        )
        for cpd in compounds
    }

    for cpd1, cpd2 in combinations(fps.keys(), 2):
        fp1 = fps[cpd1]
        fp2 = fps[cpd2]
        diff = DataStructs.TanimotoSimilarity(fp1, fp2)
        H.add_edge(
            f"CPD{compounds.index(cpd1)}",
            f"CPD{compounds.index(cpd2)}",
            similarity=diff,
        )

    for (cpd, _), pros in d.items():
        for pro in pros:
            if not pro:
                continue
            H.add_edge(
                f"CPD{compounds.index(cpd)}", f"SEQ{proteins.index(pro)}"
            )

    nx.write_gml(H, f"{input}.gml")

    # Clear up
    directory.cleanup()

if __name__ == "__main__":
    import os
    files = [i for i in os.listdir() if i.startswith("2020_05_11")]
    for f in files:
        run(f)
