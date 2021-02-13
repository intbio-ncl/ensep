from pathlib import Path
import numpy as np
import networkx as nx
from collections import defaultdict

def quantify_connections_graph(input):
    G = nx.read_gml(input)
    edges = [sorted([i,j]) for i, j in G.edges() if (i + j).count("CPD") == 1]
    # Compound connections
    cpd_connections = defaultdict(int)
    pro_connections = defaultdict(int)

    for i, j in edges:
        cpd_connections[i] += 1
        pro_connections[j] += 1

    print(f"Input: {input}")
    print("Compounds:")
    print(f"\tMax: {max(cpd_connections.values())}")
    print(f"\tMin: {min(cpd_connections.values())}")
    print(f"\tQuartiles: {np.percentile(list(cpd_connections.values()), [25, 50, 75])}")
    print("Proteins:")
    print(f"\tMax: {max(pro_connections.values())}")
    print(f"\tMin: {min(pro_connections.values())}")
    print(f"\tQuartiles: {np.percentile(list(pro_connections.values()), [25, 50, 75])}")
    print()

p = Path("./results")

graphs = ["2020_05_11_icrc.msgpack.gml", "2020_05_11_ecrc.msgpack.gml", "2020_05_11_phenol_hydrox.msgpack.gml"]
for graph in graphs:
    quantify_connections_graph(p / graph)
