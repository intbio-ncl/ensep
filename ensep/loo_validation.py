from functools import partial
import multiprocessing
import random
from itertools import product
import os
import pdb

import matplotlib.pyplot as plt
import mcpt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from tqdm import tqdm

# logger.remove()
logfile = open("validation.log", "a")
logger.add(logfile, colorize=True)


def shuffle_labels(G):
    pros = [node for node in G.nodes() if node.startswith("SEQ")]
    new_pros = pros[:]
    random.shuffle(new_pros)
    pro_map = {i: j for i, j in zip(pros, new_pros)}

    cpds = [node for node in G.nodes() if node.startswith("CPD")]
    new_cpds = cpds[:]
    random.shuffle(new_cpds)
    cpd_map = {i: j for i, j in zip(cpds, new_cpds)}

    # CPD - Pro edges
    edges = []
    for edge in G.edges():
        i, j = sorted(edge)
        if i.startswith("CPD") and j.startswith("SEQ"):
            edges.append((i, j))

    new_edges = [(cpd_map[i], pro_map[j]) for i, j in edges]

    G.remove_edges_from(edges)
    G.add_edges_from(new_edges)

    return G


def get_similarity_of_most_similar_cpd(cpd, G):
    cpd_neighbours = [
        neighbour
        for neighbour in G.neighbors(cpd)
        if neighbour.startswith("CPD") and neighbour != cpd
    ]

    if cpd_neighbours:
        closest_cpd_neighbour = max(
            cpd_neighbours, key=lambda i: G.get_edge_data(cpd, i)["similarity"]
        )
        closest_cpd_similarity = G.get_edge_data(cpd, closest_cpd_neighbour)[
            "similarity"
        ]
    else:
        closest_cpd_similarity = np.NaN

    return closest_cpd_similarity


def get_identity_of_most_similar_seq(seq, G):
    seq_neighbours = [
        neighbour
        for neighbour in G.neighbors(seq)
        if neighbour.startswith("SEQ") and neighbour != seq
    ]
    if seq_neighbours:
        closest_seq_neighbour = max(
            seq_neighbours, key=lambda i: G.get_edge_data(seq, i)["identity"]
        )
        closest_seq_similarity = G.get_edge_data(seq, closest_seq_neighbour)[
            "identity"
        ]
    else:
        closest_seq_similarity = np.NaN

    return closest_seq_similarity


def test_recovery_negative_edges(G, cpd_clusters, seq_clusters):
    logger.info("Testing recover of negative edges")

    # Get all true edges
    true_edges = [sorted(edge) for edge in G.edges]
    # Sort edges, and get those that are between compounds and sequences
    true_edges = [
        (i, j)
        for i, j in true_edges
        if i.startswith("CPD") and j.startswith("SEQ")
    ]
    logger.info(f"There are {len(true_edges)} known edges")

    # Get the compounds and protein sequences
    seqs = [node for node in G.nodes() if node.startswith("SEQ")]
    cpds = [node for node in G.nodes() if node.startswith("CPD")]

    # Create all false edges
    false_edges = [
        (cpd, seq)
        for cpd, seq in product(cpds, seqs)
        if not G.has_edge(seq, cpd)
    ]

    # If there are > 25,000 false edges, take a sample
    if len(false_edges) > 25_000:
        false_edges = random.sample(false_edges, 25_000)

    # Make a map of edge to recovery density, setting all at infinity
    # initially
    densities = {edge: np.inf for edge in false_edges}

    # For all possible combinations of clustering thresholds
    for cpd_cluster, seq_cluster in product(
        cpd_clusters.values(), seq_clusters.values()
    ):
        # Create a graph, and add the clusters as nodes.
        H = nx.Graph()
        H.add_nodes_from(set(seq_cluster.values()))
        H.add_nodes_from(set(cpd_cluster.values()))

        # For each edge in the true set, add an edge between clusters
        # containing the member
        for i, j in true_edges:
            H.add_edge(cpd_cluster[i], seq_cluster[j])

        # Calculate the density of the bipartite configuration
        density = nx.algorithms.bipartite.density(H, set(cpd_cluster.values()))

        # For each false edge, if the clusters containing the members are
        # connected, then update the recovery density if it is lower than
        # the current recorded recovery density
        for edge in false_edges:
            i, j = edge
            if H.has_edge(cpd_cluster[i], seq_cluster[j]):
                if densities[edge] > density:
                    densities[edge] = density

    densities = {k: np.log10(v) for k, v in densities.items()}
    return densities


def obtain_true_edges(G):
    true_edges = [sorted(edge) for edge in G.edges()]
    true_edges = [
        (i, j)
        for i, j in true_edges
        if i.startswith("CPD") and j.startswith("SEQ")
    ]

    return true_edges


def test_recovery_true_edge(
    true_edge, G, cpd_clusters, seq_clusters, true_edges
):
    true_cpd, true_seq = true_edge

    # Assume the lowest recovery density is infinity
    lowest_density = np.inf

    # Generate some statistics

    closest_cpd_similarity = get_similarity_of_most_similar_cpd(true_cpd, G)
    closest_seq_similarity = get_identity_of_most_similar_seq(true_seq, G)

    for cpd_clustering, seq_clustering in product(
        cpd_clusters.values(), seq_clusters.values()
    ):
        # Construct the bipartite graph
        H = nx.Graph()
        H.add_nodes_from(cpd_clustering.values())
        H.add_nodes_from(seq_clustering.values())

        for edge in G.edges():
            edge = sorted(edge)
            # Skip the edge if it's the same as the edge we're trying to
            # recover
            if tuple(edge) == true_edge:
                continue

            i, j = edge
            if i.startswith("CPD") and j.startswith("SEQ"):
                H.add_edge(cpd_clustering[i], seq_clustering[j])

        # Calculate the density of the network
        density = nx.algorithms.bipartite.density(
            H, set(cpd_clustering.values())
        )

        # If the clusters containing the target sequence and cluster are
        # connected, then update the density if it's lower than the
        # current recorded density.
        target_edge = (cpd_clustering[true_cpd], seq_clustering[true_seq])
        if density < lowest_density and H.has_edge(*target_edge):
            lowest_density = density

    if np.isinf(lowest_density):
        lowest_density = 1

    return {
        "log10 recovery density": np.log10(lowest_density)
        if not np.isnan(lowest_density)
        else lowest_density,
        "compound degree": sum(
            1 for node in G.neighbors(true_edge[0]) if node.startswith("SEQ")
        ),
        "protein degree": sum(
            1 for node in G.neighbors(true_edge[1]) if node.startswith("CPD")
        ),
        "closest sequence similarity": closest_seq_similarity,
        "closest compound similarity": closest_cpd_similarity,
    }


def test_recovery_true_edges(G, cpd_clusters, seq_clusters, true_edges=None):
    if not true_edges:
        logger.info("Obtaining a selection of known edges")
        # Get all true edges
        true_edges = obtain_true_edges(G)
        # If there are > 25,000 false edges, take a sample
        if len(true_edges) > 25_000:
            true_edges = random.sample(true_edges, 25_000)
    else:
        logger.info("Known passed, not generating.")

    # Hash map to store recovery densities for true edges
    densities = {}

    fx = partial(
        test_recovery_true_edge,
        G=G,
        cpd_clusters=cpd_clusters,
        seq_clusters=seq_clusters,
        true_edges=true_edges,
    )

    with multiprocessing.Pool(multiprocessing.cpu_count()) as P:
        for edge, result in tqdm(
            zip(true_edges, P.imap(fx, true_edges)),
            total=len(true_edges),
            leave=False,
        ):
            densities[edge] = result

    return densities


def run_wpgma(G, nodes, *, similarity_key, threshold, cluster_label):
    # Induce subgraph
    seqs_sg = G.subgraph(nodes)
    # Create similarity matrix
    similarity_mat = nx.to_numpy_matrix(
        seqs_sg, nodelist=nodes, weight=similarity_key
    )
    distance_mat = 1 - similarity_mat

    compressed = squareform(distance_mat)
    linkage_mat = linkage(compressed, method="weighted")
    cluster_assignment = fcluster(
        linkage_mat, t=threshold, criterion="distance"
    )
    clusters = {
        node: f"{cluster_label}_{cluster}"
        for node, cluster in zip(nodes, cluster_assignment)
    }

    return clusters


def main():
    files = [
        f"results/{i}"
        for i in (
            "2020_05_11_icrc.msgpack.gml",
            "2020_05_11_ecrc.msgpack.gml",
            "2020_05_11_phenol_hydrox.msgpack.gml",
            "2020_05_11_o2a.msgpack.gml",
        )
    ]
    files.sort(key=lambda i: os.path.getsize(i))

    for infile in files:
        logger.info(f"Reading GML file: {infile}")

        G = nx.read_gml(infile)
        true_edges = obtain_true_edges(G)
        print(f"There are {len(true_edges)} true edges")
        # G = shuffle_labels(G)

        cpds = [node for node in G.nodes if node.startswith("CPD")]
        seqs = [node for node in G.nodes if node.startswith("SEQ")]
        print(f"There are {len(cpds)} compounds")
        print(f"There are {len(seqs)} proteins")

        [G.add_edge(i, i, identity=1) for i in seqs]
        [G.add_edge(i, i, similarity=1) for i in cpds]

        logger.info("Generating compound clusters")
        cpd_clusters = {
            round(t, 1): run_wpgma(
                G,
                cpds,
                similarity_key="similarity",
                threshold=t,
                cluster_label="CPD_CLUSTER",
            )
            for t in np.linspace(0, 1, 11)
        }

        logger.info("Generating protein clusters")
        seq_clusters = {
            round(t, 1): run_wpgma(
                G,
                seqs,
                similarity_key="identity",
                threshold=t,
                cluster_label="SEQ_CLUSTER",
            )
            for t in np.linspace(0, 0.7, 8)
        }

        neg_dens = test_recovery_negative_edges(G, cpd_clusters, seq_clusters)
        pos_dens = test_recovery_true_edges(
            G, cpd_clusters, seq_clusters, true_edges=true_edges
        )

        min_density = min(
            list(neg_dens.values()) +
            [i["log10 recovery density"] for i in pos_dens.values()]
        )

        # Create a data frame.
        true_df = pd.DataFrame(pos_dens.values())
        logger.warning(f"There are {len(true_df[true_df['log10 recovery density'] == np.log10(2)])} irrecoverable edges.")

        # ROC curve.
        thresholds = set(neg_dens.values())
        thresholds = thresholds | set(true_df["log10 recovery density"])
        thresholds = sorted(thresholds)

        x, y = [0], [0]
        vline = -float("inf")

        for threshold in thresholds:
            # True-positive rate is TP / (TP + FN) -- the proportion of truths correctly called true.
            tpr = len(true_df[ true_df["log10 recovery density"] <= threshold]) / len(true_df)
            # False-positive rate is FP / (TN + FP) -- the proportion of falsehoods incorrectly called true.
            fpr = sum(1 for v in neg_dens.values() if v <= threshold) / len(neg_dens)
            if threshold <= -1:
                vline = fpr
            y.append(tpr)
            x.append(fpr)
        x.append(1)
        y.append(1)

        auc = np.trapz(y, x=x)
        print(auc)

        plt.plot(x, y, marker=".", label="Recovery density")
        plt.plot([0,1], [0,1], linestyle="--", label="Random")
        if vline != -float("inf"):
            plt.axvline(vline, color="gray", linestyle="--", alpha=0.3)
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(f"Reciever operating characteristic curve comparing\nrecovery density to a random search\n(AUC={auc:.2f})\n")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_ROC.png", dpi=1000)


        one_degree_protein = len(true_df[true_df["protein degree"] == 1])
        one_degree_compound = len(true_df[true_df["compound degree"] == 1])
        one_degree_both = len(true_df[(true_df["compound degree"] == 1) & (true_df["protein degree"] == 1)])
        logger.info(f"There are {one_degree_protein} edges with protein degree == 1")
        logger.info(f"There are {one_degree_compound} edges with compound degree == 1")
        logger.info(f"There are {one_degree_both} edges with both compound and protein degree == 1")

        # There is no such thing as an "irrecoverable edge", because at worst, the
        # recovery density should be 0 (when thresholds are 0, and everything is
        # lumped together in one cluster on each side).

        # Make scatter plot of recovery density vs the closest compound of a pair.
        points = plt.scatter(
            true_df["closest compound similarity"],
            true_df["closest sequence similarity"],
            c=true_df["log10 recovery density"],
            cmap="Spectral",
            edgecolor='black',
            linewidth=0.1,
            marker=".",
            s=50,
        )
        plt.xlim(-0.05, 1.05)
        plt.ylim(-0.05, 1.05)

        cbar = plt.colorbar(points)
        cbar.ax.set_ylabel("$log_{10}$ recovery density")

        plt.title(
            "Distribution of known edges, $C_n-P_n$,\n"
            "according to protein identity and compound similarity\n"
        )
        plt.xlabel("Similarity with most similar compound to $C_n$")
        plt.ylabel("Identity with most similar protein to $P_n$")
        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig1.png", dpi=1000)
        plt.close()

        # Add an extra column that marks whether the compound degree is greater than one.
        true_df["compound degree > 1"] = np.where(
            true_df["compound degree"] > 1, "> 1", "1"
        )
        true_df["protein degree > 1"] = np.where(
            true_df["protein degree"] > 1, "> 1", "1"
        )


        # Make violin plot of the recovery density.
        sns.violinplot(
            x="compound degree > 1",
            y="log10 recovery density",
            data=true_df,
            bw=0.1,
        )
        plt.title(
            "Violin plot of recovery densities by number of protein\n"
            "neighbours of compound, $C_n$, in known pairs, $C_n-P_n$\n"
        )
        plt.xlabel("Number of protein neighbours of compound, $C_n$")
        plt.ylabel("$log_{10}$ recovery density")
        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig2.png", dpi=1000)
        plt.close()

        # Make violin plot of the recovery density.
        sns.violinplot(
            x="protein degree > 1",
            y="log10 recovery density",
            data=true_df,
            bw=0.1,
        )
        plt.title(
            "Violin plot of recovery densities by number of compound\n"
            "neighbours of protein, $P_n$, in known pairs, $C_n-P_n$\n"
        )
        plt.xlabel("Number of compound neighbours of protein, $P_n$")
        plt.ylabel("$log_{10}$ recovery density")
        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig3.png", dpi=1000)
        plt.close()

        # Make scatter plot.
        coef, p = spearmanr(
            true_df["compound degree"], true_df["log10 recovery density"]
        )
        sns.scatterplot(
            x="compound degree",
            y="log10 recovery density",
            data=true_df,
            marker=".",
            s=6,
        )
        log10_str = "$log_{10}$"
        plt.xlabel("Number of protein neighbours of compound, $C_n$")
        plt.ylabel("$log_{10}$ recovery density")

        plt.title(
            f"Scatter plot of {log10_str} recovery density\n"
            "by the number of protein neighbours of compound, $C_n$,\n"
            "in known pairs, $C_n-P_n$"
        )
        props = dict(boxstyle="round", alpha=0.5, facecolor="white")
        ax = plt.gca()
        ax.text(
            0.8,
            0.95,
            f"$p$ = {coef:.2f}\np = {p:.2f}",
            bbox=props,
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=10,
        )

        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig4.png", dpi=1000)
        plt.close()

        # Distribution plots.
        sns.distplot(
            [i for i in neg_dens.values()],
            label="Unknown edges",
            kde=False, bins=np.linspace(min_density, 0, 150),
            color="red"
        )

        plt.xlabel("$log_{10}$ recovery density")
        plt.ylabel("Frequency")
        plt.title(
            "Histogram showing distributions\nof recovery density by edge type\n"
        )
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig5a.png", dpi=1000)
        plt.close()

        sns.distplot(
            [
                i["log10 recovery density"]
                for _, i in true_df.iterrows()
                if not np.isnan(i["log10 recovery density"])
                and i["compound degree"] == 1
            ],
            label="Known, compound degree = 1",
            color="blue",
            kde=False, bins=np.linspace(min_density, 0, 150),
        )

        sns.distplot(
            [
                i["log10 recovery density"]
                for _, i in true_df.iterrows()
                if not np.isnan(i["log10 recovery density"])
                and i["compound degree"] > 1
            ],
            label="Known, compound degree > 1",
            color="green",
            kde=False, bins=np.linspace(min_density, 0, 150),
        )
        plt.xlabel("$log_{10}$ recovery density")
        plt.ylabel("Frequency")
        plt.title(
            "Histogram showing distributions\nof recovery density by edge type\n"
        )
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig5b.png", dpi=1000)
        plt.close()

        logger.info("Recovery density comparisons of edges by compound degree:")

        logger.debug(f'Median recovery density for degree = 1: {np.median(true_df[true_df["compound degree"] == 1]["log10 recovery density"])}')
        logger.debug(f'Median recovery density for degree > 1: {np.median(true_df[true_df["compound degree"] > 1]["log10 recovery density"])}')
        logger.debug(f"Median recovery density for false edges: {np.median(list(neg_dens.values()))}")

        logger.debug("Degree == 1 vs negative edges")
        r = mcpt.permutation_test(
            list(neg_dens.values()),
            list(
                true_df[true_df["compound degree"] == 1]["log10 recovery density"]
            ),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")

        logger.debug("Degree > 1 vs negative edges")
        r = mcpt.permutation_test(
            list(neg_dens.values()),
            list(
                true_df[true_df["compound degree"] > 1]["log10 recovery density"]
            ),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")

        logger.debug("Degree > 1 vs Degree == 1")
        r = mcpt.permutation_test(
            list(
                true_df[true_df["compound degree"] == 1]["log10 recovery density"]
            ),
            list(
                true_df[true_df["compound degree"] > 1]["log10 recovery density"]
            ),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")

        # Make scatter plot.
        coef, p = spearmanr(
            true_df["protein degree"], true_df["log10 recovery density"]
        )
        sns.scatterplot(
            x="protein degree", y="log10 recovery density", data=true_df
        )

        log10_str = "$log_{10}$"
        plt.xlabel("Number of compound neighbours of protein, $P_n$")
        plt.ylabel("$log_{10}$ recovery density")
        plt.title(
            f"Scatter plot of {log10_str} recovery density\n"
            "by the number of compound neighbours of protein, $P_n$,\n"
            "in known pairs, $C_n-P_n$"
        )
        props = dict(boxstyle="round", alpha=0.5, facecolor="white")
        ax = plt.gca()
        ax.text(
            0.8,
            0.95,
            f"$p$ = {coef:.2f}\np = {p:.2f}",
            bbox=props,
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=10,
        )
        plt.tight_layout()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig6.png", dpi=1000)
        plt.close()

        # Distribution plots.
        sns.distplot(
            [i for i in neg_dens.values()],
            label="Unknown edges",
            kde=False, bins=np.linspace(min_density, 0, 150),
            color="red",
        )

        plt.xlabel("$log_{10}$ recovery density")
        plt.ylabel("Frequency")
        plt.title(
            "Histogram showing distributions\nof recovery density by edge type\n"
        )
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig7a.png", dpi=1000)
        plt.close()

        sns.distplot(
            [
                i["log10 recovery density"]
                for _, i in true_df.iterrows()
                if not np.isnan(i["log10 recovery density"])
                and i["protein degree"] == 1
            ],
            label="Known, protein degree = 1",
            kde=False, bins=np.linspace(min_density, 0, 150),
            color="blue",
        )

        sns.distplot(
            [
                i["log10 recovery density"]
                for _, i in true_df.iterrows()
                if not np.isnan(i["log10 recovery density"])
                and i["protein degree"] > 1
            ],
            label="Known, protein degree > 1",
            kde=False, bins=np.linspace(min_density, 0, 150),
            color="green",
        )

        plt.xlabel("$log_{10}$ recovery density")
        plt.ylabel("Frequency")
        plt.title(
            "Histogram showing distributions\nof recovery density by edge type\n"
        )
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{infile.rsplit('.', 1)[0]}_Fig7b.png", dpi=1000)
        plt.close()

        logger.info("Recovery density comparisons of edges by protein degree:")

        logger.debug(f'Median recovery density for degree = 1: {np.median(true_df[true_df["protein degree"] == 1]["log10 recovery density"])}')
        logger.debug(f'Median recovery density for degree > 1: {np.median(true_df[true_df["protein degree"] > 1]["log10 recovery density"])}')
        logger.debug(f"Median recovery density for false edges: {np.median(list(neg_dens.values()))}")

        logger.debug("Degree == 1 vs negative edges")
        r = mcpt.permutation_test(
            list(neg_dens.values()),
            list(
                true_df[true_df["protein degree"] == 1]["log10 recovery density"]
            ),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")

        logger.debug("Degree > 1 vs negative edges")
        r = mcpt.permutation_test(
            list(neg_dens.values()),
            list(true_df[true_df["protein degree"] > 1]["log10 recovery density"]),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")

        logger.debug("Degree > 1 vs Degree == 1")
        r = mcpt.permutation_test(
            list(
                true_df[true_df["protein degree"] == 1]["log10 recovery density"]
            ),
            list(true_df[true_df["protein degree"] > 1]["log10 recovery density"]),
            "median",
            side="greater",
            n=1_000,
        )
        logger.debug(f"{r}")
        logger.info("Done!")

# pdb.set_trace()
main()
