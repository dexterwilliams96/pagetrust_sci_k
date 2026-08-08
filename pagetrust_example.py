import json
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from pagetrust import pagetrust


def main():
    with open("subgraph.json", "r") as f:
        graph = json.load(f)
    origin_doi = "10.1016/j.ijantimicag.2020.105949"
    origin_node = next(node for node in graph["nodes"] if node["doi"] == origin_doi)
    G = nx.DiGraph()
    for edge in graph["edges"]:
        G.add_edge(edge[0], edge[1], weight=edge[2])
    print(f"Subgraph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    forced_zero_nodes = [
        node["doi"] for node in graph["nodes"] if node.get("retracted", False)
    ]
    scores, _ = pagetrust(
        G, alpha=0.85, M=1, beta=1.0, forced_zero_nodes=forced_zero_nodes
    )
    print("PageTrust Ranks:")
    for node, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        print(f"Node {node}: {score:.4f}")
    # Change scores variable into a ranking of scores for visualization
    scores = {
        node: rank
        for rank, (node, _) in enumerate(
            sorted(scores.items(), key=lambda item: item[1], reverse=False), start=1
        )
    }

    # Sorting keeps angles stable across runs (ancestors/descendants are sets)
    backlinks = sorted(nx.ancestors(G, origin_doi))
    forwardlinks = sorted(nx.descendants(G, origin_doi))
    print(len(backlinks), len(forwardlinks))
    pos_cache = "positions.json"
    if os.path.exists(pos_cache):
        with open(pos_cache, "r") as f:
            pos = {node: np.array(xy) for node, xy in json.load(f).items()}
    else:
        rng = np.random.default_rng(0)
        pos = nx.spring_layout(G, seed=0)
        pos[origin_doi] = [0, 0]
        for i, node in enumerate(backlinks):
            hops = nx.shortest_path_length(
                G.to_undirected(), source=origin_doi, target=node
            )
            angle = 3.14159 * (i / len(backlinks))
            radius = hops * 0.25
            radius += rng.uniform(-0.02, 0.02)
            pos[node] = [radius * np.cos(angle), radius * np.sin(angle)]
        for i, node in enumerate(forwardlinks):
            hops = nx.shortest_path_length(
                G.to_undirected(), source=origin_doi, target=node
            )
            angle = 3.14159 + 3.14159 * (i / len(forwardlinks))
            radius = hops * 0.25
            radius += rng.uniform(-0.02, 0.02)
            pos[node] = [radius * np.cos(angle), radius * np.sin(angle)]
        with open(pos_cache, "w") as f:
            json.dump(
                {node: [float(x), float(y)] for node, (x, y) in pos.items()}, f
            )

    plt.figure(figsize=(9, 7))
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[origin_doi],
        node_color="red",
        node_size=200,
        node_shape="o",
        edgecolors="black",
        linewidths=0.8,
    )
    general_nodes = [
        node
        for node in G.nodes()
        if node != origin_doi and node not in backlinks and node not in forwardlinks
    ]
    linked_nodes = backlinks + forwardlinks
    linked_nodes = [node for node in linked_nodes if node not in forced_zero_nodes]
    node_colors = [scores[node] for node in linked_nodes]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=linked_nodes,
        node_color=node_colors,
        cmap=plt.cm.viridis_r,
        node_size=50,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[
            node
            for node in forced_zero_nodes
            if node not in general_nodes and node != origin_doi
        ],
        node_color="red",
        node_size=50,
    )
    plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax = plt.gca()
    ax.text(
        0.02,
        0.98,
        "papers that cite Gautret et al. and their descendants",
        transform=ax.transAxes,
        fontstyle="italic",
        color="gray",
        ha="left",
        va="top",
    )
    ax.text(
        0.02,
        0.02,
        "papers Gautret et al. cites and their ancestors",
        transform=ax.transAxes,
        fontstyle="italic",
        color="gray",
        ha="left",
        va="bottom",
    )
    plt.axis("off")
    plt.savefig("figure.png", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
