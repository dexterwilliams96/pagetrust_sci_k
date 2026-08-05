import json
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
    scores = {node: rank for rank, (node, _) in enumerate(sorted(scores.items(), key=lambda item: item[1], reverse=False), start=1)}

    pos = nx.spring_layout(G, seed=42)
    pos[origin_doi] = [0, 0]
    node_colors = [scores.get(node, 0.0) for node in G.nodes() if node != origin_doi and node not in forced_zero_nodes]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[origin_doi],
        node_color="red",
        node_size=500,
        node_shape="s",
    )
    # Draw nodes in rings around the origin node based on how many hops away they are from the origin node
    # Nodes should be positioned in rings, color should be based on score
    # The bottom half of the ring should be formed of all backlinks and the top half of the ring should be formed of forward links
    # ALL predecessors of origin (recursive
    backlinks = list(nx.ancestors(G, origin_doi))
    # ALL successors of origin (recursive)
    forwardlinks = list(nx.descendants(G, origin_doi))
    print(len(backlinks), len(forwardlinks))
    for node in backlinks:
        hops = nx.shortest_path_length(G.to_undirected(), source=origin_doi, target=node)
        angle = 3.14159 * (list(G.nodes()).index(node) / G.number_of_nodes())
        radius = hops * 0.5
        pos[node] = [radius * np.cos(angle), radius * np.sin(angle)]
    for node in forwardlinks:
        hops = nx.shortest_path_length(G.to_undirected(), source=origin_doi, target=node)
        angle = 3.14159 + 3.14159 * (list(G.nodes()).index(node) / G.number_of_nodes())
        radius = hops * 0.5
        pos[node] = [radius * np.cos(angle), radius * np.sin(angle)]
    # For all other nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[node for node in G.nodes() if node != origin_doi and node not in forced_zero_nodes],
        node_color=node_colors,
        cmap=plt.cm.viridis_r,
        node_size=50,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=forced_zero_nodes,
        node_color="red",
        node_size=50,
    )
    # Draw a horizontal line through the center of the figure to separate the backlinks from the forward links
    plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
    plt.axis("off")
    plt.show()

if __name__ == "__main__":
    main()
