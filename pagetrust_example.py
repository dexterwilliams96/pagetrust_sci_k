import json
import matplotlib.pyplot as plt
import networkx as nx

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
    node_colors = [scores.get(node, 0.0) for node in G.nodes()]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[origin_doi],
        node_color="red",
        node_size=500,
        node_shape="s",
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=[node for node in G.nodes() if node != origin_doi],
        node_color=[scores.get(node, 0.0) for node in G.nodes() if node != origin_doi],
        cmap=plt.cm.viridis_r,
        node_size=200,
    )
    plt.axis("off")
    plt.show()

if __name__ == "__main__":
    main()
