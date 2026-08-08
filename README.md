# PageTrust for retracted papers
This repository contains code to recreate the figure from our paper "Retraction without Traction: Towards Handling Post-publication Invalidation in Open Science Infrastructure". To run the code you will need a python environment with networkx, numpy, matplotlib, and scipy installed.

`pagetrust.py` contains our implementation of the PageTrust algorithm based on de Kerchove & Van Dooren (2008); this algorithm should function with any directed networkx graph, however, it has not been tested extensively and is computationally expensive.

`subgraph.json` contains a subgraph of OpenAIRE centering on one retracted paper (Gautret et al., 2020) and a random snowball of references and citations.

`pagetrust_example.py` contains an example of how to apply the PageTrust algorithm to the subgraph. The visualization created is the figure found in our paper, written to `figure.png`. Node positions are computed once and cached in `positions.json` so that repeated runs reproduce the same layout; delete that file to generate a new one.



