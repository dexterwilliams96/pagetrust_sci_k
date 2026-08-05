# Sample code for our paper TITLE DBD
This repository contains code to recreate the figure from our paper "TITLE TBD". To run the code you will need a python environment with networkx, numpy, matplotlib, and scipy installed.

`pagetrust.py` contains our implementation of the PageTrust algorithm; this algorithm should function with any directed networkx graph, however, it has not been tested extensively and is computationally expensive.

`subgraph.json` contains a subgraph of OpenAIRE centering on one retracted paper (Gautret et al., 2020) and a random snowball of references and citations.

`pagetrust_example.py` contains an example of how to apply the PageTrust algorithm to the subgraph. The visualization created is the figure found in our paper.

# Reference
If you use this code in your research, please cite our paper:

```
TODO
```



