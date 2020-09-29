# DCAlign
# Aligning biological sequences using Direct Coupling Analysis models

# Authors

Anna Paola Muntoni, Andrea Pagnani, Martin Weigt and Francesco Zamponi

# Description

This package contains the implementation of DCAlign, a tool for aligning biological sequences to a DCA model of a seed, discussed in [arXiv:2005.08500](https://arxiv.org/abs/2005.08500). DCAlign requires, together with the set of sequences to be aligned, a model for the seed composed of:

- a DCA model;
- a set of penalties for opening and extending insertions; they can be site-dependent;
- a penalty for adding a gap at the beginning or the end of the aligned sequence ("external") or between two symbols ("internal").

We provide in the ```test``` folder, the parameters of the DCA models, the gap/insertion penalties used to align the [Pfam](https://pfam.xfam.org/) and [Rfam](https://rfam.xfam.org/) families PF00035, PF00677, PF00684, PF00763, RF00059, RF00162, RF00167, RF01734. For an arbitrary seed, we provide the sub-module [DCAbuild](https://github.com/anna-pa-m/DCAbuild) which performs the learning of all the necessary parameters (see DCAbuild page for the documentation).

The code is written in [Julia](https://julialang.org/).

# Usage

A step-by-step description of the alignment procedure of a single sequence is provided in the Jupyter notebook in ```notebook```. 
