# DCAlign: Aligning biological sequences using Direct Coupling Analysis models

# Authors

Anna Paola Muntoni, Andrea Pagnani, Martin Weigt and Francesco Zamponi

# Description

This package contains the implementation of `DCAlign`, a tool for aligning
biological sequences to a DCA model of a seed alignment, discussed in [Aligning
biological sequences by exploiting residue conservation and
coevolution](https://link.aps.org/doi/10.1103/PhysRevE.102.062409).

The latest version is described in [DCAlign v1.0: Aligning biological sequences
using co-evolution models and informed priors](https://doi.org/10.1093/bioinformatics/btad537)

## Latest release v1.0

`DCAlign` requires, together with the set of sequences to be aligned, a model for the seed composed of:

- a DCA model;
- an empirical prior modeling insertions and deletions.

The learning of the model has been integrated in `v1.0` (see `Usage` and the notebook in the folder `notebook`) and requires a few seconds of computation. A description of the new release is available at [bioRxiv](https://biorxiv.org/cgi/content/short/2022.05.18.492471v1). 


The code is written in [Julia](https://julialang.org/).

## Installation

The repository `DCAlign` has not yet been registered, but it can be installed in two ways:

1. Within the Julia PackageManager (activate it using the `]` key), type

    ```(@v1.?) pkg> add https://github.com/infernet-h2020/DCAlign #v1.0```

    To use the package in the Julia terminal use

    ```julia> using DCAlign```

2. First clone the package 

    ```git clone  https://github.com/infernet-h2020/DCAlign```

    and set the `v1.0` branch

    ```git checkout v1.0```

    Then, instantiate the package from the Julia REPL

    ```julia> using Pkg; Pkg.instantiate()```
    
This second option allows one to easily use the explicatory notebook. 


# Usage

A step-by-step description of the alignment procedure of a single sequence is provided in the Jupyter notebook in the folder `notebook`. 

As an example, we learn a seed model for Pfam PF00035, PF00684 or Rfam RF00167 whose seed sequences are present in the folder `seeds` in Stockholm format. We then try to align one of the available query sequences stored in the folder `test` and we compare our aligned sequence to the alignment performed by HMMER or Infernal.

### Algorithm parameters

The last implementation of `DCAlign` performs an annealing scheme of the update equations by increasing the value of `β` within the iterations. It stops when the approximate marginal probabilities are sufficiently concentrated. To tune the annealing scheme and the dispersion of the marginals, one can modify the increment `Δβ` (default: 0.05), the frequency of the annealing, performed each `Δt` (default: 10) iterations, and the threshold probability `thP` (default: 0.30). See [bioRxiv](https://biorxiv.org/cgi/content/short/2022.05.18.492471v1) for further details.

### Limitations

For now, due to a bug in a compiled dependency, the package is not guaranteed to
work on Windows.