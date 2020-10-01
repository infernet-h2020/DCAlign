# DCAlign: Aligning biological sequences using Direct Coupling Analysis models

# Authors

Anna Paola Muntoni, Andrea Pagnani, Martin Weigt and Francesco Zamponi

# Description

This package contains the implementation of ```DCAlign```, a tool for aligning biological sequences to a DCA model of a seed, discussed in [arXiv:2005.08500](https://arxiv.org/abs/2005.08500). ```DCAlign``` requires, together with the set of sequences to be aligned, a model for the seed composed of:

- a DCA model;
- a set of penalties for opening and extending insertions; they can be site-dependent;
- a penalty for adding a gap at the beginning or the end of the aligned sequence ("external") or between two symbols ("internal").

We provide in the ```test``` folder, the parameters of the DCA models, the gap/insertion penalties used to align the [Pfam](https://pfam.xfam.org/) and [Rfam](https://rfam.xfam.org/) families PF00035, PF00677, PF00684, PF00763, RF00059, RF00162, RF00167, RF01734. For an arbitrary seed, we provide the sub-module [DCAbuild](https://github.com/anna-pa-m/DCAbuild) which performs the learning of all the necessary parameters (see ```DCAbuild``` page for the documentation).

The code is written in [Julia](https://julialang.org/).

## Installation

The repository ```DCAlign``` can be installed within the PackageManager (typing the ] key) using

```(@v1.?) pkg> add https://github.com/anna-pa-m/DCAlign (remember to add the correct name)```

and then typing 

```julia> using DCAlign```

### Dependencies

The package uses the modules ```FastaIO``` and ```OffsetArrays``` which must be installed before adding ```DCAlign```.


# Usage

A step-by-step description of the alignment procedure of a single sequence is provided in the Jupyter notebook in the folder ```notebook```. To align a set of sequences, it is possible to use the script ```align_all.jl``` which takes as input:


+ ```q``` : the length of the alphabet. The package assumes that for protein sequences we set ```q = 21``` while for RNA sequences ```q = 5```; <br>
+ ```L``` : the number of sites (i.e. the length of the aligned sequence) of the desired multiple sequence alignment; <br>
+ ```filename_par```: name of the file where the DCA parameters are stored. By default, the package assumes that they are written using the syntax of [bmDCA](https://github.com/matteofigliuzzi/bmDCA) and [adabmDCA](https://github.com/anna-pa-m/adabmDCA) output:

>> **J** _i_ _j_ _a_ _b_ for _i_ ∈ \[0,...,L-1\], _j_ ∈ \[i+1,...,L-1\], _a_ ∈ \[0,...,q-1\] and _b_ ∈ \[0,...,q-1\] <br>
>> **h** _i_ _a_ for _i_ ∈ \[0,...,L-1\] and _a_ ∈ \[0,..., q-1\] <br> 

+ ```filename_full```: name of the file containing the unaligned sequences (in FASTA format)
+ ```inspen_file```: name of the file where the insertion penalties are stored. The script assumes to be written as a two-columns file of length L
+ ```μext``` : penalties associated with "external" gaps;
+ ```μint``` : penalties associated with "internal" gaps;

## Extra-input
This script allows for a comparison between the DCAlign alignment and a given one. In this case one should specify:
+ ```filename_align``` : file written in FASTA format containing the alignment (with gaps);
+ ```filename_ins``` : file written in FASTA format containing the aligned part of sequences as well as the gaps and insertions (written in lower case). See the Jupyter notebook for an example. <br>
It is possible to align a portion of the unaligned sequence, cutted around the hit in ```filename_align``` . The number of symbols to be included, before and after the target hit is specified in ```delta``` (note that, for RNA sequences, DCAlign consider the full length sequences as they are often already centered in the correct hits). <br>

If the DCA parameters are learned using [PlmDCA](https://github.com/pagnani/PlmDCA), one assumes that the format is
>>  **J** _a_ _b_ _i_ _j_ for _i_ ∈ \[1,...,L\], _j_ ∈ \[i+1,...,L\], _a_ ∈ \[1,...,q\] and _b_ ∈ \[1,...,q\] <br>
>> **h** _a_ _i_ for _i_ ∈ \[1,...,L\] and _a_ ∈ \[1,..., q\] <br> 

and the gap symbol translates to _q_. To use and properly read them in the alignment process it suffices to set ```typel = :plm```.

### Algorithm parameters
It is possible to set the value of the (inverse) temperature used for the global Hamiltonian using ```β```. Regarding the convergence criterium, it is possible to set the maximum number of iterations of DCAlign using ```maxiter```. 
Then:
- ```epsconv``` : the minimum error, between the marginal probabilities of two consecutive iterations, required to stop the algorithm (for ```β < ∞```); default: 1e-5
- ```mindec``` : the minimum number of consecutive iterations in which the assigned variables do not change (for ```β = ∞```); default: 50

## Ouput files






