# DCAlign: Aligning biological sequences using Direct Coupling Analysis models

# Authors

Anna Paola Muntoni, Andrea Pagnani, Martin Weigt and Francesco Zamponi

# Description

This package contains the implementation of `DCAlign`, a tool for aligning biological sequences to a DCA model of a seed, discussed in [arXiv:2005.08500](https://arxiv.org/abs/2005.08500). `DCAlign` requires, together with the set of sequences to be aligned, a model for the seed composed of:

- a DCA model;
- a set of penalties for opening and extending insertions; they can be site-dependent;
- a penalty for adding a gap at the beginning or at the end of the aligned sequence ("external") and one for inserting a gap between two symbols ("internal").

We provide in the `test` folder, the parameters of the DCA models, the gap/insertion penalties used to align the [Pfam](https://pfam.xfam.org/) and [Rfam](https://rfam.xfam.org/) families PF00035, PF00677, PF00684, PF00763, RF00059, RF00162, RF00167, RF01734 whose performances are described in [arXiv:2005.08500](https://arxiv.org/abs/2005.08500). We also add a set of synthetic sequences analyzed in the same work.
For an arbitrary seed, we provide the sub-module `DCAbuild` in order to learn of all the necessary parameters. 


The code is written in [Julia](https://julialang.org/).

## Installation

The repository `DCAlign` can be installed within the PackageManager (typing the ] key) using

```(@v1.?) pkg> add https://github.com/infernet-h2020/DCAlign```

and then typing 

```julia> using DCAlign```

### Dependencies

The package uses the modules `FastaIO` and `OffsetArrays` which must be installed before adding `DCAlign`.


# Usage

A step-by-step description of the alignment procedure of a single sequence is provided in the Jupyter notebook in the folder `notebook`. To align a set of sequences, it is possible to use the script `align_all.jl` which takes as input:


+ `q` : the length of the alphabet. The package assumes that for protein sequences we set ```q = 21``` while for RNA sequences ```q = 5```; <br>
+ `L` : the number of sites (i.e. the length of the aligned sequences) of the desired multiple sequence alignment; <br>
+ `filename_par`: name of the file where the DCA parameters are stored. By default, the package assumes that they are written using the syntax of [bmDCA](https://github.com/matteofigliuzzi/bmDCA) and [adabmDCA](https://github.com/anna-pa-m/adabmDCA) (not yet available) output: <br>
  ```J i j a b value``` <br> 
  ```h i a value``` <br>
  where _i_ and _j_ run over the columns of the MSA while _a_ and _b_ run over the symbol (`-` translates in `0`).
+ `filename_full` : name of the file containing the unaligned sequences (in FASTA format)
+ `inspen_file`: name of the file where the insertion penalties are stored. The script assumes to be written as a two-columns file of length L
+ `μext` : penalty associated with "external" gaps;
+ `μint` : penalty associated with "internal" gaps;

## Extra-input
This script allows for a comparison between the `DCAlign` alignment and a given one. In this case one should specify:
+ `filename_align` : file written in FASTA format containing the alignment (with gaps);
+ `filename_ins` : file written in FASTA format containing the aligned part of sequences as well as the gaps and insertions (written in lower case). See the Jupyter notebook for an example. <br>
It is possible to align a portion of the unaligned sequence, cut around the hit in `filename_align` . The number of symbols to be included, before and after the target hit is specified in `delta` (note that, for RNA sequences, `DCAlign` considers the full-length sequences as they are often already centered in the correct hits). <br>

If the DCA parameters are learned using [PlmDCA](https://github.com/pagnani/PlmDCA), one assumes that the format is <br>

 ```J a b i j value```              
 ```h a i value```                    

and the gap symbol translates to `q`. To use and properly read them in the alignment process it suffices to set `typel = :plm`. Note that the output parameters of `DCAbuild` will be written using this syntax.

### Algorithm parameters
It is possible to set the value of the (inverse) temperature used for the global Hamiltonian through the parameter `β`. Regarding the convergence criterium, it is possible to set: 
- `epsconv` : the minimum error, between the marginal probabilities of two consecutive iterations, required to stop the algorithm (for `β < ∞`); default: 1e-5
- `mindec` : the minimum number of consecutive iterations in which the assigned variables do not change (for `β = ∞`); default: 50

Furthermore, it is possible to modify the maximum number of iterations of `DCAlign` using `maxiter` and to print the output of `DCAlign` iterations using `verbose = true`.

## Ouput files
As default, `DCAlign` produces a file in FASTA format containing the best aligned sequence between the output of the `β < ∞` and `β = ∞` versions of the algorithm. It is possible to modify the file name using `filename_out`. <br>

When a comparison to a given MSA is performed, `DCAlign` produces an additional output file containing several information about the aligned sequences. The default name can be changed using `filename_flag`. For sake of simplicity, we report here an example:

```
>E7EVI1_HUMAN/188-252	sat: true beta: 1.00 muext: 2.5 muint: 0.0 T0: true Dist: 9 2 3 4 Time: 8.61 
test     --SLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSV-GEFSAEGEGNSKKLSKKRAATTVLQEL -150.493 -145.493 
dcalign  EISLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSVGEFSA-EGEGNSKKLSKKRAA-TVLQEL -184.124 -181.232 
ins test     --SLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSV-GEFSAEGEGNSKKLSKKRAATTVLQEL 
ins dcalign  EISLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSVGEFSA-EGEGNSKKLSKKRAA-tTVLQEL 

>E7EVI1_HUMAN/188-252	sat: true beta: 1.00 muext: 2.5 muint: 0.0 T0: false Dist: 8 1 3 4 Time: 2.12 
test     --SLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSV-GEFSAEGEGNSKKLSKKRAATTVLQEL -150.493 -145.493 
dcalign  EISLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSVGEFSA-EGEGNSKKLSKKRAATTVLQEL -186.791 -186.791 
ins test     --SLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSV-GEFSAEGEGNSKKLSKKRAATTVLQEL 
ins dcalign  EISLVFEIALKRNMPVSFEVIKESGPPHM-KSFVTRVSVGEFSA-EGEGNSKKLSKKRAATTVLQEL
```

We have aligned a domain within the protein `E7EVI1_HUMAN` belonging to the PF00035 family and we have compared it to the output of [HMMer](http://hmmer.org/), denoted here as `test`. The first line contains:
- the name of the protein and the position of the domain, `188-252`
- `sat`: if `true` the short-range constraints are all satisfied at the convergence of the algorithm
- `beta`: it refers to `β `
- `muext` and `muint` refer to the gap penalties used
- `T0`: if `true` (`false`) the aligned sequence was obtained by the `β = ∞` (`β < ∞`) implementation
- `Dist`: site-by-site comparison between the `test` and the `dcalign` sequence. The four numbers denote the Hamming distance, Gap + (i.e. how many times `dcalign` puts a `-` when `test` has a residue), Gap - (i.e. how many times `dcalign` puts a residue when `test` puts a `-`) and Mismatches (i.e. how many times both `dcalign` and `test` decide not to put a `-` but they choose different residues for the same column).
- `Time` : running time for `DCAlign`, in seconds.

The rows `test` and `dcalign` contain the aligned sequence and its energy according to the DCA model and the complete model including the gap and insertion penalties. The lines starting with `ins` report the aligned sequences complete of the insertions as lower case symbols.


## Examples
Let us align the set of full-length sequences for the RF00162 family in `test/RF00162` using `align_all`. We would like to rename the output file containing the MSA as `RF00162.ali` and have information about the iterations of `DCAlign`. To do that, one needs to run:

```
DCAlign.align_all(5, 108, "../test/RF00162/Parameters_bm_RF00162seed_potts.dat", "../test/RF00162/Test_RF00162.full",
                          "../test/RF00162/Lambda_RF00162.dat", 3.5, 3.0, filename_out = "../test/RF00162/RF00162.ali", verbose = true)
```

Let us now align the test sequences for the Pfam famility PF00677 of the `test/PF00677` folder and let us compare them to the HMMer outcomes provided in the same folder. We also would like to add `10` amino-acids before and after the hits already selected by HMMer and to rename the output files, the alignment and the comparison files, as `PF00677.ali` and `PF00677.out`. In this case we should run:

```
DCAlign.align_all(21, 87, "../test/PF00677/Parameters_bm_PF00677seed_potts.dat", "../test/PF00677/Test_PF00677.full",
                          "../test/PF00677/Lambda_PF00677.dat", 0.0, 2.0, filename_out = "../test/PF00677/PF00677.ali",
                          filename_flag = "../test/PF00677/PF00677.out", 
                          filename_align = "../test/PF00677/Test_PF00677.fasta", filename_ins = "../test/PF00677/Test_PF00677.ins", delta = 10)

```


