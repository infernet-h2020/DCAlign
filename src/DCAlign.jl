module DCAlign
using FastaIO
using Printf, Random, Logging, Distributed, DelimitedFiles
using ExtractMacro, OffsetArrays
using StatsBase, LinearAlgebra, Statistics
export Seq, palign, AllVar, align_seed_mafft, align_seed_pfam, align_seed_rfam

import Base.show
### All and Back-forw algorithm
include("types.jl")     # general types (all marginals and BP/MS msg)
include("iterate_bplc.jl")   # BPlc algoritm at chosen temperature
include("utils.jl")     # utilities
include("palign_bplc.jl")    # the main alignment function (BPlc at chosen temperature)
include("dataset.jl")   # read full length seqs or generate them
### Seed utils
include("fasta_utils.jl")
include("seed_utils.jl")
include("insertions.jl") # insertions prior
### Nucleation 
include("decimation.jl") # all decimation routines

end #end module
