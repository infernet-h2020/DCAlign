module DCAlign
using FastaIO
using Printf, Random
using ExtractMacro, OffsetArrays
export Seq, palign
import Base.show
using Distributed
using DelimitedFiles

include("types.jl")     # general types (all marginals and BP/MS msg)
include("iterate_bplc.jl")   # BPlc algoritm at chosen temperature
include("utils.jl")     # utilities
include("palign_bplc.jl")    # the main alignment function (BPlc at chosen temperature)
include("dataset.jl")   # read full length seqs or generate them
include("iterate_bplc0.jl") # BPlc update at zero temperature
include("decimation.jl") # all decimation routines
include("align_all.jl") # align a set of sequences

end #end module
