struct Output
    seq::String
    seqins::String
    start::Int
    finish::Int
    en::Float64
    sat::Bool
end

mutable struct Alg
    verbose::Bool #verbosity of output
    maxiter::Int  # maximal number of sweeps
    epsconv::Float64 # tolerance for the convergence of T > 0 alg
    mindec::Int # convergence for T = 0 alg, assignment must be consecutively repeated for at least mindec times
    damp::Float64  # damping factor: Pnew = damp*Pold + (1-damp)*Pnew
    μext::Float64 # gap external penalty    exp[-μext], nᵢ ∈ {N+1, 0}, xᵢ = 0
    λo::Vector{Float64} # penalty for opening exp[-λo]I(Δn > 0)
    λe::Vector{Float64} # penalty for extending exp[-λe(Δn - 1)I(Δn >0)]
    μint::Float64 #gap internal penalty exp[-μint], N ≥ n ≥ 1, xᵢ = 0
    nprint::Int # print info every nprint iterations
end


struct Seq
    header::String
    strseq::String
    intseq::Vector{Int}
    ctype::Symbol
    function Seq(h::String,s::String,v::Vector{Int},ctype::Symbol)
        if length(s) != length(v) 
            error("length s != length v")
        end
        for i in eachindex(s) #check alignment between strseq and inteseq
            if letter2num(s[i],ctype) != v[i] 
                error("intseq not representing strseq at index $i")
            end
        end
        return new(h::String,s::String,v::Vector{Int},ctype::Symbol)
    end

end

function Seq(h::String,s::String,ctype::Symbol)
    v = fill(-1,length(s))
    s = Vector{Char}(s)
    if ctype == :nbase
        replace!(s, 'T' => 'U')
    end
    for i in eachindex(s)
        v[i] = letter2num(s[i],ctype)
    end
    s = String(s)
    Seq(h,s,v,ctype)
end

function Base.show(io::IO, x::Seq)
    println(io,x.header)
    #print_with_color(:cyan,io,x.strseq)
    println(io,x.strseq)
end


struct Jh
    J::Array{Float64,4}  # couplings q × q × L × L (L length of Potts Model)
    h::Array{Float64,2}  # field     q × L
    function Jh(J::Array{Float64,4}, h::Array{Float64,2})
        qJ1,qJ2,NJ1,NJ2 = size(J)
        qh,Nh = size(h)
        qJ1 != qJ2 && error("not square color submatrix")
        NJ1 != NJ2 && error("not square sites submatrix")
        NJ1 != Nh  && error("inconsistent number of sites between fields and couplings")
        qJ1 != qh  && error("inconsistent number of colors between fields and couplings")
        return new(J, h)
    end
end

const Marginal = OffsetArrays.OffsetArray{Float64, 2, Array{Float64,2}}
Marginal(N::Int) = OffsetArray(rand(2,N+2), 0:1, 0:N+1)

# We adopt the following convention:
# xᵢ = 1 → match,
# xᵢ = 0 → gap,
# nᵢ = 0 if xᵢ = 0 and no match has been yet found
# N + 1 > nᵢ > 0, xᵢ = 1 for all the matched cases
# N + 1 > nᵢ > 0, xᵢ = 0 for all the internal gaps
# nᵢ = N + 1 if xᵢ = 0 for last gaps

struct PBF
    L::Int # Potts-model length
    N::Int # Sequence-to-align length
    q::Int # Alphabet (20 amino acids + "-") gap = 21.
    P::Vector{Marginal} # P[i] → p[nᵢ, xᵢ] where nᵢ∈ [0,N+1] , xᵢ ∈ [0,1] ∀ i ∈ [1..L] 
    B::Vector{Marginal}
    F::Vector{Marginal}
    newP::Marginal
    newB::Marginal
    newF::Marginal
    forward::Marginal
    central::Marginal
    bckward::Marginal
end

function PBF(L::Int, N::Int, q::Int)

    fwds = [Marginal(N) for i in 1:L]
    bcks = [Marginal(N) for i in 1:L]
    PBF(L,N,q,
        [Marginal(N) for i in 1:L],
        bcks,
        fwds,
        Marginal(N),
        Marginal(N),
        Marginal(N),
        Marginal(N),
        Marginal(N),
        Marginal(N))
end

function PBF(x::Jh,seq::String)
    q,L = size(x.h)
    N = length(seq)
    return PBF(L,N,q)
end

function PBF(x::Jh, seq::Seq)
    q,L = size(x.h)
    N = length(seq.intseq)
    return PBF(L,N,q)
end

struct AllVar
    pbf::PBF # marginal
    jh::Jh   # coupling
    seq::Seq # sequence to align
    alg::Alg # algorithm parameters.
end


