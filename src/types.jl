const Marginal = OffsetArrays.OffsetArray{Float64, 2, Array{Float64,2}}
Marginal(N::Int) = OffsetArray(rand(2,N+2), 0:1, 0:N+1)

# We adopt the following convention:
# xᵢ = 1 → match,
# xᵢ = 0 → gap,
# nᵢ = 0 if xᵢ = 0 and no match has been yet found
# N + 1 > nᵢ > 0, xᵢ = 1 for all the matched cases
# N + 1 > nᵢ > 0, xᵢ = 0 for all the internal gaps
# nᵢ = N + 1 if xᵢ = 0 for last gaps

struct Data
    J::Array{Float64,4}
    h::Array{Float64,2}
    Λ::OffsetArray{Float64,3,Array{Float64,3}}
end

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
    damp::Float64  # damping factor: Pnew = damp*Pold + (1-damp)*Pnew
    Λ::OffsetArray{Float64,3, Array{Float64,3}} #all insertions penalties (short and long range)
    nprint::Int # print info every nprint iterations
    thP::Float64
    Δβ::Float64
    Δt::Int64
    μint::Float64
    μext::Float64
    function Alg(verb::Bool, maxit::Int,  damp::Float64, Λ::OffsetArray{Float64,3, Array{Float64,3}}, 
                nprint::Int, N::Int, L::Int, pcount::Float64; 
                thP::Float64 = 0.3, Δβ::Float64 = 0.01, μint::Float64=0.0, μext::Float64=0.0, Δt::Int64=10)

            Nseed = length(view(Λ,1,2,:))
            Λnew = OffsetArray(fill(0.0, (L,L,N+2)), 1:L, 1:L, 0:N+1)

            for i = 1:L, j in 1:L
                for n in 0:N+1
                    if n < Nseed
                        Λnew[i,j,n] = (1.0 - pcount) * Λ[i,j,n] + pcount + 1e-4 * rand()
                    else
                        Λnew[i,j,n] = pcount + 1e-4 * rand() #pseudocount and small noise
                    end
                end
                Λnew[i,j,:] .= Λnew[i,j,:] ./ sum(view(Λnew,i,j,:))
            end
            return new(verb::Bool, maxit::Int, damp::Float64, Λnew::OffsetArray{Float64,3, Array{Float64,3}}, 
                nprint::Int, thP::Float64, Δβ::Float64, Δt::Int64, μint::Float64, μext::Float64)
    end
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
    data::Data # Original J, h, and Λ
end


