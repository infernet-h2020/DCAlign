function palign(lseq, J::Array{Float64,4}, h::Array{Float64,2}, Λ::OffsetArray{Float64,3,Array{Float64,3}}, ctype::Symbol;
    verbose::Bool = true,
    maxiter::Integer = 100,
    nprint::Integer = 500,
    seed::Int = 0,
    damp::Float64 = 0.0,
    pcount::Float64 = 0.0001,
    Δβ::Float64 = 0.05,
    thP::Float64 = 0.30,
    Δt::Int64 = 10,
    μext::Float64 = 0.0,
    μint::Float64 = 0.0
)

    verbose && println("Run DCAlign for RNA or protein alignment")
    Random.seed!(seed)
    en = Inf
    seq = Seq(lseq.header, lseq.strseq, lseq.intseq, ctype) # crappy way to avoid recompilation
    jh = Jh(deepcopy(J), deepcopy(h))
    pbf = PBF(jh, seq)
    @extract pbf : L N
    if L > N
        verbose && println("Looking for a fragment...more time needed")
        maxiter = max(5000, maxiter)
    end
    alg = Alg(verbose, maxiter, damp, Λ, nprint, pbf.N, pbf.L, pcount, thP = thP, Δβ = Δβ, μext = μext, μint = μint, Δt = Δt)
    aux = lseq.strseq
    data = Data(J, h, deepcopy(alg.Λ))
    allvar = AllVar(pbf, jh, seq, alg, data)
    
    if verbose
        println("L = $L: length of the Potts model")
        println("N = $N: length of the sequence A to be aligned")
        println("A: $aux")
    end
    maxiter, flagconv, en = update!(allvar)
    return maxiter, flagconv, allvar, en
end

function palign(
    sequn::String,
    J::Array{Float64,4},
    h::Array{Float64,2},
    Λ::OffsetArray{Float64,3,Array{Float64,3}},
    ctype::Symbol;
    kwdargs...
)
    seq = readunalignedseq(sequn, ctype)
    palign(seq, J, h, Λ, ctype; kwdargs...)
end
