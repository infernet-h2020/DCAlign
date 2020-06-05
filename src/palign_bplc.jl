function palign(lseq, J::Array{Float64,4}, h::Array{Float64,2}, λo::Array{Float64,1}, λe::Array{Float64,1};

                verbose::Bool=true,
                maxiter::Integer=100,
                nprint::Integer=500,
                epsconv::Real=1e-4,
                seed::Int=0,
                damp::Real=0,
                μint::Real=1.0,
                μext::Real=1.0,
		λ::Real=1.0,
	        T0::Bool=false,
                mindec::Int = 10,
                ctype::Symbol="amino"
               )

    Random.seed!(seed)
    en = Inf
    seq = Seq(lseq.header, lseq.strseq, lseq.intseq,ctype) # crappy way to avoid recompilation
    jh = Jh(J,h)
    alg = Alg(verbose, maxiter, epsconv, mindec,damp, μext, λo,λe, μint, nprint)
    pbf = PBF(jh, seq)
    aux = lseq.strseq
    allvar = AllVar(pbf, jh, seq, alg)
    @extract pbf : L N
    @extract pbf : newP 
    if T0 == false
        if verbose
            println("Run Belief Propagation, large-connectivity approximation, for proteins alignment")
            println("L $L = length of the Potts model")
            println("N $N = length of the sequence A to be aligned")
            println("A: $aux")
        end
        maxiter, flagconv, en= update!(allvar)
    else
        if verbose
            println("Run Belief Propagation, large-connectivity approximation, for proteins alignment")
            println("L $L = length of the Potts model")
            println("N $N = length of the sequence A to be aligned")
            println("A: $aux")
            println("T = 0, using β → +∞ limit")
        end
        maxiter, flagconv, en = updateT0!(allvar)
    end
    return maxiter, flagconv, allvar, en
end

function palign(filefasta::String, J::Array{Float64,4}, h::Array{Float64,2}, λo::Array{Float64,1}, λe::Array{Float64, 1}, ctype::Symbol; kwdargs...)
    seq = readunalignedfasta(filefasta, ctype)[1]
    palign(seq,J,h; kwdargs...)
end
