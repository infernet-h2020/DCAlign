
function update!(allvar::AllVar)
    @extract allvar : alg
    @extract alg : nprint verbose epsconv maxiter μext μint
    ΔP,ΔB,ΔF = 0.0,0.0,0.0
    en = Inf
    it = maxiter
    delta = 0
    initialize_all(allvar)
    for it in 1:maxiter
        ΔP,ΔB,ΔF = onesweep!(allvar)
        en = compute_en(allvar)
        if mod(it,nprint) == 0 && verbose
	    @printf("it = %d -- en = %.2f μext = %.2f μint = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, en, μext, μint,ΔP,ΔB,ΔF)
         end
         if maximum(max.(ΔP,ΔB,ΔF)) <= epsconv
            if verbose
	        	@printf("it = %d -- en = %.2f μext = %.2f μint = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, en, μext,μint, ΔP,ΔB,ΔF)
            end
            return it, :converged, en
        end
    end
    if verbose
        @printf("it = %d -- en = %.2f μext = %.2f μint = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, en, μext, μint, ΔP,ΔB,ΔF)
    end
    return maxiter, :unconverged, en
end

function initialize_all(allvar::AllVar)
    @extract allvar : pbf alg
    @extract pbf : L N P B F newP newB newF

    for marg in (newP, newB, newF)
        marg[1,0] = 0.0
        marg[1,N+1] = 0.0
        normm = sum(marg)
        marg ./= normm
    end
    for i in 1:L
        for marg in (P[i], B[i], F[i])
            marg[1,0] = 0.0
            marg[1,N+1] = 0.0
            normm = sum(marg)
            marg ./= normm
        end
    end
end

function onesweep!(allvar::AllVar)
    @extract allvar : pbf alg
    @extract pbf : L N P B F newP newB newF
    @extract alg : damp
    ΔP, ΔF, ΔB = 0.0, 0.0, 0.0
    v = randperm(L)
    for idx in 1:L # i-th marginal 1 ≤ i ≤ L
        i = v[idx]
        newPBFi!(allvar, i::Int)
        p,f,b = P[i], F[i], B[i]
        if any(isnan, newP)
            warn("node $i has NaN in marginal")
        end
        if any(isnan, newB)
            warn("node $i has NaN in backward term")
        end
        if any(isnan, newF)
            warn("node $i has NaN in forward term")
        end
        ΔP = max(ΔP, maximum(abs.(newP-p)))
        p .= damp .* p .+ (1.0 - damp) .* newP # inplace sum
        ΔF = max(ΔF, maximum(abs.(newF-f)))
        f .= damp .* f .+ (1.0 - damp) .* newF
        ΔB = max(ΔB, maximum(abs.(newB-b)))
        b .= damp .* b .+ (1.0 - damp) .* newB
    end
    return ΔP, ΔF, ΔB
end

function newPBFi!(allvar::AllVar, i::Int)
    @extract allvar : pbf
    @extract pbf : L N forward central bckward newP newB newF F B

    for fcb in (forward,central,bckward)
        fill!(fcb, 0.0)
    end

    central!(central, allvar, i)
    i > 1 &&  forward!(forward, allvar, i)
    i < L && backward!(bckward, allvar, i)

    if i==1 # case site = 1 first site
        newF .= central
        newP .= central .* bckward
        newB .= newP
    elseif 1 < i < L #internal sites
        newF .= forward.*central
        newP .= newF .* bckward
        newB .= central .* bckward
    elseif i == L# case site = L
        newF .= forward.*central
        newB .= central
        newP .= newF
    else
        error("the end of the world has come ...")
    end

    for marg in (newP,newB,newF)
        normi = sum(marg)
        if normi > 0
            marg ./= normi
        else
            warn("norm site = $i = $normi")
            fill!(marg,inv(2(N+2)))
        end
        if any(isnan,marg)
            warn("marg $i contains NaN!")
            fill!(marg,inv(2(N+2)))
        end
    end
end


function central!(central::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P L q
    @extract jh : J h
    @extract seq : intseq
    @extract alg : μext μint

    # case xᵢ = 1 (match)
    A0 = q
    @inbounds for nᵢ in 1:N
        Anᵢ = intseq[nᵢ]
        expscra = h[Anᵢ,i]
        for j in 1:i-2
            Pj = P[j]
            @simd for nⱼ = 0:nᵢ-1 # light cone constraint nⱼ < nᵢ
               Anⱼ = nⱼ == 0 ? A0 : intseq[nⱼ]
               expscra += J[Anⱼ, Anᵢ, j, i] * Pj[1,nⱼ]
               expscra += J[A0, Anᵢ, j, i] * Pj[0,nⱼ]
            end
        end
        for j in i+2:L
            Pj = P[j]
            @simd for nⱼ = nᵢ:N+1 # light cone constraint nⱼ > nᵢ
               Anⱼ = nⱼ == N+1 ? A0 : intseq[nⱼ]
               expscra += nᵢ == nⱼ ? 0 : J[Anⱼ,Anᵢ,j,i] * Pj[1,nⱼ]
               expscra += J[A0, Anᵢ,j,i] * Pj[0,nⱼ]
            end
        end
        central[1, nᵢ] = expscra
    end
    # case xᵢ = 0 (gap)
    @inbounds for nᵢ in 0:N+1
       expscra = 0.0
       for j in 1:i-2
          Pj = P[j]
          @simd for nⱼ in 0:nᵢ
            Anⱼ = (nⱼ == 0 || nⱼ == N+1) ? A0 : intseq[nⱼ]
            expscra += J[Anⱼ,A0, j,i] * Pj[1,nⱼ]
            expscra += J[A0, A0, j,i] * Pj[0,nⱼ]
         end
      end
      for j in i+2:L
         Pj = P[j]
         @simd for nⱼ in nᵢ:N+1
            Anⱼ = (nⱼ == 0 || nⱼ == N+1) ? A0 : intseq[nⱼ]
            expscra += nᵢ == nⱼ ? 0 : J[Anⱼ, A0, j, i] * Pj[1,nⱼ]
            expscra += J[A0, A0, j, i] * Pj[0,nⱼ]
         end
      end
      central[0,nᵢ] = expscra
    end
	@inbounds @simd for nᵢ = 0:N+1
	    if i == 1
	        central[0, nᵢ] = nᵢ == 0 ? central[0, nᵢ] - μext + h[A0, i] : -Inf
	    elseif i == L
	        central[0, nᵢ] = nᵢ == N + 1 ? central[0, nᵢ] - μext + h[A0, i] : -Inf
	    else
	        central[0, nᵢ] =
	            (nᵢ == 0 || nᵢ == N + 1) ? central[0, nᵢ] - μext + h[A0, i] :
	            central[0, nᵢ] - μint + h[A0, i]
	    end
	end
    #cannot match n = 0 pointer
    central[1,0] = -Inf
    central[1,N+1] = -Inf
    @inbounds @simd for l in eachindex(central)
        central[l] = exp(central[l])
    end
end

function forward!(forward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P F q L
    @extract jh : J
    @extract seq : intseq
    @extract alg : λo λe
    fwdm1 = F[i-1]

    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for nᵢ in 1:N
        scra = 0.0
        Anᵢ = intseq[nᵢ]
        @simd for nᵢm in 0:nᵢ-1
            Anᵢm = nᵢm == 0 ? A0 : intseq[nᵢm]
            Δn = nᵢ- nᵢm -1
            pen = (Δn > 0)*(nᵢm > 0)*(-λo[i] - λe[i]*(Δn - 1));
            # case xᵢm = 1 (match)
            scra += fwdm1[1,nᵢm] * exp(J[Anᵢm,Anᵢ,i-1,i] + pen)
            # case xᵢm = 0 (gap)
            scra += fwdm1[0,nᵢm] * exp(J[A0,Anᵢ,i-1,i] + pen)
        end
        forward[1, nᵢ] = scra
    end
    forward[1,0] = 0.0
    forward[1,N+1] = 0.0
    # case xᵢ = 0 gap
    @inbounds @simd for nᵢ in 0:N+1
       Anᵢ = (nᵢ == 0 || nᵢ == N+1) ? A0 : intseq[nᵢ]
       forward[0,nᵢ] = (nᵢ == 0 || nᵢ == N+1) ? 0.0 : fwdm1[1,nᵢ] * exp(J[Anᵢ,A0,i-1,i])
       forward[0,nᵢ] += fwdm1[0,nᵢ] * exp(J[A0, A0,i-1,i])
    end
    @inbounds @simd for nᵢm = 1:N
        Anᵢm = intseq[nᵢm]
        forward[0,N+1] += fwdm1[1,nᵢm] * exp(J[Anᵢm,A0,i-1,i])
    end
end

function backward!(backward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P B q
    @extract jh : J
    @extract seq : intseq
    @extract alg : λo λe

    bckp1 = B[i+1]

    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for nᵢ in 1:N
        Anᵢ = intseq[nᵢ]
        # case xᵢp = 0 (gap)
        scra = bckp1[0,nᵢ] * exp(J[A0,Anᵢ,i+1,i])
        scra += bckp1[0,N+1] * exp(J[A0,Anᵢ,i+1,i])
        # case xᵢp = 1 (match)
        @simd for nᵢp in nᵢ+1:N
            Anᵢp = intseq[nᵢp]
            Δn = nᵢp - nᵢ - 1
            pen = (Δn > 0)*(-λo[i+1] - λe[i+1]*(Δn - 1))
            scra += bckp1[1,nᵢp] * exp(J[Anᵢp, Anᵢ,i+1,i] + pen )
        end
        backward[1,nᵢ] = scra
    end
    # case xᵢ = 0 (gap)
    @inbounds for nᵢ in 0:N+1
        # case xᵢp = 0 (gap)
        scra = bckp1[0,nᵢ] * exp(J[A0, A0,i+1,i])
        # case xᵢp = 1 (match)
        @simd for nᵢp in nᵢ+1:N
            Anᵢp = intseq[nᵢp]
            Δn = nᵢp - nᵢ- 1
            pen = (Δn > 0)*(nᵢ > 0)*(-λo[i+1] -λe[i+1]*(Δn - 1))
            scra += bckp1[1,nᵢp] * exp(J[Anᵢp, A0, i+1,i] + pen)
        end
        backward[0,nᵢ] = scra
    end
    backward[1,0] = 0.0
    backward[1,N+1] = 0.0

end

function compute_en(allvar::AllVar)

    @extract allvar : seq pbf jh alg
    @extract jh : J h
    @extract pbf : L P q N
    @extract seq : intseq
    @extract alg : λo λe μext μint

    en = 0.0
    match = zeros(Int,L)
    n = zeros(Int,L)
    hpart = 0.0
    @inbounds for i in 1:L
        H = P[i]
        ai = argmax(H)
        match[i] = ai[1]; n[i] = ai[2]
        idx = n[i]
        An = match[i] == 1 ? intseq[idx] : q
        if n[i] == 0 || n[i] == N + 1
            en += match[i] == 0 ? μext -h[An,i] : -h[An,i]
            hpart += -h[An,i]
        else
            hpart += -h[An,i]
            en += match[i] == 0 ? μint -h[An,i] : -h[An,i]
        end
    end
    Jpart = 0.0

    @inbounds for i in 1:L
        ni = n[i]
        Ai = match[i] == 1 ? intseq[ni] : q
        @simd for j in i+1:L
            nj = n[j]
            Aj = match[j] == 1 ? intseq[nj] : q
            Jpart += -J[Ai,Aj,i,j]
            en += -J[Ai,Aj,i,j]
        end
    end
    @inbounds for i in 1:L-1
        ni = n[i]
        Ai = match[i] == 1 ? intseq[ni] : q
        if match[i+1] == 1
            nj = n[i+1]
            Δn = nj - ni - 1
            pen = (Δn > 0)*(ni > 0)*(λo[i+1] + λe[i+1]*(Δn - 1))
            en += pen
        end
    end
    #println("Jpart ", Jpart, " hpart ", hpart)
    return en
end
