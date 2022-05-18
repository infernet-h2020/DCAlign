function update!(allvar::AllVar)
    @extract allvar : alg pbf data
    @extract alg : nprint verbose maxiter thP Δβ Δt
    @extract pbf : L N
    ΔP,ΔB,ΔF = 0.0,0.0,0.0
    en = Inf
    it = maxiter
    β = 1.0
    initialize_all!(allvar)
    for it in 1:maxiter
        ΔP,ΔB,ΔF = onesweep!(allvar)
        en = compute_en(allvar, β = β)
        minp, sat = check_solution(allvar)
        if mod(it,nprint) == 0 && verbose
	        @printf("it = %d -- en = %.2f sat = %d min(max P) = %.2e β = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it,
                 en, sat, minp, β, ΔP,ΔB,ΔF)
        end
        if Δβ > 0.0
            if sat && minp > thP
                if verbose
	        	    @printf("it = %d -- en = %.2f sat = %d min(max P) = %.2e β = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, 
                    en, sat, minp, β, ΔP,ΔB,ΔF)
                end
                return it, :converged, en
            end
        else
            if maximum([ΔP, ΔB, ΔF]) < 1e-4
                if verbose
	        	    @printf("it = %d -- en = %.2f sat = %d min(max P) = %.2e β = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, 
                    en, sat, minp, β, ΔP,ΔB,ΔF)
                end
                return it, :converged, en
            end
        end
        if mod(it, Δt) == 0
            β += Δβ
            allvar.jh.J .= β .* data.J
            allvar.jh.h .= β .* data.h
            allvar.alg.Λ .= data.Λ .^ β
        end
        flush(stdout)
        flush(stderr)
    end
    if verbose
        @printf("it = %d -- en = %.2f (ΔP,ΔB,ΔF) = %.2e %.2e %.2e \n", it, en, ΔP,ΔB,ΔF)
    end
    
    return maxiter, :unconverged, en
end

function initialize_all!(allvar::AllVar)
    @extract allvar : pbf alg
    @extract pbf : L N P B F newP newB newF

    
    for marg in (newP, newB, newF)
        marg .= OffsetArray(rand(2,N+2), 0:1, 0:N+1)
        marg[1,0] = 0.0
        marg[1,N+1] = 0.0
        normm = sum(marg)
        marg ./= normm
    end
    for i in 1:L
        for marg in (P[i], B[i], F[i])
            marg .= OffsetArray(rand(2,N+2), 0:1, 0:N+1)
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
            @warn "node $i has NaN in marginal"
            initialize_all!(allvar)
        end
        if any(isnan, newB)
            @warn "node $i has NaN in backward term"
            initialize_all!(allvar)
        end
        if any(isnan, newF)
            initialize_all!(allvar)
            @warn "node $i has NaN in forward term"
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
            @warn "norm site = $i = $normi"
            initialize_all!(allvar)
            #fill!(marg,inv(2(N+2)))
        end
        if any(isnan,marg)
            @warn "marg $i contains NaN!"
            initialize_all!(allvar)
            #fill!(marg,inv(2(N+2)))
        end
    end
end


function central!(central::Marginal, allvar::AllVar, i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P L q
    @extract jh : J h
    @extract seq : intseq
    @extract alg : Λ μext μint

    # case xᵢ = 1 (match)
    A0 = q
    @inbounds for nᵢ = 1:N
        Anᵢ = intseq[nᵢ]
        expscra = h[Anᵢ, i]
        for j = 1:i-2
            Pj = P[j]
            for nⱼ = 0:nᵢ-1 # light cone constraint nⱼ < nᵢ
                Anⱼ = (nⱼ == 0) ? A0 : intseq[nⱼ]
                Δn = nᵢ - nⱼ
                expscra +=  J[Anⱼ, Anᵢ, j, i] * Pj[1, nⱼ] * Λ[j,i,Δn] 
                expscra +=  J[A0, Anᵢ, j, i] * Pj[0, nⱼ] * Λ[j,i,Δn] 
            end
        end
        for j = i+2:L
            Pj = P[j]
            for nⱼ = nᵢ:N+1 # light cone constraint nⱼ > nᵢ
                Anⱼ = nⱼ == N + 1 ? A0 : intseq[nⱼ]
                Δn = nⱼ - nᵢ 
                expscra += (nᵢ == nⱼ) ? 0.0 : J[Anⱼ, Anᵢ, j, i] * Pj[1, nⱼ] * Λ[i,j,Δn] 
                expscra += J[A0, Anᵢ, j, i] * Pj[0, nⱼ] * Λ[i,j,Δn] 
            end
        end
        central[1, nᵢ] = expscra
    end
    # case xᵢ = 0 (gap)
    @inbounds for nᵢ = 0:N+1
        expscra = 0.0
        for j = 1:i-2
            Pj = P[j]
            for nⱼ = 0:nᵢ
                Δn = nᵢ - nⱼ 
                Anⱼ = (nⱼ == 0 || nⱼ == N + 1) ? A0 : intseq[nⱼ]
                expscra += J[Anⱼ, A0, j, i] * Pj[1, nⱼ] * Λ[i,j,Δn] 
                expscra += J[A0, A0, j, i] * Pj[0, nⱼ] * Λ[i,j,Δn] 
            end
        end
        for j = i+2:L
            Pj = P[j]
            for nⱼ = nᵢ:N+1
                Δn = nⱼ - nᵢ 
                Anⱼ = (nⱼ == 0 || nⱼ == N + 1) ? A0 : intseq[nⱼ]
                expscra += (nᵢ == nⱼ) ? 0.0 :  J[Anⱼ, A0, j, i] * Pj[1, nⱼ] * Λ[i,j,Δn] 
                expscra += J[A0, A0, j, i] * Pj[0, nⱼ] * Λ[i,j,Δn] 
            end
        end
        central[0, nᵢ] = expscra
    end
    @inbounds for nᵢ = 0:N+1
        if i == 1
            central[0, nᵢ] = nᵢ == 0 ? central[0, nᵢ] + h[A0, i] -μext : -Inf
        elseif i == L
            central[0, nᵢ] =
                nᵢ == N + 1 ? central[0, nᵢ] + h[A0, i] -μext : -Inf
        else
            central[0, nᵢ] =
                (nᵢ == 0 || nᵢ == N + 1) ? central[0, nᵢ] + h[A0, i] -μext :
                central[0, nᵢ] + h[A0, i] -μint
        end
    end
    #cannot match n = 0 pointer
    central[1, 0] = -Inf
    central[1, N+1] = -Inf
    @inbounds for l in eachindex(central)
        central[l] = exp(central[l])
    end
end

function forward!(forward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P F q L
    @extract jh : J
    @extract seq : intseq
    @extract alg : Λ
    fwdm1 = F[i-1]

    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for nᵢ in 1:N
        scra = 0.0
        Anᵢ = intseq[nᵢ]
         for nᵢm in 0:nᵢ-1
            Anᵢm = nᵢm == 0 ? A0 : intseq[nᵢm]
            Δn = nᵢ- nᵢm 
            # case xᵢm = 1 (match)
            scra += fwdm1[1,nᵢm] * exp(J[Anᵢm,Anᵢ,i-1,i]) * Λ[i-1,i,Δn] 
            # case xᵢm = 0 (gap)
            scra += fwdm1[0,nᵢm] * exp(J[A0,Anᵢ,i-1,i]) * Λ[i-1,i,Δn] 
        end
        forward[1, nᵢ] = scra
    end
    forward[1,0] = 0.0
    forward[1,N+1] = 0.0
    # case xᵢ = 0 gap
    @inbounds  for nᵢ in 0:N+1
       Anᵢ = (nᵢ == 0 || nᵢ == N+1) ? A0 : intseq[nᵢ]
       forward[0,nᵢ] = (nᵢ == 0 || nᵢ == N+1) ? 0.0 : fwdm1[1,nᵢ] * exp(J[Anᵢ,A0,i-1,i]) * Λ[i-1,i,0] 
       forward[0,nᵢ] += fwdm1[0,nᵢ] * exp(J[A0, A0,i-1,i]) * Λ[i-1,i,0]
    end
    @inbounds  for nᵢm = 1:N
        Anᵢm = intseq[nᵢm]
        forward[0,N+1] += fwdm1[1,nᵢm] * exp(J[Anᵢm,A0,i-1,i]) * Λ[i-1,i,N+1-nᵢm]
    end
end

function backward!(backward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P B q
    @extract jh : J
    @extract seq : intseq
    @extract alg : Λ  

    bckp1 = B[i+1]

    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for nᵢ in 1:N
        Anᵢ = intseq[nᵢ]
        # case xᵢp = 0 (gap)
        scra = bckp1[0,nᵢ] * exp(J[A0,Anᵢ,i+1,i]) * Λ[i+1,i,0] 
        scra += bckp1[0,N+1] * exp(J[A0,Anᵢ,i+1,i]) * Λ[i+1,i,N+1-nᵢ] 
        # case xᵢp = 1 (match)
         for nᵢp in nᵢ+1:N
            Anᵢp = intseq[nᵢp]
            Δn = nᵢp - nᵢ 
            scra += bckp1[1,nᵢp] * exp(J[Anᵢp, Anᵢ,i+1,i]) * Λ[i+1,i,Δn] 
        end
        backward[1,nᵢ] = scra
    end
    # case xᵢ = 0 (gap)
    @inbounds for nᵢ in 0:N+1
        # case xᵢp = 0 (gap)
        scra = bckp1[0,nᵢ] * exp(J[A0, A0,i+1,i]) * Λ[i+1,i,0] 
        # case xᵢp = 1 (match)
         for nᵢp in nᵢ+1:N
            Anᵢp = intseq[nᵢp]
            Δn = nᵢp - nᵢ
            scra += bckp1[1,nᵢp] * exp(J[Anᵢp, A0, i+1,i]) * Λ[i+1,i,Δn]
        end
        backward[0,nᵢ] = scra
    end
    backward[1,0] = 0.0
    backward[1,N+1] = 0.0
end

function compute_en(allvar::AllVar; β::Float64 = 1.0)

    @extract allvar : seq pbf jh alg
    @extract jh : J h
    @extract pbf : L P q N
    @extract seq : intseq

    en = 0.0
    match = zeros(Int,L)
    n = zeros(Int,L)
    hpart = 0.0
    @inbounds for i in 1:L
        H = P[i]
        ai = argmax(H)
        match[i] = ai[1]
		n[i] = ai[2]
        idx = n[i]
        An = match[i] == 1 ? intseq[idx] : q
        if n[i] == 0 || n[i] == N + 1
            en += match[i] == 0 ? -h[An,i] : -h[An,i]
            hpart += -h[An,i]
        else
            hpart += -h[An,i]
            en += match[i] == 0 ? -h[An,i] : -h[An,i]
        end
    end
    Jpart = 0.0

    @inbounds for i in 1:L
        ni = n[i]
        Ai = match[i] == 1 ? intseq[ni] : q
         for j in i+1:L
            nj = n[j]
            Aj = match[j] == 1 ? intseq[nj] : q
            Jpart += -J[Ai,Aj,i,j]
            en += -J[Ai,Aj,i,j]
        end
    end
    return en / β      
end

function check_solution(allvar::AllVar)

    @extract allvar : pbf
    @extract pbf : L P N

    minP = Inf
    for i in 1:L
        a, _ = findmax(P[i])
        minP = min(minP, a)
    end
    @assert minP > 0.0
    sat = check_assignment(P, false, N)

    return minP, sat

end