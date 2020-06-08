function updateT0!(allvar::AllVar)
    @extract allvar : alg pbf
    @extract pbf : L
    @extract alg : verbose nprint epsconv maxiter mindec μext μint
    en = Inf
    assign = zeros(Int,L,2)
    ndec = 0
    initialize_msg0(allvar)

    for it in 1:maxiter
        onesweepT0!(allvar)
        ndec, assign = assignment!(allvar,ndec,assign)
        en = compute_en(allvar)
        if mod(it, nprint) == 0 && verbose
            @printf("it = %d -- en %.3f μext %.2f μint %.2f ndec %d\n", it, en,μext, μint,ndec)
        end
        if ndec >= mindec
            if verbose
                @printf("it = %d -- en %.3f μext %.2f μint %.2f ndec %d\n", it, en,μext, μint,ndec)
            end
            return it, :converged, en
        end
    end
    if verbose
        @printf("it = %d -- en %.3f μext %.2f μint %.2f ndec %d\n", maxiter, en,μext, μint,ndec)
    end
    return maxiter, :unconverged, en
end

function initialize_msg0(allvar::AllVar)
    @extract allvar : pbf
    @extract pbf : L N P B F newP newB newF

    for marg in (newP, newB, newF)
        marg[1,0] = -Inf
        marg[1,N+1] = -Inf
        value, idx = findmax(marg)
        marg .= marg .- marg[idx]
    end
    for i in 1:L
        for marg in (P[i], F[i], B[i])
            marg[1,0] = - Inf
            marg[1,N+1] = -Inf
            if i == 1
                for n = 1:N+1
                    marg[0,n] = -Inf
                end
            end
            if i == L
                for n = 0:N
                    marg[0,n] = -Inf
                end
            end
            value, idx = findmax(marg)
            marg .= marg .- marg[idx]
        end
    end

end

function onesweepT0!(allvar::AllVar)
    @extract allvar : pbf alg
    @extract pbf : L P B F newP newB newF
    @extract alg : damp

    v = randperm(L)
    @inbounds for idx in 1:L # i-th marginal 1 ≤ i ≤ L
        i = v[idx]
        newPBFiT0!(allvar, i::Int)
        p,f,b = P[i], F[i], B[i]
        p .= newP; f .= newF; b .= newB;
        #p .= damp .* p .+ (1.0 - damp) .* newP # inplace sum
        #f .= damp .* f .+ (1.0 - damp) .* newF
        #b .= damp .* b .+ (1.0 - damp) .* newB
    end
end

function newPBFiT0!(allvar::AllVar, i::Int)
    @extract allvar : pbf
    @extract pbf : L N forward central bckward newP newB newF

    for fcb in (forward,central,bckward)
        fill!(fcb, 0.0)
    end

    centralT0!(central, allvar, i)
    i > 1 &&  forwardT0!(forward, allvar, i)
    i < L && backwardT0!(bckward, allvar, i)


    if i==1 # case site = 1 first site
        newF .= central
        newP .= central .+ bckward
        newB .= newP
    elseif 1 < i < L #internal sites
        newF .= forward .+ central
        newP .= newF .+ bckward
        newB .= central .+ bckward
    elseif i == L# case site = L
        newF .= forward .+ central
        newB .= central
        newP .= newF
    else
        error("the end of the world has come ...")
    end
    newP[1,0] = -Inf
    newB[1,0] = -Inf
    newF[1,0] = -Inf
    newP[1,N+1] = -Inf
    newB[1,N+1] = -Inf
    newF[1,N+1] = -Inf
    if i == 1
        for n in 1:N+1
            newP[0,n] = -Inf
            newB[0,n] = -Inf
            newF[0,n] = -Inf
        end
    end
    if i == L
        for n in 0:N
            newP[0,n] = -Inf
            newB[0,n] = -Inf
            newF[0,n] = -Inf
        end
    end
    for marg in (newP,newB,newF)
        v, idx = findmax(marg)
        marg .= marg .- marg[idx]
    end
end

# function myargmax(Pj,nbeg,nend)
#     #argm = -Inf
#     val = (-1,-1,-Inf)
#     @inbounds for n = nbeg:nend
#         for x = 0:1
#             val = ifelse(Pj[x, n]>val[3],(x,n,Pj[x,n]),val)
#         end
#     end
#     @assert val[1] > -1
#     @assert val[2] > -1
#     return val[1],val[2]
# end


function myargmax(Pj,nbeg,nend)
    argm = -Inf
    xj = -1
    nj = -1
    @inbounds for n = nbeg:nend
        for x = 0:1
            if Pj[x, n] > argm
                argm = Pj[x, n]
                xj = x
                nj = n
            end
        end
    end
    @assert(xj > -1)
    @assert(nj > -1)
    return xj,nj
end


function centralT0!(central::Marginal, allvar::AllVar, i::Int)
    @extract allvar:pbf jh seq alg
    @extract pbf:N P L q
    @extract jh:J h
    @extract seq:intseq
    @extract alg:μext μint

    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for ni = 1:N
        Ani = intseq[ni]
        scra = h[Ani, i]
        for j = 1:i-2
            xj, nj = myargmax(P[j], 0, ni - 1)
            Anj = (nj == 0 || xj == 0) ? A0 : intseq[nj]
            scra += J[Ani, Anj, i, j]
        end
        for j = i+2:L
            xj, nj = myargmax(P[j], ni + 1, N + 1)
            Anj = (xj == 0 || nj == N + 1) ? A0 : intseq[nj]
            scra += J[Ani, Anj, i, j]
        end
        central[1, ni] = scra
    end
    # case xᵢ = 0 (gap)
    @inbounds for ni = 0:N+1
        scra = 0
        Ani = A0
        for j = 1:i-2
            xj, nj = myargmax(P[j], 0, ni)
            Anj = (nj == 0 || nj == N + 1 || xj == 0) ? A0 : intseq[nj]
            scra += J[Ani, Anj, i, j]
        end
        for j = i+2:L
            xj, nj = myargmax(P[j], ni, N + 1)
            Anj = (nj == 0 || nj == N + 1 || xj == 0) ? A0 : intseq[nj]
            scra += J[Ani, Anj, i, j]
        end
        central[0, ni] = scra
    end
    @inbounds for ni = 0:N+1
        central[0, ni] =
            (ni == 0 || ni == N + 1) ? central[0, ni] - μext + h[A0, i] :
            central[0, ni] - μint + h[A0, i]
    end
    #cannot match n = 0 pointer
    #central[1,0] = -Inf
    #mcentral = maximum(central)
    #central = central .- mcentral
end

function forwardT0!(forward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P F q
    @extract jh : J
    @extract seq : intseq
    @extract alg : λo λe
    fwdm1 = F[i-1]
    A0 = q
    # case xᵢ = 1 (match)
    @inbounds for ni in 1:N
        Ani = intseq[ni]
        scra = -Inf
        @simd for nim in 0:ni-1
           Anim = nim == 0 ? A0 : intseq[nim]
           Δn = ni - nim - 1
           pen = (Δn > 0)*(nim > 0)*(- λo[i] - λe[i]*(Δn - 1))
           scra = max(scra, fwdm1[0,nim] + J[Ani, A0,   i, i-1] + pen)
           scra = max(scra, fwdm1[1,nim] + J[Ani, Anim, i, i-1] + pen)
        end
        forward[1, ni] = scra
    end
    # case xᵢ = 0 gap
    @inbounds @simd for ni in 1:N
        Ani = intseq[ni]
        forward[0,ni] = max(fwdm1[1,ni]+J[A0,Ani,i,i-1],fwdm1[0,ni]+J[A0,A0,i,i-1])
    end
    forward[0,0] = fwdm1[0,0] + J[A0,A0,i,i-1]
    scra = -Inf
    @inbounds @simd for nim = 1:N
        Anim = intseq[nim]
        scra = max(scra, fwdm1[1,nim] + J[A0,Anim,i,i-1])
    end
    forward[0,N+1] = max(scra, fwdm1[0,N+1] + J[A0,A0,i,i-1])
end

function backwardT0!(backward::Marginal,allvar::AllVar,i::Int)
    @extract allvar : pbf jh seq alg
    @extract pbf : N P B q
    @extract jh : J
    @extract seq : intseq
    @extract alg : λo λe

    bckp1 = B[i+1]
    A0 = q

    # case xᵢ = 1 (match)
    @inbounds for ni in 1:N
        # case xᵢp = 0 (gap)
        Ani = intseq[ni]
        scra = bckp1[0,ni] + J[A0,Ani,i+1,i]
        scra = max(scra, bckp1[0,N+1] + J[A0,Ani,i+1,i])
        for nip in ni+1:N
            Anip = intseq[nip]
            Δn = nip - ni- 1
	        pen = (Δn > 0)*(- λo[i+1] - λe[i+1]*(Δn - 1))
            scra = max(scra, bckp1[1,nip] + J[Anip, Ani, i+1,i] + pen)
        end
        backward[1,ni] = scra
    end
    # case xᵢ = 0 (gap)
    @inbounds for ni in 0:N+1
        scra = bckp1[0,ni] + J[A0,A0,i+1,i]
        for nip in ni+1:N
            Anip = intseq[nip]
            Δn = nip - ni - 1
	        pen = (ni > 0)*(Δn > 0)*(-λo[i+1] - λe[i+1]*(Δn - 1))
            scra = max(scra, bckp1[1,nip], J[Anip,A0,i+1,i] + pen)
        end
        backward[0,ni] = scra
    end
    #mbackward = maximum(backward)
    #backward = backward .- mbackward
end

function assignment!(allvar::AllVar, ndec::Int, assign::Matrix{Int})

    @extract allvar : pbf
    @extract pbf : L N P
    count = 0
    @inbounds for i in 1:L
        Pi = P[i]
        #v, idx = findmax(Pi)
        new_a = argmax(Pi)
        #display(new_a)
        #count how many nodes change state
        if (new_a[1] != assign[i,1]) || (new_a[2] != assign[i,2])
            count += 1
        end
        assign[i,1] = new_a[1]
        assign[i,2] = new_a[2]
    end

    if count == 0
        ndec += 1
    else
        ndec = 0
    end

    return ndec, assign
end
