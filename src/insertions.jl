

function count_insertions(seed::Dict, L::Int64)

nseq = length(seed)
ins = zeros((nseq,L))
pointseq = 0

for name in keys(seed)
	s = seed[name]
	pointseq += 1
	pos = L
	for k in 0:length(s)-1
		if isuppercase(s[length(s)-k]) || s[length(s)-k] == '-'
			go = true
			j = 1
			delta = 0
			while go
				if length(s)-k-j <= 0
					go = false
				elseif  isuppercase(s[length(s)-k-j]) || s[length(s)-k-j] == '-'
					go = false
				elseif s[length(s)-k-j] == '.'
					j += 1
				elseif islowercase(s[length(s)-k-j])
					delta += 1
					j += 1
				end
			end
			ins[pointseq,pos] = delta
			pos -= 1
		end
	end
end



return ins

end

function maximum_likelihood(ins, L, nseq)

maxit = 1e8
tol = 1e-5
eta = 1e-5
lo0 = 0.0
le0 = 0.0
lambda_o = zeros(L)
lambda_e = zeros(L)


for i in 1:L
	allz = false
	f0 = sum(x -> x == 0, ins[:,i]) / nseq
	m = mean(ins[:,i])
	if f0 == 1.0 && le0 != 0.0 && lo0 != 0.0
		lambda_e[i] = le0
		lambda_o[i] = lo0
		continue
	end
	if f0 == 1.0
		allz = true
		f0 = 1.0 - 1e-3
	end
	lo = 1.0
	le = 1.0
	eps = 1.0
	it = 1
	while eps > tol && it < maxit
		dLdle = ( exp(-lo - le) *((1.0 - exp(-le))^-2) )/(1.0 + exp(-lo) * ( (1.0-exp(-le))^-1 )) - m + (1.0-f0) - (2.0/nseq)*le
		le = le + eta * dLdle
		dLdlo = ( exp(-lo) * ((1.0-exp(-le))^-1) ) / (1.0 + exp(-lo)* ( (1.0 -exp(-le))^-1 )) - (1.0 - f0) - (2.0/nseq) *lo
		lo = lo + eta * dLdlo
		eps = max(abs(dLdle), abs(dLdlo))
		it += 1
	end
	lambda_o[i] = lo
	lambda_e[i] = le
	if allz && le0 == 0.0 && lo0 == 0.0
		le0 = le
		lo0 = lo
	end
	@printf("site: %i it: %2.1e eps: %.2e f0: %.2f m: %.2f ins pen: ( %.3f, %.3f )\n", i, it, eps, f0, m,  lo, le)
end

return lambda_o, lambda_e

end

function infer_ins_pen(fileseed::String, L::Int64; ctype::Symbol=:amino)

	seed = readfull(fileseed, ctype = ctype, pos = true)
	ins = count_insertions(seed, L)
	lo, le = maximum_likelihood(ins, L, length(seed))

	return lo, le
end

function compute_dist(S::Dict, L::Int64)

	names = keys(S)

	dist = zeros(Int64,length(names), L, L)
	count = 1
	for key in names
		seq = S[key]
		idxi = 0
		ninstot = 0
		for i = 1:length(seq)
			if islowercase(seq[i])
				ninstot += 1
				continue
			end
			idxi += 1
			j = i+1
			nins = 0
			nsymb = 0
			while j <= length(seq)
				if islowercase(seq[j])
					nins += 1
				end
				if isuppercase(seq[j])
					nsymb += 1
				end
				if isuppercase(seq[j]) || seq[j] == '-'
					idxj = j - nins - ninstot
					dist[count,idxi,idxj] = convert(Int64, nsymb + nins)
					dist[count,idxj,idxi] = convert(Int64, nsymb + nins)
				end
				j += 1
			end
		end
		count += 1
	end

	return dist
end

function deltan_prior(fileseed_ins::String, L::Int)

	S = readfull(fileseed_ins)
	Mseed = length(S)
	dist = compute_dist(S,L);
	
	N = maximum(dist)

	println("Maximum Δn ", N+1)
    Λ = OffsetArray(fill(0.0, (L,L,N+2)), 1:L, 1:L, 0:N+1)
	for i in 1:L, j in 1:L
	    f = fit(Histogram, dist[:,i,j], 0:1:N+1)
	    f2 = normalize(f, mode = :probability)
	    n = length(f2.weights)
	    for a in 1:n
	        if f2.weights[a] > 0.0
	            Δn = convert(Int64,f2.edges[1][a])
				if Δn <= N
	            	Λ[i,j,Δn] = f2.weights[a]
				else
					println("Δn: ", Δn, " w: ", f2.weights[a], "(i,j): ", i, " ",j)
				end
	        end
	    end
	end


	return Λ, Mseed, dist
end
