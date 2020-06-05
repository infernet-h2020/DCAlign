
function constrain_neigh_plus!(prob, uni, N, xsol, nsol, T0::Bool)
	M = prob;
	M .= M .+ uni
	zeroprob = (T0 == true) ? -Inf : 0.0
	if xsol == 1
		for n in 0:N+1
			M[0,n] = (n == nsol || n == N + 1) ? M[0,n]  : zeroprob
			M[1,n] = (n > nsol) ? M[1,n] : zeroprob
		end
	else
		for n in 0:N+1
			M[0,n] = (n == nsol) ? M[0,n]  : zeroprob
			M[1,n] = (n > nsol) ? M[1,n] : zeroprob
		end
	end
	M[1,0] = zeroprob
	M[1,N+1] = zeroprob
	if T0
		M .= M .- maximum(M)
	else
		M .= M./sum(M)
	end
	return M
end

function constrain_neigh_minus!(prob, uni, N, xsol, nsol, T0::Bool)

	M = prob;
	M .= M .+ uni
	zeroprob = (T0 == true) ? -Inf : 0.0
	if xsol == 0
	    if nsol == 0
	       for n = 1:N+1
		  M[0,n] = zeroprob
	       end
	       for n = 0:N+1
		  M[1,n] = zeroprob
	       end
	    elseif nsol == N + 1
	       for n = 0:N
		  M[0,n] = zeroprob
	       end
	    else
	        M[0,0] = zeroprob
		M[0, N+1] = zeroprob
		for x in 0:1, n in 1:N
		     M[x,n] = (n == nsol) ? M[x,n] : zeroprob
		end
	    end
	else
		for x in 0:1, n in 0:N+1
		     M[x,n] = (n < nsol) ? M[x,n]  : zeroprob
		end
	end
	M[1,0] = zeroprob
	M[1,N+1] = zeroprob
	if T0
		M .= M .- maximum(M)
	else
		M .= M./sum(M);
	end
	return M
end

function decimate_post(allvar::DCAlign.AllVar, T0::Bool)

	@extract allvar : seq pbf jh alg
	@extract pbf : P L q N
	sol = zeros(Int64,L,2)
	fixed = zeros(Bool,L)
	can_assign = ones(Bool,L)
	Nfree = count(fixed[i] == false for i in 1:L)
	uni = 1/(2*(N+2))
	while Nfree > 0
		Hmax = -Inf
		imax = NaN; xmax = NaN; nmax = NaN;
		if Nfree == L
			for i in 1:L
				for x in 0:1, n in 0:N+1
					if(P[i][x,n] > Hmax)
						Hmax = P[i][x,n]
						xmax = x
						nmax = n
						imax = i
					end
				end
			end
		else
			idx = findall(can_assign .== true)
			if length(idx) > 2
				display(["we have a problem, indices are" length(idx)])
			elseif length(idx) == 2
				if idx[1] < idx[2]
					idxm = idx[1]
					idxM = idx[2]
				else
					idxm = idx[2]
					idxM = idx[1]
				end
				P[idxm] = constrain_neigh_minus!(P[idxm], uni, N, sol[idxm+1,1], sol[idxm+1,2], T0)
				P[idxM] = constrain_neigh_plus!(P[idxM], uni, N, sol[idxM-1,1], sol[idxM-1,2], T0)
			elseif length(idx) == 1
				if idx[1]-1 > 0
				    if fixed[idx[1]-1] == true
					P[idx[1]] = constrain_neigh_plus!(P[idx[1]], uni, N, sol[idx[1]-1,1], sol[idx[1]-1,2],T0)
				    end
				end
				if idx[1]+1 <= L
				    if fixed[idx[1]+1] == true
					P[idx[1]] = constrain_neigh_minus!(P[idx[1]], uni,N, sol[idx[1] + 1,1], sol[idx[1] + 1,2], T0)
				    end
				end
			end
			for j in 1:length(idx)	
			    for x in 0:1, n in 0:N+1
			    	if(maximum(P[idx[j]][x,n]) > Hmax)
					Hmax = P[idx[j]][x,n]
					xmax = x
					nmax = n
					imax = idx[j]
				end
			    end
			end
		end
		fixed[imax] = true
		sol[imax,1] = xmax; sol[imax,2] = nmax;
		if Nfree == L
			can_assign = zeros(Bool, L,)
			if imax == 1
				can_assign[2] = true
			elseif imax == L
				can_assign[L-1] = true
			else
				can_assign[imax-1] = true
				can_assign[imax+1] = true
			end
		else
			can_assign[imax] = false
			if imax == 1
				can_assign[2] = (fixed[2] == false) ? true : false
			elseif imax == L
				can_assign[L-1] = (fixed[L-1] == false) ? true : false
			else
				can_assign[imax-1] = (fixed[imax-1] == false) ? true : false
				can_assign[imax+1] = (fixed[imax+1] == false) ? true : false
			end
		end
		Nfree = count(fixed[i] == false for i in 1:L)
	end
	s = ""
	for i in 1:L
		x = sol[i,1]; n = sol[i,2];
		if x == 0
			s *= "-"
		else
			s *= seq.strseq[n]
		end
	end
	return s, P


end


