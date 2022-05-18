
struct DecodedPosterior
    seq::String          # aligned sequence usign argmax(P)
    score::String        # score = max(P[argmax])
    frag::String         # aligned sequence usign max(P[argmax]) or `-` if max(P[argmax]) < 5 (used when it does not converge -> fragment)
    seqins::String       # aligned sequence with insertions
    fragins::String      # aligned fragment with insertions
    fullfrag::String     # full sequence with match states in uppercase and gaps (associated with fragment)
    fullseq::String      # full sequence with match states in uppercase and gaps
    start::Int           # start and end index of the aligned part w.r.t the full sequence
    finish::Int
    frag_start::Int
    frag_finish::Int
end

function Base.show(io::IO, x::DecodedPosterior)
    println(io, x.seq)
    println(io, x.seqins)
    print_with_color(:cyan, io, x.fds, "\n")
    print(io, x.fdp)
end


function check_assignment(P, verbose, N)

    count = 0
    L = length(P)
    if verbose == true
        println("Let us check the assignment...")
    end
    a = CartesianIndex{2}
    for i = L:-1:2
        a = argmax(P[i])
        x = a[1]; n = a[2];
        a = argmax(P[i-1])
        xj = a[1]; nj = a[2]
        sat = (x == 1) ? (n > nj) : (n == nj || n == N + 1)
        if sat == false && verbose == true
            println("- $i → ($x, $n) $i-1 → ($xj, $nj)")
        end
        count += 1 - convert(Int, sat)
    end

    if count == 0
        if verbose == true
            println("The subsequence satisfies the constraints")
        end
        return true
    else
        if verbose == true
            println("There are $count short-range constraints not satisfied")
        end
        return false
    end

end


#amino = ['A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
#nbase = ['A', 'U', 'C', 'G', '-']


aminoalphabet = Dict{String,Int64}()
nbasealphabet = Dict{String,Int64}()
aminoalphabet = Dict('A' => 1, 'B' => 21, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 'G' => 6, 'H' => 7, 'I' => 8, 'J' => 21, 'K' => 9, 'L' => 10, 'M' => 11, 'N' => 12, 'O' => 21, 'P' => 13, 'Q' => 14, 'R' => 15, 'S' => 16, 'T' => 17, 'U' => 21, 'V' => 18, 'W' => 19, 'X' => 21, 'Y' => 20, 'Z' => 21, '-' => 21)
nbasealphabet = Dict('A' => 1, 'U' => 2, 'C' => 3, 'G' => 4, '-' => 5, 'T' => 2, 'N' => 5, 'R' => 5, 'X' => 5, 'V' => 5, 'H' => 5, 'D' => 5, 'B' => 5, 'M' => 5, 'W' => 5, 'S' => 5, 'Y' => 5, 'K' => 5)

function letter2num(c::Union{Char,UInt8}, ctype::Symbol)
    if ctype == :amino
        return aminoalphabet[c]
    end
    if ctype == :nbase
        return nbasealphabet[c]
    end
end


function compute_potts_en(J::Array{Float64,4}, h::Array{Float64,2}, seq::String, L::Int, ctype::Symbol)

    en = 0
    for i = 1:L
        Aᵢ = letter2num(seq[i], ctype)
        en += -h[Aᵢ, i]
        for j = i+1:L
            Aⱼ = letter2num(seq[j], ctype)

            en += -J[Aᵢ, Aⱼ, i, j]
        end
    end
    return en
end

function count_gaps_ins(v::String)

    N = length(v)
    Ng = 0
    Ni = 0
    Nb = 0
    types = zeros(N)
    for i = 1:N
        if isuppercase(v[i])
            types[i] = 1
        end
        if islowercase(v[i])
            types[i] = -1
        end
        if v[i] == '-'
            types[i] = 0
        end
    end
    for i = 1:N-1
        if v[i] == '-'
            Ng += 1
        elseif islowercase(v[i])
            Ni += 1
        end
        if types[i] != types[i+1]
            Nb += 1
        end
    end
    if v[N] == '-'
        Ng += 1
    end
    if islowercase(v[N])
        Ni += 1
    end

    return Ng, Ni, Nb
end

function hammingdist(seqaux, seqtrue)
    length(seqaux) == length(seqtrue) || error("seqs of different lengths")
    ctr = 0
    moreg = 0
    lessg = 0
    diffm = 0
    for i in eachindex(seqaux)
        if seqaux[i] !== seqtrue[i]
            ctr += 1
        end
        if seqaux[i] == '-' && seqtrue[i] != '-'
            moreg += 1
        end
        if seqaux[i] != '-' && seqtrue[i] == '-'
            #display(["lg" seqaux[i] seqtrue[i]])
            lessg += 1
        end
        if seqaux[i] != '-' && seqtrue[i] != '-' && seqaux[i] != seqtrue[i]
            diffm += 1
        end
    end
    return ctr, moreg, lessg, diffm
end



function read_parameters(filename::String, q::Int, L::Int; gap::Int = 0, typel::Symbol = :bm) # generic function if input is a file

    if typel == :plm
        @info "Assuming J a b i j and h a i format"
    else
        @info "Assuming J i j a b and h i a format"
    end
    @info "Output tersors: J[a b i j] and h[a i]"
    @info "Gap in input file $gap now in $q"
    J = zeros(q, q, L, L)
    h = zeros(q, L)

    if gap == q
        offset = 0
    else
        offset = 1
    end
    open(filename) do file
        for ln in eachline(file)
            line = split(ln, ' ')
            if occursin('J', ln)
                if typel == :bm
                    i = parse(Int64, line[2]) + offset
                    j = parse(Int64, line[3]) + offset
                    a = gap == q ? parse(Int64, line[4]) + offset :
                        parse(Int64, line[4])
                    b = gap == q ? parse(Int64, line[5]) + offset :
                        parse(Int64, line[5])
                else
                    i = parse(Int64, line[4]) + offset
                    j = parse(Int64, line[5]) + offset
                    a = gap == q ? parse(Int64, line[2]) + offset :
                        parse(Int64, line[2])
                    b = gap == q ? parse(Int64, line[3]) + offset :
                        parse(Int64, line[3])
                end
                if a == gap && gap == 0
                    a = q
                end
                if b == gap && gap == 0
                    b = q
                end
                J[a, b, i, j] = parse(Float64, line[6])
                J[b, a, j, i] = parse(Float64, line[6])
            end
            if occursin('h', ln)
                if typel == :bm
                    i = parse(Int64, line[2]) + offset
                    a = gap == q ? parse(Int64, line[3]) + offset :
                        parse(Int64, line[3])
                else
                    i = parse(Int64, line[3]) + offset
                    a = gap == q ? parse(Int64, line[2]) + offset :
                        parse(Int64, line[2])
                end
                if a == gap && gap == 0
                    a = q
                end
                h[a, i] = parse(Float64, line[4])
            end
        end
    end
    return J, h
end

function decodeposterior(P, strseq; thP::Float64=0.5)

    N = length(strseq)
    L = length(P)
    seq = ""
    score = ""
    frag = ""
    for i = 1:L
        maxP = -Inf
        idxx = []
        idxn = []
        for x = 0:1, n = 0:N+1
            if P[i][x, n] > maxP
                maxP = P[i][x, n]
                idxx = x
                idxn = n
            end
        end
        if maxP == -Inf
            println("Problem with $i")
        end
        sc = convert(Int,mod(floor(Int,maxP*10),10))
        score *= (sc == 10) ? '*' : string(sc,pad=1)
        n = idxn
        x = idxx
        if x == 0
            frag *= '-'
            seq *= '-'
        else
            frag *= (maxP > thP) ? strseq[n] : '-'
            seq *= strseq[n]
        end
    end
    seqins = ""
    fragins = ""
    fullseq = ""
    fullfrag = ""
    i = 1
    start = false
    nold = 0
    f = 0
    l = 0
    ff = 0
    ll = 0
    while i <= L
        maxP = -Inf
        idxx = []
        idxn = []
        for x = 0:1, n = 0:N+1
            if P[i][x, n] > maxP
                maxP = P[i][x, n]
                idxx = x
                idxn = n
            end
        end
        n = idxn
        x = idxx
        if x == 1 && start == true
            delta = n - nold - 1
            if delta > 0
                for d = 1:delta
                    seqins *= lowercase(strseq[nold+d])
                    fragins *= lowercase(strseq[nold+d])
                end
            end
            fragins *= (maxP > thP) ? strseq[n] : lowercase(strseq[n]) * '-'
            seqins *= strseq[n]
            nold = n
            l = n
            ll = (maxP > thP) ? n : ll
        end
        if x == 1 && start == false
            start = true
            fragins *= (maxP > thP) ? strseq[n] : lowercase(strseq[n]) * '-'
            seqins *= strseq[n]
            nold = n
            f = n
            ff = (maxP > thP) ? n : ff
        end
        if x == 0
            seqins *= '-'
            fragins *= '-'
        end
        i += 1
    end
    fullseq *= lowercase.(strseq[1:f-1])
    fullseq *= seqins
    fullseq *= lowercase.(strseq[l+1:end])
    fullfrag *= lowercase.(strseq[1:f-1])
    fullfrag *= fragins
    fullfrag *= lowercase.(strseq[l+1:end])
    #@show strseq[l+1:end], f, l, ff, ll
    DecodedPosterior(seq, score, frag, seqins, fragins, fullfrag, fullseq, f, l, ff, ll)
end
