struct SummaryAlign{T<:Real,S<:Integer}
    Pij::Matrix{T}
    Pi::Vector{T}
    Pijpc::Matrix{T}
    Pipc::Vector{T}
    Header::Vector{String}
    Seqs::Vector{String}
    Z::Matrix{S}
    W::Vector{T}
    M::Int
    N::Int
    q::Int
end

function SummaryAlign(
    Pij::Matrix{T},
    Pi::Vector{T},
    Pijpc::Matrix{T},
    Pipc::Vector{T},
    Header::Vector{String},
    Seqs::Vector{String},
    Z::Matrix{S},
    W::Vector{T},
) where {T<:Real} where {S<:Integer}
    M1, N1 = size(Pij)
    N2 = length(Pi)
    N, M = size(Z)
    if M1 != N1 || N1 != N2
        throw(DimensionMismatch("size(Pij) = $M1,$N1 length(Pi) = $N2. They should be all equal"))
    end

    M2, N2 = size(Pijpc)
    N3 = length(Pipc)

    (M1 == M2 && N1 == N2) ||
        throw(DimensionMismatch("Pij and Pijpc have different sizes"))
    N3 == N2 || throw(DimensionMismatch("Pi and Pij have diffeent sizes"))
    Nz, Mz = size(Z)
    Mw = length(W)
    Mz == Mw || throw(DimensionMismatch("W and Z have different sizes"))
    MH = length(Header)
    MH == Mz || throw(DimensionMismatch("Header and Z have different sizes"))
    MS = length(Seqs)
    MS == MH || throw(DimensionMismatch("Seqs and Z have different sizes"))
    qmin, qmax = extrema(Z)
    qmin < 0 && error("qmin = $qmin < 0")
    qmax > 21 && error("qmax = $qmax > 21")
    q = div(M1, N)
    qmin ≤ q ≤ qmax || error("q = $qmin  ∉ [$qmin,$qmax]")
    return SummaryAlign(Pij, Pi, Pijpc, Pipc, Header, Seqs, Z, W, M1, N, q + 1)
end



"""
    function summary_align(Z::Matrix{Int8};
        pseudocount::Real = 0.8,
        mode::Symbol = :nodiag,
        theta = :auto,
        max_gap_fraction::Real = 0.9,
        score::Symbol = :frob,
        min_separation::Integer = 5,
        remove_dups::Bool = false) -> Corr
Compute one and two points statistics. Returns a `Corr` structure containing
the `N(q-1) × N(q-1)` two point correlation matrix `Pij`, `Pi`, and `Z` the
`q × M` alignment data translated in `Int8`.
"""
function summary_align(
    Z::Array{Int8,2},
    Headers::Array{String,1},
    Seqs::Array{String,1};
    pseudocount::Real = 0.0,
    mode::Symbol = :nodiag,
    theta = 0.0,
    max_gap_fraction::Real = 1.0,
    score::Symbol = :frob,
    min_separation::Integer = 1,
    remove_dups::Bool = true,
)

    0 <= pseudocount <= 1 ||
        throw(DomainError("pseudocount = $pseudocount should be in the interval [0,1]"))

    println("pseudocount = $pseudocount")
    GaussDCA.use_threading(true)
    if remove_dups
        Z, _ = GaussDCA.remove_duplicate_seqs(Z)
    end
    N, M = size(Z)
    q = Int(maximum(Z))
    q > 32 && error("parameter q=$q is too big (max 32 is allowed)")

    Pi_true, Pij_true, Meff, W = GaussDCA.compute_new_frequencies(Z, q, theta)
    N == div(length(Pi_true), q - 1) ||
        throw(DimensionMismatch("something wrong here"))
    W .= W / sum(W)
    Pi, Pij = GaussDCA.add_pseudocount(Pi_true, Pij_true, pseudocount, N, q)

    return SummaryAlign(Pij_true, Pi_true, Pij, Pi, Headers, Seqs, Z, W)
end


"""
    function summary_align(nomefile::String;
        pseudocount::Real = 0.8,
        mode::Symbol = :nodiag,
        theta = :auto,
        max_gap_fraction::Real = 0.9,
        score::Symbol = :frob,
        min_separation::Integer = 5,
        remove_dups::Bool = false) -> SummaryAlign

Compute one and two points statistics of a `FASTA` msa contained in `nomefile`.
Returns a `Corr` structure containingthe `N(q-1) × N(q-1)` two point correlation
matrix `Pij`, `Pi`, and `Z` the`q × M` alignment data translated in `Int8`.
"""
summary_align(nomefile::String; args...) =
    summary_align(read_fasta_alignment(nomefile, 1.0)...; args...)

function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)

    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 ||
        error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2
    Z = Array{Int8}(undef, fseqlen, length(seqs))
    Header = Array{String}(undef, length(seqs))
    Seqs = Array{String}(undef,length(seqs))
    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        Header[seqid] = name
        Seqs[seqid] = seq
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    return Z, Header, Seqs
end

let alphabet = [
        1,
        21,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        21,
        9,
        10,
        11,
        12,
        21,
        13,
        14,
        15,
        16,
        17,
        21,
        18,
        19,
        21,
        20,
    ]
    # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end


"""
    function read_stockholm(iFile::String) -> Dict{String,String}

Read stokholm formatted `iFile` and return a `Dict{String,String}` of the alignment:
`keys -> header`, `values -> sequence`. 

"""
function read_stockholm(iFile::String)

    fasta_dict = Dict{String,String}()
    header_order = String[]
    open(iFile) do iFid
        for line in eachline(iFid)
            line = chomp(line)
            line == "" && continue
            line[1] == '#' && continue
            line == "//" && continue
            header, seq = split(line)
            if haskey(fasta_dict, header)
                fasta_dict[header] = string(fasta_dict[header], seq)
            else
                fasta_dict[header] = seq
                push!(header_order, header)
            end
        end
    end
    fasta_dict
end


function extract_ins(iFile::String)

    f = FastaReader(iFile)
    fasta_dict = Dict{String, String}()
    for (name,line) in f
        newline = replace(line, r"[.]" => s"")
        f_match = 0
        l_match = length(newline)
        i = 1
        idx_s = 0
        idx_e = 0
        count = 0
        while i <= length(newline)
            if newline[i] != '-'
                count = count + 1
            end
            if (newline[i] == '-' || isuppercase(newline[i])) && idx_s == 0
                idx_s = i
            end
            if isuppercase(newline[i]) && f_match == 0
                f_match = count
            end
            i = i + 1
        end
        i = length(newline) 
        count = count + 1
        while i >= 1
            if newline[i] != '-'
                count = count - 1
            end
            if (newline[i] == '-' || isuppercase(newline[i])) && idx_e == 0
                idx_e = i
            end
            if isuppercase(newline[i]) && l_match == length(newline)
                l_match = count
            end
            i = i - 1
        end
        if occursin("/", name)
            newname = name #pfam
        else
            newname = name * "/" * string(f_match) * "-" * string(l_match) #mafft
        end
        fasta_dict[newname] = newline[idx_s:idx_e]
    end

    return fasta_dict

end

function extract_align(f)

    fasta_dict = Dict{String, String}()
    L = 0
    for (name,line) in f
        line = replace(line, r"[.]" => s"")
        newline = ""
        for i = 1:length(line)
            if !islowercase(line[i]) 
                newline = newline * line[i]
            end
        end
        L = length(newline)
        fasta_dict[name] = newline
    end

    return fasta_dict, L


end
"""
    linsi(iFile::String, oFile::String; nthread::Int=1 eFile::Union{String,Nothing}=nothing)

Runs linsi (mafft accurate alignment strategy) on input fasta file `iFile` and writes
the output in `oFile` (aligned fasta file). Optional arguments are, `nthread::Int`
(default = 1) and `eFile::String` the file to collect errors (default = `nothing`,
i.e it writes on `stderr`)
"""
function linsi(
    iFile::String,
    oFile::String;
    nthreads::Int = 1,
    eFile::Union{String,Nothing} = nothing,
)
    if eFile === nothing
        run(pipeline(
            `linsi --thread $nthreads --anysymbol $iFile`;
            stdout = oFile,
            stderr = string(oFile, ".stderr"),
        ))
    else
        run(pipeline(
            `linsi --thread $nthreads --anysymbol $iFile`;
            stdout = oFile,
            stderr = eFile,
        ))
    end
end

function hmmbuild(iFile::String,  ohmmFile::String;  nthreads::Int=1, extra::String="")
    if length(extra) > 1
        run(pipeline(`$(HMMER_jll.hmmbuild()) $extra --cpu $nthreads $ohmmFile $iFile`), devnull)
    else
        run(pipeline(`$(HMMER_jll.hmmbuild())  --cpu $nthreads $ohmmFile $iFile`), devnull)
    end
end

function cmbuild(iFile::String, oFile::String, cmFile::String;
                    nthreads::Int=1, extra::String="")
    run(pipeline(`$(Infernal_jll.cmbuild()) -F -O $oFile $cmFile $iFile`, devnull))
end

function hmmalign(hmmFile::String,iFile::String, oFile::String)
    run(pipeline(`$(HMMER_jll.hmmalign()) -o $oFile $hmmFile $iFile`, devnull))
end
