## Pipeline file to build an alignment from a full length sequence
struct SimpleAlign
    Header::Vector{String}
    Seqs::Vector{String}
    function SimpleAlign(H::Vector{String},S::Vector{String})
        l = length(S[1])
        length(H) == length(S) || throw(DimensionMismatch("header and seqs have different lenghts"))
        for i in eachindex(S)
            length(S[i]) == l || throw(DimensionMismatch("not an alignment"))
        end
        new(H,S)
    end
end

SimpleAlign(d::SummaryAlign) = SimpleAlign(d.Header, d.Seqs)

function idx_insert(d::SummaryAlign, thr_ins)
    gapfreq =  vec(1.0 .- sum(reshape(d.Pi,(d.q-1,d.N)),dims=1))
    return findall(gapfreq .> thr_ins)
end

function insertify(d::SummaryAlign, idx)
    headers, seqs = d.Header,d.Seqs
    newseqs = similar(seqs)
    for i in eachindex(seqs)
        charseq = collect(seqs[i])
        for j in idx
            _c = charseq[j]
            charseq[j] = isletter(_c) ? lowercase(_c) : '.'
        end
        newseqs[i] = String(charseq)
    end
    SimpleAlign(headers,newseqs)
end



function writefasta(fou::String, d::Dict{String,String})
    FastaWriter(fou) do fw
        for (k,v) in d
            writeentry(fw,k,v)
        end
    end
end

function writefasta(fou::String,d::SimpleAlign)
    FastaWriter(fou) do fw
        h,s = d.Header,d.Seqs
        for i in eachindex(h)
            writeentry(fw,h[i],s[i])
        end
    end
end

"""
    function align_seed_mafft(
            fin::String,
            fou::String;
            nthreads::Integer = 8,
            thr_ins::Real = 1.0
        )
Read `fin` full-length seed in FASTA format and produces `fou.fasta` and `fou.ins` alignments.

"""
function align_seed_mafft(
    fin::String,
    fou::String;
    nthreads::Integer = 1,
    thr_ins::Real = 0.5,
)
    mktempdir() do tmpdir
        @info "Working in $tmpdir"
        linsiout = joinpath(tmpdir, "all_linsi.fasta")
        @info "### Aligning $fin with mafft-linsi. This may take a while ... "
        linsi(fin, linsiout, nthreads = nthreads)
        @info "Alignment done!"
        rescor = summary_align(linsiout)
        tmp_align = if thr_ins < 1.0 # make insertion by hand
            idxins = idx_insert(rescor, thr_ins)
            insertify(rescor, idxins)
        else
            SimpleAlign(rescor)
        end
        @info "Converting in FASTA (+ insertions) format"
        file_tmp_align = joinpath(tmpdir, "tmp_linsi.fasta")
        writefasta(file_tmp_align, tmp_align)
        clean_align = extract_ins(file_tmp_align)
        writefasta(fou * ".ins", clean_align)
        clean_fasta, L = extract_align(clean_align)
        @info "Converting in FASTA format"
        writefasta(fou * ".fasta", clean_fasta)
        @info "L = $L"
       
        @info "Done!"
    end

    return L
end


"""
    function align_seed_pfam(
            fseqs::String,
            fout::String;
            nthreads::Integer = 8,
        )
Read `fseqs` full-length seed downloaded from Pfam, in Stockholm format. It produces `fou.fasta` and `fou.ins` alignments from the output of `hmmbuild`.

"""

function align_seed_pfam(fseqs::String, fou::String; nthreads::Integer = 1, extra::String = "")

    L = 0
    mktempdir() do tmpdir
        @info "Working in $tmpdir"
        @info "### Computing the HMM and get aligned seed (with and without insertions)"
        hmmbuild(fseqs, joinpath(tmpdir, "tmp_stk.out"), joinpath(tmpdir, "tmp_stk.hmm"), nthreads = nthreads, extra = extra)
        file_tmp_stk = joinpath(tmpdir, "tmp_stk.out")
        d = read_stockholm(file_tmp_stk)
        #tmp, g = extract_align(d)
        writefasta(joinpath(tmpdir, "tmp.fasta"), d)
        clean_ins = extract_ins(joinpath(tmpdir, "tmp.fasta"))
        clean_fasta, L = extract_align(clean_ins)
        writefasta(fou * ".ins", clean_ins)
        writefasta(fou * ".fasta", clean_fasta)
        @info "L = $L"
        @info "Done!"
    end
    return L
end

function align_seed_rfam(fseqs::String, fou::String; nthreads::Integer = 1, extra::String = "")
    L = 0
    mktempdir() do tmpdir
        @info "Working in $tmpdir"
        @info "### Computing the CM and get aligned seed (with and without insertions)"
        cmbuild(fseqs, joinpath(tmpdir, "tmp_stk.out"), joinpath(tmpdir, "tmp_stk.cm"), nthreads = nthreads, extra = extra)
        file_tmp_stk = joinpath(tmpdir, "tmp_stk.out")
        d = read_stockholm(file_tmp_stk)
        #tmp, g = extract_align(d)
        writefasta(joinpath(tmpdir, "tmp.fasta"), d)
        clean_ins = extract_ins(joinpath(tmpdir, "tmp.fasta"))
        clean_fasta, L = extract_align(clean_ins)
        writefasta(fou * ".ins", clean_ins)
        writefasta(fou * ".fasta", clean_fasta)
        @info "L = $L"
        @info "Done!"
    end
    return L
end
