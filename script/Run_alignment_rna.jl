# %% [markdown]
# # How to use DCAlign
# 
# In the following, we will show:
# (i) how to get a seed model, i.e. the DCA parameters as well as the priors, from the seed (unalign or Pfam);
# (ii) the main steps needed to align a sequence to a seed model. 
# As an example, we deal with the Pfam family PF00684 
# 

# %% [markdown]
# #### 1. Let us load the DCAlign.jl package:
# 

# %%
using Revise
using Pkg
Pkg.activate("..")
using DCAlign
using DelimitedFiles
#include("../../AlignPotts/src/AlignPotts.jl")

# %% [markdown]
# #### 2. Select the family

#= %%
fam = "PF00684"        
q = 21;
ctype=:amino
#Pfam seed in Stockholm format
fileseed = "../seeds/" * fam * "_seed.txt" 
extra = ""
=#
# %%
fam = "RF01734"  
#Malign = 10000
q = 5; 
ctype=:nbase
#Pfam seed in Stockholm format
fileseed = "../seeds/" * fam * "_seed.txt" 
extra = ""
#= %%
fam = "PF00278"        
q = 21; 
ctype=:amino
#Pfam seed in Stockholm format
fileseed = "../seeds/" * fam * "_seed.txt" 
extra = "--hand"
=#
# %% [markdown]
# #### 3. Get the seed in FASTA and FASTA+insertions format

# %%
# Pay attention that to get the aligned sequences of the correct length, one has to install HMMER 3.1b2 (used in seed HMM)
out_label = "../seeds/" * fam * "_seed"
L = align_seed_pfam(fileseed, "../seeds/" * fam * "_seed", extra = extra)
seedfasta = out_label * ".fasta"
seedins = out_label * ".ins"



# %% [markdown]
# #### 4.1 Learn the DCA seed model using adabmDCA

# %%
# Download the binary files from https://github.com/anna-pa-m/adabmDCA and copy the executable in the notebook folder
# run(`./adabmDCA -f $seedfasta -d 1e-5 -P -L -e 150 -t 75 -u 0.01 -v 0.01 -m 100 -k $fam"seed"`)

# %% [markdown]
# #### 5. Learn the prior for $P\left( \Delta n \right)$

# %%
#using DelimitedFiles
Λ, Mseed, dist = DCAlign.deltan_prior(seedins, L); # i j Δn
#λo, λe = DCAlign.infer_ins_pen(seedins, L)
#μext, μint = readdlm("../test/" * fam * "/Gap_Ext_Int.dat")

# %% [markdown]
# #### 4.2 Learn the DCA seed model using PlmDCA

# %%
#J, h = DCAlign.read_parameters("../test/" * fam * "/Parameters_bm_" * fam * "seed_potts.dat", q, L, gap = 0, typel = :bm)

reg = abs(L - Mseed) < 100 ? 1e-6 : 1e-2
using PlmDCA
if ctype == :amino
    PlmOut = plmdca(seedfasta, lambdaJ = reg, lambdaH = reg)
else
    Z, W = DCAlign.get_Z_W(seedfasta, ctype)
    PlmOut = plmdca(Z,W, lambdaJ = reg, lambdaH = reg)
end
J = PlmOut.Jtensor
h = PlmOut.htensor;

# %% [markdown]
# #### 6. Load sequences
# 
# The sub-routine **DCAlign.enveloptoalign** takes as input 
# 1. a full set of non-aligned sequences alone 
# 2. a full set of non-aligned sequences together with a known alignment in two formats (a standard insertions-free MSA and a MSA where insertions are added as lower case symbols). 
# 
# All files should be written according to the FASTA layout
# 
# Ex.
# 
# PF1_full.fasta 
# 
#     >Seq01
#     SLSTAQLLQPSGGLQASVISNIVLMKGQAKGLGFSIVGGDSIYSPIGIYVRTIFAGRAAAADGRLQEGDEILELNGESMAGLTHQDALQKFK
#     QAKKGLLTLTVRTRLTAPHALGGPLSPPLSRS
# 
# PF1_align.fasta
# 
#     >Seq01/22-104
#     -IVLMKGQAKGLGFSIVGGDSIYSPIGIYVRTIFAGRAAAADGLQEGDEILELNGESMAGLTHQDALQKFKQAKKLLTLTVR
# 
# PF1_ins.fasta
# 
#     >Seq01/22-104
#     -IVLMKGQAKGLGFSIVGGDSIYSPIGIYVRTIFAGRAAAADGrLQEGDEILELNGESMAGLTHQDALQKFKQAKKgLLTLTVR
#     
# Only for RNA, the names of full-length sequences must contain the position of the hit (in agreement with the standard Rfam format).
# 
# You must also specify in _ctype_ the type of variables (amino-acids or nucleic bases) and, when the aligned sequences are added as input, the cut-off size _delta_ of the full length sequences (the final length will be _delta_ + _L_ + _delta_)
# 
# Usage
# 1. >seq = DCAlign.enveloptoalign("PF1_full.fasta", ctype=Symbol("amino"))
# 2. >seq = DCAlign.enveloptoalign("PF1_full.fasta", "PF1_align.fasta", "PF1_ins.fasta", delta = 20, ctype=Symbol("amino"))
# 
# In the fist case _seq_ is a dictionary whereas in the last case it contains for all sequences in _i_ = 1,...,M using the following format:
#     
#     seq[i][1] : name
#     seq[i][2] : full length sequence
#     seq[i][3] : aligned sequence
#     seq[i][4] : aligned sequence with insertions
# 

# %%

delta = 30;
#al = DCAlign.enveloptoalign(fam * "_missing.fasta",
al = DCAlign.enveloptoalign("../test/" * fam * "/" * fam * "_full_length_sequences.fasta",
                            "../test/" * fam * "/" * fam * "_full.fasta", 
                            "../test/" * fam * "/" * fam * "_full.ins" ,
                            delta = delta, ctype = ctype
                            )

M = length(al)
Malign = M
println("Tot. number of sequences ", M, ". To be aligned ", Malign)

#alignment_cm = Vector{Tuple{String, String}}(undef, M)
#alignment_plm =  Vector{Tuple{String, String}}(undef, M)

Threads.@threads for idx0 = 1:M
    # consider one sequence as an example
    ti = Threads.threadid()
    #println("Consider this sequence: ")
    (aux,garb) = split(al[idx0][1], "/")
    seqsol =  al[idx0][3];
    seqins = al[idx0][4];
    seq = Seq(seqsol,  al[idx0][2], ctype) 
    N = length(al[idx0][2])
    #=
    open("out_" * fam * "_cm_th_" * string(ti) * ".dat", "a") do io
        write(io, ">" * al[idx0][1] * "\n")
        write(io, seqsol * "\n")
        flush(io)
        close(io)
    end
    =#
    #=
    println(al[idx0][1])
    println("Full length: ")
    println(al[idx0][2])
    println("Aligned by HMMer (without and with insertions): ")
    println(seqsol)
    println(seqins)
    =#


# %% [markdown]
# #### 7. Run DCAlign
# 
# Run the approximate message-passing algorithm at a chosen temperature

# %%
# BP - large connectivity, (inv) temperature β, or β → + ∞ if T0 = true  (here put the parameter β = 1)

    damp = 0.0
    seed = 0
    pcount = 1.0/Mseed
    @time _, conv, res,_ = DCAlign.palign(seq, deepcopy(J), deepcopy(h), deepcopy(Λ), ctype,
        nprint = 100, maxiter = 2000, seed = seed, 
        damp = damp, verbose = true, pcount = pcount);


    marg = deepcopy(res.pbf.P)

    P = copy(res.pbf.P)
    out = DCAlign.decodeposterior(P, seq.strseq, thP = res.alg.thP);
    seqpa = out.seq # aligned seq 
    seqpo = out.seqins # aligned seq with insertions
    sat = DCAlign.check_assignment(P, true, N)

    if sat == false
    
        (seqpa, P) = DCAlign.decimate_post(res, false)
        out = DCAlign.decodeposterior(P, seq.strseq, thP = res.alg.thP)
        seqpo = out.seqins
        println("Nucleation sol (a): ", out.seq)
        println("Nucleation sol (ins): ", out.seqins)
        println("(Potts) energy nucl sol: ", DCAlign.compute_potts_en(J,h,out.seq, L,ctype))
    
    end

    Ngap, Nins, Nb = DCAlign.count_gaps_ins(seqpo)
    hdist, Gapp, Gapm, Sm = DCAlign.hammingdist(seqpa, seqsol)
    println("sol   ", seqsol)
    println("mf    ", seqpa)
    println("score ", out.score)
    println("fragm ", out.frag)
    println("dist(mf, test) = ", hdist)
    #en = DCAlign.compute_potts_en(J, h, seqpa, L, ctype)
    #ensol = DCAlign.compute_potts_en(J, h, seqsol, L, ctype)
    if conv == :converged
        open("out_" * fam * "_plm_th_" * string(ti) * ".dat", "a") do io
            write(io, ">" * al[idx0][1] * "\n")
            write(io, out.seq * "\n")
            flush(io)
            close(io)
        end
    else
        open("out_" * fam * "_plm_th_" * string(ti) * ".dat", "a") do io
            write(io, ">" * al[idx0][1] * "\n")
            write(io, out.frag * "\n")
            flush(io)
            close(io)
        end
    end

    #seq_c = seqpa
    #println("(Potts) energy of mean-field approach ", en)

    #println("(Potts) energy of the compared sequence ", ensol)
    #println("mf ins:  ", seqpo)
    #println("sol ins: ", seqins)
    #println("mf full: ", out.pf)
    res = nothing
    #=
     @time res = AlignPotts.palign(seq, J, h, λo, λe, ctype,
                                nprint = 100, verbose = true, T0 = false,
                                maxiter = 2000, μext = μext, μint = μint,
                                seed = seed, damp = 0.0, epsconv = 1e-6);
    out = AlignPotts.decodeposterior(res[3].pbf.P, seq.strseq)
    seqpa = out.pa
    seqpo = out.po
    sat = AlignPotts.check_assignment(res[3].pbf.P, true, N)
    if sat == false
        (seqpa, P) = AlignPotts.decimate_post(res[3], false)
        out = AlignPotts.decodeposterior(res[3].pbf.P, seq.strseq)
        seqpo = out.po
    end

    #en_old[idx0] = DCAlign.compute_potts_en(J, H, seqpa, L, ctype)
    #alignment_old[idx0] = (name,seqpa)
    #alignins_old[idx0] = (name,seqpo)
    open("out_" * fam * "_dcaold_th_" * string(ti) * ".dat", "a") do io
        write(io, ">" * al[idx0][1] * "\n")
        write(io, seqpa * "\n")
        flush(io)
        close(io)
    end
    res = nothing
    =#
    GC.gc()
end


#=
# %%
using Plots
using LaTeXStrings
N = length(seq.intseq)

M = zeros(L, N + 2)

for i = 1:L
    m = marg[i]
    for n = 0:N+1
        M[i, n+1] = m[1, n]
    end
end


p1 = heatmap(M, c = :blues, aspect_ratio = 2.0, xlims = (-1, N + 2), ylims = (1, L), colorbar_title = L"P_{i}(n_{i}, x_{i} = 1)")
M = zeros(L, N + 2)
for i = 1:L
    m = marg[i]
    for n = 0:N+1
        M[i, n+1] = m[0, n]
    end
end
p2 = heatmap(M, c = :blues, colorbar_title = L"P_{i}(n_{i}, x_{i} = 0)", titlelocation = :left, aspect_ratio = 2.0, xlims = (-1, N + 2), ylims = (1, L))
plot(p1, p2, xlabel = L"n_{i}", ylabel = L"i")
=#



# %%



