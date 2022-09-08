using Pkg
Pkg.activate("..")
using DCAlign
using DelimitedFiles

fam = "PF00035"        
q = 21; 
ctype=:amino
#Pfam seed in Stockholm format
fileseed = "../seeds/" * fam * "_seed.txt" 
extra = ""
# %%
# Pay attention that to get the aligned sequences of the correct length, one has to install HMMER 3.1b2 (used in seed HMM)
out_label = "../seeds/" * fam * "_seed"
L = align_seed_pfam(fileseed, "../seeds/" * fam * "_seed", extra = extra)
seedfasta = out_label * ".fasta"
seedins = out_label * ".ins"

Λ, Mseed, dist = DCAlign.deltan_prior(seedins, L); # i j Δn

using PlmDCA

reg = abs(L - Mseed) < 100 ? 1e-6 : 1e-2
using PlmDCA

Z, W = DCAlign.get_Z_W(seedfasta, ctype)
PlmOut = plmdca(Z,W, lambdaJ = reg, lambdaH = reg)

J = PlmOut.Jtensor
h = PlmOut.htensor;


# %%
delta = 30;
al = DCAlign.enveloptoalign("../test/" * fam * "/" * fam * "_full_length_sequences.fasta",
                            "../test/" * fam * "/" * fam * "_full.fasta", 
                            "../test/" * fam * "/" * fam * "_full.ins" ,
                            delta = delta, ctype = ctype
                            )
M = length(al)
println("Tot. number of sequences ", M)
Malign = min(10000, M)

Threads.@threads for idx0 = 1:Malign
    # consider one sequence as an example
    ti = Threads.threadid()
    #println("Consider this sequence: ")
    (aux,garb) = split(al[idx0][1], "/")
    seqsol =  al[idx0][3];
    seqins = al[idx0][4];
    seq = Seq(seqsol,  al[idx0][2], ctype) 
    N = length(al[idx0][2])

    damp = 0.0
    seed = 0
    pcount = 1.0/Mseed
    @time _, conv, res,_ = DCAlign.palign(seq, deepcopy(J), deepcopy(h), deepcopy(Λ), ctype,
        nprint = 100, maxiter = 2000, seed = seed, 
        damp = damp, verbose = true, pcount = pcount);

    marg = deepcopy(res.pbf.P)

    P = copy(res.pbf.P)
    out = DCAlign.decodeposterior(P, seq.strseq, thP = res.alg.thP);
    sat = DCAlign.check_assignment(P, true, N)

    if sat == false
    
        (seqpa, P) = DCAlign.decimate_post(res, false)
        out = DCAlign.decodeposterior(P, seq.strseq, thP = res.alg.thP)
        println("Nucleation sol (a): ", out.seq)
        println("Nucleation sol (ins): ", out.seqins)
        println("(Potts) energy nucl sol: ", DCAlign.compute_potts_en(J,h,out.seq, L,ctype))
    
    end

    Ngap, Nins, Nb = DCAlign.count_gaps_ins(out.seqins)
    hdist, Gapp, Gapm, Sm = DCAlign.hammingdist(out.seq, seqsol)
    println("sol   ", seqsol)
    println("mf    ", out.seq)
    println("score ", out.score)
    println("fragm ", out.frag)
    println("dist(mf, test) = ", hdist)
    en = DCAlign.compute_potts_en(J, h, out.seq, L, ctype)
    ensol = DCAlign.compute_potts_en(J, h, seqsol, L, ctype)
    open("out_" * fam * "_plm_th_" * string(ti) * ".dat", "a") do io
        write(io, ">" * al[idx0][1] * "\n")
        write(io, out.seq * "\n")
        flush(io)
        close(io)
    end
    
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

