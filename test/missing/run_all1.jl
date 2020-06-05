

using Pkg
Pkg.activate("../..")
using AlignPotts
using DelimitedFiles
m = AlignPotts

using Printf
q = 5; L = 105; 
ctype=Symbol("nbase")
typel=Symbol("bm")
fam = "RF00059"
file_par = @sprintf("Parameters_bm_%s_phmm.dat", fam);
J, H = m.read_parameters(file_par, q, L, gap=0, typel=typel);

delta = 100;
file_full = @sprintf("%smissing_phmm2.full", fam);
file_out = open(@sprintf("%smissing_phmm.fasta", fam), "a");
al = m.enveloptoalign(file_full,ctype = ctype);
M = length(al)

ins_file =  @sprintf("Lambda_%s.dat", fam)
Lambda_all = readdlm(ins_file)
lambda_o = Lambda_all[:,1];
lambda_e = Lambda_all[:,2];

μ = 0.50;
μint = 1.50;
fake = '-' ^ L
β = 1.0; μ = β * μ; μint = β * μint;
damp = 0.0
seed = 0
println("Tot number of seqs ", M)

for fname in keys(al)
    fullseq = al[fname]
    (name,hit) = split(fname, "/")
    (start, finish) = split(hit, "-")
    if start > finish
        aux = finish
        finish = start
        start = aux
    end

    seq = m.Seq(fake, fullseq,ctype)
    N = length(fullseq)
    
    
    T0 = true
    @time res=m.palign(seq, β.*J, β.*H, β.*lambda_o, β.* lambda_e, ctype = ctype ,mindec = 50, nprint=100, maxiter=500, T0 = T0,  μint =  μint,seed=seed, damp=damp, μ=μ,epsconv=1e-10 ,verbose=true); 
  
    P = res[3].pbf.P 
    out = m.decodeposterior(P, seq.strseq);
    
    seqpa0 = out.pa # aligned seq 
    sat = m.check_assignment(P,true,N)
    if sat == false
        (seqpa0, P) = m.decimate_post(res[3],T0,idx0 = 1)
        out = m.decodeposterior(P, seq.strseq)
    end
    b0 = out.start
    e0 = out.finish
    enT0 = m.compute_potts_en( J, H, seqpa0, L,ctype)

        
    T0 = false
    @time res=m.palign(seq, β.*J, β.*H, β.*lambda_o, β.* lambda_e, ctype = ctype ,mindec = 50, nprint=100, maxiter=500, T0 = T0,  μint =  μint,seed=seed, damp=damp, μ=μ,epsconv=1e-10 ,verbose=true); 
   
    P = res[3].pbf.P 
    out = m.decodeposterior(P, seq.strseq);
    seqpa = out.pa # aligned seq 
    sat = m.check_assignment(P,true,N)
    if sat == false
        (seqpa, P) = m.decimate_post(res[3],T0,idx0 = 1)
        out = m.decodeposterior(P, seq.strseq)
    end
    enT = m.compute_potts_en( J, H, seqpa, L,ctype)
    b = out.start
    e = out.finish
    println(start, finish, b0, e0, b ,e)
    if enT < enT0
        @printf(file_out, ">%s/%i-%i\n",name, b + parse(Int64,start)-1, e + parse(Int64, finish)-1)
        @printf(file_out, "%s\n",seqpa)
    else
        @printf(file_out, ">%s/%i-%i\n",name, b0 + parse(Int64,start)-1, e0 + parse(Int64, finish)-1)
        @printf(file_out, "%s\n",seqpa0)
    end

    flush(file_out)
        
end
    

