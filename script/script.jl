cd(joinpath(ENV["HOME"],"CODE/DCAlign"))
using Pkg,DelimitedFiles,BenchmarkTools, Profile
Pkg.activate(".")
using Revise,DCAlign
m = DCAlign

q=21; N=87; gap=0;
ctype=Symbol("amino")
typel=Symbol("bm")

J,H =  m.read_parameters("../AlignPotts/test/PF00677/Parameters_bm_PF00677potts.dat",q,N,gap=0)
Lambda_all = readdlm("test/PF00677/Lambda_PF00677.dat")
lambda_o = Lambda_all[:,1]
lambda_e = Lambda_all[:,2]

delta = 50;
al = m.enveloptoalign( "test/PF00677/PF00677.full", "test/PF00677/PF00677.align", "test/PF00677/PF00677.ins", delta = delta, ctype = ctype);
M = length(al)

idx0 = 1
(aux,garb) = split(al[idx0][1], "/")
seqhmm = al[idx0][3];
seqins = al[idx0][4];
seq = m.Seq(seqhmm,al[idx0][2],ctype)
N = length(al[idx0][2])
display(al[idx0][1])
display(al[idx0][2])
display(seqhmm)
display(seqins)


μext = 0.00;
μint = 0.00;
β=1; T0=false; seed=33; damp=0.0;
@benchmark resalign=m.palign($seq, $J, $H, $lambda_o, $lambda_e,
    ctype = ctype ,mindec = 50, nprint=1, maxiter=2, T0 = T0, μext=μext, μint =  μint,seed=seed,
    damp=damp, epsconv=1e-10 ,verbose=false)

@benchmark resalign=m.palign($seq, $J, $H, $lambda_o, $lambda_e,
        ctype = ctype ,mindec = 50, nprint=1, maxiter=2, T0 = true, μext=μext, μint =  μint,seed=seed,
        damp=damp, epsconv=1e-10 ,verbose=false)


T0=false; seed=33; damp=0.0;
Profile.clear()
@profiler resalign=m.palign(seq, J, H, lambda_o, lambda_e,
        ctype = ctype ,mindec = 50, nprint=1, maxiter=20,
        T0 = T0, μext=μext, μint =  μint,seed=seed,
        damp=damp, epsconv=1e-10 ,verbose=false) combine=true

resalign=m.palign(seq, J, H, lambda_o, lambda_e,
                ctype = ctype ,mindec = 50, nprint=1, maxiter=20000,
                T0 = false, μext=μext, μint =  μint,seed=seed,
                damp=damp, epsconv=1e-10 ,verbose=true)
