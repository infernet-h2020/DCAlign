function print_results(fileflag, seqpa, seqpo, J, H, L, ctype, lambda_o, lambda_e, μext, μint, seqtest, seqins, init, fin, start, fname, satmf, beta, T0, el_time)

	enmf = compute_potts_en(J, H, seqpa, L,ctype)
	enhmm = compute_potts_en(J,H, seqtest, L ,ctype)
	costmf = compute_cost_function(J,H,seqpo,L,ctype,lambda_o,lambda_e,μext,μint)
	costhmm = compute_cost_function(J,H,seqins,L,ctype,lambda_o,lambda_e,μext,μint)
	distmf, mgmf, lgmf, dmmf = hammingdist(seqpa, seqtest)

	@printf(fileflag, ">%s/%d-%d\t", fname, init + parse(Int64, start) - 1, fin + parse(Int64, start) - 1)
	@printf(fileflag, "sat: %s beta: %.2f muext: %.1f muint: %.1f T0: %s Dist: %d %d %d %d Time: %.2f\n", satmf, beta, μext, μint, T0, distmf,mgmf,lgmf,dmmf,el_time)
	@printf(fileflag, "test     %s %.3f %.3f\n", seqtest, enhmm, costhmm)
	@printf(fileflag, "dcalign  %s %.3f %.3f\n", seqpa, enmf, costmf)
	@printf(fileflag, "ins test     %s\n", seqins)
	@printf(fileflag, "ins dcalign  %s\n", seqpo)
	flush(fileflag)
end

function align_all(
    q::Int,
    L::Int,
    filename_par::String,
    filename_full::String,
    inspen_file::String,
    μext::Real,
    μint::Real;
    delta::Int = 20,
    typel::Symbol = :bm,
    filename_align::String = "",
    filename_ins::String = "",
    β::Real = 1.0,
    filename_out::String = "out.fasta",
    filename_flag::String = "flag.dat",
    seedrng::Int = 0,
    mindec::Int = 50,
    nprint::Int = 100,
    maxiter::Int = 1000,
    epsconv::Real = 1e-5,
    verbose::Bool = false,
)

    if q == 5
        #ctype = Symbol("nbase")
        ctype = :nbase
        delta = 1000
    else
        #ctype = Symbol("amino")
        ctype = :amino
    end

    if typel == :bm
        gap = 0
    else
        gap = q
    end
    println("### Reading parameters from file ###")
    J, H = read_parameters(filename_par, q, L, gap = gap, typel = typel)
    file_out = open(filename_out, "a")
    println("### Reading sequences ###")
    typedic = true
    if filename_align == "" || filename_ins == ""
        println("No delta")
        println("File full: ", filename_full)
        al = enveloptoalign(filename_full, ctype = ctype)
        names = collect(keys(al))
        fake = '-'^L
    else
        println("### Cutting the full length sequence using delta =", delta)
        file_flag = open(filename_flag, "a")
        typedic = false
        println("File full: ", filename_full)
        println("File aligned seqs: ", filename_align)
        println("File ins seqs: ", filename_ins)
        al = enveloptoalign(
            filename_full,
            filename_align,
            filename_ins,
            delta = delta,
            ctype = ctype,
        )
    end

    M = length(al)

    Lambda_all = readdlm(inspen_file)

    lambda_o = Lambda_all[:, 1]
    lambda_e = Lambda_all[:, 2]
    damp = 0.0
    println("Tot number of sequences to be aligned ", M)


    tenp = ceil(Int64, M / 10)
    println("### Starting loop ###")
    for k = 1:M
        if !typedic
            (name, hit) = split(al[k][1], "/")
            (start, finish) = split(hit, "-")
            seqtest = al[k][3]
            seqins = al[k][4]
            fullseq = al[k][2]
            seq = Seq(seqtest, fullseq, ctype)
        else
            fname = names[k]
            fullseq = al[fname]
            if occursin("/", fname)
                (name, hit) = split(fname, "/")
                (start, finish) = split(hit, "-")
            else
                name = fname
                start = "1"
                finish = @sprintf("%s", L)
            end
            seq = Seq(fake, fullseq, ctype)
        end
        if parse(Int64, start) > parse(Int64, finish) # it often happens in RNA seq
            aux = finish
            finish = start
            start = aux
        end

        N = length(fullseq)
        T0 = true
        stime = time()
        res = palign(
            seq,
            β .* J,
            β .* H,
            β .* lambda_o,
            β .* lambda_e,
            ctype = ctype,
            mindec = mindec,
            nprint = nprint,
            maxiter = maxiter,
            T0 = T0,
            μint = β * μint,
            seed = seedrng,
            damp = damp,
            μext = β * μext,
            epsconv = epsconv,
            verbose = verbose,
        )
        el_time = time() -stime
        P = res[3].pbf.P
        out = decodeposterior(P, seq.strseq)
        seqpa0 = out.pa
        seqpo0 = out.po
        sat = check_assignment(P, false, N)
        if sat == false
            (seqpa0, P) = decimate_post(res[3], T0)
            out = decodeposterior(P, seq.strseq)
            seqpo0 = out.po
        end
        b0 = out.start
        e0 = out.finish
        enT0 = compute_potts_en(J, H, seqpa0, L, ctype)
        costT0 = compute_cost_function(J,H,seqpo0,L,ctype,lambda_o,lambda_e,μext,μint)

        if !typedic
            print_results(
                file_flag,
                seqpa0,
                seqpo0,
                β .* J,
                β .* H,
                L,
                ctype,
                β * lambda_o,
                β * lambda_e,
                β * μext,
                β * μint,
                seqtest,
                seqins,
                b0,
                e0,
                start,
                name,
                sat,
                β,
                T0,
                el_time,
            )
        end

        T0 = false
        stime = time()
        res = palign(
            seq,
            β .* J,
            β .* H,
            β .* lambda_o,
            β .* lambda_e,
            ctype = ctype,
            mindec = mindec,
            nprint = nprint,
            maxiter = maxiter,
            T0 = T0,
            μint = β * μint,
            seed = seedrng,
            damp = damp,
            μext = β * μext,
            epsconv = epsconv,
            verbose = verbose,
        )
        el_time = time() - stime
        P = res[3].pbf.P
        out = decodeposterior(P, seq.strseq)
        seqpa = out.pa
        seqpo = out.po
        sat = check_assignment(P, false, N)
        if sat == false
            (seqpa, P) = decimate_post(res[3], T0)
            out = decodeposterior(P, seq.strseq)
            seqpo = out.po
        end
        enT = compute_potts_en(J, H, seqpa, L, ctype)
        costT = compute_cost_function(J,H,seqpo,L,ctype,lambda_o,lambda_e,μext,μint)
        b = out.start
        e = out.finish

        #if enT < enT0
        if costT < costT0
            @printf(
                file_out,
                ">%s/%i-%i\n",
                name,
                b + parse(Int64, start) - 1,
                e + parse(Int64, start) - 1
            )
            @printf(file_out, "%s\n", seqpa)
        else
            @printf(
                file_out,
                ">%s/%i-%i\n",
                name,
                b0 + parse(Int64, start) - 1,
                e0 + parse(Int64, start) - 1
            )
            @printf(file_out, "%s\n", seqpa0)
        end
        flush(file_out)


        if !typedic
            print_results(
                file_flag,
                seqpa,
                seqpo,
                β .* J,
                β .* H,
                L,
                ctype,
                β * lambda_o,
                β * lambda_e,
                β * μext,
                β * μint,
                seqtest,
                seqins,
                b,
                e,
                start,
                name,
                sat,
                β,
                T0,
                el_time,
            )
        end

        if mod(k, tenp) == 0
            println("Done...", k, "/", M)
        end

    end
    println("### End ###")
end
