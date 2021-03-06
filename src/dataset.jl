# need one:
# 1. FASTA  with full non redundant sequences
# 2. FASTA with envelop of domain
# return: envelopes to align

function enveloptoalign(filefull; ctype::Symbol=:amino, pos::Bool=true) 
    dful = readfull(filefull, ctype=ctype, pos=pos)
    #dful = readenvelop(filefull, ctype=ctype)
    return dful
end

function enveloptoalign(filefull,fileenv,fileins; delta = 0, ctype::Symbol=:amino) 
    denv = readenvelop(fileenv, ctype=ctype)
    #println(length(denv))
    #print(denv)
    if ctype == :amino
        dful = readfull(filefull, ctype=ctype, pos=false)
    else
        dful = readfull(filefull, ctype=ctype, pos=true)
    end
    #println(length(dful))
    dins = readenvelop(fileins, ctype= ctype)
    #println(length(dins))
    al = Vector{Tuple{String,String,String,String}}()
    for (name,seq) in dful
        #print("full name ", name)
        l = length(seq)
        if haskey(denv, name)
            envelops = denv[name]
            envinsall = dins[name]
            for env in envelops
                auxins = []
                for ins in envinsall
                    if env[1] == ins[1]
                        auxins = ins
                    end
                end
                if ctype == :amino
                    lb = max(1,env[1]-delta)
                    ub = min(l,env[2]+delta)
                else
                    lb = env[1]
                    ub = env[2]
                end
                if ctype == :amino
                    envseq = seq[lb:ub]
                else
                    envseq = seq
                end
                envhmm = env[3]
                envins = auxins[3]
		if occursin('/', name)
			envname = name
		else
                	envname = name*"/$lb-$ub"
		end
                push!(al,(envname,envseq,envhmm,envins))
            end
        end
    end
    return al
end

function readfull(filefull; ctype::Symbol=:amino, pos::Bool=true)
    ffull = FastaIO.FastaReader(filefull)
    dheader = Dict{String,String}()
    for (_name,seq) in ffull
        name=String(split(_name)[1])
        if occursin('|', name)
            (aux,name) = split(name, "|")
        end
        if occursin('/', name) && !pos
            (name,aux) = split(name,"/")
        end
        if !haskey(dheader,name)
            dheader[name] = seq
        end
    end
    close(ffull)
    dheader
end

function readenvelop(filenv; ctype::Symbol=:amino)
    fenve = FastaIO.FastaReader(filenv)
    dheader = Dict{String, Any}()
    for (name,seq) in fenve
        if occursin('|', name)
            (aux,name) = split(name, "|")
        end
        if occursin('/', name) && ctype == :amino
            name,_envelop = split(name,"/")
            envelop = parse.(Int, split(_envelop,"-"))
        else
            envelop = [0, length(seq)]
        end
        if envelop[2] < envelop[1]
            if ctype == :amino
                error("strange envelop in $name")
            end
        end
        tname = String(name)
        tname = strip(tname)
        if haskey(dheader, tname)
            push!(dheader[tname],[envelop[1],envelop[2],seq])
        else
            dheader[tname] = [[envelop[1],envelop[2],seq],]
        end
    end

    close(fenve)
    dheader
end

"""
function writenewal(fileout, al)
    fp = FastaIO.FastaWriter(fileout)
    for i in eachindex(al)
        header = al[i].name
        seq = al[i].pa        
        length(seq) > 0 && writeentry(fp,header,seq)         
    end
    close(fp)
end
"""

