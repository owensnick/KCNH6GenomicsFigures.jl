function readhomermotif(file)
    p = Vector{Float64}[]
    seq = ""
    name = ""
    data = ""
    
    for line in eachline(file)
        fields = split(line, "\t")
        if startswith(line, ">")
            data=fields
            seq = fields[1][2:end]
            name = replace(fields[2], r"\([A-Za-z0-9-\(\)\/\.]*" => "")
        else
            push!(p, parse.(Float64, fields))
        end
    end
    pwm = reduce(hcat, p)
    (name=name, seq=seq, data=data, pwm = pwm)
    
end

function loadhomer(file)
    df = CSV.read(file, DataFrame)
    rename!(df, [:Motif, :Seq, :pval, :logpval, :qval, :num_target, :percent_target, :num_background, :percent_background]) 
    df
end


function loadmotifs(dir)
    sampledirs = glob("hm*bgT", dir)
    motiffiles = mapreduce(s -> glob("*.motif", joinpath(s, "knownResults")), vcat, sampledirs)
    motifs = readhomermotif.(motiffiles)
    motdict = Dict{String, eltype(motifs)}()

    for m in motifs
        motdict[m.data[2]] = m
    end
    motdict
end

function loadhomerres(dir)

    motifdict = loadmotifs(dir)
    sampledirs = glob("hm*bgT", dir) 
    files = joinpath.(sampledirs, "knownResults.txt.gz")
    samples = basename.(sampledirs)
    fi = isfile.(files)
    if any(.!fi)
        println("Samples Missing")
        println(join(samples[.!fi], ","))
    end
    
    files = files[fi]
    samples = samples[fi]
    
    tables = loadhomer.(files);
    for (t, l) in zip(tables, samples)
        t[!, :Sample] .= l
        t.logpval = -t.logpval
        rename!(t, :logpval => :nlogpval)
    end
    res = unique(reduce(vcat, tables))
    
    sigtable = unstack(res, :Motif, :Sample, :nlogpval);
    sigtable[!, :hm_max] = vec(maximum(Matrix(sigtable[!, r"hm_"]), dims=2));
       
    res.percent_target     = parse.(Float64, replace.(res.percent_target, Ref("%" => ""))) 
    res.percent_background = parse.(Float64, replace.(res.percent_background, Ref("%" => ""))) 
    res.FC = res.percent_target./res.percent_background
    
    fctable = unstack(res, :Motif, :Sample, :FC);
    fctable[!, :fc_max] = vec(maximum(Matrix(fctable[!, r"hm_"]), dims=2));
    rename!(fctable, replace.(names(fctable), Ref(r"hm_" => "fc_")));

    sigtable = innerjoin(sigtable, fctable, on=:Motif)
    sort!(sigtable, :hm_max, rev=true)

    selection = excludefam(sigtable)
    selection.Sel = falses(size(selection, 1))
    selection.Sel[1:7] .= true
    selection.Sel[occursin.("ZNF143", selection.Motif)] .= false
    selection.Sel[occursin.(r"^Atf1", selection.Motif)] .= true

    (res=res, sigtable=selection, motifdict=motifdict)
end

function excludefam(tbl)
    exind = trues(size(tbl, 1)) 
    exind .&= .!ex_max_example_ind("ETS", tbl)
    exind .&= .!ex_max_example_ind("Sox", tbl)
    exind .&= .!ex_max_example_ind(r"Hox|HOX", tbl)
    exind .&= .!ex_max_example_ind(r"CTCF|BORIS", tbl)


    exind .&= .!ex_cand_ind(r"^Sp|^KLF|^Klf", "Sp1(Zf)/Promoter/Homer", tbl)
    exind .&= .!ex_cand_ind(r"^GFY|^Ronin", "ZNF143|STAF(Zf)/CUTLL-ZNF143-ChIP-Seq(GSE29600)/Homer", tbl)
    exind .&= .!occursin.("PRDM10", tbl.Motif)
    
    tbl[exind, :]
end

function ex_max_example_ind(pattern, tbl)
    fi = findall(occursin.(pattern, tbl.Motif))
    mi = argmax(tbl.hm_max[fi])
    ind = falses(size(tbl, 1))
    ind[fi] .= true
    ind[fi[mi]] = false
    ind
end

function ex_cand_ind(pattern, retain, tbl)
    ind = occursin.(pattern, tbl.Motif)
    ind[tbl.Motif .== retain] .= false
    ind
end