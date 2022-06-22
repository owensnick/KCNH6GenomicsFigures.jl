function loaddata()

    meta = loadmeta()    
    isotpm, isoweight, tpm = load_rsem_iso_tables(meta.Label, meta.File)
    stats, filtind = tablestats_filter(meta, tpm, lrt=0.4)
    ids = idtable()
    
    
    meta, tpm, isotpm, isoweight, stats, filtind, ids
end

function loadmeta(; datadir="rsem", projdir=getprojectdir(), sampledir = "c:\\home\\projects\\pot\\aligns-rsemstar-xt9\\")

    files = glob("*.isoforms.results.gz", joinpath(projdir, "data", datadir))

    samples = replace.(basename.(files), r"rsemstar.|.isoforms.results.gz" => "")
    fields = split.(samples, "_")

    labels = replace.(getindex.(fields, 2), Ref("-" => "_"))
    treat = first.(split.(labels, "_"))
    time = parse.(Int, last.(split.(labels, "_")))

    meta = DataFrame(Study="KCNH6", Treat=treat, Time=time, Label=labels, SampleName=samples, File=files)
    

    meta = sort(meta, [order(:Treat, by=x -> Dict("UIC" => 1, "hiK" => 2)[x]), :Time])
    meta[!, :Index] = 1:size(meta, 1)

    meta
end

function load_rsem_iso(file, verbose=false)
    verbose && println("Loading $file")
    data = CSV.read(file, DataFrame)
    data[!, :Gene]  = rsem_valid_geneids(data.gene_id)
    gns = split.(data.Gene, "|")
    @assert sort(unique(length.(gns))) == [1, 2]
    data[!, :Isoform] = string.(first.(gns), "|", data.transcript_id, ifelse.(length.(gns) .== 1, "", string.("|", last.(gns))))
    
    sort!(data[!, [:Gene, :Isoform, :TPM, :IsoPct, :pme_TPM, :IsoPct_from_pme_TPM]], [:Gene, :Isoform])
end


### There are some inconsistent geneids and genenames from origin gff, fix those
function rsem_valid_geneids(geneids, verbose=false)
    fields = split.(geneids, "_")
    geneid = first.(fields)
    genename = [join(f[2:end]) for f in fields]
    

    idn = Dict{String, String}()
    mids = false
    for (gid, gn) in zip(geneid, genename)
        if haskey(idn, gid)
            if isempty(idn[gid]) 
                idn[gid] = gn
            elseif (idn[gid] != gn) && !isempty(gn)
                println("Multiple ids: $gid:\t#$gn#\t", idn[gid], "#")
                mids = true
            end
        else
            idn[gid] = gn
        end
    end
    verbose && !mids && println("All Gene ID Name pairs consistent")

    verbose && println("[RVG]\t Total IDs         : ", length(idn))
    verbose && println("[RVG]\t Total w/Gene Name : ", sum(!isempty, values(idn)), "\t", sum(!isempty, values(idn))/length(idn))

    idname = Dict(gid => ifelse(isempty(gn), gid, string(gid, "|", gn)) for (gid, gn) in idn)
    getindex.(Ref(idname), geneid)
end

function load_rsem_iso_tables(labels, files; tpmfield = :TPM, isopctfield=:IsoPct)
    D = @showprogress map(load_rsem_iso, files)
    genes = first(D).Gene
    @assert all(d -> genes == d.Gene, D)

    isoforms = first(D).Isoform
    @assert all(d -> isoforms == d.Isoform, D)

    isotpm    = [DataFrame(Gene=genes, Isoform=isoforms) DataFrame(mapreduce(d -> d[!, tpmfield], hcat, D), labels)]
    isoweight = [DataFrame(Gene=genes, Isoform=isoforms) DataFrame(mapreduce(d -> d[!, isopctfield], hcat, D), labels)]

    G = [combine(groupby(d, :Gene), tpmfield => sum => tpmfield) for d in D]
    gs = first(G).Gene
    @assert all(d -> gs == d.Gene, G)
    tpm = [DataFrame(Gene=gs) DataFrame(mapreduce(g -> g[!, tpmfield], hcat, G), labels)]

    sort!(tpm, :Gene)
    sort!(isotpm, [:Gene, :Isoform])
    sort!(isoweight, [:Gene, :Isoform])

    @assert isotpm.Isoform == isoweight.Isoform
    @assert isotpm.Gene == isoweight.Gene
    @assert tpm.Gene == unique(isotpm.Gene)

    isotpm, isoweight, tpm
end

function tablestats_filter(meta, tpm; lrt=0.4, rl = 6, ns=2)

    T = Matrix(tpm[!, Symbol.(meta.Label)])

    minv  = vec(minimum(T, dims=2))
    maxv  = vec(maximum(T, dims=2))
    meanv = vec(mean(T, dims=2))
    minv  = vec(minimum(T, dims=2))
    
    ui_ind = meta.Treat .== "UIC"
    hk_ind = meta.Treat .== "hiK"


    lr_UI = [longest_run(c[ui_ind], lrt) for c in eachrow(T)]
    lr_HK = [longest_run(c[hk_ind], lrt) for c in eachrow(T)]

    stats = DataFrame(Gene=tpm.Gene, minv=minv, meanv=meanv, maxv=maxv, lr_UI=lr_UI, lr_HK=lr_HK, Index=1:length(minv))

    # filtind = mapreduce(l -> stats[!, l] .>= rl, +, [:lr_UI, :lr_UI, :lr_HK]) .>= ns
    filtind = mapreduce(l -> stats[!, l] .>= rl, +, [:lr_UI, :lr_HK]) .>= ns

    stats, filtind
end


#### load xenbase ids to enable searching by secondary gene names
function idtable(file="c:\\home\\resource\\Xt\\GenePageGeneralInfo_ManuallyCurated.txt")

    geneids = String[]
    altids  = String[]
    io = open(file) |> GzipDecompressorStream
    for line in eachline(file)
        fields = split(line, '\t')
        geneid = fields[2]
        altid = split(fields[5], "|")
        for ai in altid
            push!(geneids, geneid)
            push!(altids, ai)
        end
    end
    close(io)

    # println("Total gene ids: ", length(unique(geneids)))
    # println("Total alt  ids: ", length(unique(altids)), "  ", length(altids))

    df = DataFrame(GeneID=geneids, AltID=altids)
    
    df
end

### determin the longest run of time series data above threshold τ
function longest_run(x, τ = 0)
    rl = 0
    mrl  = 0
    @inbounds for v ∈ x
        rl  = ifelse(v > τ, rl + 1, 0)
        mrl = max(mrl, rl)
    end
    mrl
end


"""
    loadpromoters(isoweight, file="XENTR_9.1_Xenbase_spike.promoter.500.tsv.gz"

    Loads promoter coordinates, selects the maximally expressed isoform for each gene
"""
function loadpromoters(isoweight, file="XENTR_9.1_Xenbase_spike.promoter.500.tsv.gz", projdir=getprojectdir())
    filepath = joinpath(projdir, "data", file)
    proms = CSV.read(filepath, DataFrame)
    promiso = leftjoin(@subset(isoweight, .!occursin.(r"^ERCC", :Gene)), proms, on=[:Gene, :Isoform])
    @assert !any(ismissing, promiso.chrom)
    dropmissing!(promiso)

    promgene = combine(groupby(promiso, :Gene)) do df
        i = argmax(df.Weight) 
        df[i, Not([:Gene, :Isoform])]
    end

    promgene
end