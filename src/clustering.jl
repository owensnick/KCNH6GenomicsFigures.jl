
maxnorm(X) = X./maximum(X, dims=2)

function reordercluster(KMK, nr = 1:KMK.k) 
    KSI = KMK.KSI[nr]
    KSSI = sortperm(KSI)
    A = KSSI[KMK.KM.assignments]
    SI = sortperm(A)
    (KSI = KSI, A=A, SI=SI, k=KMK.k, KM=KMK.KM, X=KMK.X)
end

function pot_kmeans_ud(meta, tpm, filtind, kd, ku; seed=16)
    uic_labels = filter(n -> occursin("UIC_", string(n)), names(tpm))
    hik_labels = filter(n -> occursin("hiK_", string(n)), names(tpm))
    

    U = Matrix{Float64}(tpm[filtind, uic_labels])
    H = Matrix{Float64}(tpm[filtind, hik_labels])

    u_ge_h = vec(mean(U .- H, dims=2)) .≥ 0
    ind_u_ge_h = falses(length(filtind))
    ind_u_ge_h[filtind] .= u_ge_h


    UH_D = maxnorm([U H][  u_ge_h, :])
    UH_U = maxnorm([U H][.!u_ge_h, :])

    KMK_D = kmeansorder(UH_D', kd, seed)
    KMK_U = kmeansorder(UH_U', ku, seed)

    KMK_U = reordercluster(KMK_U, [1, 2, 4, 3]); ## for seed = 16, reorder so that the transcriptional target cluster is last


    #### Make cluster table
    genename = last.(split.(tpm.Gene, "|"))
    gtable = DataFrame(Gene=tpm.Gene, GeneName=genename, Named=occursin.("|", tpm.Gene) .& .!occursin.(r"^LOC|^Xetrov", genename),
                       Ind=filtind, D_ind = ind_u_ge_h, U_ind = .!ind_u_ge_h, ClusterD=zeros(Int, length(filtind)), ClusterU=zeros(Int, length(filtind)))
    
    gtable.ClusterD[filtind .&   ind_u_ge_h] = KMK_D.A
    gtable.ClusterU[filtind .& .!ind_u_ge_h] = KMK_U.A

    addtable = tpm[!, Not(r"^Gene|^UIC_|^hiK_")]
    if !isempty(addtable)
        gtable = [gtable addtable]
    end

    

    (KMK_D=KMK_D, KMK_U=KMK_U, gtable=gtable, meta=meta)

end

function cluster_silhouettes_ud(meta, tpm, filtind, ks=2:20; seed=16)
    uic_labels = filter(n -> occursin("UIC_", string(n)), names(tpm))
    hik_labels = filter(n -> occursin("hiK_", string(n)), names(tpm))
    

    U = Matrix{Float64}(tpm[filtind, uic_labels])
    H = Matrix{Float64}(tpm[filtind, hik_labels])

    u_ge_h = vec(mean(U .- H, dims=2)) .≥ 0
    ind_u_ge_h = falses(length(filtind))
    ind_u_ge_h[filtind] .= u_ge_h


    UH_D = maxnorm([U H][  u_ge_h, :])
    UH_U = maxnorm([U H][.!u_ge_h, :])

    pw_D = pairwise(Euclidean(), UH_D')
    pw_U = pairwise(Euclidean(), UH_U')

    
    KMK_D = @showprogress [kmeansorder(UH_D', k, seed) for k in ks]
    KMK_U = @showprogress [kmeansorder(UH_U', k, seed) for k in ks]

    

    silhouettes_D = [silhouettes(KMK.KM, pw_D) for KMK in KMK_D]
    silhouettes_U = [silhouettes(KMK.KM, pw_U) for KMK in KMK_U]
    
    psd = boxplot(ks, silhouettes_D, leg=false, xticks=ks, title="Repressed")
    psu = boxplot(ks, silhouettes_U, leg=false, xticks=ks, title="Activated")

    

    plot(psd, psu, layout=(1, 2), size=(600, 300))
end


"""
    viscluster_double(xp, KMK ; dy=2, kwargs...)

    Stacked cluster heatmap of U and H
"""
function viscluster_double(xp, KMK ; dy=2, kwargs...)
    X = KMK.X'
    delta = xp[2] - xp[1]
    dxp = [xp ; xp[end] .+ delta .+ xp]
    

    heatmap(dxp, 1:dy:size(X, 1), average_heatmap(X[KMK.SI, :], 1, dy), yflip=true; kwargs...)
    vline!([last(xp) .+ delta], c=:white, lab="")
    cn = KMK.KM.counts[KMK.KSI]
    hline!(cumsum(cn), c=:white, lab="")
    yticks!(cumsum(cn) .- cn/2, string.(1:KMK.k))
end


"""
    plotclustermeansd(KMK, meta; X = KMK.X, titlelabel = "", labels = string.(titlelabel, meta.Label), lt = :row2, clusterby=["U", "H"], timescale=true, stageticks=true, minplot=4, kwargs...)

    Pointwise cluster summaries mean + std
"""
function plotclustermeansd(KMK, meta; X = KMK.X, titlelabel = "", labels = string.(titlelabel, meta.Label), lt = :row2, clusterby=["U", "H"], timescale=true, stageticks=true, minplot=4, kwargs...)

    
    phs = Plots.Plot[]
    cc = [:steelblue, :orange]
    for k = 1:KMK.k
        
        M = vec(mean(X[:, KMK.A .== k], dims=2))
        S = vec( std(X[:, KMK.A .== k], dims=2))
    
        p = plot(grid=false, title=string(titlelabel, " C", k, " : ", sum(KMK.A .== k)), titlefont=font(9, "helvetica"), framestyle=:box; kwargs...)
        for (i, cl) in enumerate(clusterby)
            if cl == "U"
                ind = meta.Treat .== "UIC"
            elseif cl == "H"
                ind = meta.Treat .== "hiK"
            end
            
            t = meta.Time[ind]
            xts = 1:2:13
            if timescale && !stageticks
                t = t./2 .+ 4
                xts = 4:2:10
            end
            
            if stageticks
                xts = stageticks_short()
            end

            plot!(t, M[ind], lab=ifelse(k == 1, cl, ""), ribbon=S[ind],  fillalpha=0.21, fillcolor=cc[i], c=cc[i], xticks=xts, xlims=extrema(t))
        end
        
        push!(phs, p)
    end

    if length(phs) < minplot
        push!(phs, plot(axis=false, framestyle=nothing, grid=false, xticks=false))
    end

    [plot!(p, left_margin=-1mm, right_margin=-1mm, yformatter=y->"", yticks=false) for p in phs[2:end]]
    plot!(phs[1], right_margin=-1mm, yticks=0:1)

    plot(phs..., layout=(1, length(phs)), ylims=(0, 1.05))

end

function plotclusterud_stack(KMK_UD)

    
    phu = viscluster_double(KMK_UD.meta.Time[KMK_UD.meta.Treat .== "UIC"], KMK_UD.KMK_U, c=:viridis, colorbar=false, ylabel="Activated",  titlefont=font(9, "helvetica"), xticks=false, bottom_margin=-2mm, title="UIC               hiK")
    phd = viscluster_double(KMK_UD.meta.Time[KMK_UD.meta.Treat .== "UIC"], KMK_UD.KMK_D, c=:viridis, colorbar=false, ylabel="Repressed",  titlefont=font(9, "helvetica"), top_margin=-1mm, xticks=stageticks_short_double(), xlabel="Stage", bottom_margin=3mm)

    pcd = plotclustermeansd(KMK_UD.KMK_D, KMK_UD.meta, titlelabel="Down", xlabel="Stage", bottom_margin=3mm)
    pcu = plotclustermeansd(KMK_UD.KMK_U, KMK_UD.meta, titlelabel="Up")


    
    lta = @layout [a ; b{0.3h}]
    ltb = @layout [a{0.2w} b{0.05w} [ c; d]]
    pa = plot(phu, phd, layout=lta)

    h2 = scatter([0,0], [0,1], zcolor=[0,1], clims=(0, 1),  xlims=(1,1.1), label="", c=:viridis, framestyle=:none) ## trick to add colorbar
    plot(pa, h2, pcu, pcd, layout=ltb, size=(900, 300), xguidefont=font(8, "helvetica"))
    plot!(size=(800*1.1, 300*1.1), fontfamily="helvetica")
end

enrichr_term_selection() = [ 
        ("mTOR signaling pathway", "mTOR", "KEGG_2019_Human"), 
        ("Autophagy", "Autophagy", "KEGG_2019_Human"), 
      

      ("ubiquitin-protein transferase activity (GO:0004842)", "Ubiquitin", "GO_Molecular_Function_2018"),
    
      ("Protein processing in endoplasmic reticulum", "ER", "KEGG_2019_Human"),
      ("Spliceosome", "Splicesome", "KEGG_2019_Human"), 
      ("Wnt signaling pathway", "WNT Signalling", "KEGG_2019_Human"),
      ("Signaling pathways regulating pluripotency of stem cells", "Pluripotency Signalling", "KEGG_2019_Human"),
      ("Endoderm Differentiation WP2853", "Mesendoderm", "WikiPathways_2019_Human")] 

function plot_cluster_enrich_vert(KMK_UD, en_ud, tgs, hommot, args... ; ec=:viridis, kwargs...)
    phu = viscluster_double(KMK_UD.meta.Time[KMK_UD.meta.Treat .== "UIC"], KMK_UD.KMK_U, c=:viridis, colorbar=true, ylabel="Up", titlefont=font(10, "helvetica"),xticks=stageticks_short_double(), xlabel="Stage", bottom_margin=3mm, title="UIC               hiK")

    phh = enrich_bubble(KMK_UD.KMK_U, en_ud, tgs, colorbar=true, framestyle=:box, title="Gene Set Enrichment", ylabel="Cluster")
    
    
    pcd, phs = plotclustermeansd_vert(KMK_UD.KMK_U, KMK_UD.meta)# xguidefont=font(2, "helvetica"))
    plot!(phs[1], title="Cluster Mean")
    plot!(phs[1], bottom_margin=-1mm, xticks=false)
    [plot!(p, top_margin=-1mm, bottom_margin=-1mm,  xticks=false) for p in phs[2:end-1]]
    plot!(phs[end], top_margin=-1mm, xlabel="Stage")
    pmot = motif_bubble(hommot)
    plot!(title="Motif Enrichment", ylabel="Cluster")


    plot!(phs[end], bottom_margin=5mm)
    plot!(phu,      bottom_margin=5mm)
    plot!(phh,      bottom_margin=5mm, xrotation=0)
    plot!(pmot,      bottom_margin=5mm, xrotation=0)
    

    fgs = ["atg13", "rictor", "foxh1", "sox3",  "pou5f3.3", "ventx1.2", "smad4.1", "eomes"]
    ytf = [[0, 25], [0, 45], [0, 1900], [0, 1200], [0, 3000], [0, 500], [0, 100], [0, 450]]
    ghs = [plotgpset(g, fgenes, GPC, showparams=false, samples=["U", "H"], fontfamily="helvetica") for g in fgs]
    for (i, (p, l, yt)) in enumerate(zip(ghs, fgs, ytf))
        cl = @subset(KMK_UD.gtable, occursin.(l, :Gene)).ClusterU[1]
        plot!(p, title=string(l, ":", cl), xticks=stageticks_short(), xlabel="Stage", framestyle=:box, grid=false, bottom_margin=5mm, fontfamily="helvetica") 
        ylm = ylims()[2]
        plot!(p, right_margin=-1mm,left_margin=-1mm, ylims=(0, ylm), yticks=yt)
        (i == 1) && plot!(p, ylabel="TPM", left_margin=5mm)
    end

    lt = @layout [ a{0.425w} grid(4, 1) e{0.425w}]
    lt = @layout [ a{0.25w} grid(4, 1) e{0.425w}]
    lt = @layout [[ a{0.3w} grid(4, 1){0.1w} [d 
                                              e]]
                            grid(1, length(ghs)){0.2h}]

    
    plot(phu, phs..., phh,  pmot, ghs..., layout=lt ; kwargs...)

end

function clabel(cd, cu)
   
    if (cd == 0) && (cu == 0)
        return "N"
    elseif (cd == 0) && (cu != 0)
        return string("U", cu)
    elseif (cd != 0) && (cu == 0)
        return string("D", cd)
    else
       error("CL Assignment: $cd, $cu") 
    end
        
end


function savepromoterclusterbed(promclusters;  promlen=500, projdir=getprojectdir())
    prombed = @subset(promclusters[!, [:chrom, :PromoterStart, :PromoterStop, :CLabel, :strand, :Length, :Gene, :GeneName, :ClusterU, :ClusterD]], :Length .== promlen); ## subset just excludes genes on short scaffolds
    prombed.PromoterStart .-= 1

    resultsdir = joinpath(projdir, "results")
    mkpath(resultsdir)
    filepath = joinpath(resultsdir, "pot_xen_hik_cluster_promoters.500.bed")
    CSV.write(filepath, prombed, delim='\t', header=false)

end


function bubbleplot(C, S=C; mso=5, msm=2, sp=0.5, kwargs...)
    @assert size(C) == size(S)
    XY = [(j, i) for i = 1:size(C, 1), j=1:size(C, 2)][:]    
    scatter(first.(XY), last.(XY), ms=msm.*sqrt.(S[:]) .+ mso, zcolor=C[:], yflip=true, grid=false, xlims=(1 - sp, size(C, 2) + sp), ylims=(1 - sp, size(C, 1)+sp), marker=stroke(0); kwargs...)
end

function enrich_bubble(KMK, en_ud, tgs ; yts=[1, 2, 3, 4], c=cgrad(:PuOr_4, rev=true), kwargs...)
    
    P = mapreduce(tg -> get.(Ref(get_field_set(en_ud, first(tg), last(tg), :padj)), 1:KMK.k, 1.0), hcat, tgs)


    field = :Combined;  min_t = 30; mso = 0.85;  msm = 3.5; ladder = [30, 20, 10, 1]
    # field = :zscore;  min_t = 15; mso = 1;  msm = 5.5; ladder = [15, 10, 5, 0]

    C = mapreduce(tg -> get.(Ref(get_field_set(en_ud, first(tg), last(tg), field)), 1:KMK.k, 0.0), hcat, tgs)
    C = min.(min_t, C)
    
    PT = min.(7, abs.(-log10.(P)))
    
    PT = [PT zeros(size(PT, 1))]
    C  = [C ladder]
    

    p =  bubbleplot(PT, C, mso=mso, msm=msm, lab="", xticks=(1:length(tgs), replace.(getindex.(tgs, 2), Ref(" " => "\n"))), xrotation=45, yticks=(1:length(yts), string.(yts)), c=cgrad(:PRGn_4, rev=false); kwargs...)
    plot!(xlims=(0.5, 9.5))
    p
end

function get_field_set(en_ud, term, geneset, field ; disp =false)
    tbl = @subset(en_ud, :Term .== term, :GeneSet .== geneset, :Direction .== "Up") 
    disp && showwide(tbl)
    Dict(parse(Int, cl) => pv for (cl, pv) in zip(replace.(tbl.CLC, Ref("," => "")), tbl[!, field]))
end


function plotclustermeansd_vert(KMK, meta; X = KMK.X, titlelabel = "", labels = string.(titlelabel, meta.Label), lt = :row2, clusterby=["U", "H"], timescale=true, stageticks=true, minplot=4, kwargs...)

    # pointwise mean and sd M and SD
    phs = Plots.Plot[]
    cc = [:steelblue, :orange]
    for k = 1:KMK.k
        
        M = vec(mean(X[:, KMK.A .== k], dims=2))
        S = vec( std(X[:, KMK.A .== k], dims=2))
    
        p = plot(grid=false, ylabel=string(titlelabel, " C", k, " : ", sum(KMK.A .== k)), titlefont=font(9, "helvetica"), yticks=0:1, ylims=(0, 1.05), framestyle=:box; kwargs...)
        for (i, cl) in enumerate(clusterby)
            if cl == "U"
                ind = meta.Treat .== "UIC"
            elseif cl == "H"
                ind = meta.Treat .== "hiK"
            end
                
            t = meta.Time[ind]
            xts = 1:2:13
            if timescale && !stageticks
                t = t./2 .+ 4
                xts = 4:2:10
            end
            
            if stageticks
                xts = stageticks_short()
            end


            plot!(t, M[ind], lab=ifelse(k == 1, cl, ""), ribbon=S[ind],  fillalpha=0.21, fillcolor=cc[i], c=cc[i], xticks=xts, xlims=extrema(t))
        end
        
        push!(phs, p)
    end

    if length(phs) < minplot
        push!(phs, plot(axis=false, framestyle=nothing, grid=false, xticks=false))

    end

    plot(phs..., layout=(length(phs), 1)), phs

end


function motif_bubble(hommot; c=cgrad(:PuOr_4, rev=true)) #cgrad(:PRGn_4, rev=false))


    yp = 1:7
    ypi = 
    M = Matrix{Float64}(hommot[yp, r"hm_U_[1234]_bgT"])' 
    P = log2.(Matrix{Float64}(hommot[yp, r"fc_U_[1234]_bgT"])' )
    
    P[P .< 0] .= 0

    KMK = kmeansorder(M, 4);
    
    MP = [min.(M[:, KMK.SI], 30) ones(size(M, 1))]
    MP = [min.(M[:, KMK.SI], 25) ones(size(M, 1))]
    PP = [min.(P[:, KMK.SI], 0.75) log2.([1.6, 1.4, 1.2, 1.0])]

    bubbleplot(MP, PP, mso=2, msm=25.0, size=(400, 300) ,colorbar=true, framestyle=:box, c=c, fontfamily="helvetica", lab="")
    yls = first.(split.(hommot.Motif[yp][KMK.SI], "/"))
    yls = replace.(yls, Ref("Atf1" => "CRE/Atf1"))
    yls = replace.(yls, Ref(r"\([0-9a-zA-Z,]*\)" => ""))
    plot!(yticks=1:4, xticks=(yp, yls))

    
end

function motif_supplemental(hommots; ind=1:15)
    seltfu = hommots.sigtable[!, r"Motif|U_[1234]_bgT$|hm_max"] |> excludefam 
    X = Matrix{Float64}(seltfu[ind, r"hm_U_[1234]_bgT"]);
    
    XC  = deepcopy(X)
    vind = [argmax(r) for r in eachrow(X)] .== 4
    XC[vind, :] = XC[vind, :].*[0 0 0 1]
    # display(XC)

    osi = sortperm(DataFrame(XC, :auto), rev=true)
    heatmap(X[osi, :], yflip=true, clims=(0, 20), yticks=(ind, replace.(first.(split.(seltfu.Motif[ind][osi], r"[\|/(]")), Ref("([A-Za-z]*)" => ""))), size=(200, 310), grid=false, fontfamily="helvetica", c=:inferno)

end


function seqlogo_selection(hommots)
    motifs = hommots.sigtable.Motif[hommots.sigtable.Sel]
    for m in motifs
        seqlogo(hommots.motifdict[m], title=first(split(m, r"[\|/(]")), size=(200, 100), titlefont=font(12, "helvetica"), grid=false, framestyle=:none) |> display
     end
end

function calcium_responsive_tf_enrichment(KMK_UD, en_ud)
    clustercounts = DataFrame(Direction=repeat(["Up", "Down"], inner=4), CLC=string.([1:4 ; 1:4]), ClusterCount=[KMK_UD.KMK_U.KM.counts[KMK_UD.KMK_U.KSI] ; KMK_UD.KMK_D.KM.counts[KMK_UD.KMK_D.KSI]])
    crebm_df = innerjoin(@subset(en_ud, occursin.(r"CRE[BM]|ETS[ 0-9A-Za-z_-]*JURKAT", :Term), occursin.("ChIP-Seq", :Term), .!occursin.(",", :CLC)), clustercounts, on=[:Direction, :CLC])
    
    cc = hcat(cgrad(:viridis, 4, categorical=true)[1:end]...)
    pp = @with @subset(crebm_df, :Direction .== "Up") groupedbar(replace.(:Term, Ref(" " => "\n")), -log10.(:padj), group=:CLC, bottom_margin=10mm, leg=^(:topleft), c=cc)
    plot!(ylabel="-log10 FDR")
    pc = @with @subset(crebm_df, :Direction .== "Up") groupedbar(replace.(:Term, Ref(" " => "\n")), 100*:NumGenes./:ClusterCount, group=:CLC, bottom_margin=10mm, leg=^(:topleft), c=cc)
    plot!(ylabel="Percent genes activated", title="") #title="Intersection of genes in activated clusters found in proximity to  public CREM/CREB ChIP-seq")
    p = plot(pc, pp, fontfamily="helvetica", size=(750, 340))
end