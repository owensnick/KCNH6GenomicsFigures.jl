
maxnorm(X) = X./maximum(X, dims=2)

function pot_kmeans_ud(meta, tpm, filtind, kd, ku)
    uic_labels = filter(n -> occursin("UIC_", string(n)), names(tpm))
    hik_labels = filter(n -> occursin("hiK_", string(n)), names(tpm))
    

    U = Matrix{Float64}(tpm[filtind, uic_labels])
    H = Matrix{Float64}(tpm[filtind, hik_labels])

    u_ge_h = vec(mean(U .- H, dims=2)) .≥ 0
    ind_u_ge_h = falses(length(filtind))
    ind_u_ge_h[filtind] .= u_ge_h


    UH_D = maxnorm([U H][  u_ge_h, :])
    UH_U = maxnorm([U H][.!u_ge_h, :])

    KMK_D = kmeansorder(UH_D', kd)
    KMK_U = kmeansorder(UH_U', ku)


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

function cluster_silhouettes_ud(meta, tpm, filtind, ks=2:20)
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

    
    KMK_D = @showprogress [kmeansorder(UH_D', k) for k in ks]
    KMK_U = @showprogress [kmeansorder(UH_U', k) for k in ks]

    

    silhouettes_D = [silhouettes(KMK.KM, pw_D) for KMK in KMK_D]
    silhouettes_U = [silhouettes(KMK.KM, pw_U) for KMK in KMK_U]
    
    psd = boxplot(ks, silhouettes_D, leg=false, xticks=ks)
    psu = boxplot(ks, silhouettes_U, leg=false, xticks=ks)

    

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
end