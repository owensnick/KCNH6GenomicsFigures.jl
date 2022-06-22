
### savedisplay save svg and display current figure
function savedisplay(label, fmts=[".svg", ".png"], projdir=getprojectdir(), verbose=true)
    p = plot!(fontfamily="helvetica")
    files = joinpath.(projdir, "figures", string.(label, fmts))
    verbose && println("Writing ", first(files))
    mkpath(dirname(first(files)))
    map(savefig, files)
    
    display(p)
end


function plotspikemodel(models, tables, labels, pred_field=:Conc; pred_trans=x -> log(x + 1), plot_model=true, kwargs...)
    phs = [plotmodel(m, t, pred_field, l, plot_model=plot_model) for (m, t, l) in zip(models, tables, labels)]
    for p in phs
       plot!(p, [0,5], [0, 5], c=:black, lab="", colorbar=false, xlabel="log Spike Conc", framestyle=:box, ylabel=ifelse(plot_model, "log Corrected TPM", "log TPM"), xlims=(0, 5), ylims=(0, 5))
    end
    pl = plotmodel(models[end], tables[end], pred_field, labels[end], plot_model=plot_model)
    plot!(framestyle=:none, colorbar=true, xlims=(-200, -199), title="GC")
    lt = @layout [a b d{0.1w}]
    plot(phs..., pl, layout=lt, size=(600, 300), link=:y, top_margin=5mm, titlefont=font(10, "helvetica"), fontfamily="helvetica" ; kwargs...)
end


function plotmodel(model, data, pred_field, label; pred_trans=x->log10(x + 1), plot_model=true)

    if plot_model
        y = exp.(GLM.predict(model))
        yp = GLM.predict(model)
    else
        y = data[!, model.mf.f.rhs.terms[2].args_parsed[1].sym]
        yp = log.(y)
    end
    yt = log.(data[!, pred_field])
    # yt = data[!, pred_field]

    s = 1/log(10)

    r2e = 1 - sum((yt.*s .- yp.*s).^2)./sum((mean(yt.*s) .- yp.*s).^2)   
    
    if !plot_model || (typeof(model.model) <: GeneralizedLinearModel)
        r2s = string("R2: ", round(r2e, digits=3))
    else
        r2s = string("R2: ", round(r2(model), digits=3)) #, " | R2E: ", round(r2e, digits=3))
    end

    # title = string(label, "\nDev: ", round(deviance(model), digits=3), "  |  ", r2s)
    title = string(label, "\n", r2s)
    @with data scatter(pred_trans.(cols(pred_field)), pred_trans.(y), group=:ID, zcolor=:GC, marker=(stroke(0), 0.2), fmt=^(:png), leg=false, colorbar=true, c=^(cgrad(:blues, rev=false)), title=title, titlefont=font(11, "helvetica"))
    

end



############## plot gene 

findgene(gene, table::DataFrame; field=1) = findgene(gene, table[!, field])
function findgene(gene, genes)
    ind = findall(occursin.(gene, genes))
    isempty(ind) && error("$gene not found")
    (length(ind) > 1) && println("Multiple matches for $gene:\n$(genes[ind])")
    first(ind)
end

plotgene(args...; kwargs...) = (plot(); plotgene!(args...; kwargs...))
function plotgene!(gene, meta, tpm, altids=[] ; kwargs...)
    if isempty(altids)
        plotgene!(findgene(gene, tpm), meta, tpm; kwargs...)
    else
        # ind = findene(gene, altids)
        ind = findall(occursin.(gene, altids.AltID))

        if isempty(ind)
            plotgene!(findgene(gene, tpm), meta, tpm; kwargs...)
        else
            (length(ind) > 1) && println("Multiple matches for $gene:\n$(genes[ind])")
            plotgene!(findgene(altids.GeneID[first(ind)], tpm), meta, tpm; suptitle=string(gene, "|"), kwargs...)
        end
    end
end
function plotgene!(gene::Int, meta, tpm ; suptitle="", samples=["U", "H"], cc=[:steelblue, :orange], kwargs...)

    x = [tpm[gene, Symbol(i)] for i in meta.Label]

    for (c, s) in zip(cc, samples)
        ind = startswith.(meta.Treat, uppercase(s)) .| startswith.(meta.Treat, lowercase(s))
        plotgeneset!(meta.Time, x, ind, c, s; kwargs...)
    end
    plot!(title=string(suptitle, tpm[gene, 1]) ; kwargs...)
    ylims!(0, last(ylims()))
end
function plotgeneset!(xp, x, ind, c=:auto, l="" ; kwargs...)
    if sum(ind) > 13
        plot!(xp[ind], x[ind], c=c, lab=l ; kwargs...)
    else
        plot!(xp[ind], x[ind], marker=:auto, c=c, lab=l ; kwargs...)
    end
end

#### tick marks time to stage

function stageticks_long()
    data = [(1, 8), (2, 8.5), (4, 9), (6.5, 10), (8, 10.5), (10, 11), (12.5, 12)]
    Float64.(first.(data)), string.(last.(data))
end


function stageticks_long_double()
    data = [(1, 8), (4, 9), (6.5, 10), (10, 11), (12.5, 12),
            (14, 8), (17, 9), (6.5+13, 10), (10+13, 11), (12.5+13, 12)]
    Float64.(first.(data)), string.(last.(data))
end

function stageticks_short()
    data = [(1, 8), (4, 9), (6.5, 10), (10, 11), (12.5, 12)]
    Float64.(first.(data)), string.(last.(data))
end

function stageticks_short_double()
    data = [(1, 8), (10, 11), (14, 8), (23, 11)]
    Float64.(first.(data)), string.(last.(data))
end


### QC plots

function pairwisecorrelationplot(meta, tpm ; kwargs...)  
    M = Matrix(tpm[!, meta.Label])
    heatmap(corspearman(M) ; kwargs...)
    vline!([13.5], c=:black, lab="")
    hline!([13.5], c=:black, lab="")
    
    xticks!(stageticks_short_double()...)
    yticks!(stageticks_short_double()...)
end


function pcaplot(meta, tpm ; kwargs...)
    
    M = log10.(Matrix(tpm[!, meta.Label]) .+ 1/3)
    f = fit(PCA, M)
    TD = MultivariateStats.transform(f, M)'
    scatter(TD[:, 1], TD[:, 2], group=meta.Treat, marker=(stroke(0), :auto), zcolor=meta.Time, c=cgrad(:viridis, rev=false) ; kwargs...)
    ofs = ifelse.(meta.Treat .== "UIC", -1, 1)*4
    explained = principalvars(f)/tvar(f)
    xlabel!(string("PC 1: ", round(explained[1], digits=2)))
    ylabel!(string("PC 2: ", round(explained[2], digits=2)))
    annotate!([(TD[i, 1] + ofs[i], TD[i, 2], text(meta.Time[i], font(10, "helvetica", ifelse(ofs[i] > 0, :black, :red)))) for i = 1:size(TD, 1)])
end
