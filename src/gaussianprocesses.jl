### data transforms
dtm(m, α=1, β=1000) = x -> sqrt(β*(x/m) + α)
itm(m, α=1, β=1000) = y -> ifelse(y > sqrt(α), m*(y*y - α)/β, 0)
sqd(x) = ifelse(x > 0, x*x, 0);

function gp_set_dt(gi, meta, tpm; α=1, β=1000)
    t = Float64.(unique(meta.Time))
    @assert t == 1:13
    st = range(0, 13, length=100)
    yu = [v for v in tpm[gi, r"UIC_[0-9]+"]]::Vector{Float64}
    yh = [v for v in tpm[gi, r"hiK_[0-9]+"]]::Vector{Float64}

    tf = [t ;  t]
    yf = [yu ;  yh];

    m = maximum(yf)
    

    dt = dtm(m, α, β)
    it = itm(m, α, β)
    
    gpru = gp_reg(t, yu, st, datatrans=dt, invtrans=it);
    gprh = gp_reg(t, yh, st, datatrans=dt, invtrans=it);
    gprf = gp_reg(tf, yf, st, datatrans=dt, invtrans=it);

    (u=gpru, h=gprh, f=gprf)
end

function gp_reg(t, y, st; initparams = (1.0, 0.0, 4.0), datatrans=sqrt, invtrans=sqd, f_prior = Normal(1.4, 4.0), ℓ_prior = Normal(1.2, 1.0), n_prior=Normal(1.0, .75))
    σf, ℓ, σn = initparams

    kern = Mat52Iso(ℓ, σf)
    set_priors!(kern, [ℓ_prior, f_prior])
    gp = GP(t, datatrans.(y), MeanZero(), kern, σn)
    set_priors!(gp.logNoise, [n_prior] )
    
    @suppress_out optimize!(gp) ## suppress std out during optimisation to avoid pos def errors encountered
    μ, v = predict_y(gp, st)
    sf  = invtrans.(μ)
    cil = invtrans.(μ .- 1.96*sqrt.(v))
    ciu = invtrans.(μ .+ 1.96*sqrt.(v))

    (gp=gp, μ=μ, v=v, sf=sf, cil=cil, ciu=ciu, t=t, y=y, st=st)
end

function gptable(genes, GPC; dt = 3)
    
    st = GPC[1].u.st[1:dt:end]
    SF = mapreduce(g -> mapreduce(f -> Float64.(g[f].sf[1:dt:end]), vcat, [:u, :h]), hcat, GPC);
    metagp = DataFrame(Treat=repeat(["UIC", "hiK"], inner=length(st)), Time=repeat(st, outer=2))
    metagp.Label = string.(metagp.Treat, "_", repeat(1:length(st), outer=2))
    tpmc_gp = [DataFrame(Gene=genes) DataFrame(SF', Symbol.(metagp.Label))];
    
    
    gpc_pt = innerjoin(paramtable(genes, GPC), gp_de_table(genes, GPC), on=[:Gene, :Index])
    
    metagp, tpmc_gp, gpc_pt

end

function gp_de_table(genes, GPC)
    lr_bic = [gp_bic(gp.f, (gp.u, gp.h)) for gp in GPC]
    DataFrame(Gene=genes, Index=1:length(genes), LR=first.(lr_bic), BIC=last.(lr_bic), CD=maxcohensd.(GPC))
end

function paramtable(genes, G, fields=[:u, :h, :f])
    σf = [[sqrt.(g[f].gp.kernel.σ2)      for g in G] for f in fields]
    ℓ  = [[g[f].gp.kernel.ℓ       for g in G] for f in fields]
    σn = [[exp.(g[f].gp.logNoise.value) for g in G] for f in fields]
    SNR = [f./n for (f, n) in zip(σf, σn)]


    df = DataFrame(σf, Symbol.(string.(:σf, "_", fields)))
    dl = DataFrame(ℓ, Symbol.(string.(:ℓ, "_", fields)))
    dn = DataFrame(σn, Symbol.(string.(:σn, "_", fields)))
    ds = DataFrame(SNR, Symbol.(string.(:SNR, "_", fields)))
    [DataFrame(Gene=genes, Index=1:length(genes)) df dl dn ds]
end
### gaussianproces plotting

### plot DE summaries

function lrhistogram(gpc_pt)
    bins = range(-20, 20, length=150)

    ph = plot()
    stephist!(gpc_pt.LR, fill=0, fillalpha=0.2, bins=bins, c=:steelblue, lab=string("LR > 0 : ", sum(gpc_pt.LR .> 0)))
    vline!([0], c=:black, lab="", ls=:dash)
    plot!(xlabel="LR", title="Histogram of LR")


    pe = plot(leg=:topright, xlabel="LR", title="Inverse cumulative LR")
    plot!(bins, size(gpc_pt, 1)*(1 .- ecdf(gpc_pt.LR)(bins)), c=:black, lab=string("LR > 0 : ", sum(gpc_pt.LR .> 0)))
    vline!([0], c=:black, lab="", ls=:dash)


    hline!([sum(gpc_pt.LR .> 0)], c=:black, lab="")
    annotate!([(20, 750 + sum(gpc_pt.LR .> 0), text(string(sum(gpc_pt.LR .> 0), " genes"), font(:black, :right)))])

    plot(ph, pe, size=(375*2*1.00, 300*1.00), bottom_margin=10mm, ylabel="# Genes", left_margin=5mm, framestyle=:box, fontfamily="helvetica")
end


function lr_cd_hist(gpc_pt)
    histogram2d(gpc_pt.LR, gpc_pt.CD, bins=250, ylims=(-20, 20), framestyle=:zerolines, size=(375*1*1.00, 300*1.00),
    xlabel="UIC vs HIK Likelihood Ratio", ylabel="Max Divergence (z-score)", fontfamily="helvetica", title="Differentially Expressed Genes")
end

### plot gene

plotgpset(args...; kwargs...) = (plot(); plotgpset!(args...; kwargs...))
function plotgpset!(gene, genes, gpr, altids=[]; plotall=false, kwargs...)
    
    if isempty(altids)
        
        ind = findall(occursin.(gene, genes))
        
        isempty(ind) && error("$gene not found")
        if length(ind) > 1
            println("Multiple matches for $gene:\n$(genes[ind])")
            phs = [plotgpset(gpr[i]; suptitle=genes[i], kwargs...) for i in ind]
            plot(phs..., size=(1200, 500); kwargs...)
        else
            plotgpset!(gpr[first(ind)] ; suptitle=genes[first(ind)], kwargs...)
        end
        
    else
        # ind = findene(gene, altids)
        ind = findall(occursin.(gene, altids.AltID))

        if isempty(ind)
            ind = findgene(gene, genes)
            plotgpset!(gpr[first(ind)] ; suptitle=genes[first(ind)], kwargs...)
        else
            candgenes = unique(altids.GeneID[ind])
            (length(candgenes) > 1) && println("Multiple matches for $gene:\n$candgenes")
            # @show plotall
            if plotall &&  (length(candgenes) > 1)
                # phs = [plotgpset(c, genes, gpr, altids; kwargs...) for c in candgenes]
                # plot(phs..., size=(1200, 500); kwargs...)
            else
                gind = findgene(first(candgenes), genes)
                plotgpset!(gpr[first(gind)] ; suptitle=string(gene, "|", genes[first(gind)]), kwargs...)
            end
        end
    end

end



function plotgpset!(gpr; samples=["U", "H"], showparams=false, suptitle="", kwargs...)
    p = plot(; kwargs...)
    ("U" ∈ samples) && plotgp!(gpr.u, c=:steelblue, lab="U", showparams=showparams; kwargs...)
    ("H" ∈ samples) && plotgp!(gpr.h, c=:orange, lab="H", showparams=showparams; kwargs...)
    ("F" ∈ samples) && plotgp!(gpr.f, c=:lightgrey, lab="F", showparams=showparams; kwargs...)

    if showparams
        lml_strs = String[]
        ("U" ∈ samples) && push!(lml_strs, string("lml(u) = ", tt(gpr.u.gp.mll)))
        ("H" ∈ samples) && push!(lml_strs, string("lml(h) = ", tt(gpr.h.gp.mll)))

        param_strs = String[]
        ("U" ∈ samples) && push!(param_strs, gpparamstring(gpr, :u))
        ("H" ∈ samples) && push!(param_strs, gpparamstring(gpr, :h))

        s = string(" ", join(lml_strs, " | "), "\n",
                   "lr, bic = ", tt.(gp_bic(gpr.f, [gpr.u, gpr.h])), "\n",
                   join(param_strs, "\n"), "\n")
        plot!(titlefont=font(8), top_margin=10mm)
    else
        lr, bic = gp_bic(gpr.f, [gpr.u, gpr.h])
        s = string("  LR: ", tt(lr))#, ifelse(bic > 0, " ***", ""))
        plot!(titlefont=font(10))
    end
    plot!(title=string(suptitle, s))#, titlefont=font(8, :white), top_margin=10mm)
    p
end


plotgp(args...; kwargs...) = (plot(); plotgp!(args...; kwargs...))
tt(x) = round(x, digits=2)
function plotgp!(gpr; c=:auto, marker=:auto, plotribbon=true, showparams=false, lab="", kwargs...)

    if length(gpr.t) == 13
        plot!(gpr.t, gpr.y, marker=(:circle), c=c, lab=""; kwargs..., alpha=0.75)
    else
        scatter!(gpr.t[1:13], gpr.y[1:13], marker=(:steelblue), c=c, lab="U data"; kwargs...)
        scatter!(gpr.t[14:end], gpr.y[14:end], marker=(:orange), c=c, lab="H data"; kwargs...)
    end

     if plotribbon
         p = plot!(gpr.st, gpr.sf, ribbon=[gpr.sf - gpr.cil gpr.ciu - gpr.sf], c=c, lab=lab; kwargs...)
     else
         p = plot!(gpr.st, gpr.sf,  c=c, lab=lab, kwargs...)
     end

     if showparams
        s = string("LML = ", tt(gpr.gp.mll), "\np = (",
         tt(gpr.gp.kernel.σ2), ", ",
         tt(gpr.gp.kernel.ℓ), ", ",
         tt(exp(2*gpr.gp.logNoise.value)), ")")
        plot!(title=s, titlefont=font(8, :white))
     end
     p
end

function gpparamstring(gpr, field)
    string("p(", field, ") = (",
     tt(gpr[field].gp.kernel.σ2), ", ",
     tt(gpr[field].gp.kernel.ℓ), ", ",
     tt(exp(2*gpr[field].gp.logNoise.value)), ")")
end


# gp stats
function gp_bic(gprf, gprs)
    f_mll = gprf.gp.mll
    s_mll = mapreduce(gpr -> gpr.gp.mll, +, gprs)
    lr = s_mll - f_mll
    bic = lr - 3*log(gprf.gp.nobs)/2
    lr, bic
end

function maxcohensd(gpr, fieldA=:u, fieldB=:h)
    md = 0.0
    mdi = 0
    
    for (μA, μB, vA, vB) in zip(gpr[fieldA].μ, gpr[fieldB].μ, gpr[fieldA].v, gpr[fieldB].v)
        v = abs(μA - μB)/sqrt(0.5*(vA + vB))
        # v = abs(μA - μB)/sqrt(1.0*(vA + vB))
        if v > md
            md = v
            if μA > μB
                mdi = -1
            elseif μA < μB
                mdi = 1
            else
                mdi = 0
            end
        end
    end
    md*mdi
end



function plot_gp_selection(fgs, fgenes, GPC; kwargs...)
    phs = [plotgpset(g, fgenes, GPC, showparams=false, samples=["U", "H"], fontfamily="helvetica") for g in fgs]
    for (i, (p, l)) in enumerate(zip(phs, fgs))
        plot!(p, title=l) 
        ylm = ylims()[2]
        l = max(10^round(log10(ylm) - 1), 2)
        yt = l*round(ylm/l)
        plot!(p, yticks=[0, yt], ylims=(0, max(ylm, yt)))
        plot!(p, right_margin=-1mm)
        (i == 1) && plot!(p, ylabel="TPM", left_margin=5mm)
    end
    plot(phs..., layout=(1, length(phs)), size=(162.5*length(fgs), 200), xticks=stageticks_short(), xlabel="Stage", framestyle=:box, grid=false, bottom_margin=10mm, fontfamily="helvetica"; kwargs...)
end