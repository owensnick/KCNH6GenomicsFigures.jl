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