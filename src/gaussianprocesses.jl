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
    
    optimize!(gp)
    μ, v = predict_y(gp, st)
    sf  = invtrans.(μ)
    cil = invtrans.(μ .- 1.96*sqrt.(v))
    ciu = invtrans.(μ .+ 1.96*sqrt.(v))

    (gp=gp, μ=μ, v=v, sf=sf, cil=cil, ciu=ciu, t=t, y=y, st=st)
end