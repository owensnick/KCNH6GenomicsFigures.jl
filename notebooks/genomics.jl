include(joinpath("..", "src", "project.jl"));
ENV["GKSwstype"] = "100" ### this needed for GR when running on server without a display terminal, comment out should you want GR to display figures rather than saving

meta, tpm, isotpm, isoweight, stats, filtind, ids = loaddata();

@show sum(filtind), mean(filtind);

@time sgm = spike_gc_model(isotpm, 2, meta, pseudocount=2, verbose=true); ### build spike model
@time isomodel = isogc_mrna(isotpm, 2, sgm.models, meta, pseudocount=2, verbose=true); ## apply correction
tpmc = isomodel.tpmc; ### spike gc corrected tpm

isoc_weight = select(groupby(combine(groupby(isomodel.igc, [:Gene, :Isoform]), 
            :ModelPredict => mean => :ModelPredict), :Gene), 
            :Isoform, :ModelPredict => (x -> replace(x/sum(x), NaN => 1/length(x))) => :Weight) # calculate reweighted informs
promgene = loadpromoters(isoc_weight); ## load promoter of the maximally expressed isoform.

plotspikemodel(sgm.models, sgm.sgc_tables, ["UIC", "HIK"], plot_model=false, bottom_margin=10mm); savedisplay("sup_fig9b_errcc_uncorrected");
plotspikemodel(sgm.models, sgm.sgc_tables, ["UIC", "HIK"], plot_model=true, bottom_margin=10mm); savedisplay("sup_fig9b_errcc_corrected");


cc = pairwisecorrelationplot(meta, tpmc[filtind, :]) #, clims=(0.3, 1))
plot!(size=(375, 300), title="Pairwise Spearman Correlation", xlabel="Stage", ylabel="Stage", fontfamily="helvetica")
savedisplay("sup_fig9c_pairwise_spearman");

p = pcaplot(meta,  tpmc[filtind, :], grid=true, framestyle=:box, title="Principal Components Analysis", fontfamily="helvetica")
plot!(fontfamily="helvetica", size=(375*1.0, 300*1.0), xlims=(-60, 55))
savedisplay("sup_fig9d_pca");

filtfc(fc, l=-2, u=4) = .!isnan.(fc) .& .!isinf.(fc) .& (l .≤ fc .≤ u)

fct_tpm  = (l=-2.5, u=4.5)
fct_pair = (l=-3, u=3)
exind_fc_tpm  = @with isomodel.modelstats_gene filtfc(:mu_fc_uic, fct_tpm...) .& filtfc(:mu_fc_hik, fct_tpm...)
exind_fc_pair = @with isomodel.modelstats_gene filtfc(:fc_H_U, fct_pair...)


spikeind = occursin.(r"^ERCC", tpmc.Gene)
miss_ind = filtind .& .!spikeind .& (.!exind_fc_tpm .| .!exind_fc_tpm);

@show sum(exind_fc_tpm), sum(filtind), sum(exind_fc_pair)
@show sum(miss_ind)

combine(groupby(DataFrame(LongestRunCondition = filtind, FC_Correction=exind_fc_tpm, FC_condition=exind_fc_pair), [:LongestRunCondition, :FC_Correction, :FC_condition], sort=true), nrow => :count)

exfiltind = filtind .& exind_fc_tpm .& exind_fc_pair .& .!spikeind;
println("Final selected: ", sum(exfiltind), " genes")

tpmc[!, :filtind] = filtind;
f_ind = findall(exfiltind);
fgenes = tpmc.Gene[f_ind];


GPC = @showprogress "Gaussian Process Differential Expression: " [gp_set_dt(i, meta, tpmc) for i ∈ f_ind];

metagp, tpmc_gp, gpc_pt = gptable(fgenes, GPC);

lrhistogram(gpc_pt); savedisplay("lr_histogram")
lr_cd_hist(gpc_pt); savedisplay("sup_fig9e_lr_divergence")

cl_ind = (gpc_pt.LR .> 0)
cluster_silhouettes_ud(metagp, tpmc_gp, cl_ind, 2:15, seed=16)

KMK_UD = pot_kmeans_ud(metagp, tpmc_gp, cl_ind, 4, 4, seed=16);
plotclusterud_stack(KMK_UD)
savedisplay("sup_fig9fg_cluster_overview")

promclusters = leftjoin(KMK_UD.gtable, promgene, on=:Gene) 
promclusters.CLabel = clabel.(promclusters.ClusterD, promclusters.ClusterU);
promclusters.loc = @with promclusters string.(:chrom, ":", :PromoterStart, "-", :PromoterStop);
prombed = @subset(promclusters[!, [:chrom, :PromoterStart, :PromoterStop, :CLabel, :strand, :Length, :Gene, :GeneName, :ClusterU, :ClusterD]], :Length .== 500);
savepromoterclusterbed(promclusters);

hommots = loadhomerres(joinpath(getprojectdir(), "results",  "homer_vertebrate_500"));

@time en_ud = enrichr_ud(KMK_UD);

enc = @subset(en_ud, .!occursin.(",", :CLC))
enc.Genes = join.(enc.Genes, ", ")
rename!(enc, :GeneSet => :GeneSetLibrary, :CLC => :Cluster)
enc = enc[!, [:Direction, :Cluster, :GeneSetLibrary, :Rank, :Term, :pvalue, :zscore, :Combined, :padj, :NumGenes, :Genes]]
enc.Direction = getindex.(Ref(Dict("Up" => "Activated", "Down" => "Repressed")), enc.Direction)
mkpath(joinpath("..", "results"))
CSV.write(joinpath("..", "results", "pot_degenes_enrichr_genesets.tsv.gz"), enc, delim='\t', compress=true);

plot_cluster_enrich_vert(KMK_UD, en_ud, enrichr_term_selection(), hommots.sigtable[hommots.sigtable.Sel, :], size=(900*1.35,500*1.35), fontfamily="helvetica")
savedisplay("fig3")

seqlogo_selection(hommots)

motif_supplemental(hommots, ind=1:16)
savedisplay("sup_fig9h")

fgs = ["ets1" , "crem", "creb1", "atf1"]
plot_gp_selection(fgs, fgenes, GPC)
savedisplay("sup_fig9i_calcium_responsive_tfs")

calcium_responsive_tf_enrichment(KMK_UD, en_ud)
savedisplay("sup_fig9j_calcium_responsive_tfs_enrichments")
