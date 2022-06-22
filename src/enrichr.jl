#### functions to use Enrichr API

function human_terms(gts)
    gts = replace.(gts, Ref("POU5F3" => "POU5F1"))
    gts = replace.(gts, Ref("MIX1" => "MIXL1"))
    gts = replace.(gts, Ref("DPPA2" => "DPPA4"))
    gts = replace.(gts, Ref("LEFTY" => "LEFTY2"))
    gts = replace.(gts, Ref(r"VENTX[0-9]" => "NANOG"))
    gts = replace.(gts, Ref("MESPB" => "MESP1"))
    gts = replace.(gts, Ref("SOX17B" => "SOX17"))
    gts = replace.(gts, Ref("SOX17A" => "SOX17"))
    gts = unique(gts)
end

function xenopus_terms(gts)
    gts = replace.(gts, Ref("POU5F1" => "POU5F3"))
    gts = replace.(gts, Ref("MIXL1" => "MIXL1"))
    gts = replace.(gts, Ref("DPPA4" => "DPPA2"))
    gts = replace.(gts, Ref("LEFTY2" => "LEFTY"))
    gts = replace.(gts, Ref("NANOG" => "VENTX"))
    gts = replace.(gts, Ref("MESP1" => "MESPB"))
    gts = replace.(gts, Ref("SOX17" => "SOX17A"))
    gts = unique(gts)
end

function get_genesets()
    gs = ["KEGG_2019_Human",
            "WikiPathways_2019_Human", 
            "BioPlanet_2019", 
            "GO_Biological_Process_2018", 
            "GO_Molecular_Function_2018", 
            "GO_Cellular_Component_2018",
            "ChEA_2016"]
end

clustercombos() = [[1], [2], [3], [4], [1, 2], [2, 3], [3, 4], [1, 2, 3], [2, 3, 4], [1, 2, 4], [1, 2, 3, 4]];
function enrichr_ud(KMK_UD, genesets=get_genesets();  u_clsets = clustercombos(), d_clsets = clustercombos(), humanise=true)

    
    dfs = DataFrame[]

    for (cls, f, label) in zip([u_clsets, d_clsets], [:ClusterU, :ClusterD], ["Up", "Down"])
        @showprogress for g in cls
            ind = KMK_UD.gtable[!, f] .∈ Ref(g)
            uptbl = KMK_UD.gtable[ind, :]
            engenes = filter(f -> !occursin(r"LOC|gene[0-9]*", f), uptbl.GeneName) .|> (x -> replace(x, r"\.[0-9]$" => "")) .|> uppercase
            if humanise
                engenes = human_terms(engenes)
            end
            config = load_genelist(engenes)
            for gs in genesets
                et = enrichment(config, gs)
                et[:, :Direction] .= label
                et[:, :CLC] .= join(g, ",")
                push!(dfs, et)
            end
        end
    end
    en = reduce(vcat, dfs)   
    en.NumGenes = length.(en.Genes)

    en
end