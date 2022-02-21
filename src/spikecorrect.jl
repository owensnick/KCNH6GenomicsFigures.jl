### functions to correct GC bias in sequencing using ERCC spikes



getspikefasta(projdir=getprojectdir()) = joinpath(projdir, "data", "ERCC92.fa")
getspikequantfile(projdir=getprojectdir()) = joinpath(projdir, "data", "spikemeta.tsv")


"""
    spike_gc_model(isotpm, k, meta;  spikefile=getspikefasta(), spikequantfile=getspikequantfile(), showplot=false, pseudocount=1, top_spike_n=20)

    Calculate the full spike correction
"""
function spike_gc_model(isotpm, k, meta;  spikefile=getspikefasta(), spikequantfile=getspikequantfile(), showplot=false, pseudocount=1, top_spike_n=20)

    ### setup labels
    samplelabels = Symbol.(meta.Label)
    label = string("sg_k", k, "_s", pseudocount)

    ### 0. kmer labels
    KML = Symbol.(string.("K_", getkmers(k)))
    println("[SGC]\tk = $k\t|kmers| = ", length(KML))

    ### 1. Spike annotation
    println("[SGC]\tAnnotating spikes...")
    spikemeta = fastacomposition(spikefile, max(k, 1), end_filter=25, pseudocount=pseudocount)
    spikequant = CSV.read(spikequantfile, DataFrame)
    rename!(spikequant, [:NID, :ID, :Group, :Conc, :Mix2, :EFR, :log2_EFR])
    spd = leftjoin(spikemeta, spikequant[!, [:ID, :Conc]], on=:ID)
    @assert iszero(sum(ismissing.(spd.Conc)))
    spd[!, :Conc] = coalesce.(spd.Conc, 0)


    ### 4. Separate spike TPM
    println("[SGC]\tStacking spikes...")
    spike_ind = occursin.(r"ERCC-", isotpm.Gene)
    ssi = sortperm(maximum(Matrix(isotpm[spike_ind, samplelabels]), dims=1) |> vec, rev=true)

    spike_stack = stack(isotpm[spike_ind, [:Gene ; samplelabels]], samplelabels, :Gene, variable_name=:Sample, value_name=:TPM)[!, [:Gene, :Sample, :TPM]]
    rename!(spike_stack, :Gene => :ID)
    sfields = split.(string.(spike_stack.Sample), "_")
    spike_stack[!, :SampleType] = first.(sfields);
    spike_stack[!, :Time]       = parse.(Int, last.(sfields));

    # join stacked spike onto spike data
    spikegc = innerjoin(spike_stack, spd[!, [[:ID, :GC, :Conc] ; KML]], on=:ID)
    @assert size(spikegc, 1) == size(spike_stack, 1)


    ### 6. calc model
    println("[SGC]\tCalculating model...")
    sgc, models, sgc_tables = spikegcmodel(spikegc, k)


    (label=label, sgc=sgc, sgc_tables=sgc_tables, models=models, spd=spd)
end


"""
    fastacomposition(file, k=2; end_filter=0, pseudocount=1)

    Kmer composition of a FASTA file

"""
function fastacomposition(file, k=2; end_filter=0, pseudocount=1)
    reader = open(FASTA.Reader, file)

    seqmeta = DataFrame(ID=String[], A=Int[], C=Int[], G=Int[], T=Int[], len=Int[], GC=Float64[])
    kmers = DNAMer{k}.(UInt.(0:(4^k - 1)))
    KML = Symbol.(string.("K_", kmers))
    kmerfreq = DataFrame([Float64[] for i = 1:length(kmers)], KML)

    record = FASTA.Record()
    while !eof(reader)
        read!(reader, record)
        seq = FASTA.sequence(record)[1:(end - end_filter)]
        comp = composition(seq)
        push!(seqmeta, (FASTA.identifier(record), comp[DNA_A], comp[DNA_C], comp[DNA_G], comp[DNA_T], length(seq), (comp[DNA_C] + comp[DNA_G])/length(seq)))
        CK = composition(each(DNAMer{k}, seq))
        KC = [getindex(CK, k) + pseudocount for k in kmers]
        push!(kmerfreq, KC./sum(KC))
    end
    [seqmeta kmerfreq]
end



"""
    spikegcmodel(spikegc, k; τ=1)

    Setup spike KMer linear model for UIC and hiK independently
"""
function spikegcmodel(spikegc, k; τ=1)
    println("[SGM]\tMaking tables...")
    
    sgc_uic = @subset(spikegc, :TPM .> τ, :ID .!= "ERCC-00116", :SampleType .== "UIC") ## ERCC-00116 performs poorly see https://doi.org/10.1016/j.celrep.2015.12.050
    sgc_hik = @subset(spikegc, :TPM .> τ, :ID .!= "ERCC-00116", :SampleType .== "hiK") ## ERCC-00116 performs poorly see https://doi.org/10.1016/j.celrep.2015.12.050

    sgc_tables = [sgc_uic, sgc_hik]
    
    F = get_kGC_formula(k)
    println("[SGM]\tFormula: ", F)
    println("[SGM]\tFitting models... ")
    models = [lm(F, t) for t in sgc_tables]
    println("[SGM]\tModel fit complete, predicting...")
    for (m, st) in zip(models, sgc_tables)
        st[!, :ModelPredict] = exp.(GLM.predict(m))
        st[!, :FC] = log2.((st.TPM .+ 0.01)./(st.ModelPredict .+ 0.01))
    end

    println("[SGM]\tBuilding tables... ")
    sgc = reduce(vcat, sgc_tables)
    println("[SGM]\tComplete")
    sgc, models, sgc_tables
end



##### kmer functions

getkmers(k) = DNAMer{k}.(UInt.(0:(4^k - 1)))
kformstring(k) = join(string.("log(K_", getkmers(4), ")"), " + ")


### current formulas are hardcode, up to k = 4, TODO automatic formula generation
function get_kGC_formula(k)
    if k == 0
        return @formula log(Conc) ~ log(TPM) + log(GC)
    elseif k == 1
        return @formula log(Conc) ~ log(TPM) + log(K_A) + log(K_C) + log(K_G) + log(K_T)
    elseif k == 2
        return @formula log(Conc) ~ log(TPM) + log(K_AA) + log(K_AC) + log(K_AG) + log(K_AT) + log(K_CA) + log(K_CC) + log(K_CG) + log(K_CT) + log(K_GA) + log(K_GC) + log(K_GG) + log(K_GT) + log(K_TA) + log(K_TC) + log(K_TG) + log(K_TT)
    elseif k == 3
        return @formula log(Conc) ~ log(TPM) + log(K_AAA) + log(K_AAC) + log(K_AAG) + log(K_AAT) + log(K_ACA) + log(K_ACC) + log(K_ACG) + log(K_ACT) + log(K_AGA) + log(K_AGC) + log(K_AGG) + log(K_AGT) + log(K_ATA) + log(K_ATC) + log(K_ATG) + log(K_ATT) + log(K_CAA) + log(K_CAC) + log(K_CAG) + log(K_CAT) + log(K_CCA) + log(K_CCC) + log(K_CCG) + log(K_CCT) + log(K_CGA) + log(K_CGC) + log(K_CGG) + log(K_CGT) + log(K_CTA) + log(K_CTC) + log(K_CTG) + log(K_CTT) + log(K_GAA) + log(K_GAC) + log(K_GAG) + log(K_GAT) + log(K_GCA) + log(K_GCC) + log(K_GCG) + log(K_GCT) + log(K_GGA) + log(K_GGC) + log(K_GGG) + log(K_GGT) + log(K_GTA) + log(K_GTC) + log(K_GTG) + log(K_GTT) + log(K_TAA) + log(K_TAC) + log(K_TAG) + log(K_TAT) + log(K_TCA) + log(K_TCC) + log(K_TCG) + log(K_TCT) + log(K_TGA) + log(K_TGC) + log(K_TGG) + log(K_TGT) + log(K_TTA) + log(K_TTC) + log(K_TTG) + log(K_TTT)
    elseif k == 4
        return @formula log(Conc) ~ log(TPM) + log(K_AAAA) + log(K_AAAC) + log(K_AAAG) + log(K_AAAT) + log(K_AACA) + log(K_AACC) + log(K_AACG) + log(K_AACT) + log(K_AAGA) + log(K_AAGC) + log(K_AAGG) + log(K_AAGT) + log(K_AATA) + log(K_AATC) + log(K_AATG) + log(K_AATT) + log(K_ACAA) + log(K_ACAC) + log(K_ACAG) + log(K_ACAT) + log(K_ACCA) + log(K_ACCC) + log(K_ACCG) + log(K_ACCT) + log(K_ACGA) + log(K_ACGC) + log(K_ACGG) + log(K_ACGT) + log(K_ACTA) + log(K_ACTC) + log(K_ACTG) + log(K_ACTT) + log(K_AGAA) + log(K_AGAC) + log(K_AGAG) + log(K_AGAT) + log(K_AGCA) + log(K_AGCC) + log(K_AGCG) + log(K_AGCT) + log(K_AGGA) + log(K_AGGC) + log(K_AGGG) + log(K_AGGT) + log(K_AGTA) + log(K_AGTC) + log(K_AGTG) + log(K_AGTT) + log(K_ATAA) + log(K_ATAC) + log(K_ATAG) + log(K_ATAT) + log(K_ATCA) + log(K_ATCC) + log(K_ATCG) + log(K_ATCT) + log(K_ATGA) + log(K_ATGC) + log(K_ATGG) + log(K_ATGT) + log(K_ATTA) + log(K_ATTC) + log(K_ATTG) + log(K_ATTT) + log(K_CAAA) + log(K_CAAC) + log(K_CAAG) + log(K_CAAT) + log(K_CACA) + log(K_CACC) + log(K_CACG) + log(K_CACT) + log(K_CAGA) + log(K_CAGC) + log(K_CAGG) + log(K_CAGT) + log(K_CATA) + log(K_CATC) + log(K_CATG) + log(K_CATT) + log(K_CCAA) + log(K_CCAC) + log(K_CCAG) + log(K_CCAT) + log(K_CCCA) + log(K_CCCC) + log(K_CCCG) + log(K_CCCT) + log(K_CCGA) + log(K_CCGC) + log(K_CCGG) + log(K_CCGT) + log(K_CCTA) + log(K_CCTC) + log(K_CCTG) + log(K_CCTT) + log(K_CGAA) + log(K_CGAC) + log(K_CGAG) + log(K_CGAT) + log(K_CGCA) + log(K_CGCC) + log(K_CGCG) + log(K_CGCT) + log(K_CGGA) + log(K_CGGC) + log(K_CGGG) + log(K_CGGT) + log(K_CGTA) + log(K_CGTC) + log(K_CGTG) + log(K_CGTT) + log(K_CTAA) + log(K_CTAC) + log(K_CTAG) + log(K_CTAT) + log(K_CTCA) + log(K_CTCC) + log(K_CTCG) + log(K_CTCT) + log(K_CTGA) + log(K_CTGC) + log(K_CTGG) + log(K_CTGT) + log(K_CTTA) + log(K_CTTC) + log(K_CTTG) + log(K_CTTT) + log(K_GAAA) + log(K_GAAC) + log(K_GAAG) + log(K_GAAT) + log(K_GACA) + log(K_GACC) + log(K_GACG) + log(K_GACT) + log(K_GAGA) + log(K_GAGC) + log(K_GAGG) + log(K_GAGT) + log(K_GATA) + log(K_GATC) + log(K_GATG) + log(K_GATT) + log(K_GCAA) + log(K_GCAC) + log(K_GCAG) + log(K_GCAT) + log(K_GCCA) + log(K_GCCC) + log(K_GCCG) + log(K_GCCT) + log(K_GCGA) + log(K_GCGC) + log(K_GCGG) + log(K_GCGT) + log(K_GCTA) + log(K_GCTC) + log(K_GCTG) + log(K_GCTT) + log(K_GGAA) + log(K_GGAC) + log(K_GGAG) + log(K_GGAT) + log(K_GGCA) + log(K_GGCC) + log(K_GGCG) + log(K_GGCT) + log(K_GGGA) + log(K_GGGC) + log(K_GGGG) + log(K_GGGT) + log(K_GGTA) + log(K_GGTC) + log(K_GGTG) + log(K_GGTT) + log(K_GTAA) + log(K_GTAC) + log(K_GTAG) + log(K_GTAT) + log(K_GTCA) + log(K_GTCC) + log(K_GTCG) + log(K_GTCT) + log(K_GTGA) + log(K_GTGC) + log(K_GTGG) + log(K_GTGT) + log(K_GTTA) + log(K_GTTC) + log(K_GTTG) + log(K_GTTT) + log(K_TAAA) + log(K_TAAC) + log(K_TAAG) + log(K_TAAT) + log(K_TACA) + log(K_TACC) + log(K_TACG) + log(K_TACT) + log(K_TAGA) + log(K_TAGC) + log(K_TAGG) + log(K_TAGT) + log(K_TATA) + log(K_TATC) + log(K_TATG) + log(K_TATT) + log(K_TCAA) + log(K_TCAC) + log(K_TCAG) + log(K_TCAT) + log(K_TCCA) + log(K_TCCC) + log(K_TCCG) + log(K_TCCT) + log(K_TCGA) + log(K_TCGC) + log(K_TCGG) + log(K_TCGT) + log(K_TCTA) + log(K_TCTC) + log(K_TCTG) + log(K_TCTT) + log(K_TGAA) + log(K_TGAC) + log(K_TGAG) + log(K_TGAT) + log(K_TGCA) + log(K_TGCC) + log(K_TGCG) + log(K_TGCT) + log(K_TGGA) + log(K_TGGC) + log(K_TGGG) + log(K_TGGT) + log(K_TGTA) + log(K_TGTC) + log(K_TGTG) + log(K_TGTT) + log(K_TTAA) + log(K_TTAC) + log(K_TTAG) + log(K_TTAT) + log(K_TTCA) + log(K_TTCC) + log(K_TTCG) + log(K_TTCT) + log(K_TTGA) + log(K_TTGC) + log(K_TTGG) + log(K_TTGT) + log(K_TTTA) + log(K_TTTC) + log(K_TTTG) + log(K_TTTT)
    else
        error("1 <= K <= 4")
    end
end

