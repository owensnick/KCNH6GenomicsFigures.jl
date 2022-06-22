using DataFrames, DataFramesMeta, CSV, Glob
using CodecZlib, ProgressMeter
using Statistics, StatsBase, StatsModels, GLM, MultivariateStats
using BioSequences, FASTX
using Plots, StatsPlots, Measures
using GaussianProcesses, Distributions, Suppressor
using Clustering, ClusterOrderTools, Distances

using Enrichr
using MotifScanner

theme(:wong2)

include("loaddata.jl")
include("spikecorrect.jl")
include("plots.jl")
include("gaussianprocesses.jl")
include("clustering.jl")
include("enrichr.jl")
include("homer.jl")

function showwide(table, nc=10000)
    c = ENV["COLUMNS"]
    ENV["COLUMNS"] = string(nc)
    display(table)
    ENV["COLUMNS"] = c;
    nothing;
end

function getprojectdir()
    d = pwd()
    if basename(d) == "notebooks"
        return dirname(d)
    else
        return d
    end
end
