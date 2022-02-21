using DataFrames, DataFramesMeta, CSV, Glob
using CodecZlib, ProgressMeter
using Statistics, StatsModels, GLM
using BioSequences, FASTX


include("loaddata.jl")
include("spikecorrect.jl")

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