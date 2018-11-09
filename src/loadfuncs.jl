using Distributed

@everywhere using SharedArrays
@everywhere using LinearAlgebra
@everywhere using Distributions
@everywhere using RCall
@everywhere using LightGraphs

@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/nichemodelweb.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/cascade.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/plotweb.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/smallwebs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/src/trophic.jl")
