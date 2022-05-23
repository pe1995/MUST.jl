module MUST

#= Julia modules =#
using PyCall
using DataFrames
using FortranFiles: FortranFile, read, readlines
using Printf
using Statistics
using HDF5
using Random
using Interpolations
using Glob
using Distributed
using ArgParse
using Interpolations
using DelimitedFiles
import Base.filter, Base.filter!
import Base.getindex
import Base.length
import Base.keys
import Base.read!, Base.write

#= Python modules =#
const numpy             = PyNULL()
const scipy_interpolate = PyNULL()

__init__() = begin 
    copy!(scipy_interpolate,pyimport("scipy.interpolate"))
    copy!(numpy,pyimport("numpy"))
end

#= MUST interface =#
export import_dispatch, in_dispatch
export Space, spacebox, Box, add!
#export read!, write

#= Julia code files =#
include("_argparse.jl")
include("_parallel.jl")
include("_help.jl")
include("_scaling.jl")
include("_stagger.jl")
include("_dispatch.jl")
include("_namelist.jl")
include("_atmos.jl")
include("_mustgrid.jl")
include("_running.jl")

end
