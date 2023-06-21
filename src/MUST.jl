module MUST

#= Julia modules =#
using PyCall
using DataFrames
using CSV
using FortranFiles
using Printf
using Statistics
using HDF5
using Random
using Interpolations
using Glob
using Distributed
#using NumericalIntegration
using ArgParse
using DelimitedFiles
using TimerOutputs
using NetCDF
import Base.filter, Base.filter!
import Base.getindex
import Base.length
import Base.keys
import Base.read!, Base.write
import Base.size, Base.axes
import Base.Broadcast.broadcastable


#= Python modules =#
const numpy             = PyNULL()
const scipy_interpolate = PyNULL()

__init__() = begin 
    copy!(scipy_interpolate,pyimport("scipy.interpolate"))
    copy!(numpy,pyimport("numpy"))
end


#= Abstract types =#
abstract type AbstractSpace end
abstract type AbstractInitialModel end

"""
Abstract MUSTGrid. Depending on the type of grid different methods are available.
    All concrete types are required to have an info field:
    name :: String
    info :: DataFrame
    So that a namelist can be constructed from the grid.
"""
abstract type AbstractMUSTGrid end

# Gridding
"""
Abstract Grid of a MUST.Box object.
"""
abstract type AbstractBoxGrid
end

"""
Abstract Axis of a MUST.Box object.
"""
abstract type AbstractBoxAxis
end



#= MUST interface =#
export import_dispatch, in_dispatch
export Space, spacebox, Box, add!
export pick_snapshot
export ginterpolate, gevaluate, gevaluate!
#export read!, write


#= Julia code files =#
include("_argparse.jl")
include("_grids.jl")
include("_grid_interpolation.jl")
include("_eos.jl")
include("_parallel.jl")
include("_help.jl")
include("_scaling.jl")
include("_stagger.jl")
include("_dispatch.jl")
include("_namelist.jl")
include("_atmos.jl")
include("_grid_construction.jl")
include("_mustgrid.jl")
include("_stagger_grid.jl")
include("_initial_conditions.jl")
include("_running.jl")
include("_atmos2legacy.jl")
include("_atmos2multi.jl")

end
