module MUST

#= Julia modules =#
using CondaPkg
using PythonCall
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
using ArgParse
using DelimitedFiles
using TimerOutputs
using NetCDF
using Integrals
using Mmap
using Downloads
import Base.filter, Base.filter!
import Base.getindex
import Base.getproperty
import Base.length
import Base.keys
import Base.read!, Base.write
import Base.size, Base.axes
import Base.Broadcast.broadcastable
import Base.:+
import Base.dirname


#= Python modules =#
const numpy = PythonCall.pynew()
const f90nml = PythonCall.pynew()
const scipy_interpolate = PythonCall.pynew()

__init__() = begin 
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
    PythonCall.pycopy!(numpy, pyimport("numpy"))
    PythonCall.pycopy!(f90nml, pyimport("f90nml"))
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
export @import_dispatch, @in_dispatch, @import_m3dis, @in_m3dis
export Space, spacebox, Box, add!, profile, time_average_profile, flip!
export multiBox, snapshotBox
export pick_snapshot, list_of_snapshots, converted_snapshots, list_snapshots
export ginterpolate, gevaluate, gevaluate!, gresample, Grid
export plane_statistic 

# Multi
export whole_spectrum, spectrum, Teff, flux, M3DISRun, window

# Marcs
export marcsBox

export axis, closest
export ingredients

# monitoring
export defaultWatchDog, monitor




#= Julia code files =#
include("_constants.jl")
include("_argparse.jl")
include("_grids.jl")
include("_grid_interpolation.jl")
include("_eos.jl")
include("_parallel.jl")
include("_help.jl")
include("_dispatch.jl")
include("_scaling.jl")
include("_stagger.jl")
include("_namelist.jl")
include("_spectra.jl")
include("_atmos.jl")
include("_grid_construction.jl")
include("_mustgrid.jl")
include("_stagger_grid.jl")
include("_initial_conditions.jl")
include("_running.jl")
include("_multi.jl")
include("_atmos2legacy.jl")
include("_atmos2multi.jl")
include("_read_marcs.jl")
include("_convert.jl")

# co5bold reader
#include("_co5bold.jl")

# physical utilities
include("_physical_quantities.jl")

# watchdog
include("_watchdog.jl")

# timers
const timers = [
    boxingTime,              # time the collection of patches
    opticalDepthTime,        # time the optical depth computation
    heightScaleTime,         # time the optical depth scale interpolation
    multiTime                # time the creation of M3D cubes
]

const detailedTimers = [
    timers...,               # General timers
    detailedBoxingTimers...  # more detail for converting patches to Boxes
]

end
