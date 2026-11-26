module MUST

# ============================================================================
# Julia modules
# ============================================================================
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
using NCDatasets
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

# ============================================================================
# Python modules
# ============================================================================
const numpy = PythonCall.pynew()
const f90nml = PythonCall.pynew()
const scipy_interpolate = PythonCall.pynew()

__init__() = begin
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
    PythonCall.pycopy!(numpy, pyimport("numpy"))
    PythonCall.pycopy!(f90nml, pyimport("f90nml"))
end

# ============================================================================
# Abstract types
# ============================================================================
abstract type AbstractSpace end
abstract type AbstractInitialModel end

abstract type AbstractMUSTGrid end

abstract type AbstractBoxGrid
end

abstract type AbstractBoxAxis
end

# ============================================================================
# MUST interface
# ============================================================================
export @import_dispatch, @in_dispatch, @import_m3dis, @in_m3dis, @import_tumult, @in_tumult, @ingredient
export Space, spacebox, Box, add!, profile, time_average_profile, flip!
export multiBox, snapshotBox
export pick_snapshot, list_of_snapshots, converted_snapshots, list_snapshots, find_snapshots, BoxRef, availableData, switch, switch!, migrateBoxes
export ginterpolate, gevaluate, gevaluate!, gresample, Grid
export plane_statistic
export height_scale

export save_must1d, save_must1d_all, read_must1d, must1d_columns, averageBox

export synthesize, spectrum, composition_settings, Teff, flux, M3DISRun, window, TUMULTRun

export marcsBox

export mafagsBox

export axis, closest
export ingredients

export defaultWatchDog, monitor

export Spectrum, spectra, load_spectra, snapshot_spectra, average_spectra, μ_z, intensity, composition_string, curve_of_growth

# ============================================================================
# Julia code files
# ============================================================================
include("_constants.jl")
include("_heightscales.jl")
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
include("_atmos.jl")
include("_grid_construction.jl")
include("_mustgrid.jl")
include("_model_information.jl")
include("_stagger_grid.jl")
include("_initial_conditions.jl")
include("_running.jl")
include("_multi.jl")
include("_spectra.jl")
include("_atmos2legacy.jl")
include("_atmos2multi.jl")
include("_must1d.jl")
include("_read_marcs.jl")
include("_mafags.jl")
include("_convert.jl")
include("_restart.jl")

include("_muram.jl")

include("_physical_quantities.jl")

include("_watchdog.jl")

const timers = [
    boxingTime,
    opticalDepthTime,
    heightScaleTime,
    multiTime
]

const detailedTimers = [
    timers...,
    detailedBoxingTimers...
]

end
