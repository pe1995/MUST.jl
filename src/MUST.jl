module MUST

#= Julia modules =#
using PyCall
using DataFrames
using FortranFiles: FortranFile, read, readlines
using Printf
import Base.filter, Base.filter!
import Base.getindex

#= Python modules =#
const numpy             = PyNULL()
const scipy_interpolate = PyNULL()

__init__() = begin 
    copy!(scipy_interpolate,pyimport("scipy.interpolate"))
    copy!(numpy,pyimport("numpy"))
end

#= MUST interface =#
export import_dispatch, in_dispatch
export Space, spacebox, Box

#= Julia code files =#
include("_dispatch.jl")
include("_atmos.jl")

end
