using Pkg; Pkg.activate("."); 
using PyPlot
using LaTeXStrings
using NetCDF
using MUST
using StatsBase
using TSO
using DelimitedFiles

model_path = ARGS[1]
model_name = split(last(split(model_path, "/")),".nc") |> first
extension  = length(ARGS) > 1 ? "_$(ARGS[2])" : ""

model_name = "$(model_name)$(extension)"

@show model_name

U  = ncread(model_path, "U")
W  = ncread(model_path, "W")
T  = ncread(model_path, "T")
P  = ncread(model_path, "P")
By = ncread(model_path, "By")
V  = ncread(model_path, "V")
Bx = ncread(model_path, "Bx")
R  = ncread(model_path, "R")
E  = ncread(model_path, "E")
Bz = ncread(model_path, "Bz")

s = size(U)

dx = ncgetatt(model_path, "Global", "dx")
dy = ncgetatt(model_path, "Global", "dy")
dz = ncgetatt(model_path, "Global", "dz")

time  = ncgetatt(model_path, "Global", "time")
paras = MUST.AtmosphericParameters() 
paras.time = time

x =  range(-(s[1]-1) * dx/2,   step=dx, length=s[1])
y =  range(-(s[2]-1) * dy/2,   step=dy, length=s[2])
z =  range(-(s[3]-1) * dz*1/2, step=dz, length=s[3])

xx, yy, zz = MUST.meshgrid(x, y, z)

data = Dict{Symbol, Array{eltype(U), 3}}(   :ux  => U ,
                                            :uz  => W ,
                                            :T   => T ,
                                            :P   => P ,
                                            :By  => By,
                                            :uy  => V ,
                                            :Bx  => Bx,
                                            :d   => R ,
                                            :E   => E ,
                                            :Bz  => Bz )

b = MUST.Box(xx, yy, zz, data, paras)


#======================================= optical depth ===#

eos_path = "/u/peitner/Turbospectrum/opacity_tables/tests/DIS_MARCS_v0.4.4/"
eos      = reload(SqEoS, joinpath(eos_path, "eos.hdf5"))

b.data[:kr] = exp.(lookup(eos, :lnRoss, log.(b[:d]), log.(b[:T])))

τ = MUST.optical_depth(b, opacity=:kr, density=:d)
b.data[:τ_ross] = τ


#================================================ Save ===#

MUST.save(b, name="box_MURaM_$(model_name)", folder="")
