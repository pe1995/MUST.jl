#= Running dispatch on a Grid =#

#========== Type Definitions ==========#

mutable struct RangeMUSTGrid{DF<:DataFrame,F<:AbstractFloat} <: AbstractMUSTGrid
    seeds  :: Dict{Symbol, Dict{String, F}}
    limits :: Dict{Symbol, Dict{String, F}}
    ngrid  :: Int
    grid   :: Dict{Symbol, Dict{String, Vector{F}}}
    name   :: String
    info   :: DF
end

mutable struct IdentityMUSTGrid{DF<:DataFrame} <: AbstractMUSTGrid
    name :: String
    info :: DF
end

"""Restart Grid. Shall be constructed from previous grid."""
mutable struct RestartMUSTGrid{DF<:DataFrame} <: AbstractMUSTGrid
    grid        :: Dict{Symbol, Dict{String, Vector{Any}}}
    from_name   :: String
    name        :: String
    info        :: DF
end



#========== Constructors ==========#

RangeMUSTGrid(seeds::T, limits::T; phase="phase1") where {F<:AbstractFloat, T<:Dict{Symbol, Dict{String, F}}} = begin
    RangeMUSTGrid(seeds, limits,  0, Dict{Symbol, Dict{String, Vector{F}}}(), phase, DataFrame())
end

RangeMUSTGrid(info::DataFrame; t::Type=Float64, phase="phase1") = RangeMUSTGrid(Dict{Symbol, Dict{String, t}}(), Dict{Symbol, Dict{String, t}}(), 0, Dict{Symbol, Dict{String, Vector{t}}}(), phase, info)

IdentityMUSTGrid(; phase="phase_identity") = IdentityMUSTGrid(phase, DataFrame())

RestartMUSTGrid(from::AbstractMUSTGrid; grid::T=T(), phase="phase2") where {T<:Dict{Symbol, Dict{String, Vector{<:Any}}}} = begin
    from_phase = from.name
    grid_in = Dict{Symbol, Dict{String, Vector{Any}}}()
    grid    = deepcopy(grid)

    for k1 in keys(grid)
        grid_in[k1] = Dict{String, Vector{Any}}()
        for k2 in keys(grid[k1])
            grid_in[k1][k2] = grid[k1][k2]
        end
    end

    !(:restart_params in keys(grid_in)) ? grid_in[:restart_params] = Dict{String, Vector{Any}}() : nothing
    grid_in[:restart_params]["run"] = ["'"*String(split(n, ".nml")[1])*"'" for n in from.info[!,"$(from_phase)_name"]]

    RestartMUSTGrid(grid_in, from_phase, phase, deepcopy(from.info))
end



#========== Functionality of the types ==========#

string_from_keys(inkeys...; dlim="!") = join(String.(inkeys), dlim)
keys_from_string(instr ; dlim="!")    = split(instr, dlim)

randomgrid!(grid::RangeMUSTGrid{DF, F}, ngrid::Int) where {DF, F} = begin
    grid.ngrid = ngrid
    summary    = Dict{String, Vector{Any}}()

    for field in keys(grid.seeds)
        grid.grid[field] = Dict{String, Vector{F}}()

        for p in keys(grid.seeds[field])
            lims = [ grid.seeds[field][p] - abs(grid.seeds[field][p]*grid.limits[field][p]), 
                     grid.seeds[field][p] + abs(grid.seeds[field][p]*grid.limits[field][p]) ]

            grid.grid[field][p] = MUST.randrange(lims..., F, ngrid) 
            
            summary[string_from_keys(field, p)] = grid.grid[field][p]
        end
    end
    
    grid.info  = DataFrame(summary)
    nothing
end

sync_nml_grid!(nml::AbstractNamelist, grid::IdentityMUSTGrid, args...; kwargs...) = nothing
sync_nml_grid!(nml::AbstractNamelist, grid::AbstractMUSTGrid, idx::Int) = begin
    for field in keys(grid.grid)
        data = getfield(nml, field)
        for parameter in keys(grid.grid[field])
            data[lowercase(parameter)] = grid.grid[field][parameter][idx]
        end
    end
    nothing
end

modify!(f::Function, grid::AbstractMUSTGrid, field::Symbol, parameter::String, args...; kwargs...) = begin
    grid.grid[field][parameter] = f.(grid.grid[field][parameter], args...; kwargs...)
    grid.info[!,string_from_keys(field, parameter)] = grid.grid[field][parameter]
    nothing
end




#========== Grid <-> Namelist creation ==========#

"""
Create ngrid namelists from existing namelist.
"""
function create_namelists!(grid::AbstractMUSTGrid; 
                            default_nml::Union{Nothing, StellarNamelist}=nothing, start_index::Int=1)
    ngrid = nrow(grid.info)
    phase = grid.name

    dummy_nml = isnothing(default_nml) ? 
                MUST.StellarNamelist(MUST.@in_dispatch("stellar.nml")) : default_nml
    
    summary = Dict("$(phase)_name"=>[])

    for i in 1:ngrid
        # name of namelists
        id = "grid$(start_index-1+i)"
        name = "$(id)_$(phase).nml"

        # copy the dummy namelist
        nml = deepcopy(dummy_nml)

        # Set the values
        sync_nml_grid!(nml, grid, i)
        
        # write the namelists
        MUST.write(nml, MUST.@in_dispatch(name))

        # save the summary
        append!(summary["$(phase)_name"], [name])
    end

    for key in keys(summary)
        grid.info[!,key] = summary[key]
    end

    grid.name = phase
end
