#========================================================== Type definitions =#

"""
Information of a part of the Stagger grid. Contains location of
snapshots and averages.
"""
struct Atmos1DGrid{DF<:DataFrame} <:AbstractMUSTGrid
    name ::String
	path ::String
    info ::DF
end




#============================================================== Constructors =#

"""
	Atmos1DGrid(path)

Load a already computed average grid. 
The save location and layout of the Stagger grid may
vary from case to case, so only the reading function is 
provided here. For an example case please go to the 
examples folder and look at the `initial_models/from_stagger.ipynb` notebook.

Methods for the interaction with the `TSO.jl` module are located in the same 
folder for now to keep the depencies separate. At a later stage, `TSO.jl` might
be added as a dependecy here.

__NOTE__: The path of the grid is essential, because paths to the average models 
within are assumed to be relative to that path. This is needed, because the 
same grid file should be usable whereever the code is used, however it should still
be flexible enough to put your own favorite grid whereever you want and execute the code
from somewhere else. 
"""
Atmos1DGrid(path) = begin
    info = CSV.read(path, DataFrame)
    name = basename(path)[1:findlast('.', basename(path))-1]
	path = abspath(path)

    Atmos1DGrid(name, path, info)
end

StaggerGrid(args...; kwargs...) = Atmos1DGrid(args...; kwargs...)


#==================================================================== Saving =#

save(grid::Atmos1DGrid, path) = CSV.write(path, grid.info)
save(grid::Atmos1DGrid) = CSV.write(grid.path, grid.info)



#=============================================================== Convenience =#

Base.getindex(g::Atmos1DGrid, k::String, i=!) = g.info[i, k]
Base.:+(g1::Atmos1DGrid, g2::Atmos1DGrid) = Atmos1DGrid(g1.name, g1.path, vcat(g1.info, g2.info))
Base.dirname(g::Atmos1DGrid) = dirname(g.path)




#===================================================================== Paths =#

"""
	absolute_path!(grid; from=@__FILE__)

Convert grid paths from relative to absolute. Assume they are relative to `from`
"""
absolute_path!(grid; from=@__FILE__) = begin
    if "av_path" in names(grid.info) 
		grid.info[!, "av_path"] .= abspath.(joinpath.(dirname(from), grid.info[!, "av_path"]))
	end
	if "avo_path" in names(grid.info) 
    	grid.info[!, "avo_path"] .= abspath.(joinpath.(dirname(from), grid.info[!, "avo_path"]))
    end
	if "eos_root" in names(grid.info) 
		grid.info[!, "eos_root"] .= abspath.(joinpath.(dirname(from), grid.info[!, "eos_root"]))
	end
	if "matching_eos" in names(grid.info) 
		grid.info[!, "matching_eos"] .= abspath.(joinpath.(dirname(from), grid.info[!, "matching_eos"]))
	end
	if "matching_opa" in names(grid.info) 
		grid.info[!, "matching_opa"] .= abspath.(joinpath.(dirname(from), grid.info[!, "matching_opa"]))
	end
	if "matching_sopa" in names(grid.info) 
		grid.info[!, "matching_sopa"] .= abspath.(joinpath.(dirname(from), grid.info[!, "matching_sopa"]))
	end
end

relative_path(from, to) = begin
	fstag = basename(from)
	feos = basename(to)
	stag_path = split(dirname(from), "/", keepempty=false)
	eos_path = split(dirname(to), "/", keepempty=false)

	i = findfirst(x->!(x in eos_path), stag_path)
	ieos = findfirst(eos_path .== stag_path[i-1]) + 1

	path_difference = stag_path[i:end]
	eos_path_difference = eos_path[ieos:end]
	
	joinpath([".." for _ in path_difference]..., eos_path_difference..., feos)
end

"""
	relative_path!(grid; from=@__FILE__)

Convert grid paths from absolute to relative to `from`.
"""
relative_path!(grid; to=@__FILE__) = begin
    if "av_path" in names(grid.info) 
		grid.info[!, "av_path"] .= relative_path.(to, grid.info[!, "av_path"])
	end
	if "avo_path" in names(grid.info) 
    	grid.info[!, "avo_path"] .= relative_path.(to, grid.info[!, "avo_path"])
    end
	if "eos_root" in names(grid.info) 
		grid.info[!, "eos_root"] .= abspath.(joinpath.(dirname(from), grid.info[!, "eos_root"]))
	end
	if "matching_eos" in names(grid.info) 
		grid.info[!, "matching_eos"] .= relative_path.(to, grid.info[!, "matching_eos"])
	end
	if "matching_opa" in names(grid.info) 
		grid.info[!, "matching_opa"] .= relative_path.(to, grid.info[!, "matching_opa"])
	end
	if "matching_sopa" in names(grid.info) 
		grid.info[!, "matching_sopa"] .= relative_path.(to, grid.info[!, "matching_sopa"])
	end
end


#============================================================== Replace Path =#

"""
	replace_relative_path!(grid; kwargs...)

Replace e.g. the EoS of a grid that contains relative paths with a new array 
of paths (that are absolute!).
"""
function replace_relative_path!(grid; kwargs...)
	# make paths absolute
	absolute_path!(grid)

	# replace paths with new ones
	replace_absolute_path!(grid; kwargs...)

	# make paths relative again
	relative_path!(grid)

	grid
end

"""
	replace_absolute_path!(grid; kwargs...)

Replace e.g. the EoS of a grid that contains absolute paths with a new array 
of paths (that are absolute!).
"""
function replace_absolute_path!(grid; kwargs...)
	# paths are already absolute
	for (k, v) in kwargs
		grid.info[k] = v
	end
	
	grid
end




#============================================================= Interpolation =#
 
"""
	interpolate_quantity(grid::Atmos1DGrid, what; teff, logg, feh, method="linear")

Interpolate quantity `what` from within the grid to the new teff, logg and feh values.
Uses scipy griddata for scattered interpolation. Convert output to julia.

# Examples

```julia
vmin = interpolate_quantity(grid, "vmin"; teff=teff, logg=logg, feh=feh)
```
"""
function interpolate_quantity(grid::Atmos1DGrid, what; teff, logg, feh, method="linear")
	feh_gr  = grid.info[!, "feh"]
	femask = feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]

	logg_gr = grid.info[femask, "logg"]
	teff_gr = grid.info[femask, "teff"]
	what_gr = grid.info[femask, what]

	a = pyconvert(
		Any,
		first(
			scipy_interpolate.griddata(
				(logg_gr, teff_gr), what_gr, ([logg], [teff]), 
				method=method
			)
		), 
	)

	if isnan(a)
		pyconvert(
			Any,
			first(
				scipy_interpolate.griddata(
					(logg_gr, teff_gr), what_gr, ([logg], [teff]), 
					method="nearest"
				)
			), 
		)
	else
		a
	end
end




#=================================================================== running =#

allowed_namelists(grid::Atmos1DGrid) = grid.info[!,"namelist_name"]

stage_namelists(grid::Atmos1DGrid, folder="run_grid"; clean_logs=true, clean_namelists=false) = begin
	nml_names = basename.(grid["namelist_name"])
	nml_path = grid["namelist_name"]
	folder_run = @in_dispatch(folder)

	if !isdir(folder_run)
		mkdir(folder_run)
	end

	# clean up dispatch from logs and namelists
	glob("*", folder_run) .|> rm

	for (i, n) in enumerate(nml_names)
		@info "copy $(nml_path[i]) to $(joinpath(folder_run, n))."
		cp(nml_path[i], joinpath(folder_run, n))
	end

	if clean_namelists
		glob("grid_*.nml", @in_dispatch("")) .|> rm
		for (i, n) in enumerate(nml_names)
			cp(joinpath(folder_run, n), joinpath(@in_dispatch(""), n))
		end
	end
	if clean_logs
		glob("grid_*.log", @in_dispatch("")) .|> rm
		glob("grid_*.err", @in_dispatch("")) .|> rm
	end
end


#=============================================================================#