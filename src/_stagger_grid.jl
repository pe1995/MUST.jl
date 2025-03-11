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

"""
	ModelInformation

Information from model headers in grid naming convention.
"""
Base.@kwdef mutable struct ModelInformation
	category = ""
	teff = ""
	logg = ""
	feh = ""
	extension = ""
	wrong_format = false
	original_name = ""
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

Convert grid paths from relative to absolute. Assume they are relative to `from`.
First check if the path is not already absolute.
"""
absolute_path!(grid; from=@__FILE__) = begin
    if "av_path" in names(grid.info) 
		grid.info[!, "av_path"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "av_path"]]
	end
	if "avo_path" in names(grid.info) 
    	grid.info[!, "avo_path"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "avo_path"]]
    end
	if "eos_root" in names(grid.info) 
		grid.info[!, "eos_root"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "eos_root"]]
	end
	if "matching_eos" in names(grid.info) 
		grid.info[!, "matching_eos"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_eos"]]
	end
	if "matching_opa" in names(grid.info) 
		grid.info[!, "matching_opa"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_opa"]]
	end
	if "matching_sopa" in names(grid.info) 
		grid.info[!, "matching_sopa"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_sopa"]]
	end
end

relative_path(relative_to, abs_path) = begin
	filename = basename(abs_path)

	relative_to_components = split(dirname(relative_to), "/", keepempty=false)
	abs_path_components = split(dirname(abs_path), "/", keepempty=false)

	n_components = length(abs_path_components)
	n_new_components = length(relative_to_components)
	i = findfirst(x->!(x in relative_to_components), abs_path_components)
	i = isnothing(i) ? n_components + 1 : i
	different_components = i>n_components ? [] : abs_path_components[i:end]
	from_node_to_new_path_components = i>n_new_components ? [] : relative_to_components[i:end]

	joinpath([".." for i in from_node_to_new_path_components]..., different_components..., filename)
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
		grid.info[!, "eos_root"] .= relative_path.(to, grid.info[!, "eos_root"])
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


#====================================== Information of grid models from Name =#

"""
	_find_datastring(parts)

Find the datastring within the name of a grid model in naming convention.
"""
_find_datastring(parts) = begin
	i_part = Ref(-1)
	for (i, x) in enumerate(parts)
		if occursin("t", x)&occursin("g", x)&occursin("m", x)
			p = try 
				_parameters_from_string(x)
			catch
				nothing
			end
			if isnothing(p)
				continue
			else
				i_part[] = i
				break
			end
		else
			continue
		end
	end

	i_part[]
end

"""
	_parameters_from_string(para)

Extract stellar parameters from name of the model in naming convention.
"""
_parameters_from_string(para) = begin
	t1 = findfirst(x->x=='t', para) + 1
	t2 = findfirst(x->x=='g', para) - 1
	t3 = findfirst(x->x=='m', para) - 1

	t = parse(Float64, para[t1:t2]) .* 100
	g = parse(Float64, para[t2+2:t3]) ./10
	m = parse(Float64, para[t3+2:end])

	d = Dict("teff"=>t, "logg"=>g, "feh"=>m)
end





#= Comparison of model information =#

"""
	same_parameters(a, b)

Check if a and b have the same stellar parameters.
"""
same_parameters(a::ModelInformation, b::ModelInformation) = (a.teff≈b.teff) & (a.logg≈b.logg) & (a.feh≈b.feh)

"""
	same_id(a, b, field)

Check if a and b have the same `field` entry. `False` if entries are empty.
"""
same_id(a::ModelInformation, b::ModelInformation, field) = if (length(getfield(a, field)) > 0) & (length(getfield(b, field)) > 0)
	getfield(a, field) == getfield(b, field)
else
	false
end

"""
	same_id(a, b, field)

Check if a and b have the same `category` and `extension` entries.
"""
same_class(a::ModelInformation, b::ModelInformation) = same_id(a, b, :category) & same_id(a, b, :extension)





#= building new name string from information =#

"""
	name_string(minfo::ModelInformation; add_E=false)

Build a new name string from a model information.
"""
name_string(minfo::ModelInformation; add_E=false) = begin
	tname = MUST.@sprintf "%.2f" minfo.teff/100
	gname = MUST.@sprintf "%.2f" minfo.logg*10
	fname = MUST.@sprintf "%.3f" minfo.feh
	name = "t$(tname)g$(gname)m$(fname)"
	
	components = if add_E
		[minfo.category, "E", name, minfo.extension]
	else
		[minfo.category, name, minfo.extension]
	end

	components = [c for c in components if length(c)>0]
	join(components, '_')
end




"""
	ModelInformation(name)

Extract category, teff, logg, feh and version from standard m3dis format names.
The format is `category_teff/1000_logg*10_±feh_extension`.
"""
ModelInformation(name) = begin
	name_parts = split(name |> strip, "_")

	# if there is an E in the parameter name we remove it
	name_parts = [p for p in name_parts if !(p=="E")]

	try
		# find the data string
		idatastring = _find_datastring(name_parts)

		if idatastring == -1
			# this has the wrong format
			ModelInformation(wrong_format=true, original_name=name)
		else
			# everyting before the datastring is prefix, the rest is extension
			prefix = try
				join(name_parts[1:idatastring-1], '_')
			catch
				""
			end

			extension = try
				join(name_parts[idatastring+1:end], '_')
			catch
				""
			end
	
			p = _parameters_from_string(name_parts[idatastring])
			ModelInformation(
				teff=p["teff"],
				logg=p["logg"],
				feh=p["feh"],
				category=prefix,
				extension=extension,
				original_name=name
			)
		end
	catch
		ModelInformation(wrong_format=true, original_name=name)
	end
end

#=============================================================================#