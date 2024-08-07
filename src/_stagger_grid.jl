#========================================================== Type definitions =#

"""
Information of a part of the Stagger grid. Contains location of
snapshots and averages.
"""
struct StaggerGrid{DF<:DataFrame} <:AbstractMUSTGrid
    name :: String
    info :: DF
end




#============================================================== Constructors =#

"""
    StaggerGrid(path)

Load a already computed average Stagger grid. 
The save location and layout of the Stagger grid may
vary from case to case, so only the reading function is 
provided here. For an example case please go to the 
examples folder and look at the `initial_models/from_stagger.ipynb` notebook.

Methods for the interaction with the `TSO.jl` module are located in the same 
folder for now to keep the depencies separate. At a later stage, `TSO.jl` might
be added as a dependecy here.
"""
StaggerGrid(path) = begin
    info = CSV.read(path, DataFrame)
    name = basename(path)[1:findlast('.', basename(path))-1]

    StaggerGrid(name, info)
end




#==================================================================== Saving =#

save(grid::StaggerGrid, path) = CSV.write(path, grid.info)




#=============================================================== Convenience =#

Base.getindex(g::StaggerGrid, k::String, i=!) = g.info[i, k]
Base.:+(g1::StaggerGrid, g2::StaggerGrid) = StaggerGrid(g1.name, vcat(g1.info, g2.info))




#============================================================= Interpolation =#
 
"""
	interpolate_quantity(grid::StaggerGrid, what; teff, logg, feh, method="linear")

Interpolate quantity `what` from within the grid to the new teff, logg and feh values.
Uses scipy griddata for scattered interpolation. Convert output to julia.

# Examples

```julia
vmin = interpolate_quantity(grid, "vmin"; teff=teff, logg=logg, feh=feh)
```
"""
function interpolate_quantity(grid::StaggerGrid, what; teff, logg, feh, method="linear")
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

allowed_namelists(grid::StaggerGrid) = grid.info[!,"namelist_name"]

stage_namelists(grid::StaggerGrid, folder="run_grid"; clean_logs=true, clean_namelists=false) = begin
	nml_names = basename.(grid["namelist_name"])
	nml_path = grid["namelist_name"]
	folder_run = @in_dispatch(folder)

	if !isdir(folder_run)
		mkdir(folder_run)
	end

	# clean up dispatch from logs and namelists
	glob("*", folder_run) .|> rm

	for (i, n) in enumerate(nml_names)
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