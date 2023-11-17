using MUST 
using TSO



#==================================================================== Models =#

function initial_model(grid::MUST.AbstractMUSTGrid, point::Int; kwargs...)
	try
		eos, opa = equation_of_state(grid, point; kwargs...)
		@optical Average3D(eos, grid.info[point, "av_path"]; logg=grid.info[point, "logg"]) eos opa
	catch
		Average3D(grid.info[point, "av_path"]; logg=grid.info[point, "logg"]) 
	end
end

function initial_model(grid::MUST.AbstractMUSTGrid, point::String; kwargs...)
	i = findfirst(grid.info[!, "name"] .== point)
	initial_model(grid, i; kwargs...)
end



#======================================================================= EoS =#

function equation_of_state(grid::MUST.AbstractMUSTGrid, point::Int; energy=false)
	f = if energy
		grid.info[1, "binned_E_tables"]
	else
		grid.info[1, "binned_tables"]
	end
	
	eos = reload(SqEoS, 
		joinpath(f, "eos.hdf5"))
	opa = reload(SqOpacity, 
		joinpath(f, "binned_opacities.hdf5"))
	eos, opa
end

function equation_of_state(grid::MUST.AbstractMUSTGrid, point::String; kwargs...)
	i = findfirst(grid.info[!, "name"] .== point)
	equation_of_state(grid, i; kwargs...)
end

find_closest_eos(grid::MUST.AbstractMUSTGrid, teff, feh) = begin
	mother_tables = grid["eos_root"]
	path = if all(mother_tables .== mother_tables[1])
		mother_tables[1]
	else
		# check if there are models with the same metallicity
		fehall = grid["feh"]
		tempall = grid["teff"]

		matchingtables = feh .== fehall
		if all(.!matchingtables)
			# closest metallicity
			imin = argmin(abs.(fehall .- feh))
			fehmin = fehall[imin]

			# all tables with this metallicity
			fehminall = fehmin .≈ fehall

			# pick the one where the temp is closest
			tempminall = tempall[fehminall]
			mothertablesminall = mother_tables[fehminall]
			imin = argmin(abs.(tempminall .- teff))

			mothertablesminall[imin]
		else
			# pick the one where the temp is closest
			tempminall = tempall[matchingtables]
			mothertablesminall = mother_tables[fehminall]
			imin = argmin(abs.(tempminall .- teff))

			mothertablesminall[imin]
		end
	end

	reload(
			SqEoS,
			joinpath(path, "combined_eos.hdf5")
		)
end



#================================================ interpolate average models =#

function interpolate_average(grid::MUST.AbstractMUSTGrid; teff, logg, feh, common_size=1000, kwargs...)
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	
	# read all average models and interpolate them to the same number of points
	models = [initial_model(grid, String(name))
					for name in grid.info[!, "name"]]
	r_models = [minimum(abs.(diff(m.z))) for m in models]
	
	models = [TSO.upsample(m, common_size) for m in models]

	# now we interpolate all points to one common point, for every point
	scatter_int(v, x, y, z) = MUST.pyconvert(typeof(x),
		first(MUST.scipy_interpolate.griddata((logg_gr, teff_gr, feh_gr), v, ([x], [y], [z]), method="linear"))
	)
	
	points = zeros(eltype(models[1].z), length(models), 3)
	z = zeros(eltype(models[1].z), common_size)
	t = zeros(eltype(models[1].z), common_size)
	d = zeros(eltype(models[1].z), common_size)
	
	for i in 1:common_size
		for j in eachindex(models)
			points[j, 1] = models[j].z[i]
			points[j, 2] = models[j].lnT[i]
			points[j, 3] = models[j].lnρ[i]
		end
		z[i] = scatter_int(points[:, 1], logg, teff, feh)
		t[i] = scatter_int(points[:, 2], logg, teff, feh)
		d[i] = scatter_int(points[:, 3], logg, teff, feh)
	end

	r_target = scatter_int(r_models, logg, teff, feh)
	npoints  = ceil(Int, abs(maximum(z) - minimum(z)) / r_target)

	TSO.upsample(Model1D(z=z, lnT=t, lnρ=d, logg=logg), npoints)
end



#================================================================== Adiabats =#

function adiabat(mother_table, av_path, logg, common_size; save=false)
	eos = reload(
		SqEoS,
		joinpath(mother_table, "combined_eos.hdf5")
	)

	data = flip(Average3D(av_path, logg=logg))
	start_point = TSO.pick_point(data, 1)
	end_point = TSO.pick_point(data, length(data.z))
	
	a = TSO.upsample(TSO.adiabat(start_point, end_point, eos; kwargs...), common_size)
	a = TSO.flip(a, depth=true)
	if save
		open(av_path, "w") do f
			MUST.writedlm(f, [a.z, a.lnT, a.lnρ])
		end
	end
end


#===========================================================interpolate grid =#

function interpolate_from_grid(grid::MUST.AbstractMUSTGrid, teff::F, logg::F, feh::F) where {F<:AbstractFloat}
	# create the iniitial model from interpolating the average snapshots
	model = interpolate_average(grid, teff=teff, logg=logg, feh=feh)

	folder = "interpolated"
	mesh   = "interpolated"

	tname = MUST.@sprintf "%i" teff/100
	gname = MUST.@sprintf "%.2f" logg*10
	fname = MUST.@sprintf "%.3f" feh
	name = "t$(tname)g$(gname)m$(fname)"

	snapshot = "interpolated"
	ma_x = MUST.interpolate_quantity(grid, "ma_x"; teff=teff, logg=logg, feh=feh)
	ma_y = MUST.interpolate_quantity(grid, "ma_y"; teff=teff, logg=logg, feh=feh)
	mi_z = minimum(model.z)
	ma_z = maximum(model.z)
	vmin = MUST.interpolate_quantity(grid, "vmin"; teff=teff, logg=logg, feh=feh)
	vmax = MUST.interpolate_quantity(grid, "vmax"; teff=teff, logg=logg, feh=feh)
	
	tscale = round(
		MUST.interpolate_quantity(grid, "tscale"; teff=teff, logg=logg, feh=feh),
		sigdigits=1
	)
	
	hres = MUST.interpolate_quantity(grid, "hres"; teff=teff, logg=logg, feh=feh)

	av_path = abspath("$(name)_99999_av.dat")
	open(av_path, "w") do f
		MUST.writedlm(f, [-reverse(model.z) reverse(exp.(model.lnT)) reverse(model.lnρ)])
	end
	
	MUST.StaggerGrid(
        "interpolated",
        MUST.DataFrame(
		Dict(
			"name"     => name,
            "folder"   => folder,
            "mesh"     => mesh,
			"snapshot" => snapshot,
			"ma_x"     => ma_x,
			"mi_x"     => 0.0,
			"ma_y"     => ma_y,
			"mi_y"     => 0.0,
			"ma_z"     => ma_z,
			"mi_z"     => mi_z,
			"vmin"     => vmin,
			"vmax"     => vmax,
			"tscale"   => tscale,
			"hres"     => hres,
			"av_path"  => av_path,
			"teff"     => teff,
			"logg"     => logg,
			"feh"      => feh
		    )
	    )
    )
end

function interpolate_from_grid(grid::MUST.AbstractMUSTGrid, teff::F, logg::F, feh::F) where {F<:AbstractArray}
	subgrids = []

	for i in eachindex(teff)
		append!(subgrids, [
			interpolate_from_grid(grid, teff[i], logg[i], feh[i])
		])
	end

	g_new = deepcopy(subgrids[1])
	for i in eachindex(subgrids)
		i == 1 && continue

		g_new = g_new + subgrids[i]
	end

	g_new
end

"""
	interpolate_from_grid(grid::MUST.AbstractMUSTGrid; teff, logg, feh) = interpolate_from_grid(grid, teff, logg, feh)

Interpolate to teff, logg, feh within the given grid. Return a new grid with only the new entries. Average 1D models are created
with an interpolated resolution and resamples accordingly.

# Examples

```julia
grid = MUST.StaggerGrid("dispatch_grid.mgrid")
new_gridpoints = interpolate_from_grid(grid, teff=[4000, 5000], logg=[4.0, 4.5], feh=[0.5, 0.0])
```
"""
interpolate_from_grid(grid::MUST.AbstractMUSTGrid; teff, logg, feh) = interpolate_from_grid(grid, teff, logg, feh)

"""
	interpolate_from_grid!(grid::MUST.AbstractMUSTGrid; teff, logg, feh) = interpolate_from_grid(grid, teff, logg, feh)

Interpolate to teff, logg, feh within the given grid. Return a new grid with only the new entries. Average 1D models are created
with an interpolated resolution and resamples accordingly. Also include old grid.

# Examples

```julia
grid = MUST.StaggerGrid("dispatch_grid.mgrid")
new_and_old_gridpoints = interpolate_from_grid!(grid, teff=[4000, 5000], logg=[4.0, 4.5], feh=[0.5, 0.0])
```
"""
interpolate_from_grid!(grid::MUST.AbstractMUSTGrid; teff, logg, feh) = begin
	g = interpolate_from_grid(grid, teff, logg, feh)

	grid + g
end

