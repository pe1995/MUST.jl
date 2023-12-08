using MUST 
using TSO
using LazySets
using Polyhedra


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
	femask = feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]
	
	# read all average models and interpolate them to the same number of points
	models = [initial_model(grid, String(name))
					for name in grid.info[!, "name"]]
	r_models = [minimum(abs.(diff(m.z))) for m in models]
	
	models = [TSO.upsample(m, common_size) for m in models]

	# now we interpolate all points to one common point, for every point
	scatter_int(v, x, y) = MUST.pyconvert(typeof(x),
		first(MUST.scipy_interpolate.griddata((teff_gr[femask], logg_gr[femask]), v, ([x], [y]), method="linear"))
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
		z[i] = scatter_int(points[femask, 1], teff, logg)
		t[i] = scatter_int(points[femask, 2], teff, logg)
		d[i] = scatter_int(points[femask, 3], teff, logg)
	end

	r_target = scatter_int(r_models[femask], teff, logg)
	npoints  = ceil(Int, abs(maximum(z) - minimum(z)) / r_target)

	TSO.upsample(Model1D(z=z, lnT=t, lnρ=d, logg=logg), npoints)
end




#========================================================= Select Parameters =#

"""
	random_paramters(grid, N; [teff, logg, feh])

Random set of parameters within the convex hull of the grid.
Ranges for teff, logg and feh can be given as kwargs.
Random points will be within the give grid, so that a scatter interpolation
(e.g. scipy.griddata) can be used on them.

# Examples

```jldoctest
julia> random_paramters(grid, 10, teff=[4500, 6500], logg=[4.0, 4.5], feh=[0.0, 0.0])
10×3 Matrix{Float64}:
 5381.59  4.08404  0.0
 4639.25  4.44625  0.0
 6281.92  4.32518  0.0
 4661.24  4.47237  0.0
 5935.55  4.00568  0.0
 5893.84  4.32614  0.0
 4688.96  4.48417  0.0
 5157.83  4.21643  0.0
 6224.38  4.35535  0.0
 5161.6   4.35618  0.0
```
"""
function random_paramters(grid, N; 
	teff=[minimum(grid["teff"]), maximum(grid["teff"])], 
	logg=[minimum(grid["logg"]), maximum(grid["logg"])], 
	feh=[minimum(grid["feh"]), maximum(grid["feh"])])
	
	points = [[t, l, f] for (t, l, f) in zip(grid["teff"], grid["logg"], grid["feh"])]
	hull = convex_hull(points)
	P = VPolytope(hull)
	
	rand_points = zeros(N, 3)
	found = 0
	while found < N
		t = MUST.randrange(teff...)
		l = MUST.randrange(logg...) 
		f = MUST.randrange(feh...)
	
		if [t, l, f] ∈ P
			found += 1
			rand_points[found, 1] = t
			rand_points[found, 2] = l
			rand_points[found, 3] = f
		end
	end

	rand_points
end




#===========================================================interpolate grid =#

function interpolate_from_grid(grid::MUST.AbstractMUSTGrid, teff::F, logg::F, feh::F) where {F<:AbstractFloat}
	# create the iniitial model from interpolating the average snapshots
	model = interpolate_average(grid, teff=teff, logg=logg, feh=feh)

	folder = "interpolated"
	mesh   = "interpolated"

	tname = MUST.@sprintf "%.2f" teff/100
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

