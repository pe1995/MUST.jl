using MUST 
using TSO
using LazySets
using Polyhedra
using NaturalNeighbours

interpolationMethod="triangle"

# read the grid of interpolation coefficients
ip_coeffs_stagger = MUST.CSV.read(abspath(joinpath(dirname(@__FILE__), "..", "initial_grids", "Stagger", "interpolation_coeffs.csv")), MUST.DataFrame)


#==================================================================== Models =#

function initial_model(grid::MUST.AbstractMUSTGrid, point::Int; kwargs...)
	try
		eos, opa = equation_of_state(grid, point; kwargs...)
		a = Average3D(eos, grid.info[point, "av_path"]; logg=grid.info[point, "logg"])
		@optical a eos opa
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

find_closest_eos(grid::MUST.AbstractMUSTGrid, teff, feh; from="eos_root") = begin
	mother_tables = grid[from]
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
			mothertablesminall = mother_tables[matchingtables]
			imin = argmin(abs.(tempminall .- teff))

			mothertablesminall[imin]
		end
	end

	path
end





#================================================ interpolate average models =#

read_models_from_grid(grid::MUST.AbstractMUSTGrid; adiabats=nothing, eos=nothing, opa=nothing, common_size=1000, τbottom=8, τtop=-6, kwargs...) = begin
	# read all average models and interpolate them to the same number of points
	models = Any[initial_model(grid, String(name))
					for name in grid.info[!, "name"]]

	models = if !isnothing(adiabats)
		Any[TSO.flip!(Average3D(grid.info[i, "ad_path"], logg=grid.info[i, "logg"]) for i in 1:MUST.nrow(grid.info))]
	else
		models
	end

	eos = if ((!isnothing(eos)) & ("matching_eos" in names(grid.info)))
		@info "EoS as given in the grid is loaded for grid nodes opacity."
		[reload(
			SqEoS, grid.info[i, "matching_eos"]
		) for i in 1:MUST.nrow(grid.info)]
	else
		[eos for _ in 1:MUST.nrow(grid.info)]
	end

	models_t = if !isnothing(eos)
		ltscale = range(τbottom, τtop, length=common_size) |> collect
		for i in eachindex(models)
			models[i] = if isnothing(models[i])
				nothing
			else
				# check if there is a model that has been averaged on the optical depth scale
				if "avo_path" in names(grid.info)
					m = Average3D(grid.info[i, "avo_path"], logg=grid.info[i, "logg"])
					m.τ .= TSO.rosseland_optical_depth(eos[i], m)

					m
				else
					# If this is a non-optical model we need to add the opacity from the given EoS
					# this is not ideal! Ideally the models already contain an EoS
					if typeof(models[i]) <: TSO.Model1D
						if isnothing(opa)
							@optical(models[i], eos[i])
						else
							@optical(models[i], eos[i], opa[i])
						end
					end
				end
			end
		end
		[isnothing(m) ? nothing : TSO.interpolate_to(m; τ=ltscale, in_log=true) for m in models]
	else
		for i in eachindex(models)
			models[i] = if isnothing(models[i])
				nothing
			else
				# check if there is a model that has been averaged on the optical depth scale
				TSO.upsample(models[i], common_size)
			end
		end
	end

	models, models_t
end

add_vres_from_grid!(grid::MUST.AbstractMUSTGrid) = begin
	models = Any[initial_model(grid, String(name)) for name in grid.info[!, "name"]]
	r_models = [minimum(abs.(diff(m.z))) for m in models]
	grid.info[!, "vres"] = r_models
end




function interpolate_average(grid::MUST.AbstractMUSTGrid; teff, logg, feh, common_size=1000, models=nothing, models_mod=nothing, kwargs...)
	@info "Interpolating $teff, $logg, $feh point by point."

	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	femask = feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]
	
	# read all average models and interpolate them to the same number of points
	models, models_mod = if isnothing(models)
		read_models_from_grid(
			grid; 
			adiabats=adiabats, eos=eos, opa=opa, 
			common_size=common_size, 
			kwargs...
		) 
	else
		models, models_mod
	end

	r_models = [minimum(abs.(diff(m.z))) for m in models]

	# now we interpolate all points to one common point, for every point
	norm_teff(t) = (t - minimum(teff_gr[femask])) ./ (maximum(teff_gr[femask]) - minimum(teff_gr[femask]))
	norm_logg(t) = (t - minimum(logg_gr[femask])) ./ (maximum(logg_gr[femask]) - minimum(logg_gr[femask]))

	scatter_int = if interpolationMethod=="python"
		(v, x, y) -> MUST.pyconvert(typeof(x),
			first(MUST.scipy_interpolate.griddata((norm_teff.(teff_gr[femask]), norm_logg.(logg_gr[femask])), v, ([x], [y]), method="linear"))
		)
	else
		nodes_x = zeros(count(femask))
		nodes_y = zeros(count(femask))
		nodes_x .= norm_teff.(teff_gr[femask])
		nodes_y .= norm_logg.(logg_gr[femask])
		(v, x, y) -> NaturalNeighbours.interpolate(nodes_x, nodes_y, v, method=Triangle())(x, y)
	end

	points = zeros(eltype(models_mod[1].z), length(models_mod), 3)
	z = zeros(eltype(models_mod[1].z), common_size)
	t = zeros(eltype(models_mod[1].z), common_size)
	d = zeros(eltype(models_mod[1].z), common_size)
	
	for i in 1:common_size
		for j in eachindex(models)
			points[j, 1] = models_mod[j].z[i]
			points[j, 2] = models_mod[j].lnT[i]
			points[j, 3] = models_mod[j].lnρ[i]
		end
		z[i] = scatter_int(points[femask, 1], norm_teff(teff), norm_logg(logg))
		t[i] = scatter_int(points[femask, 2], norm_teff(teff), norm_logg(logg))
		d[i] = scatter_int(points[femask, 3], norm_teff(teff), norm_logg(logg))
	end

	#r_target = scatter_int(r_models[femask], teff, logg)
	r_target = MUST.interpolate_quantity(grid, "vres"; teff=teff, logg=logg, feh=feh)
	npoints  = ceil(Int, abs(maximum(z) - minimum(z)) / r_target)

	TSO.upsample(Model1D(z=z, lnT=t, lnρ=d, logg=logg), npoints)
end

function interpolate_average(grid::MUST.AbstractMUSTGrid, eos::SqEoS, opa=nothing; 
	teff, logg, feh, common_size=1000, τbottom=8, τtop=-6, from_snapshot="", adiabats=nothing, models=nothing, models_mod=nothing, τ_extrapolate=nothing, kwargs...)
	
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	femask = feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]
	
	# read all average models and interpolate them to the same number of points
	models, models_mod = if isnothing(models)
		read_models_from_grid(
			grid; 
			adiabats=adiabats, eos=eos, opa=opa, 
			common_size=common_size, 
			τbottom=τbottom, τtop=τtop, kwargs...
		) 
	else
		models, models_mod
	end

	r_models = [minimum(abs.(diff(m.z))) for m in models]
	ltscale = log10.(first(models_mod).τ)

	norm_teff(t) = (t - minimum(teff_gr[femask])) ./ (maximum(teff_gr[femask]) - minimum(teff_gr[femask]))
	norm_logg(t) = (t - minimum(logg_gr[femask])) ./ (maximum(logg_gr[femask]) - minimum(logg_gr[femask]))

	scatter_int = if interpolationMethod=="python"
		(v, x, y) -> MUST.pyconvert(typeof(x),
			first(MUST.scipy_interpolate.griddata((norm_teff.(teff_gr[femask]), norm_logg.(logg_gr[femask])), v, ([x], [y]), method="linear"))
		)
	else
		nodes_x = zeros(count(femask))
		nodes_y = zeros(count(femask))
		nodes_x .= norm_teff.(teff_gr[femask])
		nodes_y .= norm_logg.(logg_gr[femask])
		(v, x, y) -> NaturalNeighbours.interpolate(nodes_x, nodes_y, v, method=Triangle())(x, y)
	end
	
	mnew = if length(from_snapshot)==0
		@info "Interpolating $teff, $logg, $feh on optical depth scale."

		points = zeros(eltype(models_mod[1].z), length(models), 3)
		z = zeros(eltype(models_mod[1].z), common_size)
		t = zeros(eltype(models_mod[1].z), common_size)
		d = zeros(eltype(models_mod[1].z), common_size)
		
		for i in 1:common_size
			for j in eachindex(models_mod)
				isnothing(models_mod[j]) && continue
				points[j, 1] = models_mod[j].z[i]
				points[j, 2] = models_mod[j].lnT[i]
				points[j, 3] = models_mod[j].lnρ[i]
			end
			z[i] = scatter_int(points[femask, 1], norm_teff(teff), norm_logg(logg))
			t[i] = scatter_int(points[femask, 2], norm_teff(teff), norm_logg(logg))
			d[i] = scatter_int(points[femask, 3], norm_teff(teff), norm_logg(logg))
		end
		
		# interpolate to equally spaced in z
		# for this we need to construct a new z scale 
		TSO.flip(Model1D(z=z, lnT=t, lnρ=d, logg=logg, τ=exp10.(ltscale)), depth=true)
	else
		@info "Averaging snapshot from $(from_snapshot)."

		# load the snapshot
		b = MUST.Box(from_snapshot)
		lt, lnrho = MUST.profile(MUST.mean, b, :log10τ_ross, :logd)
		_, lnT = MUST.profile(MUST.mean, b, :log10τ_ross, :logT)
		_, z = MUST.profile(MUST.mean, b, :log10τ_ross, :z)
		m = TSO.flip(TSO.Model1D(z=z, τ=exp10.(lt), lnρ=lnrho, lnT=lnT, logg=Float32(logg)), depth=true)
		m = TSO.flip(TSO.interpolate_to(m, τ=range(τbottom, τtop, length=length(ltscale))|>collect, in_log=true), depth=true)

		m
	end

	mnew = if (!isnothing(τ_extrapolate)) && (τ_extrapolate>τtop)
		@info "Extrapolating to log10τ=$(τ_extrapolate)"
		TSO.interpolate_to(mnew, τ=range(maximum(log10.(mnew.τ)), τ_extrapolate, length=length(mnew.τ))|>collect, in_log=true)
	else
		mnew
	end

	# recompute z scale from opacity
	znew = TSO.rosseland_depth(eos, mnew)
	mnew.z .= znew

	# find the optical surface
	TSO.optical_surface!(mnew)
	TSO.flip!(mnew)

	# upsample to match the desired interpolated resolution
	#r_target = scatter_int(r_models[femask], teff, logg)
	r_target = MUST.interpolate_quantity(grid, "vres"; teff=teff, logg=logg, feh=feh) 
	npoints  = 500 #ceil(Int, abs(maximum(mnew.z) - minimum(mnew.z)) / r_target)

	# We now can interpolate everthing to this new z scale
	mnew = TSO.interpolate_to(
		mnew, 
		z=collect(range(maximum(mnew.z), minimum(mnew.z), length=npoints))
	)
	mnew
end

adiabat_from_model(model, eos) = begin
	m = TSO.flip(model)
	s = TSO.pick_point(m, 1)
	e = TSO.pick_point(m, length(m.z))
	TSO.adiabat(
		s, e, eos
	)
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
	
	feh_lim = if (minimum(feh) < minimum(grid["feh"])) | (maximum(feh) > maximum(grid["feh"]))
		@warn("Your chosen metallicity limits are outside the grid of initial model!")
		feh_min = min(max(minimum(feh), minimum(grid["feh"])), maximum(grid["feh"]))
		feh_max = max(min(maximum(feh), maximum(grid["feh"])), minimum(grid["feh"]))
		[feh_min, feh_max]
	else
		feh
	end

	points = [[t, l, f] for (t, l, f) in zip(grid["teff"], grid["logg"], grid["feh"])]
	hull = convex_hull(points)
	P = VPolytope(hull)
	
	rand_points = zeros(N, 3)
	found = 0
	while found < N
		t = MUST.randrange(teff...)
		l = MUST.randrange(logg...) 
		f = MUST.randrange(feh_lim...)
	
		if [t, l, f] ∈ P
			found += 1
			rand_points[found, 1] = t
			rand_points[found, 2] = l
			rand_points[found, 3] = f
		end
	end

	# assume that all models have average metallicity, if we could not sample
	# it due to the limits of the grid
	if (minimum(feh) < minimum(grid["feh"])) | (maximum(feh) > maximum(grid["feh"]))
		rand_points[:, 3] .= sum(feh) / 2.0
	end
	rand_points
end

"""
	random_paramters(grid, N, lowerlimit, upperlimit; [teff, logg, feh])

Random set of parameters within the convex hull of the grid.
Ranges for teff, logg and feh can be given as kwargs.
Random points will be within the give grid, so that a scatter interpolation
(e.g. scipy.griddata) can be used on them. Limit the selection between the functions
`lowerlimit` and `upperlimit`, which will are functions teff->logg setting the selection domain.
"""
function random_paramters(grid, N, lowerlimit, upperlimit; 
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
			# check if it is below the lower or above the upper limit
			lgLowLim = min(lowerlimit(t), upperlimit(t))
			lgUpLim = max(lowerlimit(t), upperlimit(t))

			if (l>=lgLowLim) & (l<=lgUpLim)
				found += 1
				rand_points[found, 1] = t
				rand_points[found, 2] = l
				rand_points[found, 3] = f
			end
		end
	end

	rand_points
end

lineLimit(teff, x1, x2) = begin
	m = (x2[2] - x1[2]) / (x2[1] - x1[1]) 
	m* teff + x2[2] - m*x2[1]
end

checkparameters(grid; teff, logg, feh) = begin
	points = [[t, l, f] for (t, l, f) in zip(grid["teff"], grid["logg"], grid["feh"])]
	hull = convex_hull(points)
	P = VPolytope(hull)
	[p ∈ P for p in points]
end

checkparameters(grid, paras) = checkparameters(grid, teff=paras[:, 1], logg=paras[:, 2], feh=paras[:, 3])



#========================================== interpolation of bottom boundary =#

get_ip_coeffs(feh, grid=ip_coeffs_stagger; offsets=nothing) = begin
	feh_grid = grid[!, :metallicity]
	feh_pick = feh_grid[argmin(abs.(feh_grid .- feh))]
	c = grid[grid[!, :metallicity] .== feh_pick, :]

	d = Dict(
		k=>parse.(Float64, split((c[!, k] |> first)[2:end-1], ','))
		for k in names(grid) if !(k=="metallicity")
	)

	if !isnothing(offsets)
		@info "Applying offsets: $(offsets)"
		d["rho_bottom_logg"][2] += offsets[1]
		d["T_bottom_teff"][2] += offsets[2]
		d["rho_top_logg"][2] += offsets[3]
		d["T_top_teff"][2] += offsets[4]
	end

	d
end

predict_bottom_boundary(teff, logg, coefficients) = begin
	rho_predicted = coefficients["rho_bottom_logg"][1] * logg + coefficients["rho_bottom_logg"][2] 
	rho_correction = coefficients["rho_bottom_correction_teff"][1] * teff + coefficients["rho_bottom_correction_teff"][2] 

	T_predicted = coefficients["T_bottom_teff"][1] * teff + coefficients["T_bottom_teff"][2] 
	T_correction = coefficients["T_bottom_correction_logg"][1] * logg + coefficients["T_bottom_correction_logg"][2] 

	rho_predicted + rho_correction, T_predicted + T_correction
end

predict_top_boundary(teff, logg, coefficients) = begin
	rho_predicted = coefficients["rho_top_logg"][1] * logg + coefficients["rho_top_logg"][2] 
	rho_correction = coefficients["rho_top_correction_teff"][1] * teff + coefficients["rho_top_correction_teff"][2] 

	T_predicted = coefficients["T_top_teff"][1] * teff + coefficients["T_top_teff"][2] 
	T_correction = coefficients["T_top_correction_logg"][1] * logg + coefficients["T_top_correction_logg"][2] 

	rho_predicted + rho_correction, T_predicted + T_correction
end


#===========================================================interpolate grid =#

function interpolate_from_grid(grid::MUST.AbstractMUSTGrid, teff::F, logg::F, feh::F; 
	eos=nothing, opa=nothing, adiabats=nothing, models=nothing, models_mod=nothing, folder="", 
	τ_extrapolate=nothing, adiabatic_extrapolation=false, boundary_extrapolation=false, extrapolation_offsets=nothing, τbottom=8, τtop=-6, τbottom_extrapolate=nothing, kwargs...) where {F<:AbstractFloat}

	folder_name = "interpolated"
	mesh   = "interpolated"

	tname = MUST.@sprintf "%.2f" teff/100
	gname = MUST.@sprintf "%.2f" logg*10
	fname = MUST.@sprintf "%.3f" feh
	name = "t$(tname)g$(gname)m$(fname)"

	snapshot = "interpolated"
	ma_x = MUST.interpolate_quantity(grid, "ma_x"; teff=teff, logg=logg, feh=feh)
	ma_y = MUST.interpolate_quantity(grid, "ma_y"; teff=teff, logg=logg, feh=feh)
	vmin = MUST.interpolate_quantity(grid, "vmin"; teff=teff, logg=logg, feh=feh)
	vmax = MUST.interpolate_quantity(grid, "vmax"; teff=teff, logg=logg, feh=feh)
	hres = MUST.interpolate_quantity(grid, "hres"; teff=teff, logg=logg, feh=feh)
	tscale = round(
		MUST.interpolate_quantity(grid, "tscale"; teff=teff, logg=logg, feh=feh),
		sigdigits=1
	)

	if !("vres" in names(grid.info))
		add_vres_from_grid!(grid)
	end
	vres = MUST.interpolate_quantity(grid, "vres"; teff=teff, logg=logg, feh=feh)


	# create the initial model from interpolating the average snapshots
	model = if isnothing(models)
		isnothing(eos) ? 
		interpolate_average(grid; teff=teff, logg=logg, feh=feh, kwargs...) : 
		interpolate_average(grid, eos, opa; teff=teff, logg=logg, feh=feh, adiabats=adiabat, τbottom=τbottom, τtop=τtop, kwargs...)
	else
		isnothing(eos) ? 
		interpolate_average(grid; teff=teff, logg=logg, feh=feh, models=models, models_mod=models_mod, kwargs...) : 
		interpolate_average(grid, eos, opa; teff=teff, logg=logg, feh=feh, adiabats=adiabats, models=models, models_mod=models_mod, τ_extrapolate=τ_extrapolate, τbottom=τbottom, τtop=τtop, kwargs...)
	end

	model = if boundary_extrapolation && (!isnothing(eos))
		# pick the interpolation coefficients that are closest to the given metallicity
		coeffs = get_ip_coeffs(feh; offsets=extrapolation_offsets)
		rho_bot, t_bot = predict_bottom_boundary(teff, logg, coeffs)
		rho_top, t_top = predict_top_boundary(teff, logg, coeffs)
		 
		mc = @optical(TSO.flip(deepcopy(model), depth=true), eos)

		# shift down the MARCS model based on the 
		if extrapolation_offsets[3] > -98.0
			offset_T = exp10(t_top) .- exp.(minimum(mc.lnT))
			offset_ρ = exp10(rho_top) .- exp.(minimum(mc.lnρ))
			mc.lnT .= log.(exp.(mc.lnT) .+ offset_T)
			mc.lnρ .= log.(exp.(mc.lnρ) .+ offset_ρ)
		else
			@info "Skipping top boundary shifting."
		end

		model_extra = TSO.reverse_adiabatic_extrapolation(
			mc, exp10(rho_bot), exp10(t_bot), eos; kwargs...
		) |> TSO.monotonic

		uniform_z = range(
			minimum(model_extra.z), maximum(model_extra.z), length=500
		) |> collect

		TSO.flip(TSO.interpolate_to(model_extra, in_log=false, z=uniform_z), depth=true)
	else
		TSO.flip(model, depth=true)
	end

	# optionally extrapolate adiabatically on this new zscale
	model = if adiabatic_extrapolation && (!isnothing(eos)) && (!isnothing(τbottom_extrapolate))
		# extrapolate until τbottom_extrapolate is reached
		model_extra = TSO.adiabatic_extrapolation(
			model, eos; τ_target=τbottom_extrapolate
		)

		model_extra = if (extrapolation_offsets[1] > 0) | (extrapolation_offsets[2] > 0) 
			offset_T = exp(extrapolation_offsets[2].+ maximum(model_extra.lnT))
			offset_ρ = exp(extrapolation_offsets[1] .+ maximum(model_extra.lnρ))
			@info "shifting bottom boundary to T=$(exp.(maximum(model_extra.lnT)))->$(offset_T), ρ=$(exp.(maximum(model_extra.lnρ)))->$(offset_ρ)"

			me = TSO.reverse_adiabatic_extrapolation(
				model_extra, offset_ρ, offset_T, eos; kwargs...
			) |> TSO.monotonic

			# now extrapolate the model back in case that has changed
			TSO.adiabatic_extrapolation(
				me, eos; τ_target=τbottom_extrapolate
			)
		else
			model_extra
		end

		uniform_z = range(
			minimum(model_extra.z), maximum(model_extra.z), length=500
		) |> collect

		TSO.flip(TSO.interpolate_to(model_extra, in_log=false, z=uniform_z), depth=true)
	elseif adiabatic_extrapolation
		mi_z = MUST.interpolate_quantity(grid, "mi_z"; teff=teff, logg=logg, feh=feh)
		ma_z = MUST.interpolate_quantity(grid, "ma_z"; teff=teff, logg=logg, feh=feh)
		model_extra = TSO.adiabatic_extrapolation(model, eos, abs(ma_z - mi_z)*1.1)
		uniform_z = range(
			minimum(model_extra.z), maximum(model_extra.z), length=500 #step=abs(vres)
		) |> collect

		TSO.flip(TSO.interpolate_to(model_extra, in_log=false, z=uniform_z), depth=true)
	else
		model
	end

	model = if !isnothing(eos)
		# interpolate to the final grid in optical depth (also extrapolate if needed)
		mnew = @optical model eos
		ltau = range(maximum(log10.(mnew.τ)), τtop, length=500) |> collect
		TSO.flip(TSO.interpolate_to(mnew, τ=ltau, in_log=true), depth=true)
	else
		model
	end
	
	mi_z = minimum(model.z)
	ma_z = maximum(model.z)
	
	av_path = abspath(joinpath(folder, "$(name)_99999_av.dat"))
	TSO.flip!(model)
	open(av_path, "w") do f
		MUST.writedlm(f, [model.z exp.(model.lnT) model.lnρ])
	end

	# Also save in the M3D format
	av_m3d_path = abspath(joinpath(folder, "$(name)_99999_m3d.dat"))
	TSO.save_text_m3d(model, av_m3d_path, header=name)
	
	MUST.Atmos1DGrid(
        "interpolated",
		abspath(joinpath(pwd(), "interpolated.mgrid")),
        MUST.DataFrame(
		Dict(
			"name"     => name,
            "folder"   => folder_name,
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
			"abs_av_path" => av_path,
			"teff"     => teff,
			"logg"     => logg,
			"feh"      => feh
		    )
	    )
    )
end

function interpolate_from_grid(grid_in::MUST.AbstractMUSTGrid, teff::F, logg::F, feh::F; 
	eos::F2=[nothing for _ in teff], opa::F3=[nothing for _ in teff], adiabats=nothing, folder="", τ_extrapolate=nothing, kwargs...) where {F<:AbstractArray, F2<:AbstractArray, F3<:AbstractArray}
	subgrids = []

	grid = deepcopy(grid_in)
	models, models_mod = read_models_from_grid(grid; eos=first(eos), opa=first(opa), adiabats=adiabats, kwargs...)
	
	if !("vres" in names(grid.info))
		@info "Vertical resolution is computed from spacing in average models of the grid."
		add_vres_from_grid!(grid)
	end

	for i in eachindex(teff)
		append!(subgrids, [
			interpolate_from_grid(
				grid, teff[i], logg[i], feh[i]; 
				eos=eos[i], opa=opa[i], 
				adiabats=adiabats, 
				models=models, models_mod=models_mod, folder=folder, τ_extrapolate=τ_extrapolate, kwargs...
			)
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
grid = MUST.Atmos1DGrid("dispatch_grid.mgrid")
new_gridpoints = interpolate_from_grid(grid, teff=[4000, 5000], logg=[4.0, 4.5], feh=[0.5, 0.0])
```
"""
interpolate_from_grid(grid::MUST.AbstractMUSTGrid; teff, logg, feh, kwargs...) = interpolate_from_grid(grid, teff, logg, feh; kwargs...)

"""
	interpolate_from_grid!(grid::MUST.AbstractMUSTGrid; teff, logg, feh) = interpolate_from_grid(grid, teff, logg, feh)

Interpolate to teff, logg, feh within the given grid. Return a new grid with only the new entries. Average 1D models are created
with an interpolated resolution and resamples accordingly. Also include old grid.

# Examples

```julia
grid = MUST.Atmos1DGrid("dispatch_grid.mgrid")
new_and_old_gridpoints = interpolate_from_grid!(grid, teff=[4000, 5000], logg=[4.0, 4.5], feh=[0.5, 0.0])
```
"""
interpolate_from_grid!(grid::MUST.AbstractMUSTGrid; teff, logg, feh, kwargs...) = begin
	g = interpolate_from_grid(grid, teff, logg, feh; kwargs...)

	grid + g
end