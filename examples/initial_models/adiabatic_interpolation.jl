### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ de2144f0-0f5d-11ee-3765-7b42d60dc151
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using Plots
end

# ╔═╡ 42366c64-3395-4cdb-8c3c-d1f4a8706df6
md"# Adiabatic Interpolation
One approach to compute initial conditions from the avarege Stagger grid is to interpolate all quantities in Teff, logg (metallicity) to a new point. One could inteprolate the 1D models directly. However, one can also just build the intermediate adiabats based on the adiabatic relation given for the other models."

# ╔═╡ 548c36a3-3054-465f-90b4-e7ae058d5b02
md"## Stagger average grid"

# ╔═╡ b09f2438-d4f0-46c9-987f-77c0a25830a8
grid = MUST.StaggerGrid("dispatch_grid.mgrid");

# ╔═╡ ca27fb93-0f95-4240-a56e-ec6770185b51
@info(grid)

# ╔═╡ 776246d4-6b64-4338-8383-c9a4fc190bb1
md"## Adiabats from this grid
We can construct adiabats for each of the grid points by converting them to an optical model, and assume that e.g. at optical depth 5 the atmosphere is fully adiabatic"

# ╔═╡ 8628054c-7d59-40b9-b853-faa875513dd4
names(grid.info)

# ╔═╡ 574afe0c-3923-47e8-8fbb-a3c44a67c591
begin
	function equation_of_state(grid, point::Int; energy=false)
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
	
	function equation_of_state(grid, point::String; kwargs...)
		i = findfirst(grid.info[!, "name"] .== point)
		equation_of_state(grid, i; kwargs...)
	end
end

# ╔═╡ 4d2583c8-e5d0-4b76-b8ad-2303b4bd8e5a
begin
	function initial_model(grid, point::Int; kwargs...)
		eos, opa = equation_of_state(grid, point; kwargs...)
		@optical Average3D(eos, grid.info[point, "av_path"]; logg=grid.info[point, "logg"]) eos opa
	end
	
	function initial_model(grid, point::String; kwargs...)
		i = findfirst(grid.info[!, "name"] .== point)
		initial_model(grid, i; kwargs...)
	end
end

# ╔═╡ d6c53eed-2c03-4103-8b52-1f420fce7ed2
eos_test, opa_test = equation_of_state(grid, "t5777g44m00")

# ╔═╡ 1d1fddd2-a6bd-40c6-b808-a8b2bbc96cd5
model_test = initial_model(grid, "t45g20m00")

# ╔═╡ c0a9af80-f52f-4021-8ef3-bf66a1b8d301
begin
	function flip!(model)
		for f in fieldnames(typeof(model))
			v = getfield(model, f)
			if typeof(v) <: AbstractArray
				reverse!(v)
			end
		end
	end
	
	function flip(model)
		m = deepcopy(model)
		flip!(m)

		m
	end
end

# ╔═╡ c2c22e2d-e040-4698-83d5-741ec823e29b
function pick_point(model, i)
	args = Dict()
	for f in fieldnames(typeof(model))
		v = getfield(model, f)
		if typeof(v) <: AbstractArray
			args[f] = [v[i]]
		else
			args[f] = v
		end
	end

	Model1D(;args...)
end

# ╔═╡ 5b509c37-bd8d-4f92-8206-4b85184adeef
begin
	"""
		adiabat(model; kwargs...)
	
	Construct a 1D adiabatic model by integrating from the bottom to the top. 
	A new height scale is constructed automatically, it is stopped when the original 
	model reaches optical depth τtop. All is done in cgs units.
	"""
	function adiabat(model_in, eos_in::SqEoS; 
						nz=10000, n_iter=30, dlnd=0.05, padding=0.05)
		model = deepcopy(model_in)
		eos = @axed eos_in
		
		# First make sure the z scale is oriented in the correct direction
		if first(model.lnT) < last(model.lnT)
			#@warn "Flipping z axis..."
			flip!(model)
		end
	
		if first(model.z) > last(model.z)
			model.z .= - model.z
		end
		
		# Find zbottom and ztop
		zbottom = 1
		ztop    = length(model.z) 
		
		# Find the initial energy from the EoS. Assert that the EoS is in the 
		# correct units for this
		@assert is_internal_energy(eos)
		ee_min,ee_max,_,_ = exp.(limits(eos))
	
		# initial point
		t0  = exp(model.lnT[zbottom])
		ρ0  = exp(model.lnρ[zbottom])
		ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
		g0  = exp10(model.logg)
		p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))
		z0  = model.z[zbottom]

		# storage arrays
		z  = zeros(2*nz)
		d  = zeros(2*nz)
		ee = zeros(2*nz)
		t  = zeros(2*nz)
		p  = zeros(2*nz)
	
		# integrate from the starting point upwards, until z is larger than ztop
		# or the maximum number of points is reached
		dee   = 0
		z[nz]  = z0
		d[nz]  = ρ0  
		ee[nz] = ee0 
		t[nz]  = t0  
		p[nz]  = p0  
		i_end = 2*nz

		δz_sponge = abs(maximum(model.z) - minimum(model.z)) * padding
		
		# upwards
		for i in nz+1:2*nz
			ee[i] = max(ee[i-1]-dee,ee_min)
			
			for iter in 1:n_iter
		        d[i] = d[i-1]*exp(-dlnd)
		        t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
				p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
	        
		        dp = p[i-1] - p[i]
		        da = (d[i-1] + d[i]) /2.
		        dz = dp / (g0*da)
				
		        z[i]  = z[i-1] + dz
		        pd    = (p[i-1] + p[i]) / (d[i-1] + d[i])
		        dee   = pd*dlnd
		        ee[i] = max(ee[i-1]-dee, ee_min)
			end
	
			if z[i] > last(model.z) + δz_sponge
				i_end = i
				break
			end
		end

		# downwards
		i_start=1
		for i in nz:-1:1
			ee[i] = min(ee[i+1]+dee, ee_max)
			
			for iter in 1:n_iter
		        d[i] = d[i+1]*exp(dlnd)
		        t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
				p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
	        
		        dp = p[i] - p[i+1]
		        da = (d[i+1] + d[i]) /2.
		        dz = dp / (g0*da)
				
		        z[i]  = z[i+1] - dz
		        pd    = (p[i+1] + p[i]) / (d[i+1] + d[i])
		        dee   = pd*dlnd
		        ee[i] = min(ee[i+1]+dee, ee_max)
			end
	
			if z[i] < first(model.z) - δz_sponge
				i_start = i
				break
			end
		end
		
		z = z[i_start:i_end]
		d = d[i_start:i_end]
		ee = ee[i_start:i_end]
		t = t[i_start:i_end] 
		p = p[i_start:i_end]  

		TSO.Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=model.logg)
	end


	"""
		adiabat(model; kwargs...)
	
	Construct a 1D adiabat just from start and end points.
	"""
	function adiabat(star_point, end_point, eos_in::SqEoS; 
						nz=10000, n_iter=30, dlnd=0.05, padding=0.05)
		eos = @axed eos_in
		
		zbottom = first(star_point.z)
		ztop = first(end_point.z)

		# Find the initial energy from the EoS. Assert that the EoS is in the 
		# correct units for this
		@assert is_internal_energy(eos)
		ee_min,ee_max,_,_ = exp.(limits(eos))
	
		# initial point
		t0  = exp(first(star_point.lnT))
		ρ0  = exp(first(star_point.lnρ))
		ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
		g0  = exp10(star_point.logg)
		p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))
		z0  = first(star_point.z)

		# storage arrays
		z  = zeros(2*nz)
		d  = zeros(2*nz)
		ee = zeros(2*nz)
		t  = zeros(2*nz)
		p  = zeros(2*nz)
	
		# integrate from the starting point upwards, until z is larger than ztop
		# or the maximum number of points is reached
		dee    = 0
		z[nz]  = z0
		d[nz]  = ρ0  
		ee[nz] = ee0 
		t[nz]  = t0  
		p[nz]  = p0  
		i_end = 2*nz

		δz_sponge = abs(first(end_point.z) - first(star_point.z)) * padding
		
		# upwards
		for i in nz+1:2*nz
			ee[i] = max(ee[i-1]-dee,ee_min)
			
			for iter in 1:n_iter
		        d[i] = d[i-1]*exp(-dlnd)
		        t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
				p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
	        
		        dp = p[i-1] - p[i]
		        da = (d[i-1] + d[i]) /2.
		        dz = dp / (g0*da)
				
		        z[i]  = z[i-1] + dz
		        pd    = (p[i-1] + p[i]) / (d[i-1] + d[i])
		        dee   = pd*dlnd
		        ee[i] = max(ee[i-1]-dee, ee_min)
			end
	
			if z[i] > ztop + δz_sponge
				i_end = i
				break
			end
		end

		# downwards
		i_start=1
		for i in nz:-1:1
			ee[i] = min(ee[i+1]+dee, ee_max)
			
			for iter in 1:n_iter
		        d[i] = d[i+1]*exp(dlnd)
		        t[i] = exp(lookup(eos, :lnT,  log(d[i]), log(ee[i])))
				p[i] = exp(lookup(eos, :lnPg, log(d[i]), log(ee[i])))
	        
		        dp = p[i] - p[i+1]
		        da = (d[i+1] + d[i]) /2.
		        dz = dp / (g0*da)
				
		        z[i]  = z[i+1] - dz
		        pd    = (p[i+1] + p[i]) / (d[i+1] + d[i])
		        dee   = pd*dlnd
		        ee[i] = min(ee[i+1]+dee, ee_max)
			end
	
			if z[i] < zbottom - δz_sponge
				i_start = i
				break
			end
		end
		
		z = z[i_start:i_end]
		d = d[i_start:i_end]
		ee = ee[i_start:i_end]
		t = t[i_start:i_end] 
		p = p[i_start:i_end]  

		TSO.Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee),
						logg=star_point.logg)
	end

	adiabat(grid, name::String; kwargs...) = begin
		eos, opa = equation_of_state(grid, name; energy=true)
		model = initial_model(grid, name)

		# First make sure the z scale is oriented in the correct direction
		if first(model.lnT) < last(model.lnT)
			@warn "Flipping z axis..."
			flip!(model)
		end
	
		if first(model.z) > last(model.z)
			model.z .= - model.z
		end

		start_point = pick_point(model, 1)
		end_point = pick_point(model, length(model.z))
		
		adiabat(start_point, end_point, eos; kwargs...)
	end
end

# ╔═╡ 7a77e8bf-0654-4b0f-88c9-4b269de4a4b8
a = adiabat(grid, "t45g20m00")

# ╔═╡ a67ad91a-4026-4a9d-b86a-1287a333c491
begin
	plot(framestyle=:box, legend=false, grid=false)

	plot!(-model_test.z, model_test.lnT)
	plot!(a.z, a.lnT)
	
	
	plot!(xlabel="τ-ross", ylabel="ln T")
end

# ╔═╡ dc086eb1-e9db-457c-b616-21d306ff6855
md"## Interpolation from Grid
The most straight forward would be to pick the next 2 logg and Teff neighbors, and interpolate first to the right logg value, and then to the right Teff value. This requires that there actually are models with those quantities."

# ╔═╡ c618d758-f28d-409e-a817-92a9946cfc7a
TSO.load_scipy_interpolate!();

# ╔═╡ 80e46962-0451-4234-9f71-56464f7616ec
sci = TSO.scipy_interpolate

# ╔═╡ 327ce40c-559a-42c9-9d8b-28e56b7abd1b
function interpolate_grid(grid, what; teff, logg)
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	v    = grid.info[!, what]

	ip = first(sci.griddata((logg_gr, teff_gr), v, ([logg], [teff])))
end

# ╔═╡ fe36a795-c6d8-403a-81ae-32fbaca80d3a
vi = interpolate_grid(grid, "mi_z", teff=6000.0, logg=4.3)

# ╔═╡ e17d4794-6873-4316-a3ba-df096ff0f0b3
begin
	plot(framestyle=:box, legend=false, grid=false)

	scatter!(grid.info[!, "logg"], grid.info[!, "teff"],
		marker_z=grid.info[!, "mi_z"], 
		vmin=minimum(grid.info[!, "mi_z"]), vmax=maximum(grid.info[!, "mi_z"]),
		markersize=5, marker=:square)
	
	scatter!([4.3], [6000.0], marker_z=[vi],
		vmin=minimum(grid.info[!, "mi_z"]), vmax=maximum(grid.info[!, "mi_z"]),
		markersize=5, marker=:circle)
end

# ╔═╡ fd5ab461-b4e7-426e-95e1-c95ffb560ee0
md"What about interpolating the start and end of the adiabats?"

# ╔═╡ b631becb-eb0d-47c5-b0c4-693864f78938
function model_edges(grid)
	start_points = []
	end_points = []
	
	for name in grid.info[!, "name"][:]
		name = String(name)
		eos, opa = equation_of_state(grid, name; energy=true)
		model = initial_model(grid, name)
	
		# First make sure the z scale is oriented in the correct direction
		if first(model.lnT) < last(model.lnT)
			@warn "Flipping z axis..."
			flip!(model)
		end
	
		if first(model.z) > last(model.z)
			model.z .= - model.z
		end
	
		start_point = pick_point(model, 1)
		end_point = pick_point(model, length(model.z))
		append!(start_points, [start_point])
		append!(end_points, [end_point])
	end

	start_points, end_points
end

# ╔═╡ 4d41c3d8-91f3-4f57-8e86-c9fcbd6adacd
function interpolate_adiabat(grid; teff, logg, feh, kwargs...)
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	

	# get an eos from the closest model in chemistry and maybe size,
	# the same for all of them for now
	eos, opa = equation_of_state(grid, String(grid.info[!, "name"][1]); energy=true)

	s, e = model_edges(grid)
	zs = [first(m.z) for m in s]
	ze = [first(m.z) for m in e]

	ds = [first(m.lnρ) for m in s]
	de = [first(m.lnρ) for m in e]

	ts = [first(m.lnT) for m in s]
	te = [first(m.lnT) for m in e]
	
	scatter_int(v, x, y, z) = first(sci.griddata((logg_gr, teff_gr, feh_gr), v, ([x], [y], [z]), method="linear"))

	
	zs_m = scatter_int(zs, logg, teff, feh)
	ze_m = scatter_int(ze, logg, teff, feh)

	ds_m = scatter_int(ds, logg, teff, feh)
	de_m = scatter_int(de, logg, teff, feh)

	ts_m = scatter_int(ts, logg, teff, feh)
	te_m = scatter_int(te, logg, teff, feh)

	sm = TSO.Model1D(z=[zs_m], lnρ=[ds_m], lnT=[ts_m])
	em = TSO.Model1D(z=[ze_m], lnρ=[de_m], lnT=[te_m])

	am = adiabat(sm, em, eos; kwargs...)

	sm, em, am
end

# ╔═╡ d1feff77-c70d-4605-893c-0d45c4e8fd8f
sm, em, am = interpolate_adiabat(grid, teff=5800.0, logg=4.5, feh=0.0, nz=400)

# ╔═╡ 4456c7ee-6ab8-485e-bcfa-fb406a2d1f45
begin
	plot(framestyle=:box, grid=false, legendposition=:bottom)

	for (i,name) in enumerate(grid.info[!, "name"])
		teff = grid.info[i, "teff"]
		logg = grid.info[i, "logg"]
		if logg < 4
			continue
		end
		mi = adiabat(grid, String(name), nz=400)
		plot!(mi.z, mi.lnT, label="teff: $(teff), logg: $(logg)")
	end
	
	scatter!(am.z, am.lnT, color=:black, lw=3)
	#hline!(sm.lnT, color=:red)
	#hline!(em.lnT, color=:red)

	#vline!(sm.z, color=:red)
	#vline!(em.z, color=:red)

	plot!(xlim=(-3e8, 1.5e8))
	plot!(ylim=(7,10.2))
	
	plot!(xlabel="z", ylabel="ln T")
end

# ╔═╡ a46da4ce-38aa-41e9-b624-15d33efc2860
md"This interpolated model can be used as a starting point for DISPATCH."

# ╔═╡ 4e6979cf-1d4f-4bf0-bfdb-8cb500817462
md"## Alternative: Interpolate average Stagger models point-wise and direct
One alternative would be to inteprolate the 1D models directly. The culprit: Also the size of the boxes change, so the models can not be interpolated on a common z scale because one will always be out-of-bounce. One idea would to interpolate both to the same number of uniformly spaced points, and then interpolate between those points regardless of their location in the atmosphere. Lets see if this returns something sensible."

# ╔═╡ 77fcef55-5922-4907-ad34-3441f6cb73fe
function interpolate_average(grid; teff, logg, feh, common_size=1000, kwargs...)
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	
	# get an eos from the closest model in chemistry and maybe size,
	# the same for all of them for now
	eos, opa = equation_of_state(grid, String(grid.info[!, "name"][1]); energy=true)

	# read all average models and interpolate them to the same number of points
	models = [initial_model(grid, String(name))
					for name in grid.info[!, "name"]]
	models = [TSO.upsample(m, common_size) for m in models]

	# now we interpolate all points to one common point, for every point
	scatter_int(v, x, y, z) = first(sci.griddata((logg_gr, teff_gr, feh_gr), v, ([x], [y], [z]), method="linear"))
	
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

	Model1D(z=z, lnT=t, lnρ=d, logg=logg)
end

# ╔═╡ fa36e059-a5da-48aa-ac1e-79c5f5d098bc
im = interpolate_average(grid, teff=6000.0, logg=4.5, feh=0.0)

# ╔═╡ 064a1d54-d2ba-4115-82b0-bb00f2cac238
begin
	plot(framestyle=:box, grid=false, 
		legendposition=:bottomleft, 
		legendforegroundcolor=nothing,
		legendbackgroundcolor=nothing, legendcolumns=2)

	for (i,name) in enumerate(grid.info[!, "name"])
		teff = grid.info[i, "teff"]
		logg = grid.info[i, "logg"]
		mi = initial_model(grid, String(name))
		plot!(-mi.z, mi.lnT, label="teff: $(teff), logg: $(logg)")
	end
	
	scatter!(-im.z, im.lnT, color=:black, lw=3)
	#hline!(sm.lnT, color=:red)
	#hline!(em.lnT, color=:red)

	#vline!(sm.z, color=:red)
	#vline!(em.z, color=:red)

	plot!(xlim=(-3e8, 1.5e8))
	plot!(ylim=(7,10.2))
	
	plot!(xlabel="z", ylabel="ln T")
end

# ╔═╡ 342755c5-7f69-455d-8c1e-bd6a491063a8
md"The big advantage of this approach is that one can directly use this to create a new opacity binning based on this model. One can furthermore direclty read off the z extend of the new model."

# ╔═╡ 37fea989-8640-4318-9da9-8d0342994ea3
md"## End-to-end generation of models on the HRD
The goal is to generate any model on the HRD by using the grid and a teff + logg. This function should create the initial adiabat for the interpolated boundaries, and then interpolate every single quantity in the table. The `mi_z`, `ma_z`, etc. might need a little adjustment in order to fit into the intial adiabat. After everything is interpolated, it should be saved as a namelist and copy the EoS over. This is the most important issue I see with this approach: It makes it impossible to re-compute the binning, one has to reply that the next closest model is reasonably close."

# ╔═╡ 6dce1923-3826-4366-ba94-fd56a663a6b1
md"An alternative is to also interpolate the average model, such that every point, including the z axis, is interpolated."

# ╔═╡ Cell order:
# ╟─42366c64-3395-4cdb-8c3c-d1f4a8706df6
# ╠═de2144f0-0f5d-11ee-3765-7b42d60dc151
# ╟─548c36a3-3054-465f-90b4-e7ae058d5b02
# ╠═b09f2438-d4f0-46c9-987f-77c0a25830a8
# ╠═ca27fb93-0f95-4240-a56e-ec6770185b51
# ╟─776246d4-6b64-4338-8383-c9a4fc190bb1
# ╠═8628054c-7d59-40b9-b853-faa875513dd4
# ╟─574afe0c-3923-47e8-8fbb-a3c44a67c591
# ╟─4d2583c8-e5d0-4b76-b8ad-2303b4bd8e5a
# ╠═d6c53eed-2c03-4103-8b52-1f420fce7ed2
# ╠═1d1fddd2-a6bd-40c6-b808-a8b2bbc96cd5
# ╟─c0a9af80-f52f-4021-8ef3-bf66a1b8d301
# ╟─c2c22e2d-e040-4698-83d5-741ec823e29b
# ╟─5b509c37-bd8d-4f92-8206-4b85184adeef
# ╠═7a77e8bf-0654-4b0f-88c9-4b269de4a4b8
# ╟─a67ad91a-4026-4a9d-b86a-1287a333c491
# ╟─dc086eb1-e9db-457c-b616-21d306ff6855
# ╠═c618d758-f28d-409e-a817-92a9946cfc7a
# ╠═80e46962-0451-4234-9f71-56464f7616ec
# ╟─327ce40c-559a-42c9-9d8b-28e56b7abd1b
# ╠═fe36a795-c6d8-403a-81ae-32fbaca80d3a
# ╟─e17d4794-6873-4316-a3ba-df096ff0f0b3
# ╟─fd5ab461-b4e7-426e-95e1-c95ffb560ee0
# ╟─b631becb-eb0d-47c5-b0c4-693864f78938
# ╟─4d41c3d8-91f3-4f57-8e86-c9fcbd6adacd
# ╠═d1feff77-c70d-4605-893c-0d45c4e8fd8f
# ╟─4456c7ee-6ab8-485e-bcfa-fb406a2d1f45
# ╟─a46da4ce-38aa-41e9-b624-15d33efc2860
# ╟─4e6979cf-1d4f-4bf0-bfdb-8cb500817462
# ╟─77fcef55-5922-4907-ad34-3441f6cb73fe
# ╠═fa36e059-a5da-48aa-ac1e-79c5f5d098bc
# ╟─064a1d54-d2ba-4115-82b0-bb00f2cac238
# ╟─342755c5-7f69-455d-8c1e-bd6a491063a8
# ╟─37fea989-8640-4318-9da9-8d0342994ea3
# ╟─6dce1923-3826-4366-ba94-fd56a663a6b1
