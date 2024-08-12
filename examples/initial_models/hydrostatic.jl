### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ b04a8d9a-5307-11ef-1268-95d5b8233977
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using DataFrames  
end

# ╔═╡ 273c40e1-f4ef-4f24-93dc-b9bb95f3abeb
plt = matplotlib.pyplot;

# ╔═╡ 77e9b680-d0e4-431f-9a59-a4a081ffd65e
md"Problem: The initial models for hot stars seem to look terrible. They crash almost immediatelly. The idea now is to look at the MARCS models and see how far off the interpolted models are. If they are far off, it might be an idea to extend the MARCS models deeper inside the star and see how that goes."

# ╔═╡ 63ac3293-0d3f-4533-8fc4-ae3695e40cc6
#mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"

# ╔═╡ 62d17d2d-6727-4632-b086-30cc6e828360
mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m10_a4_c3_vmic2_v5.0"

# ╔═╡ 83648aa8-11dd-425b-bc92-938a360e03a0
#mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m0_a0_vmic1_v5.0"

# ╔═╡ ce43f8fe-6274-4da3-b427-1f362352f1a4
#eos_mother_path = "ross_combined_eos_magg_m0_a0.hdf5"

# ╔═╡ e607ddf1-fb54-4c6c-842c-40003f2604e2
eos_mother_path = "combined_eos_magg_m10_a4_c3_vmic2.hdf5"

# ╔═╡ fc9633c5-fbd8-458c-a7de-38875987eb3b
#eos_mother_path = "combined_eos_magg_m0_a0_vmic1.hdf5"

# ╔═╡ 4fe7acc7-d74d-4c28-84b4-336630221922
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ 8f106a64-49a5-4333-865d-de2c22a3fba8


# ╔═╡ 215c05aa-4513-43ad-81c8-dd669f2c7a8f
oproot = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/"

# ╔═╡ 7a530655-c3c2-4d64-acda-cf5b58626e90
eos_name(ext) = joinpath(oproot, "TSO_M3D_$(ext)_v5.0/combined_eos_$(ext).hdf5")

# ╔═╡ 73426bbd-e476-4523-b5fc-0d0a80aadd7e
allEoS = Dict(
	0=>eos_name("magg_m0_a0_vmic1"),
	-1=>eos_name("magg_m1_a0_vmic1"),
	-2=>eos_name("magg_m2_a0_vmic1"),
	-3=>eos_name("magg_m3_a0_vmic1"),
	-4=>eos_name("magg_m4_a0_vmic1")
)

# ╔═╡ a97947c4-31bc-4fe7-a3d5-57fad759ab79
closest_eos(feh) = begin
	k = sort(keys(allEoS) |> collect)
	ik = k[argmin(abs.(k .- feh))]
	allEoS[ik]
end

# ╔═╡ a8ee335f-cd4b-479f-8c22-b54ff3a3a1e3


# ╔═╡ 158f1700-a14d-4414-95eb-e4551fb148b8
begin
	function adiabatDown(start_point, eos_in::SqEoS; kwargs...)
	    # initial point
	    t0  = exp(first(start_point.lnT))
	    ρ0  = exp(first(start_point.lnρ))
	    g0  = exp10(start_point.logg)
	    z0  = first(start_point.z)
	
		adiabatDown(t0, ρ0, z0, start_point.logg, eos_in; kwargs...)
	end
	
	function adiabatDown(t0, ρ0, z0, logg, eos_in::SqEoS; 
	                    nz=2000, n_iter=30, dlnd=0.02, padding=0.05, ee_min=nothing, z_end=-Inf)
	    eos = @axed eos_in
	
	    # Find the initial energy from the EoS. Assert that the EoS is in the 
	    # correct units for this
	    @assert is_internal_energy(eos)
	    eemin,ee_max,_,_ = exp.(limits(eos))
	
	    ee_min = if isnothing(ee_min)
	        eemin
	    else
	        ee_min
	    end
	
	    # initial point
	    ee0 = exp(lookup(eos, :lnEi, log(ρ0), log(t0)))
	    g0  = exp10(logg)
	    p0  = exp(lookup(eos, :lnPg, log(ρ0), log(ee0)))
	
	    # storage arrays
	    z  = zeros(nz)
	    d  = zeros(nz)
	    ee = zeros(nz)
	    t  = zeros(nz)
	    p  = zeros(nz)
	
	    # integrate from the starting point upwards, until z is larger than ztop
	    # or the maximum number of points is reached
	    dee    = 0
	    z[nz]  = z0
	    d[nz]  = ρ0  
	    ee[nz] = ee0 
	    t[nz]  = t0  
	    p[nz]  = p0  
	
	    # downwards
	    i_start=1
	    for i in nz-1:-1:1
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

			# check if there is an end 
			if z[i] <= z_end
				i_start = i
				break
			end
	    end

		z = z[i_start:end]
		d = d[i_start:end]
		t = t[i_start:end]
		ee = ee[i_start:end]
	    
	    TSO.Model1D(z=z, lnρ=log.(d), lnT=log.(t), lnEi=log.(ee), logg=logg)
	end
end

# ╔═╡ 6e257449-5cac-4f00-9303-6a25219f9c97
function stitchedMARCS(marcs_model, logg, eos_E; z_end=-Inf, resolution=nothing)
	i_top_marcs = (
		last(marcs_model.structure["T"]),
		last(marcs_model.structure["Density"]),
		-last(marcs_model.structure["Depth"]),
		logg
	)
	ad = adiabatDown(i_top_marcs..., eos_E, z_end=z_end)
	ad_opt = @optical TSO.flip(ad, depth=true) eos

	ee_marcs = lookup(
		eos_E, :lnEi, 
		log.(marcs_model.structure["Density"]), 
		log.(marcs_model.structure["T"])
	)

	newModel = TSO.Model1D(
		z=[marcs_model.structure["Depth"]..., ad_opt.z...], 
		lnρ=[log.(marcs_model.structure["Density"])..., ad_opt.lnρ...], 
		lnT=[log.(marcs_model.structure["T"])..., ad_opt.lnT...], 
		lnEi=[ee_marcs..., ad_opt.lnEi...], 
		logg=logg
	)
	
	newModel = if !isnothing(resolution)
		uniform_z = range(
			minimum(newModel.z), maximum(newModel.z), step=abs(resolution)
		) |> collect
		TSO.flip(TSO.interpolate_to(newModel, in_log=false, z=uniform_z), depth=true)
	else
		newModel
	end

	@optical newModel eos_E
end

# ╔═╡ 773f0cb9-10bd-4440-bdb2-dff6c969cf65


# ╔═╡ 4f1b436a-21bb-49df-b84a-750ef9e697f9
md"""
# Stitched MARCS grid
For testing purposes we create a Grid similar to the Stagger grid that contains stitched MARCS models. We utilise the Stagger grid to find the size of the simulation domain, velocities and resolution.
"""

# ╔═╡ c5325ecd-6290-443f-9384-84e3371d89a8
marcs_path = "/mnt/beegfs/gemini/groups/bergemann/users/shared-storage/gerber/marcs/marcs_standard_comp"  

# ╔═╡ 5e253ecd-ddca-4ab4-ac26-d42080b1005e
marcs_parameters_from_name(marcs_path) = begin
	name = last(split(marcs_path, "/", keepempty=false))
	name = first(split(name, ".mod"))
	name = split(name, "_", keepempty=false)
	T = parse(Float64, name[1][2:end])
	g = parse(Float64, name[2][2:end])
	feh = parse(Float64, name[6][2:end])
	turb = parse(Float64, name[4][2:end])
	alpha = parse(Float64, name[7][2:end])
	
	T, g, feh, alpha, turb
end

# ╔═╡ 790dcf69-24c7-4b35-95e2-6d5f9e6ad02a
list_marcs_models(folder) = begin
	allfiles = MUST.glob("p*", folder)
	allParameters = marcs_parameters_from_name.(allfiles)
	MUST.DataFrame(
		MUST.Dict(
			:path=>allfiles,
			:T=>getindex.(allParameters, 1),
			:logg=>getindex.(allParameters, 2),
			:feh=>getindex.(allParameters, 3),
			:alpha=>getindex.(allParameters, 4), 
			:turb=>getindex.(allParameters, 5)
		)
	)
end

# ╔═╡ 08fd37b2-2d68-4ce5-afce-a92cf6d89142
df = list_marcs_models(marcs_path)   

# ╔═╡ 23814821-fae0-402e-b8c9-696d8c9b01b2
row_selector(g) = begin
	T, logg, feh = g[1, :T], g[1, :logg], g[1, :feh]

	ipick = if MUST.nrow(g) > 1
		alphas = sort(g[!, :alpha])
		alpha = if feh < -1
			last(alphas)
		else
			first(alphas)
		end

		iT = minimum(g[g[!, :turb] .>= 1.0, :turb])

		mask = (g[!, :turb].==iT) .&& (g[!, :alpha].==alpha)
		findfirst(mask)
	else
		1
	end

	g[ipick, :]
end

# ╔═╡ 50060c7b-f021-424c-8733-ecf7db170e95
marcs_grid = combine(groupby(df, [:T, :logg, :feh]), row_selector)  

# ╔═╡ 3c2ee7a4-f85f-4505-ae48-248e9ec11fac
deleteat!(marcs_grid, (marcs_grid[!, :feh] .> 0.0) .| (marcs_grid[!, :feh] .< -4.0))

# ╔═╡ 549bbdfb-834d-4043-aecf-ef33391ef5f5
deleteat!(marcs_grid, (marcs_grid[!, :T] .> 7000.0) .| (marcs_grid[!, :T] .< 3500.0))

# ╔═╡ 646a5bfb-0fc8-4256-894b-fffc3f10525c


# ╔═╡ 7e63e0d2-7baa-4523-8340-3c709cab3366
md"Now we have a MARCS grid, that contains models with the lowest available microturbulance (greater than 0), and α-enhancement of 0.4 below metallicity of -1. This makes sure that the grid is unique in Teff, logg and FeH. We can now grab auxiliary 3D information from the Stagger grid. Those MARCS models can then be used to interpolate new models. The new models should then be adiabatically extrapolated."

# ╔═╡ 061fcfbd-45b6-4f6b-b191-64b672ac4a08
begin
	folder_marcs = "MARCS"
	if !isdir(folder_marcs) 
		mkdir(folder_marcs)
		mkdir(joinpath(folder_marcs, "av_models"))
	end
	if !isdir(joinpath(folder_marcs, "av_models"))
		mkdir(joinpath(folder_marcs, "av_models"))
	end
end

# ╔═╡ a9f10b54-4c53-497a-9586-0ac48c430276
begin
	parasmarcs = zeros(MUST.nrow(marcs_grid), 3)
	parasmarcs[:, 1] = marcs_grid[!, :T]
	parasmarcs[:, 2] = marcs_grid[!, :logg]
	parasmarcs[:, 3] = marcs_grid[!, :feh]
	folder = joinpath(folder_marcs, "av_models")
	old_path = marcs_grid[!, :path]

	eos_marcs_path = [closest_eos(feh) for feh in parasmarcs[:, 3]]
	eos_marcs = [
		reload(SqEoS, eos_marcs_path[i])
		for i in 1:size(parasmarcs, 1)
	]
end

# ╔═╡ 450f1d97-5f7f-41e8-9445-3b435e67e4ab
function simpleMarcs(oldpath, teff, logg, feh; folder, eos)
	tname = MUST.@sprintf "%.2f" teff/100
	gname = MUST.@sprintf "%.2f" logg*10
	fname = MUST.@sprintf "%.3f" feh
	name = "t$(tname)g$(gname)m$(fname)"
	av_path = abspath(joinpath(folder, "$(name)_99999_av.dat"))
	avo_path = abspath(joinpath(folder, "$(name)_99999_avo.dat"))
	
	model = MUST.MARCSModel(oldpath)
	model = TSO.Model1D(
		z=model.structure["Depth"], 
		lnρ=log.(model.structure["Density"]), 
		lnT=log.(model.structure["T"]), 
		lnEi=fill!(similar(model.structure["Depth"]), 0.0),
		logg=logg
	)
	hres = minimum(abs.(diff(model.z)))
	model = @optical model eos
	znew = TSO.rosseland_depth(eos, model)
	model.z .= znew
	TSO.optical_surface!(model)

	model = TSO.flip(TSO.interpolate_to(
		model, in_log=false,
		z=range(minimum(model.z), maximum(model.z), step=hres)
	))
	

	open(av_path, "w") do f
		MUST.writedlm(f, [model.z exp.(model.lnT) model.lnρ])
	end
	open(avo_path, "w") do f
		MUST.writedlm(f, [model.z exp.(model.lnT) model.lnρ log.(model.τ)])
	end

	model
end

# ╔═╡ ebd655db-c8aa-4fb2-9b9c-5e9ea9ffb391
function interpolate_grid(grid, teff, logg, feh; folder, oldpath)
	folder_name = "interpolated"
	mesh   = "interpolated"

	tname = MUST.@sprintf "%.2f" teff/100
	gname = MUST.@sprintf "%.2f" logg*10
	fname = MUST.@sprintf "%.3f" feh
	name = "t$(tname)g$(gname)m$(fname)"
	av_path = abspath(joinpath(folder, "$(name)_99999_av.dat"))

	snapshot = "interpolated"
	ma_x = MUST.interpolate_quantity(grid, "ma_x"; teff=teff, logg=logg, feh=feh)
	ma_y = MUST.interpolate_quantity(grid, "ma_y"; teff=teff, logg=logg, feh=feh)
	mi_z = MUST.interpolate_quantity(grid, "mi_z"; teff=teff, logg=logg, feh=feh)
	ma_z = MUST.interpolate_quantity(grid, "ma_z"; teff=teff, logg=logg, feh=feh)
	vmin = MUST.interpolate_quantity(grid, "vmin"; teff=teff, logg=logg, feh=feh)
	vmax = MUST.interpolate_quantity(grid, "vmax"; teff=teff, logg=logg, feh=feh)
	
	tscale = round(
		MUST.interpolate_quantity(grid, "tscale"; teff=teff, logg=logg, feh=feh),
		sigdigits=1
	)
	
	hres = MUST.interpolate_quantity(grid, "hres"; teff=teff, logg=logg, feh=feh)

	simpleMarcs(oldpath, teff, logg, feh; folder=folder)
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
			"av_path"  => joinpath(folder, "$(name)_99999_av.dat"),
			"abs_av_path" => av_path,
			"teff"     => teff,
			"logg"     => logg,
			"feh"      => feh
		)
	)
end

# ╔═╡ 38f1be62-784f-4efd-9a6b-4219aaac1b0f
function adiabatic_extrapolation(model, eos, dz=nothing; kwargs...)
    model_flip = TSO.flip(model)
    i_top = TSO.pick_point(model_flip, 1)

    ad_opt = if isnothing(dz)
        a = TSO.adiabatDown(i_top, eos; kwargs...)
        @optical a eos
    else
        z_end = first(model_flip) - dz
        a = TSO.adiabatDown(i_top, eos, z_end=z_end; kwargs...)
        @optical a eos
    end

	TSO.flip!(ad_opt, depth=false)
	model_flip.lnEi .= lookup(eos, :lnEi, model_flip.lnρ, model_flip.lnT)
	
    model = TSO.flip(
        TSO.Model1D(
            z=[ad_opt.z..., model_flip.z...], 
            lnρ=[ad_opt.lnρ..., model_flip.lnρ...], 
            lnT=[ad_opt.lnT..., model_flip.lnT...], 
            lnEi=[ad_opt.lnEi..., model_flip.lnEi...], 
            logg=model_flip.logg
        ), 
        depth=true
    )
	uniform_z = range(minimum(model.z), maximum(model.z), length=length(model.z)) 
	TSO.flip(
		TSO.interpolate_to(model, in_log=false, z=uniform_z |> collect),
		depth=true
	)
end

# ╔═╡ 07ac5d26-c496-40e1-8e8a-573abe05af10
function construct_grid(teff, logg, feh; folder, oldpath, eos)
	@info teff,logg,feh
	folder_name = "constructed"
	mesh   = "constructed"

	tname = MUST.@sprintf "%.2f" teff/100
	gname = MUST.@sprintf "%.2f" logg*10
	fname = MUST.@sprintf "%.3f" feh
	name = "t$(tname)g$(gname)m$(fname)"
	av_path = abspath(joinpath(folder, "$(name)_99999_av.dat"))

	snapshot = "interpolated"
	model = simpleMarcs(oldpath, teff, logg, feh; folder=folder, eos=eos)
	hres = minimum(abs.(diff(model.z)))
	
	model_extended = adiabatic_extrapolation(model, eos, nz=1000, dlnd=0.02)
	model_extended = @optical model_extended eos

	# Interpolate the extended model to be until logτ=7, find the optical surface and construct a z scale
	τ_min = max(minimum(log10.(model_extended.τ)), -8)
	τ_max = min(maximum(log10.(model_extended.τ)), 7)
	@info τ_min,τ_max
	model_clipped = TSO.flip!(TSO.interpolate_to(
		model_extended, 
		in_log=true, 
		τ=range(τ_min, τ_max, length=length(model_extended.z))
	), depth=true)

	# recompute z scale from opacity
	znew = TSO.rosseland_depth(eos, model_clipped)
	model_clipped.z .= znew

	# find the optical surface
	TSO.optical_surface!(model_clipped)
	TSO.flip!(model_clipped)

	mi_z = minimum(model_clipped.z)
	ma_z = maximum(model_clipped.z)
	ma_x = 2.0*abs(ma_z - mi_z)
	ma_y = 2.0*abs(ma_z - mi_z)

	optical_surface = TSO.interpolate_to(
		model_clipped, 
		in_log=true, 
		τ=[log10(2.0/3)]
	)
	E = is_internal_energy(@axed(eos)) ? optical_surface.lnEi : optical_surface.lnT
	pg = exp.(lookup(eos, :lnPg, optical_surface.lnρ, E))
	c = sqrt.(1.2 .* pg ./ exp.(optical_surface.lnρ)) / 7.0 |> first
	Hp = pg ./ (exp.(optical_surface.lnρ) .* exp10(logg)) |> first

	timescale = Hp / c
	lengthscale = abs(ma_z - mi_z) / 3.0
	velocityscale = lengthscale / timescale

	model_clipped
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
			"vmin"     => velocityscale,
			"vmax"     => velocityscale,
			"tscale"   => timescale,
			"hres"     => hres,
			"av_path"  => joinpath(folder, "$(name)_99999_av.dat"),
			"avo_path"  => joinpath(folder, "$(name)_99999_avo.dat"),
			"abs_av_path" => av_path,
			"teff"     => teff,
			"logg"     => logg,
			"feh"      => feh
		)
	)
end

# ╔═╡ 0b99d9d8-07d1-4884-b335-9051baa96cf4
marcs_construct = [
	construct_grid(
		parasmarcs[i, 1], parasmarcs[i, 2], parasmarcs[i, 3], 
		folder=folder, 
		oldpath=marcs_grid[i, :path],
		eos=eos_marcs[i]
	)
	for i in axes(parasmarcs, 1)
]

# ╔═╡ 2c40e250-b0c4-4a34-b626-b052e1143b1e
begin
	marcs_constructed_grid = marcs_construct[1]
	for (i, m) in enumerate(marcs_construct)
		i==1 && continue
		append!(marcs_constructed_grid, m)
	end
end

# ╔═╡ d0f927ab-b49a-4a51-ab15-f089caffb1e9
marcs_constructed_grid

# ╔═╡ 62597df8-0ffa-4a63-9d79-a5111535525f
function relative_path(from, to)
	fstag = basename(from)
	feos = basename(to)
	stag_path = split(dirname(from), "/", keepempty=false)
	eos_path = split(dirname(to), "/", keepempty=false)

	i = findfirst(x->!(x in eos_path), stag_path)
	ieos = findfirst(eos_path .== stag_path[i-1]) + 1

	path_difference = stag_path[i:end]
	eos_path_difference = eos_path[ieos:end]
	new_path = joinpath([".." for _ in path_difference]..., eos_path_difference..., feos)
end

# ╔═╡ d1998950-051d-4b26-b219-0424f92b4778
marcs_constructed_grid[!, "matching_eos"] = relative_path.(
	abspath(joinpath(folder_marcs, "marcs_grid.mgrid")),
	eos_marcs_path
)

# ╔═╡ 44fa2659-7369-42fa-a52c-996808656e9f
#marcs_constructed_grid[!, "avo_path"] = marcs_constructed_grid[!, "av_path"]

# ╔═╡ fe0512ee-38e7-4e75-aa73-e517c522aef7
MUST.save(
	MUST.Atmos1DGrid(
		"MARCS constructed", 
		joinpath(folder_marcs,"marcs_grid.mgrid"),
		marcs_constructed_grid
	)
)

# ╔═╡ Cell order:
# ╠═b04a8d9a-5307-11ef-1268-95d5b8233977
# ╟─273c40e1-f4ef-4f24-93dc-b9bb95f3abeb
# ╟─77e9b680-d0e4-431f-9a59-a4a081ffd65e
# ╠═63ac3293-0d3f-4533-8fc4-ae3695e40cc6
# ╠═62d17d2d-6727-4632-b086-30cc6e828360
# ╠═83648aa8-11dd-425b-bc92-938a360e03a0
# ╠═ce43f8fe-6274-4da3-b427-1f362352f1a4
# ╠═e607ddf1-fb54-4c6c-842c-40003f2604e2
# ╠═fc9633c5-fbd8-458c-a7de-38875987eb3b
# ╠═4fe7acc7-d74d-4c28-84b4-336630221922
# ╟─8f106a64-49a5-4333-865d-de2c22a3fba8
# ╠═215c05aa-4513-43ad-81c8-dd669f2c7a8f
# ╠═7a530655-c3c2-4d64-acda-cf5b58626e90
# ╠═73426bbd-e476-4523-b5fc-0d0a80aadd7e
# ╠═a97947c4-31bc-4fe7-a3d5-57fad759ab79
# ╟─a8ee335f-cd4b-479f-8c22-b54ff3a3a1e3
# ╟─158f1700-a14d-4414-95eb-e4551fb148b8
# ╟─6e257449-5cac-4f00-9303-6a25219f9c97
# ╟─773f0cb9-10bd-4440-bdb2-dff6c969cf65
# ╟─4f1b436a-21bb-49df-b84a-750ef9e697f9
# ╠═c5325ecd-6290-443f-9384-84e3371d89a8
# ╟─5e253ecd-ddca-4ab4-ac26-d42080b1005e
# ╟─790dcf69-24c7-4b35-95e2-6d5f9e6ad02a
# ╠═08fd37b2-2d68-4ce5-afce-a92cf6d89142
# ╟─23814821-fae0-402e-b8c9-696d8c9b01b2
# ╠═50060c7b-f021-424c-8733-ecf7db170e95
# ╠═3c2ee7a4-f85f-4505-ae48-248e9ec11fac
# ╠═549bbdfb-834d-4043-aecf-ef33391ef5f5
# ╟─646a5bfb-0fc8-4256-894b-fffc3f10525c
# ╟─7e63e0d2-7baa-4523-8340-3c709cab3366
# ╠═061fcfbd-45b6-4f6b-b191-64b672ac4a08
# ╠═a9f10b54-4c53-497a-9586-0ac48c430276
# ╟─450f1d97-5f7f-41e8-9445-3b435e67e4ab
# ╟─ebd655db-c8aa-4fb2-9b9c-5e9ea9ffb391
# ╟─38f1be62-784f-4efd-9a6b-4219aaac1b0f
# ╟─07ac5d26-c496-40e1-8e8a-573abe05af10
# ╠═0b99d9d8-07d1-4884-b335-9051baa96cf4
# ╠═2c40e250-b0c4-4a34-b626-b052e1143b1e
# ╠═d0f927ab-b49a-4a51-ab15-f089caffb1e9
# ╠═62597df8-0ffa-4a63-9d79-a5111535525f
# ╠═d1998950-051d-4b26-b219-0424f92b4778
# ╠═44fa2659-7369-42fa-a52c-996808656e9f
# ╠═fe0512ee-38e7-4e75-aa73-e517c522aef7
