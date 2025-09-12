### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1c7ee45a-8d60-11f0-29f2-abf8c4f707d2
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("..")
	using TSO
	using MUST
	using PythonPlot
	using LaTeXStrings 
	using PlutoUI
	using KernelDensity
	using PlutoUI: combine
end

# ╔═╡ ef8e7bda-8c97-418d-864f-2e7e04c6cca9
begin
	using GaussianMixtures
	using Distributions
	using LinearAlgebra
end

# ╔═╡ 29b1890e-7361-40cc-a1f2-3b874faba043
md"# Code Setup"

# ╔═╡ f1473678-24b0-4a97-967a-4f6070de46d7
begin
	@import_dispatch "../../../../dispatch2"
	@import_m3dis "../../../../Multi3D"
end

# ╔═╡ f7e7ae1d-49dd-4c76-a9ff-cec1013c86aa
begin
	auxfuncs = MUST.pyimport("m3dis.auxfuncs")
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Hoppe2024.mplstyle"))
	matplotlib.rcParams.update(Dict(
		"font.size"=> 14,
		"axes.linewidth"=> 2,
		"xtick.major.width"=> 2,
		"xtick.minor.width"=> 1.8,
		"ytick.major.width"=> 2,
		"ytick.minor.width"=> 1.8,
		"xtick.major.size"=> 6.5,
		"ytick.major.size"=> 6.5,
		"xtick.minor.size"=> 4,
		"ytick.minor.size"=> 4,
		
	))
	
	availableRuns(path; box_name="box") = begin
		box_name = box_name=="box" ? nothing : box_name
		runs = MUST.glob("*/", path)
		runs = if isnothing(box_name)
			runs
		else
			[r for r in runs if (length(MUST.glob("$(box_name)*.hdf5", r)) > 0)]
		end
		last.(split.(runs, "/", keepempty=false))
	end

	function select_for_each(runs, entries, defaults)
		return combine() do Child
			namesChild = [
				Child(runs[i], Select(entries[runs[i]], default=defaults[runs[i]]))
				for i in eachindex(runs)
			]
			inputs = [
				md""" $(name): $(
					namesChild[i]
				)"""
				for (i, name) in enumerate(runs)
			]
			
			md"""
			## Selection for each model
			$(inputs)
			"""
		end
	end

	pretty_from_name(name) = begin
		mi = MUST.ModelInformation(name)
		L"\rm T_{eff}="*"$(round(Int, mi.teff))"*L"\rm\ K,\ log(g)="*"$(round( mi.logg, digits=2))"*L"\rm,\ [Fe/H] ="*"$(round(mi.feh, digits=2))"
	end
	
	pretty_from_name_short(name) = begin
		mi = MUST.ModelInformation(name)
		"$(round(Int, mi.teff)) K, "*"$(round( mi.logg, digits=2)), "*"$(round(mi.feh, digits=2))"
	end
	pretty_from_name_no_feh(name) = begin
		mi = MUST.ModelInformation(name)
		L"\rm T_{eff}="*"$(round(Int, mi.teff))"*L"\rm\ K,\ log(g)="*"$(round( mi.logg, digits=2))"
	end
	TableOfContents(title="CEMP paper")
	
end

# ╔═╡ f7fa9c00-1bc2-4246-be6a-ac30ee3d6457
md"# Model Selection"

# ╔═╡ cdedf1ac-7de7-4e61-83f7-03282f747445
md"## 3D models"

# ╔═╡ a7c4ab87-e3eb-456c-a41b-0a5d72501da4
datadir = @in_dispatch "CEMP_models3/"

# ╔═╡ ed0947be-f0bd-4054-86f7-afd8f36641bf
snapshot_id = -1

# ╔═╡ eff9f59a-faaf-4f3c-ab1d-67c2bee0bdc0
begin
	models = String.(availableRuns(datadir))
	models_info = Dict(
		String(m) => MUST.ModelInformation(m) for m in models
	)
	cemp_models = [String(m) for m in models if MUST.ModelInformation(m).category=="CEMP"]
	scaledSolar_models = [String(m) for m in models if MUST.ModelInformation(m).category=="ScaledSolar"]
end

# ╔═╡ 9b27ec05-9e7d-40e8-b004-62d45a0e8f25


# ╔═╡ ff4b6af5-11c5-48db-95f5-48dc055e12dc
md"## MARCS models"

# ╔═╡ 895721c1-dc9d-47d4-b5fb-4e189b5184c0
md"For each model, information about the spectra is stored in the models `m3dis` cubes. Each model also contains `MARCS` models that have spectra associated."

# ╔═╡ 16f38c33-b243-45d0-8bf3-2a0b388a8d39
begin
	listOfMARCS(name) = sort(MUST.glob("*.mod", joinpath(datadir, name)))
	MARCSNames(name) = split(name, '/', keepempty=false) |> last
	
	marcsmodelspaths = Dict()
	for m in models
		l = listOfMARCS(m)
		mm = []
		for (i, p) in enumerate(l)
			(!isfile(p)) && continue
			mo = try
				MUST.marcsBox(p)
			catch
				nothing
			end
	
			if !isnothing(mo)
				append!(mm, [p])
			end
		end
		marcsmodelspaths[m] = String.(mm)
	end

	mip_models = Dict(m=>String(basename.(first(MUST.glob("marcsIP*", joinpath(datadir, m))))) for m in models)
	marcsmodels = Dict(k=>[basename.(m)..., mip_models[k]] for (k, m) in marcsmodelspaths)
	
end

# ╔═╡ 550affe7-12a4-4b3f-8edb-4ff0c8c8a41b


# ╔═╡ 06b5cd98-c60f-4d61-a723-5d2b0ac6460a
md"For each model, we select a best matching MARCS model from the list. By default, the interpolated model is used that matches the effective temperature of the 3D model."

# ╔═╡ c8469b97-df26-4753-b5be-e0b2b8a02b82


# ╔═╡ 9a3a15f2-d167-47e5-85bf-327221203325


# ╔═╡ 4f163103-84ce-4e36-9b69-2b6d3e880ca1
md"## Utility functions"

# ╔═╡ 82ce1cc4-350c-4ca3-9c01-21321414a240
begin
	get_snapshot(snapid::Int, folder, tau_name="tau500") = begin
		pick_snapshot(folder, snapid, box_name="box", tau_name=tau_name)
	end
	get_snapshot(snapid::String, folder, tau_name="tau500") = begin
		@assert tau_name=="tau500"
		b,b2 = pick_snapshot(folder, snapid, box_name="box_m3dis")

		if ndims(b[:τ500]) == 1
			b[:τ500] = reshape(b[:τ500], 1, 1, length(b[:τ500]))
		end
		b, b2
	end
end

# ╔═╡ c3f55a57-397d-4b2c-87f8-0ebc133a9ab0
get_spectra(snapid, folder) = begin
	snapid_m3dis = MUST.snapshotid(folder |> MUST.converted_snapshots, snapid)
	pick_snapshot(folder, snapid_m3dis, box_name="box_m3dis") |> first
end

# ╔═╡ 3c6ef3eb-5ef7-43b3-946f-73049a7d9515
get_m3d_spectra(snapid, folder; extension="lam_4297-4303_contr") = begin
	i = if typeof(snapid) <: Int && (snapid<0)
		"m3dis_$(MUST.list_snapshots(MUST.converted_snapshots(folder), numbered_only=true)[end+snapid+1])"
	elseif typeof(snapid) <: Int
		"m3dis_$(snapid)"
	else
		snapid
	end
	try
		MUST.M3DISRun(joinpath(folder, "spectra_$(i)_$(extension)"))
	catch
		@show joinpath(folder, "spectra_$(i)_lam_4295-4325_GBand_DC_contr")
		MUST.M3DISRun(joinpath(folder, "spectra_$(i)_lam_4295-4325_GBand_DC_contr"))
	end
end

# ╔═╡ c508746d-879e-4e76-baed-136fa36d716c
same_scaled_solar(cemp_model_name) = scaledSolar_models[
	findfirst(
		MUST.same_parameters(
			MUST.ModelInformation(a), 
			MUST.ModelInformation(cemp_model_name),
		) for a in scaledSolar_models
	)
]

# ╔═╡ 2fa849f5-bf68-4548-9d0e-234ba6edf104


# ╔═╡ d5e20abd-a153-4dc6-9895-d8142fcc6acc
md"# Structure comparison"

# ╔═╡ ba0a854b-5fa4-468d-91c9-380523e385cd
"""
	plot_profile(s, q; ax, kwargs...)

Plot profile function of optical depth cube (500).
"""
plot_profile(b, bt, q; xvar=:τ500, ax=nothing,cmap="GnBu", alpha=1.0, bins=200, kwargs...) = begin
	qs, logq = MUST.is_log(q)

	xdat = if xvar==:τ500
		log10.(b[:τ500])
	else
		b[xvar]
	end
	
	# Background distribution
	h, x, y = MUST.numpy.histogram2d(
		reshape(xdat, :), reshape(logq.(b[qs]), :), bins=bins
	)
	if !isnothing(ax)
		ax.imshow(
			h.T, 
			origin="lower",
			interpolation = "bicubic", 
			extent=[minimum(x), maximum(x), minimum(y), maximum(y)],
			cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=1), aspect="auto",
			rasterized=true,alpha=alpha
		)
	end
	
	#ax.hist2d(reshape(log10.(b[:τ500]), :), reshape(b[:T], :), rasterized=true, cmap="YlGnBu", bins=300, norm=matplotlib.colors.LogNorm())

	x, y = profile(MUST.mean, bt, :log10τ500, q)
end

# ╔═╡ 6d2305b8-b452-47c9-a0a2-b097029ade80


# ╔═╡ f4e30299-95e0-4861-919e-6a767adf1194
md"## 3D vs. 1D"

# ╔═╡ cfb5b8d3-4a10-44ba-9cc7-a566339a833f
md"""
3D model: $(@bind structure_3D_select Select(models, default=first(models)))
"""

# ╔═╡ df166dde-2837-4440-a8dc-0c2a8194866e
md"""
1D model: $(@bind structure_1D_select Select(marcsmodels[structure_3D_select], default=marcsmodels[structure_3D_select][2])))
"""

# ╔═╡ 9834da27-c370-4795-a91e-89a382f7b81d
begin
	structure_3D1D_figsave_txt = "$(structure_3D_select)_structure"
	md"Save figure at: $(@bind structure_3D1D_figsave confirm(TextField(length(structure_3D1D_figsave_txt)+1, default=structure_3D1D_figsave_txt)))"
end

# ╔═╡ 53aa57c6-2b5c-4e04-94ab-c732ef7f2ec3
md"""
Formation height to overplot (3D): $(@bind formation_height_3D confirm(TextField(15)))
"""

# ╔═╡ 64ee5874-29e7-49a3-b8b3-5f89c7156bb9
md"""
Formation height to overplot (1D): $(@bind formation_height_1D confirm(TextField(15)))
"""

# ╔═╡ 73ebcfac-d1b0-496f-8f2d-df5ea3aa7830


# ╔═╡ 7fb3939a-1709-4a77-be21-27924e2072cb
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))
	
	# tau500 snapshot
	b, bt = get_snapshot(snapshot_id, joinpath(datadir, structure_3D_select))
	plot_profile(b, bt, :T, ax=ax)
	#ax.plot(profile(mean, bt, :log10τ500, :T)..., color="steelblue", lw=5)
	
	# 1D model
	bmarcs,_ = get_snapshot(structure_1D_select, joinpath(datadir, structure_3D_select))
	@show size(bmarcs[:log10τ500]) size(bmarcs[:T])
	ax.plot(profile(mean, bmarcs, :log10τ500, :T)..., color="tomato", lw=5)


	# inset
	left, bottom, width, height = [0.15, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	vmin = 4100
	vmax = 7600
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal",
		vmin=vmin,
		vmax=vmax
	)
	#ax2.set_xlim(-2.5, 2.5)
	#ax2.set_ylim(-2.5, 2.5)
	@show minimum(Tplane) maximum(Tplane)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")
	ax2.set_xlabel(L"\rm x\ [Mm]")

	llim = 2200
	# formation height bars
	if length(formation_height_3D) > 0
		ax.hlines(llim+700, parse(Float64, formation_height_3D), 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
		ax.text(0.1, llim+700, "3D", ha="left", va="center", color="steelblue")
	end
	if length(formation_height_1D) > 0
		ax.hlines(llim+1050, parse(Float64, formation_height_1D), 0.0, ls="-", color="tomato", alpha=0.8, lw=6)
		ax.text(0.1, llim+1050, "1D", ha="left", va="center", color="tomato")
	end	

	if (length(formation_height_3D) > 0) & (length(formation_height_1D) > 0)
		xpos = 0.1
		ax.text(xpos, llim+500, L"\mathrm{line\ formation\ depth\ of\ CH\ G-band}", ha="right", va="top", color="k", alpha=1, fontsize="x-small")
	end
	
	ax.set_xlabel(L"\rm optical\ depth\ [\log_{10} \tau_{500}]")
	ax.set_ylabel(L"\rm temperature\ [K]")

	ax.set_xlim(-4.2, 1.7)
	ax.set_ylim(llim, 11000)

	ax.set_title(pretty_from_name(structure_3D_select))

	@info "Model resolution (Nx, Ny, Nz): $(size(b))"
	@info "Model resolution (km): $((maximum(b.z)-minimum(b.z)) / size(b, 3) ./1e5)"

	f.savefig(structure_3D1D_figsave_txt*"_3D1D.pdf")
	f
end

# ╔═╡ 7188fdf5-2366-47d6-b082-b476fc30ebf3
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 6), layout="tight")
	
	# tau500 snapshot
	b, bt = get_snapshot(snapshot_id, joinpath(datadir, structure_3D_select))
	plot_profile(b, bt, :T, ax=ax)
	
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	vmin = 4100
	vmax = 6600
	im = ax.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal",
		vmin=vmin,
		vmax=vmax
	)
	
	cbar = f.colorbar(im, fraction=0.046, pad=0.04)
	ax.set_title(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm x\ [Mm]")
	
	ax.set_xlabel(L"\rm X\ [Mm]")
	ax.set_ylabel(L"\rm Y\ [Mm]")

	ax.xaxis.set_major_locator(plt.MaxNLocator(7))
	ax.yaxis.set_major_locator(plt.MaxNLocator(7))

	ax.set_title(pretty_from_name(structure_3D_select))

	@info "Model resolution (Nx, Ny, Nz): $(size(b))"
	@info "Model resolution (km): $((maximum(b.z)-minimum(b.z)) / size(b, 3) ./1e5)"

	f.savefig(structure_3D1D_figsave_txt*"_opticalsurface.pdf")
	f
end

# ╔═╡ 6d6772d6-f956-4216-8731-de2fc4b32f3e


# ╔═╡ 6174ea3c-c2e5-4c24-8299-0c7e74afa2ef
md"## Two Model comparison"

# ╔═╡ 68e0ac58-ed87-442d-a7b5-13053c9c1cdf
md"""
3D model (A): $(@bind structure_3D_select_A Select(cemp_models, default=first(cemp_models)))
"""

# ╔═╡ 52a7d86a-a6e0-4b9e-8f63-adc0cb9a0359
md"""
1D model (A): $(@bind structure_1D_select_A Select(marcsmodels[structure_3D_select_A], default=marcsmodels[structure_3D_select_A][2])))
"""

# ╔═╡ e4c53afc-f220-420e-8f61-adcca2594767
md"""
3D model (B): $(@bind structure_3D_select_B Select(cemp_models, default=last(cemp_models)))
"""

# ╔═╡ ae363cc0-a747-4e26-808b-97df38543144
md"""
1D model (B): $(@bind structure_1D_select_B Select(marcsmodels[structure_3D_select_B], default=marcsmodels[structure_3D_select_B][2]))
"""

# ╔═╡ 2275bda1-cd48-42c7-9462-89aefd5e3f8c
begin
	structure_3D3D_figsave_txt = "$(structure_3D_select_A)_vs_$(structure_3D_select_B)"
	md"Save figure at: $(@bind structure_3D3D_figsave confirm(TextField(length(structure_3D3D_figsave_txt)+1, default=structure_3D3D_figsave_txt)))"
end

# ╔═╡ bf6c434a-0fd9-420b-ac06-bb8e215e9fdd
begin
	test_model_marcs = MUST.readdlm("../marcs_cemp_5800_4.5_-5.txt", skipstart=2)
end;

# ╔═╡ 951cf3cb-8f2f-4952-96fc-0318e866724e
let
	plt.close()
	f, ax = plt.subplots(2, 2, figsize=(8, 11), sharey=true)
	plt.subplots_adjust(wspace=0.02, hspace=0.15)
	#plt.subplots_adjust(wspace=0.02, hspace=0.02)
	
	xlim_tau = [-4.1, 0.2]

	plot_type = :all
	fig_label = if plot_type == :marcs
		 "_marcs.pdf"
	elseif plot_type == :carbon
		"_carbon.pdf"
	elseif  plot_type == :carbonmarcs
		"_carbonmarcs.pdf"
	elseif plot_type == :all
		"_all.pdf"
	end

	# overplot the example model
	minfo_B = MUST.ModelInformation(structure_3D_select_B)
	minfo_A = MUST.ModelInformation(structure_3D_select_A)
	ax_index = if (minfo_A.teff==5750.0) && (minfo_A.logg==4.5) && (minfo_A.feh==-5.0)
		0
	elseif (minfo_B.teff==5750.0) && (minfo_B.logg==4.5) && (minfo_B.feh==-5.0)
		1
	else
		nothing
	end

	upsample(x, y; N=30) = begin
		ma = sortperm(x)
		ip = MUST.linear_interpolation(x[ma], y[ma])
		x_new = range(minimum(x), maximum(x), length=N) |> collect
		x_new, ip.(x_new)
	end

	panel_data(b1, b2, bmarcs, x, y; average_on=:z) = begin
		_, y1 = profile(mean, b1, average_on, y)
		_, x1 = profile(mean, b1, average_on, x)
		z, tau = profile(mean, b1, average_on, :log10τ500)
		
		maskl = sortperm(tau)
		tau_ip = MUST.linear_interpolation(tau[maskl], x1[maskl], extrapolation_bc=MUST.Flat())
		xl = tau_ip.(xlim_tau)
		
		# tau500 snapshot 2
		_, y2 = profile(mean, b2, average_on, y)
		_, x2 = profile(mean, b2, average_on, x)
		
		# 1D model
		_, y3 = profile(mean, bmarcs, average_on, y)
		_, x3 = profile(mean, bmarcs, average_on, x)

		x1, y1, x2, y2, x3, y3, xl
	end
	
	plot_panel!(ax, d, T, d2, T2, d3, T3, xlim; add_test=nothing) = begin
		if (plot_type == :marcs) | (plot_type == :carbonmarcs) | (plot_type == :all) 
			mask = xlim[1] .< d3 .< xlim[2]
			ax.plot(d3[mask], T3[mask], color="tomato", lw=3, label=L"\rm 1D\  \text{scaled-solar}")

			if !isnothing(add_test)
				x_test = log10.(test_model_marcs[:, add_test])
				mask = xlim[1] .< x_test .< xlim[2]
				ax.plot(upsample(x_test[mask], test_model_marcs[mask, 3] .- 50)..., lw=1, color="tomato", marker="x", ls="", label="1D CEMP", markersize=7)
			else
				ax.plot([], [], lw=1, color="tomato", marker="x", ls="", label="1D CEMP", markersize=7)
			end
		end

		if (plot_type == :marcs) | (plot_type == :carbon) | (plot_type == :all) 
			mask = xlim[1] .< d2 .< xlim[2]
			ax.plot(upsample(d2[mask], T2[mask])..., color="steelblue", lw=3.5, label=L"\rm 3D\ \text{scaled-solar}", ls="-")
		end

		if (plot_type == :carbonmarcs) | (plot_type == :carbon) | (plot_type == :all) 
			mask = xlim[1] .< d .< xlim[2]
			ax.plot(upsample(d[mask], T[mask])..., color="steelblue", lw=1, label=L"\rm 3D\ CEMP", ls="", marker="x", markersize=7)
		end
	end
	
	# top left
	begin
		b, bt = get_snapshot(
			-1, joinpath(datadir, structure_3D_select_A)
		)
		b2, bt2 = get_snapshot(
			-1, joinpath(datadir, same_scaled_solar(structure_3D_select_A))
		)
		bmarcs,_ = get_snapshot(
			structure_1D_select_A, joinpath(datadir, structure_3D_select_A)
		)
		d1, T1, d2, T2, d3, T3, xlim1 = panel_data(
			b, b2, bmarcs, :log10d, :T, average_on=:z
		)
	end

	# top right
	begin
		b, bt = get_snapshot(
			-1, joinpath(datadir, structure_3D_select_B)
		)
		b2, bt2 = get_snapshot(
			-1, joinpath(datadir, same_scaled_solar(structure_3D_select_B))
		)
		bmarcs,_ = get_snapshot(
			structure_1D_select_B, joinpath(datadir, structure_3D_select_B)
		)
		d1B, T1B, d2B, T2B, d3B, T3B, xlim1B = panel_data(
			b, b2, bmarcs, :log10d, :T, average_on=:z
		)
	end

	# bottom left
	begin
		b, bt = get_snapshot(
			-1, joinpath(datadir, structure_3D_select_A)
		)
		b2, bt2 = get_snapshot(
			-1, joinpath(datadir, same_scaled_solar(structure_3D_select_A))
		)
		bmarcs,_ = get_snapshot(
			structure_1D_select_A, joinpath(datadir, structure_3D_select_A)
		)
		tau1, tT1, tau2, tT2, tau3, tT3, xlim2 = panel_data(
			bt, bt2, bmarcs, :log10τ500, :T, average_on=:log10τ500
		)
	end

	# bottom right
	begin
		b, bt = get_snapshot(
			-1, joinpath(datadir, structure_3D_select_B)
		)
		b2, bt2 = get_snapshot(
			-1, joinpath(datadir, same_scaled_solar(structure_3D_select_B))
		)
		bmarcs, _ = get_snapshot(
			structure_1D_select_B, joinpath(datadir, structure_3D_select_B)
		)
		tau1B, tT1B, tau2B, tT2B, tau3B, tT3B, xlim2B = panel_data(
			bt, bt2, bmarcs, :log10τ500, :T, average_on=:log10τ500
		)
	end

	# collect all panels and plot
	begin
		xlim_top = [min(xlim1[1], xlim1B[1]), max(xlim1[2], xlim1B[2])]
		xlim_bottom = [min(xlim2[1], xlim2B[1]), max(xlim2[2], xlim2B[2])]
		plot_panel!(ax[0,0], d1, T1, d2, T2, d3, T3, xlim1, add_test=ax_index==0 ? 5 : nothing)
		plot_panel!(ax[0,1], d1B, T1B, d2B, T2B, d3B, T3B, xlim1B, add_test=ax_index==1 ? 5 : nothing)
		plot_panel!(ax[1,0], tau1, tT1, tau2, tT2, tau3, tT3, xlim_bottom, add_test=ax_index==0 ? 2 : nothing)
		plot_panel!(ax[1,1], tau1B, tT1B, tau2B, tT2B, tau3B, tT3B, xlim_bottom, add_test=ax_index==1 ? 2 : nothing)

		ax[0,0].set_xlim(xlim1...)
		ax[0,1].set_xlim(xlim1B...)
		ax[1,0].set_xlim(xlim_bottom...)
		ax[1,1].set_xlim(xlim_bottom...)

		ax[0,0].set_xlabel(L"\rm density\ [g\ cm^{-3}]")
		ax[0,0].set_ylabel(L"\rm temperature\ [K]")

		
		ax[0,0].legend(loc="center", bbox_to_anchor=(0.4, 0.74), labelspacing=0.1)
		ax[1,0].legend(loc="center", bbox_to_anchor=(0.4, 0.74), labelspacing=0.1)
		

		ax[0,1].set_xlabel(L"\rm density\ [g\ cm^{-3}]")
		#ax[0,0].set_ylabel(L"\rm temperature\ [K]")
		#ax[0,0].legend()

		ax[1,0].set_xlabel(L"\rm optical\ depth\ [\log_{10} \tau_{500}]")
		ax[1,0].set_ylabel(L"\rm temperature\ [K]")

		ax[1,1].set_xlabel(L"\rm optical\ depth\ [\log_{10} \tau_{500}]")
		#ax[0,0].set_ylabel(L"\rm temperature\ [K]")
		#ax[0,0].legend()
	end

	#ax[0,1].tick_params(labelright=false, labelleft=false, labeltop=true, labelbottom=false)
	#ax[0,0].tick_params(labeltop=true, labelbottom=false)
	#ax[0,1].yaxis.set_label_position("right")
	#ax[0,0].xaxis.set_label_position("top")
	#ax[0,1].xaxis.set_label_position("top")

	ax[1,1].tick_params(labelright=false, labelleft=false)
	ax[1,1].yaxis.set_label_position("right")
	
	ax[0,0].xaxis.set_major_locator(plt.MaxNLocator(5))
	ax[0,1].xaxis.set_major_locator(plt.MaxNLocator(5))
	ax[1,0].xaxis.set_major_locator(plt.MaxNLocator(5))
	ax[1,1].xaxis.set_major_locator(plt.MaxNLocator(5))

	# parameter labels
	ax[0,0].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_A),
		ha="center",va="top", transform=ax[0,0].transAxes, 
		fontsize="medium", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)
	ax[0,1].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_B),
		ha="center",va="top", transform=ax[0,1].transAxes, 
		fontsize="medium", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)
	ax[1,0].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_A),
		ha="center",va="top", transform=ax[1,0].transAxes, 
		fontsize="medium", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)
	ax[1,1].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_B),
		ha="center",va="top", transform=ax[1,1].transAxes, 
		fontsize="medium", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)

	# panel labels
	ax[0,0].text(
		0.95,0.95,"A",
		ha="right",va="top", transform=ax[0,0].transAxes, fontweight="bold"
	)
	ax[0,1].text(
		0.95,0.95,"B",
		ha="right",va="top", transform=ax[0,1].transAxes, fontweight="bold"
	)
	ax[1,0].text(
		0.95,0.95,"C",
		ha="right",va="top", transform=ax[1,0].transAxes, fontweight="bold"
	)
	ax[1,1].text(
		0.95,0.95,"D",
		ha="right",va="top", transform=ax[1,1].transAxes, fontweight="bold"
	)

	f.savefig(structure_3D3D_figsave*fig_label, bbox_inches="tight")
	
	f
end

# ╔═╡ d883c7c9-ef6d-49cf-b1e3-46ba8b92f2ac


# ╔═╡ b0653e46-ae33-4725-98aa-443a2198c3a0
md"# Abundance Corrections"

# ╔═╡ e0d62ec9-04fc-4b9b-b215-a7160ba20205
begin
	@bind cemp_marcs_selected confirm(select_for_each(cemp_models, marcsmodels, mip_models))
end

# ╔═╡ 47b98cc4-0ecc-4d95-a257-f09076283b22


# ╔═╡ 75efbc2c-8be3-427d-bf62-16c9a64d2e6a
md"## Compute corrections"

# ╔═╡ bedda209-be73-49d4-8bc0-b1a902a1aac3
md"Specify the spectrum tag: $(@bind tag_sel confirm(TextField(35, default=\"GBand_DC\")))"

# ╔═╡ 79da8fae-65eb-4ddf-9826-93df48bb6a2f
md"Specify the element: $(@bind ele_sel confirm(TextField(10, default=\"[C/Fe]\")))"

# ╔═╡ 80b78da6-881b-4fed-9dd0-aa7d613b6d71
get_mean_spectra(run; onedimensional="") = begin
	time_average_spectra = MUST.glob("$(tag_sel)*.csv", joinpath(datadir, run))
	time_average_spectra_mask = [MUST.is_mean_spectrum(p, tag_sel, kind="flux", norm=true) for p in time_average_spectra]

	modelpath = time_average_spectra[time_average_spectra_mask]
	modelname = basename.(modelpath)
	
	# check if the model is actually 3D or 1D
	irest = [findfirst("mean_flux_norm.csv", mn) for mn in modelname]
	relevant_part = [mn[length(tag_sel)+2:first(irest[i])-2] for (i,mn) in enumerate(modelname)]
	accepted_model = if length(onedimensional)==0
		# we want a 3D model
		accs = Ref(1)
		for (i, mn) in enumerate(relevant_part)
			a = try 
				ss = split(mn, '-')
				parse(Float64, ss[1]), parse(Float64, ss[2])
				true
			catch
				false
			end
			if a
				accs[] = i
				break
			end
		end
		accs[]
	else
		findfirst(x->occursin(onedimensional, x), modelname)
	end
	#@show accepted_model modelpath[accepted_model]
	isnothing(accepted_model) ? nothing : MUST.MeanSpectra(modelpath[accepted_model])
end

# ╔═╡ 7c2b57e4-1c06-4adc-9017-a176fd36fe82
function abundance_correction_time(run, reference_abundance; λs=nothing, λe=nothing)
	corrections = Dict()
	eqws = Dict()
	abs = Dict()

	# read the time averaged spectra for this run with this tag
	sp3D = get_mean_spectra(run)
	nanmask = [all(.!isnan.(s.spectrum)) for s in sp3D]
	sp3D = sp3D[nanmask]

	if occursin("CEMP", run)
		ss_run = replace(run, "CEMP"=>"ScaledSolar")
		sp3D_ss = get_mean_spectra(ss_run)
		nanmask = [all(.!isnan.(s.spectrum)) for s in sp3D_ss]
		sp3D_ss = sp3D_ss[nanmask]

		sp3D_ref = MUST.interpolate(sp3D_ss, (ele_sel=>reference_abundance))
		
		ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
			ele_sel,
			sp3D_ref,
			sp3D,
			λs=λs,
			λe=λe,
		)
		corrections[ss_run] = Δab
	end

	snap1 = cemp_marcs_selected[Symbol(replace(run, "ScaledSolar"=>"CEMP"))]
	sp1D = get_mean_spectra(replace(run, "ScaledSolar"=>"CEMP"), onedimensional=snap1)
	sp1D_ref = MUST.interpolate(sp1D, (ele_sel=>reference_abundance))
	ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
		ele_sel,
		sp1D_ref,
		sp3D,
		λs=λs,
		λe=λe,
	)
	corrections[snap1] = Δab
	eqws[run] = ew
	abs[run] = ab

	ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
		ele_sel,
		sp1D_ref,
		sp1D,
		λs=λs,
		λe=λe,
	)
	eqws[snap1] = ew
	abs[snap1] = ab

	corrections, eqws, abs
end

# ╔═╡ 3826dca1-246c-4c61-95ee-8e497f42aa72
begin
	feh_scaling = [-6, -5, -4, -3, -2]
	cfe_scaling = [4.0, 3.0, 2.0, 1.0, 1.0]
	reference_composition(feh) = begin
		i = argmin(abs.(feh.-feh_scaling))
		cfe_scaling[i]
	end
	
	corrections_time = Dict()
	eqws_time = Dict()
	abs_time = Dict()
	for run in models
		mi = MUST.ModelInformation(run)
		reference_abundance = reference_composition(mi.feh)
		c = abundance_correction_time(run, reference_abundance, λs=4297.0, λe=4303.0)
		corrections_time[run] = c[1]
		eqws_time[run] = c[2]
		abs_time[run] = c[3]
	end
	@info "abundance corrections computed."
end

# ╔═╡ 3d760565-3e87-4c13-9e78-ed6394baecaf
corrections_time

# ╔═╡ 41af66b7-2811-4bd5-972d-5260f67fbfd7


# ╔═╡ d6226d59-d4e3-4652-a7a9-dbe41666511f
md"## Select models to show"

# ╔═╡ b606dca3-665e-4e3c-9218-eaff570551d0
begin
	metallcities = [-6.0, -5.0, -4.0, -3.0, -2.0]
	# hot MS star
	cemp_models_hot = [
		"CEMP_t62.50g40.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]
	scaledsolar_models_hot = [
		"ScaledSolar_t62.50g40.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]

	# cool sub-giant
	cemp_models_ms = [
		"CEMP_t57.50g45.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]
	scaledsolar_models_ms = [
		"ScaledSolar_t57.50g45.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]

	# hot sub-giant
	cemp_models_sg = [
		"CEMP_t52.50g30.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]
	scaledsolar_models_sg = [
		"ScaledSolar_t52.50g30.00m$(MUST.@sprintf("%.3f",m))_v1.0" 
		for m in metallcities
	]
end

# ╔═╡ b292f8a3-9288-4851-baae-31eb00e7f571


# ╔═╡ e2b8a4f0-e180-4f2a-8e15-793bbd320282
let
	plt.close()
	f, ax = plt.subplots(1, 3, figsize=(10, 4), sharey=true, sharex=true, layout="tight")
	
	plot_model(cemp_models, scaledsolar_models; kwargs...) = begin
		marcsmodel_name = [cemp_marcs_selected[Symbol(m)] for m in cemp_models]
		ax[0].plot(metallcities, [corrections_time[m][marcsmodel_name[i]] for (i,m) in enumerate(cemp_models)]; marker="s", markersize=11, kwargs...)
	
		ax[1].plot(metallcities, [corrections_time[m][scaledsolar_models[i]] for (i,m) in enumerate(cemp_models)]; marker="s", markersize=11, kwargs...)
	
		ax[2].plot(metallcities, [corrections_time[m][marcsmodel_name[i]] for (i,m) in enumerate(scaledsolar_models)]; marker="s", markersize=11, kwargs...)	
	end

	plot_model(cemp_models_hot, scaledsolar_models_hot; color="0.6", label=L"\rm 6250\, K,\ 4.0", ls=":", lw=2.5)
	plot_model(cemp_models_ms, scaledsolar_models_ms; color="steelblue", label=L"\rm 5750\, K,\ 4.5", ls="--", lw=2.5)
	plot_model(cemp_models_sg, scaledsolar_models_sg; color="tomato", label=L"\rm 5250\, K,\ 3.0", ls="-", lw=2.5)

	#ax[0].scatter([-5.0], [corrections_t57g30_1D[end]], color="magenta", marker="X", s=100, label=L"\rm 5750\, K,\ 3.0")
	#ax[0].scatter([-5.0], [corrections_t52g45_1D[end]], color="cyan", marker="X", s=100, label=L"\rm 5250\, K,\ 4.5")
	

	ax[0].set_title(L"\rm 3D\ CEMP\ -\ 1D")
	ax[1].set_title(L"\rm 3D\ CEMP\ -\ 3D\ scaled \text{-} solar")
	ax[2].set_title(L"\rm 3D\ scaled \text{-} solar - 1D")

	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[2].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm \Delta\, A(C)")

	ax[1].legend(labelspacing=0.1)
	#ax[0].legend(labelspacing=0.1, loc="upper right")
	ax[0].axhline(0.0, color="k", ls=":", alpha=0.5)
	ax[1].axhline(0.0, color="k", ls=":", alpha=0.5)
	ax[2].axhline(0.0, color="k", ls=":", alpha=0.5)

	ax[0].set_ylim(-1.1, 0.37)
	f
end

# ╔═╡ 387c093e-49cf-40cb-a757-b5a8d50f50fc


# ╔═╡ 01aab765-3346-4212-b28f-c8bc5b51174f
corrections_time

# ╔═╡ fc4cf849-36da-4006-85c4-418e85d5e23a


# ╔═╡ 66af10ce-e791-4cf0-88ef-bceabffe0d20
models_stellarParameters = MUST.ModelInformation.(sort(keys(corrections_time)|>collect))

# ╔═╡ 59ac0194-f85c-400a-8940-20fc1f4eef03
models_teff = [m.teff for m in models_stellarParameters]

# ╔═╡ 5c7a35ae-86ee-4eb8-affa-7ff0d5e99f45
models_logg = [m.logg for m in models_stellarParameters]

# ╔═╡ 85ad8731-b3a3-4960-9ebb-732cc179525e
models_feh = [m.feh for m in models_stellarParameters]

# ╔═╡ 6aad86e1-0e90-4a1f-9ffd-fb29d9c50686


# ╔═╡ c1b44ff5-924e-40b5-9914-0cf87825db4b
"""
	Finds the closest point in a set of points to a given test point.

Args:
	points: An array of points, where each point is a tuple (teff, logg).
    test_point: The test point, as a tuple (teff, logg).

Returns:
    The closest point to the test point.
"""
function closest_point(points, test_point)
  # Normalize the points to account for different magnitudes
  teff_max = maximum(p[1] for p in points)
  logg_max = maximum(p[2] for p in points)
  normalized_points = [(p[1]/teff_max, p[2]/logg_max) for p in points]
  normalized_test_point = (test_point[1]/teff_max, test_point[2]/logg_max)

  # Calculate the distances between the normalized points and the normalized test point
  distances = [sqrt((p[1] - normalized_test_point[1])^2 + (p[2] - normalized_test_point[2])^2) for p in normalized_points]

  # Find the index of the closest point
  closest_index = argmin(distances)

  # Return the original (unnormalized) closest point
  return closest_index
end

# ╔═╡ 22a6add8-76ac-4c25-8acf-c869620505a5
md"A choice needs to be made in what quantities one wants to interpolate. It seems sensible to interpolate in logg and feh. However, the difference in the correction is not solely due to Teff or logg, but a combined effect. This means it should be not be considered lineary in any of them. Maybe it is actually better to just take the closest point."

# ╔═╡ 9f9c5bc2-5926-4d5e-881c-1a403c931272
begin
	cemp_models_corr = [cemp_models_hot, cemp_models_ms, cemp_models_sg]
	ss_models_corr = [scaledsolar_models_hot, scaledsolar_models_ms, scaledsolar_models_sg]
	#cemp_models_corr = [cemp_models_sg, cemp_models_ms]
	#ss_models_corr = [scaledsolar_models_sg, scaledsolar_models_ms]
end

# ╔═╡ 4ffb9e96-89f5-458e-bb36-5abaed2ba7fa
begin
	model_parameters_corr = [MUST.ModelInformation.(m) for m in cemp_models_corr]
	model_teff_corr = [m[1].teff for m in model_parameters_corr]
	model_logg_corr = [m[1].logg for m in model_parameters_corr]
	sortmask_logg_corr = sortperm(model_logg_corr)

	get_feh_1Dcorrections(m; is_ss=false) = corrections_time[m][cemp_marcs_selected[is_ss ? Symbol(replace(m, "ScaledSolar"=>"CEMP")) : Symbol(m)]]
	get_feh_3Dcorrections(m) = corrections_time[m][replace(m, "CEMP"=>"ScaledSolar")]
	

	# feh interpolators for each model
	feh_interpolators_CEMP3D1D = [
		MUST.linear_interpolation(
			metallcities, 
			get_feh_1Dcorrections.(m),
			extrapolation_bc=MUST.Flat()
		) for (i, m) in enumerate(cemp_models_corr)
	]
	feh_interpolators_ScaledSolar3D1D = [
		MUST.linear_interpolation(
			metallcities, 
			get_feh_1Dcorrections.(m, is_ss=true),
			extrapolation_bc=MUST.Flat()
		) for (i, m) in enumerate(ss_models_corr)
	]
	feh_interpolators_CEMPScaledSolar3D = [
		MUST.linear_interpolation(
			metallcities, 
			get_feh_3Dcorrections.(m),
			extrapolation_bc=MUST.Flat()
		) for (i, m) in enumerate(cemp_models_corr)
	]

	# contruct the final interpolator
	get_correction(teff, logg, feh, interpolator) = begin
		ΔA = [f(feh) for f in interpolator]
		#MUST.linear_interpolation(
		#	model_logg_corr[sortmask_logg_corr], ΔA[sortmask_logg_corr],
		#	extrapolation_bc=MUST.Flat()
		#)(logg)
		ΔA[closest_point(zip(model_teff_corr, model_logg_corr), (teff, logg))]
	end
	get_CEMP3D1D_correction(teff, logg, feh) = get_correction(teff, logg, feh, feh_interpolators_CEMP3D1D)
	get_ScaledSolar3D1D_correction(teff, logg, feh) = get_correction(teff, logg, feh, feh_interpolators_ScaledSolar3D1D)
	get_CEMPScaledSolar3D_correction(teff, logg, feh) = get_correction(teff, logg, feh, feh_interpolators_CEMPScaledSolar3D)
end

# ╔═╡ 925ce950-46ca-45fe-8580-b2db83230a4a


# ╔═╡ 06ee2904-3e62-4d8f-a873-068384d3bfe6
md"# Data comparison"

# ╔═╡ 6cc7f3f0-f4f9-4f7a-aceb-c79e7a278154
md"## SAGA data"

# ╔═╡ 1f66327a-5581-4250-ba81-a286e0a04181
saga_all_path = "../cgisess_32aefcdc405d0a0c9fb7936f6f1324aa_data.tsv"

# ╔═╡ dcf005c3-98cd-490a-874d-f922f8d49f3f
saga_co_path = "../cgisess_co.tsv"

# ╔═╡ 9effae34-66bb-47a4-a1a9-c4a8e12475a3
begin
	saga_all_data = MUST.readdlm(saga_all_path, '\t', String, skipstart=1)
	saga_co_data = MUST.readdlm(saga_co_path, '\t', String, skipstart=1)

	teff_all, logg_all, z_all = saga_all_data[:, 5], saga_all_data[:, 6], saga_all_data[:, 7]
	teff_co, logg_co, z_co = saga_co_data[:, 5], saga_co_data[:, 6], saga_co_data[:, 7]

	cfe_all, bafe_all, eufe_all, lafe_all, feh_all = saga_all_data[:, 9], saga_all_data[:, 10], saga_all_data[:, 11], saga_all_data[:, 12], saga_all_data[:, 13]
	cfe_co, o_co, feh_co = saga_co_data[:, 10], saga_co_data[:, 9], saga_co_data[:, 8]

	object_all, reference_all = saga_all_data[:, 1], saga_all_data[:, 2]
	

	binarity_all = saga_all_data[:, 12]

	parameters_all = Dict{Any, Any}(
		"teff" => teff_all, 
		"logg" => logg_all, 
		"z"    => z_all, 
		"cfe"  => cfe_all, 
		"eufe" => eufe_all, 
		"bafe" => bafe_all, 
		"feh"  => feh_all,
		"lafe" => lafe_all,
		"b"    => binarity_all
	)

	parameters_co = Dict{Any, Any}(
		"teff" => teff_co, 
		"logg" => logg_co, 
		"z"    => z_co, 
		"feh"  => feh_co,
		"cfe"  => cfe_co, 
		"ofe" => o_co,
	)

	for (i, p) in parameters_all
		d = zeros(length(p))

		for (j, pj) in enumerate(p)
			d[j] = if pj == ""
				NaN
			elseif (pj[1] == '<') |  (pj[1] == '>')
				parse(Float64, pj[2:end])
			else
				try
					parse(Float64, pj)
				catch
					NaN
				end
			end
		end
		
		parameters_all[i] = d
	end

	for (i, p) in parameters_co
		d = zeros(length(p))

		for (j, pj) in enumerate(p)
			d[j] = if pj == ""
				NaN
			elseif (pj[1] == '<') |  (pj[1] == '>')
				#parse(Float64, pj[2:end])
				NaN
			else
				try
					parse(Float64, pj)
				catch
					NaN
				end
			end
		end
		
		parameters_co[i] = d
	end

	parameters_all["object_id"] = object_all
	parameters_all["reference"] = reference_all

	mask = (parameters_co["feh"].>0.0) 
	parameters_co["feh"][mask] .= NaN
	
	mask = (parameters_all["feh"].>0.0) 
	parameters_all["feh"][mask] .= NaN

	parameters_all
end

# ╔═╡ 8be627f4-3fe6-481a-a322-8ed2d6a015aa
md"To isolate the giant stars, one has to define a limiting log(g)."

# ╔═╡ addf8b42-5994-4073-a595-33df2f3b3fd4
logg_limit = 2.7

# ╔═╡ 481caf39-bc3f-42f5-8edd-94ea4ca6f409
md"We can furthermore define multiple sub-limits for different CEMP groups"

# ╔═╡ 54261a34-ee49-434b-8cdc-6e0e9ffe3368
# For CEMP in general
cfe_limit = 0.7

# ╔═╡ 09b8ee6d-e9f5-4cf1-ab91-305952b3972d
# For CEMP-s and CEMP-r/s
ba_limit = 1.0

# ╔═╡ 1a860a43-3fd2-4b5c-8c8f-bf4341cea882
# For CEMP-r
eu_limit = 1.0

# ╔═╡ 57b503a1-e7b2-464a-a49a-359576688e58
md"Apply r/s cut from Yoon et al. (2016) at A(C)=7.1"

# ╔═╡ dd18395c-e88f-4cb1-af1d-a5a347380464
c_limit = 7.1

# ╔═╡ d325fc52-afd9-474d-8e4b-05455fc340d9
parameters_all["r/s"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit

# ╔═╡ 5aa85e3d-8e21-47b8-a4cc-8897dc49b756
parameters_all["no"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit

# ╔═╡ f4c3987b-e002-4db5-b166-d80bf67ddf75
selection_mask = parameters_all["logg"] .>= logg_limit

# ╔═╡ 0e30cd20-184c-4015-b780-1e94be606e00
let
	teff = parameters_all["teff"][selection_mask]
	logg = parameters_all["logg"][selection_mask]

	mod_teff = model_teff_corr
	mod_logg = model_logg_corr

	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6.5, 7), layout="tight")

	# inset sub-giant
	left, bottom, width, height = [0.18, 0.57, 0.32, 0.32]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	b, bt = get_snapshot(snapshot_id, joinpath(datadir, "CEMP_t52.50g30.00m-5.000_v1.0"))
	plot_profile(b, bt, :T, ax=ax)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	vmin = 4100
	vmax = 6600
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal",
		vmin=vmin,
		vmax=vmax
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.05, location="top", orientation="horizontal")
	cbar.set_label(L"\rm temperature\ [K]", loc="left")
	ax2.set_xlabel(L"\rm [10^{3}\, km]", loc="left")
	ax2.xaxis.set_major_locator(plt.MaxNLocator(3))
	ax2.yaxis.set_major_locator(plt.MaxNLocator(3))

	# Create a Rectangle patch
	rect = matplotlib.patches.Rectangle((40, -40), 5, 5, linewidth=3, edgecolor="steelblue", facecolor="none")
	
	# Add the patch to the Axes
	ax2.add_patch(rect)
	#ax2.set_title(pretty_from_name_short("CEMP_t52.50g30.00m-5.000_v1.0"))


	# inset MS
	left, bottom, width, height = [0.71, 0.18, 0.23, 0.23]
	ax3 = f.add_axes([left, bottom, width, height])
	ax3.tick_params(labelsize=12)
	b, bt = get_snapshot(snapshot_id, joinpath(datadir, "CEMP_t57.50g45.00m-5.000_v1.0"))
	plot_profile(b, bt, :T, ax=ax)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	vmin = 4100
	vmax = 6600
	im = ax3.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal",
		vmin=vmin,
		vmax=vmax
	)
	cbar = f.colorbar(im, ax=ax3, fraction=0.046, pad=0.04, location="top", orientation="horizontal")
	cbar.set_label(L"\rm temperature\ [K]", loc="left")
	ax3.set_xlabel(L"\rm [10^{3}\, km]", loc="left")
	ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
	ax3.yaxis.set_major_locator(plt.MaxNLocator(3))
	#ax2.set_title(pretty_from_name_short("CEMP_t52.50g30.00m-5.000_v1.0"))
	

	#ax.hist2d(teff, logg, rasterized=true, cmap="Greys", bins=20)
	mean_and_nan(x) = length(x) > 0 ? mean(x) : NaN
	std_and_nan(x) = length(x) > 0 ? MUST.std(x) : NaN
	bins = range(2.5, 4.7, length=25)
	bincenters = (bins[2:end] .+ bins[1:end-1]) ./ 2
	meanT = [
		mean_and_nan(teff[bins[i-1] .< logg .< bins[i]]) for i in 2:length(bins)
	]
	sigmaT = [
		std_and_nan(teff[bins[i-1] .< logg .< bins[i]]) for i in 2:length(bins)
	]
	
	#=ax.scatter(teff, logg, s=55, rasterized=true, color="0.3", alpha=0.1, label=L"\rm SAGA", marker="o", zorder=1)
	ax.scatter(parameters_all["teff"][.!selection_mask], parameters_all["logg"][.!selection_mask], s=35, rasterized=true, color="0.8", alpha=0.15, label=L"\rm SAGA", marker="o", zorder=1)=#
	nanmask = (.!isnan.(parameters_all["teff"])) .& (.!isnan.(parameters_all["logg"]))



	
	H, xedges, yedges = MUST.numpy.histogram2d(
		parameters_all["teff"][nanmask], parameters_all["logg"][nanmask], bins=40,
		range=[[2500, 8500], [2.7, 5.0]],
	)
	#H = H + 1
	i = ax.imshow(
		H.T, interpolation="nearest", origin="lower", aspect="auto",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
		#rasterized=true,
		norm=matplotlib.colors.LogNorm(),
		cmap="binary",
		alpha=0.9
	)
	H, xedges, yedges = MUST.numpy.histogram2d(
		parameters_all["teff"][nanmask], parameters_all["logg"][nanmask], bins=40,
		range=[[2500, 8500], [1.0, 2.7]],
	)
	#H = H + 1
	ax.imshow(
		H.T, interpolation="nearest", origin="lower", aspect="auto",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
		#rasterized=true,
		norm=matplotlib.colors.LogNorm(),
		cmap="binary",
		alpha=0.45
	)
	ax.axhline(logg_limit, color="k", ls="--", alpha=0.5)

	

	
	ax.plot(meanT, bincenters, color="k", lw=4, zorder=2)
	ax.fill_betweenx(bincenters, meanT.-sigmaT, meanT.+sigmaT, color="0.5", alpha=0.3)


	

	
	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(mod_teff[1], mod_logg[1]),   # Point on the main plot (data coordinates)
	    xyB=(1, 0),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax,
	    axesB=ax2,                # Use the new axis name
	    color="tomato",
	    linestyle="-",
	    linewidth=3.,
	    alpha=1
	)
	f.add_artist(con1)
	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(mod_teff[1], mod_logg[1]),   # Point on the main plot (data coordinates)
	    xyB=(1, 1),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax,
	    axesB=ax2,                # Use the new axis name
	    color="tomato",
	    linestyle="-",
	    linewidth=3.,
	    alpha=1
	)
	f.add_artist(con1)

	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(mod_teff[2], mod_logg[2]),   # Point on the main plot (data coordinates)
	    xyB=(0, 0),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax,
	    axesB=ax3,                # Use the new axis name
	    color="steelblue",
	    linestyle="-",
	    linewidth=3.,
	    alpha=1
	)
	f.add_artist(con1)
	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(mod_teff[2], mod_logg[2]),   # Point on the main plot (data coordinates)
	    xyB=(0, 1),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax,
	    axesB=ax3,                # Use the new axis name
	    color="steelblue",
	    linestyle="-",
	    linewidth=3.,
	    alpha=1
	)
	f.add_artist(con1)

	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(40,-40),   # Point on the main plot (data coordinates)
	    xyB=(0, 0),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax2,
	    axesB=ax3,                # Use the new axis name
	    color="steelblue",
	    linestyle="--",
	    linewidth=2.5,
	    alpha=0.4
	)
	f.add_artist(con1)
	con1 = matplotlib.patches.ConnectionPatch(
	    xyA=(45,-40),   # Point on the main plot (data coordinates)
	    xyB=(0, 1),                       # Point on the inset plot (axes fraction)
	    coordsA="data",
	    coordsB="axes fraction",
	    axesA=ax2,
	    axesB=ax3,                # Use the new axis name
	    color="steelblue",
	    linestyle="--",
	    linewidth=2.5,
	    alpha=0.4
	)
	f.add_artist(con1)



	
	
	ax.scatter([mod_teff[1]], [mod_logg[1]], color="tomato", marker="X", s=600, label=L"\rm 3D\ RHD\ models", edgecolor="k", zorder=100)
	ax.scatter([mod_teff[2]], [mod_logg[2]], color="steelblue", marker="X", s=600, label=L"\rm 3D\ RHD\ models", edgecolor="k", zorder=100)
	ax.scatter([6250, 6000, 5750], [4.0,4.0,4.0], color="yellowgreen", marker="X", s=100, zorder=100)


	
	# add text
	labelnice(teff, logg) = L"\rm T_{eff}=\ "*"$(round(Int, teff))"*L"\rm \, K "*"\n"*L"\rm \log g =\ "*"$(round(logg, digits=3))"
	ax.text(mod_teff[2]+2350, mod_logg[2]+0.15, L"\mathbf{main-sequence}"*"\n"*labelnice(mod_teff[2], mod_logg[2]), ha="left", va="center", color="steelblue", fontsize=15)

	ax.text(mod_teff[1]-700, mod_logg[1]-0.04, L"\mathbf{sub-giant}"*"\n"*labelnice(mod_teff[1], mod_logg[1]), ha="left", va="center", color="tomato", fontsize=15)


	ax.set_ylim(ax.get_ylim()[1], ax.get_ylim()[0])
	ax.set_xlim(ax.get_xlim()[1], ax.get_xlim()[0])

	ax.set_ylabel(L"\rm \log g")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	#ax.legend(labelspacing=0.1, loc="lower center", ncol=2, columnspacing=0.5)

	ax.set_ylim(4.95, 1.45)
	ax.set_xlim(8200, 3200)

	f
end

# ╔═╡ b377a70f-c7cf-43d8-b6c9-cd8d42fa2835


# ╔═╡ dfe8046d-a75e-424c-bc51-a2e429b97f70
begin
	corrections_saga_all = fill!(similar(parameters_all["cfe"]), 0.0)
	
	for i in eachindex(parameters_all["cfe"])
		logg = parameters_all["logg"][i]
		teff = parameters_all["teff"][i]
		feh = parameters_all["feh"][i]
		corrections_saga_all[i] = get_CEMP3D1D_correction(teff, logg, feh)
	end

	@info "3D corrections applied."
end

# ╔═╡ 92228544-f92c-4c50-ad50-3cb818bb4523


# ╔═╡ a7707e78-18f9-4c96-b5ff-c51cb9311652
md"## Galactic CEMP Distributions"

# ╔═╡ dea25815-eebf-4308-a06f-09a7a49862c3
md"Predictions from Hartwig et al. 2018, extracted from their Figure:"

# ╔═╡ 80ba8489-3a1a-4d3b-88fe-56954397c7df
begin
	hartwig18_fiducial = MUST.readdlm("../hartwig18_fiducial.csv", ',')
	hartwig18_faint20 = MUST.readdlm("../hartwig18_faint20.csv", ',')
end;

# ╔═╡ c0a0d7cd-6b6f-47a2-a666-959e2cb68b7c
#metallicity_bin_edges = [-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]
metallicity_bin_edges = [-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]

# ╔═╡ 5709600f-4817-43bd-b34b-404e173efb43
metallicity_bin_centers = (metallicity_bin_edges[2:end] .+ metallicity_bin_edges[1:end-1]) ./ 2

# ╔═╡ c2cae081-827b-4606-91c6-a5610b3f72b9
function bin_parameter(bins, general_mask=trues(size(data, 1)); data=saga_all_data, para=parameters_all, selection=selection_mask, name="feh")
	bin_centers = (bins[2:end] .+ bins[1:end-1]) ./ 2
	bin_masks = falses(size(data, 1), length(bin_centers))
	for i in eachindex(bin_centers)
		bin_masks[:, i] .= (bins[i] .<= para[name]) .&
		(para[name] .< bins[i+1]) .& (.!isnan.(para[name]))

		bin_masks[:, i] .= bin_masks[:, i] .& general_mask .& selection
	end

	bin_masks
end

# ╔═╡ ff2843c9-a864-466c-b8fc-87b95e1da2f2
begin
	bins_general = bin_parameter(
		metallicity_bin_edges, 
		.!isnan.(parameters_all["cfe"])
	)
	count_bins_general = reshape(count(bins_general, dims=1), :)

	# binarity
	bins_binary = bin_parameter(
		metallicity_bin_edges, 
		abs.(parameters_all["b"]) .>= 0.5
	)
	count_bins_binary = reshape(count(bins_binary, dims=1), :)

	bins_nonbinary = bin_parameter(
		metallicity_bin_edges, 
		abs.(parameters_all["b"]) .< 0.5
	)
	count_bins_nonbinary = reshape(count(bins_nonbinary, dims=1), :)
	
	# CEMP
	bins_CEMP = bin_parameter(
		metallicity_bin_edges, 
		parameters_all["cfe"] .>= cfe_limit
	)
	count_bins_CEMP = reshape(count(bins_CEMP, dims=1), :)

	# CEMP (corrected)
	bins_CEMP_corr = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
	)
	count_bins_CEMP_corr = reshape(count(bins_CEMP_corr, dims=1), :)
	

	# CEMP-r/s
	bins_CEMP_rs = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		#(parameters_all["bafe"] .>= ba_limit) .& 
		#(parameters_all["eufe"] .>= eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit)
	)
	count_bins_CEMP_rs = reshape(count(bins_CEMP_rs, dims=1), :)

	# CEMP-r/s (corrected)
	bins_CEMP_rs_corr = bin_parameter(
		metallicity_bin_edges, 
		((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& 
		#(parameters_all["bafe"] .>= ba_limit) .& 
		#(parameters_all["eufe"] .>= eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .> c_limit)
	)
	count_bins_CEMP_rs_corr = reshape(count(bins_CEMP_rs_corr, dims=1), :)
	

	# CEMP-s
	bins_CEMP_s = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		(parameters_all["bafe"] .>= ba_limit) .& 
		(parameters_all["eufe"] .< eu_limit)
	)
	count_bins_CEMP_s = reshape(count(bins_CEMP_s, dims=1), :)
	

	# CEMP-no
	bins_CEMP_no = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		#(parameters_all["bafe"] .< ba_limit) .& 
		#(parameters_all["eufe"] .< eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit)
	)
	count_bins_CEMP_no = reshape(count(bins_CEMP_no, dims=1), :)

	# CEMP-no (corrected)
	bins_CEMP_no_corr = bin_parameter(
		metallicity_bin_edges, 
		((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& 
		#(parameters_all["bafe"] .< ba_limit) .& 
		#(parameters_all["eufe"] .< eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .<= c_limit)
	)
	count_bins_CEMP_no_corr = reshape(count(bins_CEMP_no_corr, dims=1), :)
	

	# CEMP-unclassified
	bins_CEMP_nothing = bins_CEMP .& (.!bins_CEMP_rs) .& (.!bins_CEMP_s) .& (.!bins_CEMP_no) 
	count_bins_CEMP_nothing = reshape(count(bins_CEMP_nothing, dims=1), :)
	

	@info "Bin counts (all): $(count_bins_general)"
	@info "Bin counts (CEMP): $(count_bins_CEMP)"
	@info "Bin counts (CEMP corrected): $(count_bins_CEMP_corr)"
	@info "Bin counts (CEMP-r/s): $(count_bins_CEMP_rs)"
	@info "Bin counts (CEMP-r/s corrected): $(count_bins_CEMP_rs_corr)"
	@info "Bin counts (CEMP-s): $(count_bins_CEMP_s)"
	@info "Bin counts (CEMP-no): $(count_bins_CEMP_no)"
	@info "Bin counts (CEMP-no corrected): $(count_bins_CEMP_no_corr)"
	@info "Bin counts (CEMP unclassified): $(count_bins_CEMP_nothing)"
	@info "Bin counts (binarity): $(count_bins_binary)"
	@info "Bin counts (non-binarity): $(count_bins_binary)"
end

# ╔═╡ 2404c7e6-b3cb-4c92-9261-bdef297e727f
get_n_below(lower_than_feh, feh, mask) = begin
	below_limit = (feh .<= lower_than_feh) .& selection_mask .& mask
	count(below_limit)
end

# ╔═╡ b9b2059c-c098-48d8-94bb-587715552774
get_n_below_frac(lower_than_feh, feh, mask) = begin
	below_limit = (feh .<= lower_than_feh) .& selection_mask
	below_limit_mask = (feh .<= lower_than_feh) .& selection_mask .& mask
	count(below_limit_mask) / count(below_limit)
end

# ╔═╡ ea48b2d1-c752-4047-8523-474aaecc82b0
begin
	@info "cumsum of those bins:"
	@info cumsum(count_bins_general)
	@info cumsum(count_bins_CEMP)
	@info cumsum(count_bins_CEMP_corr)
end

# ╔═╡ 67ddc416-4cbb-4fb5-9d3c-b83ce557d141
@info "Stars with log(g) higher than limit: $(count(selection_mask .& (parameters_all["feh"] .< 0.0)))"

# ╔═╡ 64b264bd-d556-4b08-8c64-b77c51849be2
@info "Stars with log(g) higher than limit + [C/Fe]: $(count(selection_mask .& (.!isnan.(parameters_all["cfe"])) .& (parameters_all["feh"] .< 0.0)))"

# ╔═╡ d67f968c-7405-46ff-8588-da811fc79dbc
@info "CEMP stars: $(count(selection_mask .& (.!isnan.(parameters_all["cfe"])) .& (parameters_all["feh"] .< 0.0) .& (parameters_all["cfe"] .> 0.7)))"

# ╔═╡ bd9e94c2-1b09-4cb7-82f9-df86ee21a5aa
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	#ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 1D", marker="s", markersize=10)
	#ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="steelblue", lw=2.0, label=L"\rm 3D", marker="s", markersize=10)

	#metallicities = [-6, -5, -4.5, -4, -3.75, -3.5, -3., -2.5, -2]
	metallicities = [-6, -5, -4.5, (range(-4, -1.5, step=0.25) |> collect)...]
	allm = trues(length(parameters_all["cfe"]))
	
	n1DCEMP = [
		get_n_below(
		z, parameters_all["feh"], parameters_all["cfe"] .>= cfe_limit
		) for z in metallicities
	]
	n1D = [
		get_n_below(
		z, parameters_all["feh"], allm
		) for z in metallicities
	]
	
	n3DCEMP = [
		get_n_below(
		z, parameters_all["feh"], (parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
		) for z in metallicities
	]
	n3D = [
		get_n_below(
		z, parameters_all["feh"], allm
		) for z in metallicities
	]
	ax.plot(metallicities, n1DCEMP ./ n1D*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 1D", marker="s", markersize=10)
	ax.plot(metallicities, n3DCEMP ./ n3D*100, zorder=10, color="steelblue", lw=2.0, label=L"\rm 3D", marker="s", markersize=10)
	#ax.plot(metallicity_bin_centers, count_bins_CEMP ./ count_bins_general*100, zorder=10, color="k", lw=2.0, label=L"\rm 1D", marker="s")
	#ax.plot(metallicity_bin_centers, count_bins_CEMP_corr ./ count_bins_general*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 3D", marker="s")

	ax.set_zorder(10)
	ax.patch.set_visible(false)
	
	ax.set_xlabel(L"\rm [Fe/H]")
	ax.set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax.legend(loc="lower left")

	ax.set_xlim(-6.2, -1.8)
	ax.set_ylim(0.17*100, 1.04*100)

	open("CEMP_cumsum_fractions.txt", "w") do f
		cf = n1DCEMP ./ n1D
		cf_3D = n3DCEMP ./ n3D
		write(f, "# [Fe/H] 1D 3D\n")
		for i in eachindex(metallicities)
			line = MUST.@sprintf("%.2f    %.2f    %.2f\n", metallicities[i], cf[i], cf_3D[i])
			write(f, line)
		end
	end

	f.savefig("cummulative_CEMP_fraction.pdf")
	
	f
end

# ╔═╡ 98316a96-6b40-4f04-8817-7bb8c9b992fc
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	#ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm 1D", marker="s")
	#ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 3D", marker="s")

	ax.plot(metallicity_bin_centers, count_bins_CEMP ./ count_bins_general*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 1D", marker="s", markersize=10)
	ax.plot(metallicity_bin_centers, count_bins_CEMP_corr ./ count_bins_general*100, zorder=10, color="steelblue", lw=2.0, label=L"\rm 3D", marker="s", markersize=10)

	ax.set_zorder(10)
	ax.patch.set_visible(false)
	
	ax.set_xlabel(L"\rm [Fe/H]")
	ax.set_ylabel(L"\rm N_{[Fe/H], CEMP}\ /\ N_{[Fe/H]}\ [\%]")
	ax.legend(loc="lower left")

	ax.set_xlim(-6.2, -1.8)
	ax.set_ylim(0.1*100, 1.04*100)

	open("CEMP_fractions.txt", "w") do f
		cf = count_bins_CEMP ./ count_bins_general
		cf_3D = count_bins_CEMP_corr ./ count_bins_general
		write(f, "# [Fe/H] 1D 3D\n")
		for i in eachindex(metallicity_bin_centers)
			line = MUST.@sprintf("%.2f    %.2f    %.2f\n", metallicity_bin_centers[i], cf[i], cf_3D[i])
			write(f, line)
		end
	end

	f.savefig("CEMP_fraction.pdf")
	
	f
end

# ╔═╡ 37c90a88-2615-461b-9e52-6a5df9e9d76a
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(6.5, 7), sharex=true, sharey=true)

	plt.subplots_adjust(wspace=0)

	lw = 3
	ms = 10

	#=ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=lw, marker="X", alpha=1, ls=":", markersize=ms)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="steelblue", lw=lw, marker="X", alpha=1, ls=":", markersize=ms)
	
	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=lw, marker="s", markersize=ms)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="steelblue", lw=lw, marker="s", markersize=ms)

	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=lw, marker="o", ls="--", markeredgecolor="tomato", markerfacecolor="w", markersize=ms)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="steelblue", lw=lw, marker="o", ls="--", markeredgecolor="steelblue", markerfacecolor="w", markersize=ms)
	=#


	metallicities = [-5.7, -5, -4, -3.5, -3., -2.5, -2, -1.5]
	#metallicities = [-5.7, -5, (range(-4, -1.5, length=50) |> collect)...]
	#metallicities = [-5.7, -5, (range(-4, -1.5, length=50) |> collect)...]
	allm = trues(length(parameters_all["cfe"]))
	
	rs_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .> c_limit)
	rs_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit)
	
	no_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .<= c_limit)
	no_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit)

	all_3D = (parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
	all_1D = (parameters_all["cfe"]) .>= cfe_limit

	# 1D 
	n1D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_1D) for z in metallicities
	]
	n1D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_1D) for z in metallicities
	]
	n1D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_1D) for z in metallicities
	]

	# 3D 
	n3D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_3D) for z in metallicities
	]
	n3D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_3D) for z in metallicities
	]
	n3D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_3D) for z in metallicities
	]

	ax[0].plot(metallicities, n1D_all*100, zorder=10, color="tomato", lw=lw, marker="X", alpha=1, ls=":", markersize=ms)
	ax[1].plot(metallicities, n3D_all*100, zorder=10, color="steelblue", lw=lw, marker="X", alpha=1, ls=":", markersize=ms)
	
	ax[0].plot(metallicities, n1D_no*100, zorder=10, color="tomato", lw=lw, marker="s", markersize=ms, markerfacecolor="w")
	ax[1].plot(metallicities, n3D_no*100, zorder=10, color="steelblue", lw=lw, marker="s", markersize=ms, markerfacecolor="w")

	ax[0].plot(metallicities, n1D_rs*100, zorder=10, color="tomato", lw=lw, marker="o", ls="--", markeredgecolor="tomato", markersize=ms)
	ax[1].plot(metallicities, n3D_rs*100, zorder=10, color="steelblue", lw=lw, marker="o", ls="--", markeredgecolor="steelblue", markersize=ms)

	
	ax[0].plot([], [], zorder=10, color="k", lw=2, label=L"\rm \text{CEMP-no}", marker="s", markersize=ms, markerfacecolor="w")
	ax[0].plot([], [], zorder=10, color="k", lw=2, label=L"\rm \text{CEMP-r/s}", marker="o", ls="--", markeredgecolor="k", markersize=ms)
	ax[0].plot([],[], zorder=10, color="k", lw=2, marker="X", alpha=1, ls=":", markersize=ms, label=L"\rm CEMP\ (all)")

	#=ax[0].plot(hartwig18_fiducial[:, 1], hartwig18_fiducial[:, 2].*100, color="k")
	ax[0].plot(hartwig18_faint20[:, 1], hartwig18_faint20[:, 2].*100, color="cyan")

	ax[1].plot(hartwig18_fiducial[:, 1], hartwig18_fiducial[:, 2].*100, color="k")
	ax[1].plot(hartwig18_faint20[:, 1], hartwig18_faint20[:, 2].*100, color="cyan")
	=#
	
	ax[0].set_xlabel(L"\rm [Fe/H]", fontsize="x-large")
	ax[1].set_xlabel(L"\rm [Fe/H]", fontsize="x-large")
	ax[0].set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]", fontsize="x-large")
	ax[0].legend(loc="lower center", ncol=3, bbox_to_anchor=(1.0, 0.95), fontsize="large", columnspacing=0.7)
	#ax[1].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.0))

	ax[0].text(0.95, 0.97, L"\rm \mathbf{1D}", ha="right", va="top", color="tomato", fontsize=22, transform=ax[0].transAxes)
	ax[1].text(0.95, 0.97, L"\rm \mathbf{3D}", ha="right", va="top", color="steelblue", fontsize=22, transform=ax[1].transAxes)
	
	ax[0].set_xlim(-6.1, -1.3)
	ax[0].set_ylim(-7, 107)

	f.savefig("cummulative_CEMP_nors_fraction.pdf")

	
	f
end

# ╔═╡ d6ee50f7-d83a-40ff-b6d8-d24eb55d758f
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 7))

	plt.subplots_adjust(wspace=0)

	lw = 0
	ms = 10


	#metallicities = [-5.7, -5, -4, -3.5, -3., -2.5, -2, -1.5]
	#metallicities = [-5.7, -5, (range(-4, -1.5, length=50) |> collect)...]
	metallicities = range(-6.0, -1.5, length=300) |> collect
	allm = trues(length(parameters_all["cfe"]))
	
	rs_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .> c_limit)
	rs_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit)
	
	no_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .<= c_limit)
	no_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit)

	all_3D = (parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
	all_1D = (parameters_all["cfe"]) .>= cfe_limit

	# 1D 
	n1D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_1D) for z in metallicities
	]
	n1D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_1D) for z in metallicities
	]
	n1D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_1D) for z in metallicities
	]

	# 3D 
	n3D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_3D) for z in metallicities
	]
	n3D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_3D) for z in metallicities
	]
	n3D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_3D) for z in metallicities
	]

	ax.plot(hartwig18_fiducial[:, 1], hartwig18_fiducial[:, 2].*100, color="tomato", lw=10, alpha=0.4, label="Hartwig et al. (2018), fiducial")
	ax.plot(hartwig18_faint20[:, 1], hartwig18_faint20[:, 2].*100, color="steelblue", lw=10, alpha=0.4, label="Hartwig et al. (2018), 20% faint SNe")

	#ax.plot(metallicities, n1D_all*100, color="tomato", marker="X", alpha=1, ls="", markersize=ms)
	#ax.plot(metallicities, n3D_all*100, color="steelblue", marker="X", alpha=1, ls="", markersize=ms)

	metal_low = metallicities .<-3.75
	ax.plot(metallicities[.!metal_low], n1D_no[.!metal_low]*100, color="tomato", marker="x", markersize=ms, ls="", lw=5)
	ax.plot(metallicities[.!metal_low], n3D_no[.!metal_low]*100, color="steelblue", marker="x", markersize=ms, ls="", lw=5)

	ax.plot(metallicities[metal_low], n1D_all[metal_low]*100, color="tomato", marker="o", markersize=ms*1.3, ls="", markerfacecolor="None", lw=5, markeredgewidth=1.5)
	ax.plot(metallicities[metal_low], n3D_all[metal_low]*100, color="steelblue", marker="o", markersize=ms*1.3, ls="", markerfacecolor="None", lw=5, markeredgewidth=1.5)

	#ax[0].plot(metallicities, n1D_rs*100, zorder=10, color="tomato", lw=lw, marker="o", ls="--", markeredgecolor="tomato", markersize=ms)
	#ax[1].plot(metallicities, n3D_rs*100, zorder=10, color="steelblue", lw=lw, marker="o", ls="--", markeredgecolor="steelblue", markersize=ms)

	
	ax.plot([], [], zorder=10, color="k", lw=0, label=L"\rm CEMP-no", marker="x", markersize=ms, ls="")
	#ax[1].plot([], [], zorder=10, color="k", lw=2, label=L"\rm CEMP-r/s", marker="o", ls="--", markeredgecolor="k", markersize=ms)
	ax.plot([],[], zorder=10, color="k", lw=0, marker="o", alpha=1, ls="", markersize=ms, label=L"\rm CEMP\ (all)", markerfacecolor="w", markeredgewidth=2)

	#ax[1].plot(hartwig18_fiducial[:, 1], hartwig18_fiducial[:, 2].*100, color="k")
	#ax[1].plot(hartwig18_faint20[:, 1], hartwig18_faint20[:, 2].*100, color="cyan")

	
	
	#ax[0].set_xlabel(L"\rm [Fe/H]")
	ax.set_xlabel(L"\rm [Fe/H]")
	ax.set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax.legend(loc="lower left", labelspacing=0.1)
	#ax[1].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.0))

	ax.text(-3.65, 75, L"\rm \mathbf{1D}", ha="left", va="bottom", color="tomato", fontsize=22)
	ax.text(-3.9, 51., L"\rm \mathbf{3D}", ha="right", va="top", color="steelblue", fontsize=22)
	
	ax.set_xlim(-4.15, -2.1)
	ax.set_ylim(-9, 109)

	f.savefig("cummulative_CEMP_hartwig.pdf")

	
	f
end

# ╔═╡ 96ca3162-7ccf-4027-978f-238613e5413e
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(7, 8), sharex=true, sharey=true)

	plt.subplots_adjust(wspace=0, hspace=0)

	lw = 0
	ms = 10


	#metallicities = [-5.7, -5, -4, -3.5, -3., -2.5, -2, -1.5]
	#metallicities = [-5.7, -5, (range(-4, -1.5, length=50) |> collect)...]
	metallicities = range(-6.0, -1.5, length=300) |> collect
	allm = trues(length(parameters_all["cfe"]))
	
	rs_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .> c_limit)
	rs_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit)
	
	no_3D = ((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .<= c_limit)
	no_1D = (parameters_all["cfe"] .>= cfe_limit) .& ((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit)

	all_3D = (parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
	all_1D = (parameters_all["cfe"]) .>= cfe_limit

	# 1D 
	n1D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_1D) for z in metallicities
	]
	n1D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_1D) for z in metallicities
	]
	n1D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_1D) for z in metallicities
	]

	# 3D 
	n3D_all = [
		get_n_below_frac(z, parameters_all["feh"], all_3D) for z in metallicities
	]
	n3D_rs = [
		get_n_below_frac(z, parameters_all["feh"], rs_3D) for z in metallicities
	]
	n3D_no = [
		get_n_below_frac(z, parameters_all["feh"], no_3D) for z in metallicities
	]

	ax[0].plot(hartwig18_fiducial[:, 1], hartwig18_fiducial[:, 2].*100, color="tomato", lw=10, alpha=0.4)
	ax[1].plot(hartwig18_faint20[:, 1], hartwig18_faint20[:, 2].*100, color="steelblue", lw=10, alpha=0.4)

	metal_low = metallicities .<-3.75
	ax[0].plot(metallicities, n1D_no*100, color="tomato", marker="x", markersize=ms, ls="", lw=5)
	ax[1].plot(metallicities, n3D_no*100, color="steelblue", marker="x", markersize=ms, ls="", lw=5)

	ax[0].plot(metallicities[metal_low][1:3:end], n1D_all[metal_low][1:3:end]*100, color="tomato", marker="o", markersize=ms*1.0, ls="", markerfacecolor="w", lw=5, markeredgewidth=1.9)
	ax[1].plot(metallicities[metal_low][1:3:end], n3D_all[metal_low][1:3:end]*100, color="steelblue", marker="o", markersize=ms*1.0, ls="", markerfacecolor="w", lw=5, markeredgewidth=1.9)

	
	# for the labels
	ax[0].plot([],[], color="tomato", lw=10, alpha=0.4, label=L"\rm 40 \%\ faint + 60 \%\ normal\ CCSNe")
	ax[0].plot([],[], color="steelblue", lw=10, alpha=0.4, label=L"\rm 20 \%\ faint + 80 \%\ normal\ CCSNe")
	ax[0].plot([], [], zorder=10, color="k", lw=0, label=L"\rm \text{CEMP-no}", marker="x", markersize=ms, ls="")
	ax[0].plot([],[], zorder=10, color="k", lw=0, marker="o", alpha=1, ls="", markersize=ms, label=L"\rm CEMP\ (all)", markerfacecolor="w", markeredgewidth=2)

	
	
	ax[1].set_xlabel(L"\rm [Fe/H]", fontsize="x-large")
	ax[0].set_xlabel(L"\rm [Fe/H]", fontsize="x-large")
	ax[0].set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]", fontsize="x-large")
	#ax[1].set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax[0].legend(loc="lower center", labelspacing=0.1, bbox_to_anchor=(1.0, 0.97), ncol=2, columnspacing=0.25, fontsize=15)
	#ax[1].legend(loc="lower left", labelspacing=0.1)
	#ax[1].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.0))

	ax[0].text(0.95, 0.95, L"\rm \mathbf{1D}", ha="right", va="top", color="tomato", fontsize=25, transform=ax[0].transAxes)
	ax[1].text(0.95, 0.95, L"\rm \mathbf{3D}", ha="right", va="top", color="steelblue", fontsize=25, transform=ax[1].transAxes)
	
	ax[0].set_xlim(-4.25, -2.1)
	ax[0].set_ylim(5, 105)

	f.savefig("cummulative_CEMP_hartwig_split.pdf")

	
	f
end

# ╔═╡ 556f2faa-a86e-4988-b697-af42f5b9be84
draw_ellipse_axes(ellipse, ax; N=1, skip_vertical=false, skip_horizontal=false, kwargs...) = begin
	    # Extract parameters from the Ellipse object
	    x0, y0 = MUST.pyconvert(Array, ellipse.center)
	    width = MUST.pyconvert(Float64, ellipse.width)
	    height = MUST.pyconvert(Float64, ellipse.height)
	    angle = MUST.pyconvert(Float64, ellipse.angle)
	
	    # Convert angle to radians
	    theta = MUST.pyconvert(Float64, MUST.numpy.deg2rad(angle))
	
	    # Major axis vector
	    dx_major = 0.5 * width * cos(theta)
	    dy_major = 0.5 * width * sin(theta)
	
	    # Minor axis vector
	    dx_minor = 0.5 * height * sin(theta)
	    dy_minor = -0.5 * height * cos(theta)

		@show x0,y0,angle
	    # Draw major axis line
		if !skip_horizontal
	    	ax.plot([x0 - N*dx_major, x0 + N*dx_major], [y0 - N*dy_major, y0 + N*dy_major]; kwargs...)
		end

		if !skip_vertical
		    # Draw minor axis line
		    ax.plot([x0 - N*dx_minor, x0 + N*dx_minor], [y0 - N*dy_minor, y0 + N*dy_minor]; kwargs...)
		end
	end

# ╔═╡ f7605825-8a21-4b20-87c9-cdf54f2f04de
begin
	g1_dat = MUST.readdlm("G1.csv", ',')
	g2_dat = MUST.readdlm("G2.csv", ',')
	g3_dat = MUST.readdlm("G3.csv", ',')

	
	g1_yoon = GMM(1, g1_dat, nIter=1500, nInit=1500, kind=:full)
	g2_yoon = GMM(1, g2_dat, nIter=1500, nInit=1500, kind=:full)
	g3_yoon = GMM(1, g3_dat, nIter=1500, nInit=1500, kind=:full)
end

# ╔═╡ f9032b6e-2a5c-43d5-9ff9-1cf4c442b858
begin
	transform_to_ellipse(sigma) = begin
		# 2. Calculate eigenvalues and eigenvectors
		F = eigen(sigma)
		eigenvalues, eigenvectors = F.values, F.vectors
		
		# 3. Get the geometric properties for the ellipse
		major_vec = eigenvectors[:, argmax(eigenvalues)]
		major_len = sqrt(maximum(eigenvalues))
		minor_len = sqrt(minimum(eigenvalues))
		
		# Angle in degrees for matplotlib
		angle_deg = rad2deg(atan(major_vec[2], major_vec[1]))
		
		# Width and height are the full axis lengths (2 * semi-axis)
		width = 2 * major_len
		height = 2 * minor_len

		width, height, angle_deg
	end

	lowfeh_mask_yoon = selection_mask .& (parameters_all["feh"] .<= -1)
	CEMP_mask_yoon = (parameters_all["cfe"] .>= cfe_limit) .& lowfeh_mask_yoon
	CEMP_mask_corr_yoon = (parameters_all["cfe"] .+ corrections_saga_all .>= cfe_limit) .& lowfeh_mask_yoon
	
	feh_yoon = parameters_all["feh"][CEMP_mask_yoon]
	cfe_yoon = parameters_all["cfe"][CEMP_mask_yoon]
	c_yoon = cfe_yoon .+ feh_yoon .+ 8.560
	x_1D_yoon = zeros(length(feh_yoon), 2)
	x_3D_yoon = zeros(length(feh_yoon), 2)
	x_1D_yoon[:, 1] .= feh_yoon
	x_1D_yoon[:, 2] .= c_yoon
	x_3D_yoon[:, 1] .= feh_yoon
	x_3D_yoon[:, 2] .= c_yoon .+ corrections_saga_all[CEMP_mask_yoon]
	x_1D_tuple_yoon = [[a, b] for (a,b) in zip(x_1D_yoon[:,1], x_1D_yoon[:,2])]
	x_3D_tuple_yoon = [[a, b] for (a,b) in zip(x_3D_yoon[:,1], x_3D_yoon[:,2])]

	el1_yoon = transform_to_ellipse(covars(g1_yoon)[1])
	el2_yoon = transform_to_ellipse(covars(g2_yoon)[1])
	el3_yoon = transform_to_ellipse(covars(g3_yoon)[1])
	
	n1_yoon = MvNormal(means(g1_yoon)[1,:], covars(g1_yoon)[1])
	n2_yoon = MvNormal(means(g2_yoon)[1,:], covars(g2_yoon)[1])
	n3_yoon = MvNormal(means(g3_yoon)[1,:], covars(g3_yoon)[1])

	# use these ellipses to judge the point association
	pdf1_1D = Distributions.logpdf(n1_yoon, x_1D_tuple_yoon)
	pdf2_1D = Distributions.logpdf(n2_yoon, x_1D_tuple_yoon)
	pdf3_1D = Distributions.logpdf(n3_yoon, x_1D_tuple_yoon)
	g1_mask_1D = (pdf1_1D .>= pdf2_1D) .&& (pdf1_1D .>= pdf3_1D)
	g2_mask_1D = (pdf2_1D .>= pdf1_1D) .&& (pdf2_1D .>= pdf3_1D)
	g3_mask_1D = (pdf3_1D .>= pdf1_1D) .&& (pdf3_1D .>= pdf2_1D)

	# then compute the new GMM after shifting those groups and then plot line
	g1_yoon_3D = GMM(1, x_3D_yoon[g1_mask_1D,:], nIter=1500, nInit=1500, kind=:full)
	g2_yoon_3D = GMM(1, x_3D_yoon[g2_mask_1D,:], nIter=1500, nInit=1500, kind=:full)
	g3_yoon_3D = GMM(1, x_3D_yoon[g3_mask_1D,:], nIter=1500, nInit=1500, kind=:full)
	g123_yoon_3D = GMM(3, x_3D_yoon[:,:], nIter=1500, nInit=1500, kind=:full)
	el1_yoon_3D = transform_to_ellipse(covars(g1_yoon_3D)[1])
	el2_yoon_3D = transform_to_ellipse(covars(g2_yoon_3D)[1])
	el3_yoon_3D = transform_to_ellipse(covars(g3_yoon_3D)[1])
	el1_yoon123_3D = transform_to_ellipse(covars(g123_yoon_3D)[1])
	el2_yoon123_3D = transform_to_ellipse(covars(g123_yoon_3D)[2])
	el3_yoon123_3D = transform_to_ellipse(covars(g123_yoon_3D)[3])
end

# ╔═╡ 6d3bb0c2-b17d-49e8-8141-e22131f65f33
let
	plt.close()
	f = plt.figure(figsize=(8, 8))

	gs = matplotlib.gridspec.GridSpec(
		3, 3, 
		width_ratios=[3, 1, 0.2], 
		height_ratios=[1, 3, 0.2],
        wspace=0.00, hspace=0.00
	)
	
	# Create main scatter plot
	ax = plt.subplot(gs[1, 0])
	
	# Create top histogram
	ax_xhist = plt.subplot(gs[0, 0], sharex=ax)
	
	# Create right histogram
	ax_yhist = plt.subplot(gs[1, 1], sharey=ax)


	
	ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.8), -4.5, 0.8, color="orange", alpha=0.25, ec="0.5"
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.7, 5.9), -3, 0.8, color="green", alpha=0.25, angle=45, ec="0.5"
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.6, 7.9), -2, 3, color="blue", alpha=0.25, angle=-30, ec="0.5"
	)
	ax.add_patch(ellipse1)
	ax.add_patch(ellipse2)
	ax.add_patch(ellipse3)

	ax.text(-6.7, 7.4, L"\rm \mathbf{Group\ III}", color="0.3", fontsize="small")
	ax.text(-3.2, 5.3, L"\rm \mathbf{Group\ II}", color="0.3", fontsize="small")
	ax.text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="0.3", fontsize="small", ha="right")
	
	ax.text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax.text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")

	text_group3 = """
	$(L"\rm \mathbf{CEMP-no}")
	Faint SNe
	Pop III
	"""
	ax.text(-6.7, 7.5, text_group3, color="0.5", fontsize="x-small", ha="left", va="bottom")

	
	text_group2 = """
	$(L"\rm \mathbf{CEMP-no}")
	Standard CCSNe
	Pop II
	"""
	ax.text(-3.2, 5.15, text_group2, color="0.5", fontsize="x-small", ha="left", va="top")


	text_group1 = """
	$(L"\rm \mathbf{CEMP-s\ +\ CEMP-r/s}")
	MRSNe, NS-NS or NS-BH merger
	AGB contribution
	"""
	ax.text(-1.4, 9.5, text_group1, color="0.5", fontsize="x-small", ha="right", va="bottom")
	

	
	CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe .+ feh .+ 8.560
	#ax.scatter(feh, c, s=20, marker="s", c="None", alpha=0.2, linewidth=1, edgecolor="k")
	ax.scatter(feh, c+corrections_saga_all[CEMP_mask], s=15, marker="o", c="None", alpha=0.6, edgecolor="steelblue", rasterized=true)

	
	# add histograms
	ax_xhist.hist(feh, bins=30, color="tomato")
	ax_yhist.hist(c, orientation="horizontal", bins=30, label="1D", color="tomato")
	ax_yhist.hist(c+corrections_saga_all[CEMP_mask], orientation="horizontal", bins=30, label="3D", histtype="step", lw=3, color="steelblue")
	ax_yhist.legend()
	

	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	#ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="0.5", alpha=0.3, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="tomato", lw=3, label=L"\rm 1D")

	
	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP_corr[:, i]] .+ 
		parameters_all["feh"][bins_CEMP_corr[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP_corr[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP_corr[:, i]] .+ 
		parameters_all["feh"][bins_CEMP_corr[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP_corr[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	#ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="cyan", alpha=0.3, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="steelblue", label=L"\rm 3D", lw=3)
	
	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls="--", lw=3)

	

	ax.set_xlim(-6.9, -1.1)
	ax.set_ylim(3.6, 10.9)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(loc="upper left")

	# Hide labels on shared axes
	plt.setp(ax_xhist.get_xticklabels(), visible=false)
	plt.setp(ax_yhist.get_yticklabels(), visible=false)

	f.savefig("yoon_beers_detailed.pdf")

	f
end

# ╔═╡ 716d2ed3-817e-4909-92ad-a83d7b79c171
let
	plt.close()
	f, ax = plt.subplots(2, 1, figsize=(5, 9), sharex=true, sharey=true)
	plt.subplots_adjust(wspace=0, hspace=0)

	ellipse_alpha = 0.15
	scatter_alpha = 0.6
	
	# 1D
	ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.9), -4.5, 0.9, color="orange", alpha=ellipse_alpha, ec="0.5"
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.75, 5.95), -3.5, 0.9, color="green", alpha=ellipse_alpha, angle=45, ec="0.5"
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.52, 7.9), -2, 3, color="blue", alpha=ellipse_alpha, angle=-55, ec="0.5"
	)
	ax[0].add_patch(ellipse1)
	ax[0].add_patch(ellipse2)
	ax[0].add_patch(ellipse3)
	draw_ellipse_axes(ellipse1, ax[0], N=1, skip_horizontal=true, color="blue", ls="--", alpha=0.75)
	draw_ellipse_axes(ellipse2, ax[0], N=1, skip_vertical=true, color="green", ls="--", alpha=0.75)
	draw_ellipse_axes(ellipse3, ax[0], N=1, skip_vertical=true, color="orange", ls="--", alpha=0.75)

	
	# 3D
	ellipse3 = matplotlib.patches.Ellipse(
		(-5.3, 6.3), -4.5, 0.9, color="orange", alpha=ellipse_alpha, ec="0.5", angle=11
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.75, 5.55), -3.5, 0.9, color="green", alpha=ellipse_alpha, angle=45, ec="0.5"
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.48, 7.75), -2, 3, color="blue", alpha=ellipse_alpha, angle=-55, ec="0.5"
	)
	ax[1].add_patch(ellipse1)
	ax[1].add_patch(ellipse2)
	ax[1].add_patch(ellipse3)
	draw_ellipse_axes(ellipse1, ax[1], N=1, skip_horizontal=true, color="blue", ls="--", alpha=0.75)
	draw_ellipse_axes(ellipse2, ax[1], N=1, skip_vertical=true, color="green", ls="--", alpha=0.75)
	draw_ellipse_axes(ellipse3, ax[1], N=1, skip_vertical=true, color="orange", ls="--", alpha=0.75)
	

	#=ax[0].axhline(6.9, color="orange", ls=":", alpha=0.5)
	ax[1].axhline(6.4, color="orange", ls=":", alpha=0.5)

	ax[0].axhline(5.95, color="green", ls=":", alpha=0.5)
	ax[1].axhline(5.55, color="green", ls=":", alpha=0.5)

	ax[0].axhline(7.9, color="blue", ls=":", alpha=0.5)
	ax[1].axhline(7.75, color="blue", ls=":", alpha=0.5)
	=#

	ax[0].text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="orange", fontsize="small")
	ax[0].text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="green", fontsize="small")
	ax[0].text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="blue", fontsize="small", ha="right")
	ax[1].text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="orange", fontsize="small")
	ax[1].text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="green", fontsize="small")
	ax[1].text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="blue", fontsize="small", ha="right")
	
	ax[0].text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax[0].text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")
	ax[1].text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax[1].text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")

	text_group3 = """
	$(L"\rm \mathbf{CEMP-no}")
	Faint SNe
	Pop III
	"""
	#ax[0].text(-6.7, 8.2, text_group3, color="0.5", fontsize="x-small", ha="left", va="bottom")
	#ax[1].text(-6.7, 8.2, text_group3, color="0.5", fontsize="x-small", ha="left", va="bottom")

	
	text_group2 = """
	$(L"\rm \mathbf{CEMP-no}")
	Standard CCSNe
	Pop II
	"""
	#ax[0].text(-3.1, 5.15, text_group2, color="0.5", fontsize="x-small", ha="left", va="top")
	#ax[1].text(-3.1, 5.15, text_group2, color="0.5", fontsize="x-small", ha="left", va="top")


	text_group1 = """
	$(L"\rm \mathbf{CEMP-s\ +\ CEMP-r/s}")
	MRSNe, NS-NS or NS-BH merger
	AGB contribution
	"""
	#ax[0].text(-1.4, 9.5, text_group1, color="0.5", fontsize="x-small", ha="right", va="bottom")
	#ax[1].text(-1.4, 9.5, text_group1, color="0.5", fontsize="x-small", ha="right", va="bottom")
	

	
	CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe .+ feh .+ 8.560
	#ax.scatter(feh, c, s=20, marker="s", c="None", alpha=0.2, linewidth=1, edgecolor="k")
	ax[0].scatter(feh, c, s=15, marker="o", c="None", alpha=scatter_alpha, edgecolor="tomato", rasterized=true, label="1D")
	ax[1].scatter(feh, c+corrections_saga_all[CEMP_mask], s=15, marker="o", c="None", alpha=scatter_alpha, edgecolor="steelblue", rasterized=true, label="3D")

	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax[0].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)
	ax[1].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	#ax[0].legend()
	#ax[1].legend()
	ax[0].text(0.05, 0.87, L"\mathbf{1D}", color="tomato", transform=ax[0].transAxes, fontsize="x-large")
	ax[1].text(0.05, 0.87, L"\mathbf{3D}", color="steelblue", transform=ax[1].transAxes, fontsize="x-large")
	ax[0].set_xlim(-6.9, -1.1)
	ax[0].set_ylim(3.6, 10.1)
	ax[1].set_xlabel("[Fe/H]")
	#ax[0].set_xlabel("[Fe/H]")
	ax[0].set_ylabel("A(C)")
	ax[1].set_ylabel("A(C)")
	#ax[0].legend(loc="upper left")

	
	f.savefig("yoon_beers_split.pdf")

	f
end

# ╔═╡ a1b28a53-a3b5-4e87-bc2d-5ce112c1f280
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5), layout="tight")
	plt.subplots_adjust(wspace=0, hspace=0)

	ellipse_alpha = 0.12
	scatter_alpha = 0.6

	ellipse1_test = matplotlib.patches.Ellipse(
		means(g1_yoon)[1,:], el1_yoon[1]*1.2, el1_yoon[2]*1.2, angle= el1_yoon[3], color="steelblue", alpha=ellipse_alpha*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse1_test)

	ellipse2_test = matplotlib.patches.Ellipse(
		means(g2_yoon)[1,:], el2_yoon[1], el2_yoon[2], angle= el2_yoon[3], color="yellowgreen", alpha=ellipse_alpha*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse2_test)

	ellipse3_test = matplotlib.patches.Ellipse(
		means(g3_yoon)[1,:], el3_yoon[1], el3_yoon[2], angle= el3_yoon[3], color="tomato", alpha=ellipse_alpha*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse3_test)
	draw_ellipse_axes(ellipse1_test, ax, N=1, skip_vertical=true, color="steelblue", ls="--", alpha=1, lw=4)
	draw_ellipse_axes(ellipse2_test, ax, N=1, skip_vertical=true, color="yellowgreen", ls="--", alpha=1, lw=4)
	draw_ellipse_axes(ellipse3_test, ax, N=1, skip_vertical=true, color="tomato", ls="--", alpha=1, lw=4)

	#=ax.scatter(x_1D_yoon[g1_mask_1D, 1], x_1D_yoon[g1_mask_1D, 2], marker="x", color="blue", s=10)
	ax.scatter(x_1D_yoon[g2_mask_1D, 1], x_1D_yoon[g2_mask_1D, 2], marker="s", color="green", s=10)
	ax.scatter(x_1D_yoon[g3_mask_1D, 1], x_1D_yoon[g3_mask_1D, 2], marker="v", color="orange", s=30)=#

	
	
	
	ellipse1_test_3D = matplotlib.patches.Ellipse(
		means(g1_yoon_3D)[1,:], el1_yoon_3D[1]*1.2, el1_yoon_3D[2]*1.2, angle= el1_yoon_3D[3], color="steelblue", alpha=ellipse_alpha, ec="0.5"
	)
	#ax.add_patch(ellipse1_test_3D)

	ellipse2_test_3D = matplotlib.patches.Ellipse(
		means(g2_yoon_3D)[1,:], el2_yoon_3D[1], el2_yoon_3D[2], angle= el2_yoon_3D[3], color="yellowgreen", alpha=ellipse_alpha, ec="0.5"
	)
	#ax.add_patch(ellipse2_test_3D)

	ellipse3_test_3D = matplotlib.patches.Ellipse(
		means(g3_yoon_3D)[1,:], el3_yoon_3D[1], el3_yoon_3D[2], angle= el3_yoon_3D[3], color="tomato", alpha=ellipse_alpha, ec="0.5"
	)
	#ax.add_patch(ellipse3_test_3D)
	draw_ellipse_axes(ellipse1_test_3D, ax, N=1, skip_vertical=true, color="steelblue", ls="-", alpha=1, lw=4)
	draw_ellipse_axes(ellipse2_test_3D, ax, N=1, skip_vertical=true, color="yellowgreen", ls="-", alpha=1, lw=4)
	draw_ellipse_axes(ellipse3_test_3D, ax, N=1, skip_vertical=true, color="tomato", ls="-", alpha=1, lw=4)

	#=ax.scatter(x_3D_yoon[g1_mask_1D, 1], x_3D_yoon[g1_mask_1D, 2], marker=".", color="blue", s=5, alpha=0.2)
	ax.scatter(x_3D_yoon[g2_mask_1D, 1], x_3D_yoon[g2_mask_1D, 2], marker=".", color="green", s=5, alpha=0.2)
	ax.scatter(x_3D_yoon[g3_mask_1D, 1], x_3D_yoon[g3_mask_1D, 2], marker=".", color="tomato", s=5, alpha=0.2)=#

	#ax.scatter(x_3D_yoon[:, 1], x_3D_yoon[:, 2], marker=".", color="k", s=5, alpha=0.2)
	#ax.hist2d(x_3D_yoon[:, 1], x_3D_yoon[:, 2], bins=30, cmap="Greys", vmin=1)
	
	my_cmap = plt.cm.Blues
	my_cmap.set_under("w",1)
	ax.hist2d(x_3D_yoon[g1_mask_1D, 1], x_3D_yoon[g1_mask_1D, 2], bins=15, cmap=my_cmap, vmin=1, alpha=0.4, zorder=0)
	my_cmap = plt.cm.Greens
	my_cmap.set_under("w",1)
	ax.hist2d(x_3D_yoon[g2_mask_1D, 1], x_3D_yoon[g2_mask_1D, 2], bins=15, cmap=my_cmap, vmin=1, alpha=0.4, zorder=0)
	my_cmap = plt.cm.Reds
	my_cmap.set_under("w",1)
	ax.hist2d(x_3D_yoon[g3_mask_1D, 1], x_3D_yoon[g3_mask_1D, 2], bins=7, cmap=my_cmap, vmin=1, alpha=0.5, zorder=0)


	#=ellipse1_test_3D = matplotlib.patches.Ellipse(
		means(g123_yoon_3D)[1,:], el1_yoon123_3D[1]*1.2, el1_yoon123_3D[2]*1.2, angle= el1_yoon123_3D[3], color="tomato", alpha=ellipse_alpha, ec="0.5"
	)
	ax.add_patch(ellipse1_test_3D)

	ellipse2_test_3D = matplotlib.patches.Ellipse(
		means(g123_yoon_3D)[2,:], el2_yoon123_3D[1], el2_yoon123_3D[2], angle= el2_yoon123_3D[3], color="steelblue", alpha=ellipse_alpha, ec="0.5"
	)
	ax.add_patch(ellipse2_test_3D)

	ellipse3_test_3D = matplotlib.patches.Ellipse(
		means(g123_yoon_3D)[3,:], el3_yoon123_3D[1], el3_yoon123_3D[2], angle= el3_yoon123_3D[3], color="yellowgreen", alpha=ellipse_alpha, ec="0.5"
	)
	ax.add_patch(ellipse3_test_3D)
	draw_ellipse_axes(ellipse1_test_3D, ax, N=1, skip_vertical=true, color="tomato", ls="-", alpha=1)
	draw_ellipse_axes(ellipse2_test_3D, ax, N=1, skip_vertical=true, color="steelblue", ls="-", alpha=1)
	draw_ellipse_axes(ellipse3_test_3D, ax, N=1, skip_vertical=true, color="yellowgreen", ls="-", alpha=1)=#
	
	

	ax.text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="tomato", fontsize="small")
	ax.text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="yellowgreen", fontsize="small")
	ax.text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="steelblue", fontsize="small", ha="right")
	ax.text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="tomato", fontsize="small")
	ax.text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="yellowgreen", fontsize="small")
	ax.text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="steelblue", fontsize="small", ha="right")
	
	ax.text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax.text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")
	ax.text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax.text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")

	text_group3 = """
	$(L"\rm \mathbf{CEMP-no}")
	Faint SNe
	Pop III
	"""
	
	text_group2 = """
	$(L"\rm \mathbf{CEMP-no}")
	Standard CCSNe
	Pop II
	"""

	text_group1 = """
	$(L"\rm \mathbf{CEMP-s\ +\ CEMP-r/s}")
	MRSNe, NS-NS or NS-BH merger
	AGB contribution
	"""
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax.plot([],[],ls="--", color="k",label=L"\rm 1D")
	ax.plot([],[],ls="-", color="k",label=L"\rm 3D")

	#ax.text(0.05, 0.87, L"\mathbf{1D}", color="tomato", transform=ax.transAxes, fontsize="x-large")
	#ax.text(0.05, 0.87, L"\mathbf{3D}", color="steelblue", transform=ax.transAxes, fontsize="x-large")
	ax.set_xlim(-6.9, -1.1)
	ax.set_ylim(3.6, 10.1)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(labelspacing=0)
	
	#f.savefig("yoon_beers_split.pdf")

	f
end

# ╔═╡ fac3be33-656b-4c3c-9fb4-eae4dc47d71b
begin
	lowfeh_mask = selection_mask .& (parameters_all["feh"] .<= -1)
	CEMP_mask_gmm = (parameters_all["cfe"] .>= cfe_limit) .& lowfeh_mask
	CEMP_mask_corr_gmm = (parameters_all["cfe"] .+ corrections_saga_all .>= cfe_limit) .& lowfeh_mask
	#CEMP_mask_corr_gmm = CEMP_mask_gmm
	
	feh_gmm = parameters_all["feh"][CEMP_mask_gmm]
	cfe_gmm = parameters_all["cfe"][CEMP_mask_gmm]
	c_gmm = cfe_gmm .+ feh_gmm .+ 8.560
	x_1D_gmm = zeros(length(feh_gmm), 2)

	x_1D_gmm[:, 1] .= feh_gmm
	x_1D_gmm[:, 2] .= c_gmm 


	feh_corr_gmm = parameters_all["feh"][CEMP_mask_corr_gmm]
	cfe_corr_gmm = parameters_all["cfe"][CEMP_mask_corr_gmm]
	c_corr_gmm = cfe_corr_gmm .+ feh_corr_gmm .+ 8.560
	x_3D_gmm = zeros(length(feh_corr_gmm), 2)
	x_3D_gmm[:, 1] .= feh_corr_gmm 
	x_3D_gmm[:, 2] .= c_corr_gmm .+ corrections_saga_all[CEMP_mask_corr_gmm] 
	
	gmm3_1D = GMM(3, x_1D_gmm, nIter=2500, kind=:full)
	gmm3_3D = GMM(3, x_3D_gmm, nIter=2500, kind=:full)

	# test with 1D fit only
	gmm2_1D = GMM(2, x_1D_gmm, nIter=2500, kind=:full)
	gmm2_3D = GMM(2, x_3D_gmm, nIter=2500, kind=:full)

	gmm2_1D_lin = GMM(2, x_1D_gmm[:, 2], nIter=2500, nInit=2500)
	gmm2_3D_lin = GMM(2, x_3D_gmm[:, 2], nIter=2500, nInit=2500)
end

# ╔═╡ b6020aa9-c570-4cef-b41d-629eff19f542
function confidence_ellipse!(μ, Σ, ax, n_std=1.0; kwargs...)
    cov = Σ
    pearson = cov[1, 2]/sqrt(cov[1, 1] * cov[2, 2])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = sqrt.(1 + pearson)
    ell_radius_y = sqrt.(1 - pearson)
    ellipse =  matplotlib.patches.Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2; kwargs...)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = sqrt(cov[1, 1]) * n_std
    mean_x = μ[1]

    # calculating the standard deviation of y ...
    scale_y = sqrt(cov[2, 2]) * n_std
    mean_y = μ[2]

    transf = matplotlib.transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    ax.add_patch(ellipse)
end

# ╔═╡ 0bd75590-be29-48da-ab79-92a9be005d29
function confidence_ellipse_data!(x::AbstractVector, y::AbstractVector, ax, n_std::Real=1.0; kwargs...)
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))

    # Means
    μx = mean(x)
    μy = mean(y)

    # 2×2 covariance matrix (sample covariance)
    cov_mat = cov(hcat(x, y))     # Statistics.cov on an n×2 matrix -> 2×2 matrix

    # Eigen-decomposition of symmetric covariance
    eig = eigen(Symmetric(cov_mat))
    vals = eig.values
    vecs = eig.vectors

    # Re-order so largest eigenvalue is first (major axis)
    order = sortperm(vals; rev=true)
    vals = vals[order]
    vecs = vecs[:, order]

    # Width/height (full lengths) and rotation angle (degrees)
    width  = 2 * n_std * sqrt(vals[1])
    height = 2 * n_std * sqrt(vals[2])
    # angle: arctan2(y, x) of the principal eigenvector -> convert to degrees
    angle = atan(vecs[2, 1], vecs[1, 1]) * (180 / π)

    # Create and add the ellipse via matplotlib.patches.Ellipse
    ellipse = matplotlib.patches.Ellipse((Float64(μx), Float64(μy)),
                               Float64(width),
                               Float64(height),
                               angle=Float64(angle);
                               kwargs...)
    ax.add_patch(ellipse)

    return ellipse
end

# ╔═╡ 38efffde-0ebb-4beb-97c9-8a639761a6fd
function rescale_range(arr, new_min, new_max)
    min_val = minimum(arr)
    max_val = maximum(arr)

    # Handle the case where all elements are the same (to avoid division by zero)
    if min_val == max_val
        # Return an array where every element is the new minimum
        return fill(Float64(new_min), size(arr))
    end

    # Apply the generalized min-max scaling formula
    # Formula: new_min + ( (value - min_val) * (new_max - new_min) ) / (max_val - min_val)
    return new_min .+ ((arr .- min_val) .* (new_max - new_min)) ./ (max_val - min_val)
end

# ╔═╡ 6fe595f7-e8d0-42b9-aab4-752151128bcf
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=true, sharey=false, layout="tight")

	sigdigs=3
	offset = 0.05

	c_highC = "steelblue"
	c_lowC = "tomato"
	gmm2_colors = ["k", "k"]
	
	gmm2_1D_means = means(gmm2_1D_lin)
	gmm2_1D_cov = covars(gmm2_1D_lin)
	gmm2_1D_ll = gmmposterior(gmm2_1D_lin, x_1D_gmm[:, 2:2])[1]
	gmm2_1D_ass = [argmax(gmm2_1D_ll[i, :]) for i in axes(gmm2_1D_ll, 1)]
	i_highC = argmax(gmm2_1D_means[:, 1])
	i_lowC = argmin(gmm2_1D_means[:, 1])
	gmm2_colors[i_highC] = c_highC
	gmm2_colors[i_lowC] = c_lowC
	gmm2_1D_ll_diff = 1 ./ abs.(exp.(gmm2_1D_ll[:, 2]) .- exp.(gmm2_1D_ll[:, 1])) .* 10
	ax[0].scatter(x_1D_gmm[:, 1], x_1D_gmm[:, 2], s=25, marker="o", c=[gmm2_colors[a] for a in gmm2_1D_ass], alpha=0.5, rasterized=true, label="1D")

	#=confidence_ellipse_data!(x_1D_gmm[gmm2_1D_ass .== i_highC, 1], x_1D_gmm[gmm2_1D_ass .== i_highC, 2], ax[0], 3, facecolor="none", ec=gmm2_colors[i_highC], lw=3)
	confidence_ellipse_data!(x_1D_gmm[gmm2_1D_ass .== i_lowC, 1], x_1D_gmm[gmm2_1D_ass .== i_lowC, 2], ax[0], 3, facecolor="none", ec=gmm2_colors[i_lowC], lw=3)=#

	gmm2_3D_means = means(gmm2_3D_lin)
	gmm2_3D_cov = covars(gmm2_3D_lin)
	gmm2_3D_ll = gmmposterior(gmm2_3D_lin, x_3D_gmm[:, 2:2])[1]
	gmm2_3D_ass = [argmax(gmm2_3D_ll[i, :]) for i in axes(gmm2_3D_ll, 1)]
	i_highC = argmax(gmm2_3D_means[:, 1])
	i_lowC = argmin(gmm2_3D_means[:, 1])
	gmm2_colors[i_highC] = c_highC
	gmm2_colors[i_lowC] = c_lowC
	gmm2_3D_ll_diff = 1 ./ abs.(exp.(gmm2_3D_ll[:, 2]) .- exp.(gmm2_3D_ll[:, 1])) .* 10
	ax[1].scatter(x_3D_gmm[:, 1], x_3D_gmm[:, 2], s=25, marker="o", c=[gmm2_colors[a] for a in gmm2_3D_ass], alpha=0.5, rasterized=true, label="3D")
	
	#=confidence_ellipse_data!(x_3D_gmm[gmm2_3D_ass .== i_highC, 1], x_3D_gmm[gmm2_3D_ass .== i_highC, 2], ax[1], 3, facecolor="none", ec=gmm2_colors[i_highC], lw=3)
	confidence_ellipse_data!(x_3D_gmm[gmm2_3D_ass .== i_lowC, 1], x_3D_gmm[gmm2_3D_ass .== i_lowC, 2], ax[1], 3, facecolor="none", ec=gmm2_colors[i_lowC], lw=3)=#
	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-6.9, -0.9, length=100)|>collect

	gmm2_1D_lin_means = means(gmm2_1D_lin)
	gmm2_3D_lin_means = means(gmm2_3D_lin)
	i_highC = argmax(gmm2_1D_means[:, 1])
	i_lowC = argmin(gmm2_1D_means[:, 1])
	gmm2_colors[i_highC] = c_highC
	gmm2_colors[i_lowC] = c_lowC
	ax[0].fill_between(
		x, 
		gmm2_1D_lin_means[i_highC,1]-sqrt(gmm2_1D_cov[i_highC,1]),
		gmm2_1D_lin_means[i_highC,1]+sqrt(gmm2_1D_cov[i_highC,1]),
		color=gmm2_colors[i_highC], alpha=0.1
	)
	ax[0].fill_between(
		x, 
		gmm2_1D_lin_means[i_lowC,1]-sqrt(gmm2_1D_cov[i_lowC,1]),
		gmm2_1D_lin_means[i_lowC,1]+sqrt(gmm2_1D_cov[i_lowC,1]),
		color=gmm2_colors[i_lowC], alpha=0.1
	)
	ax[0].axhline(gmm2_1D_lin_means[i_highC,1], ls="--", color=gmm2_colors[i_highC], lw=2)
	ax[0].axhline(gmm2_1D_lin_means[i_lowC,1], ls="--", color=gmm2_colors[i_lowC], lw=2)
	#ax[0].axhline((gmm2_1D_lin_means[1,1]+gmm2_1D_lin_means[2,1])/2, ls="-", lw=2, alpha=0.4, color="tomato")
	lb = L"\rm A(C) = "*"$(round(gmm2_1D_lin_means[i_lowC,1], sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(sqrt(gmm2_1D_cov[i_lowC,1]), digits=2)))"
	ax[0].text(-6.7, gmm2_1D_lin_means[i_lowC,1] -0.1, lb, ha="left", va="top", color="tomato", fontsize="small")
	lb = L"\rm A(C) = "*"$(round(gmm2_1D_lin_means[i_highC,1], sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(sqrt(gmm2_1D_cov[i_highC,1]), digits=2)))"
	ax[0].text(-6.7, gmm2_1D_lin_means[i_highC,1] +0.1, lb, ha="left", va="bottom", color="steelblue", fontsize="small")
	
	
	i_highC = argmax(gmm2_3D_means[:, 1])
	i_lowC = argmin(gmm2_3D_means[:, 1])
	gmm2_colors[i_highC] = c_highC
	gmm2_colors[i_lowC] = c_lowC
	ax[1].fill_between(
		x, 
		gmm2_3D_lin_means[i_highC,1]-sqrt(gmm2_3D_cov[i_highC,1]),
		gmm2_3D_lin_means[i_highC,1]+sqrt(gmm2_3D_cov[i_highC,1]),
		color=gmm2_colors[i_highC], alpha=0.1
	)
	ax[1].fill_between(
		x, 
		gmm2_3D_lin_means[i_lowC,1]-sqrt(gmm2_3D_cov[i_lowC,1]),
		gmm2_3D_lin_means[i_lowC,1]+sqrt(gmm2_3D_cov[i_lowC,1]),
		color=gmm2_colors[i_lowC], alpha=0.1
	)

	ax[1].axhline(gmm2_3D_lin_means[i_highC,1], ls="--", color=gmm2_colors[i_highC], lw=2)
	ax[1].axhline(gmm2_3D_lin_means[i_lowC,1], ls="--", color=gmm2_colors[i_lowC], lw=2)
	#ax[1].axhline((gmm2_3D_lin_means[1,1]+gmm2_3D_lin_means[2,1])/2, ls="-", lw=2, alpha=0.4, color="tomato")
	lb = L"\rm A(C) = "*"$(round(gmm2_3D_lin_means[i_lowC,1], sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(sqrt(gmm2_1D_cov[i_lowC,1]), digits=2)))"
	ax[1].text(-6.7, gmm2_3D_lin_means[i_lowC,1] -0.25, lb, ha="left", va="top", color="tomato", fontsize="small")
	lb = L"\rm A(C) = "*"$(round(gmm2_3D_lin_means[i_highC,1], sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(sqrt(gmm2_1D_cov[i_highC,1]), digits=2)))"
	ax[1].text(-6.7, gmm2_3D_lin_means[i_highC,1] +0.1, lb, ha="left", va="bottom", color="steelblue", fontsize="small")
	

	ax[0].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)
	ax[1].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax[0].text(0.05, 0.87, L"\mathbf{1D}", color="k", transform=ax[0].transAxes, fontsize="x-large")
	ax[1].text(0.05, 0.87, L"\mathbf{3D}", color="k", transform=ax[1].transAxes, fontsize="x-large")

	ax[0].set_xlim(x[1], -0.9)
	ax[0].set_ylim(4.2, 9.5)
	ax[1].set_xlim(x[1], -0.9)
	ax[1].set_ylim(4.2, 9.5)
	
	ax[0].set_ylabel(L"\rm A(C)")
	ax[1].set_ylabel(L"\rm A(C)")
	#ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_xlabel(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ bfb45894-1667-4582-896a-63597c6e82ab
md"Modified Assignment: With the assumption that the high-C band does not extend to low metallicities"

# ╔═╡ 2065ce06-76de-40c7-bd8f-e9beaa4cc02d
begin
	#mask_1D_highm = (x_1D_gmm[:,1] .> -3.5)
	#mask_3D_highm = (x_3D_gmm[:,1] .> -3.5)
	mask_1D_highm = (x_1D_gmm[:,1] .> -4.5)
	mask_3D_highm = (x_3D_gmm[:,1] .> -4.5) .| (x_3D_gmm[:,1] .> 6)

	Niter = 2500
	
	# split the high metallicity part in 2
	gmm_1D_group1 = GMM(2, x_1D_gmm[mask_1D_highm, 2], nIter=Niter, nInit=Niter)
	gmm_3D_group1 = GMM(2, x_3D_gmm[mask_3D_highm, 2], nIter=Niter, nInit=Niter)

	# select goup 1 targets based on this assignment
	imax_1D_group1 = argmax(means(gmm_1D_group1)[:, 1])
	prob_1D_group1 = gmmposterior(gmm_1D_group1, x_1D_gmm[:, 2:2])[1]
	mask_1D_group1 = [(argmax(prob_1D_group1[i,:]) == imax_1D_group1) .& mask_1D_highm[i]  for i in axes(prob_1D_group1, 1)]

	imax_3D_group1 = argmax(means(gmm_3D_group1)[:, 1])
	prob_3D_group1 = gmmposterior(gmm_3D_group1, x_3D_gmm[:, 2:2])[1]
	mask_3D_group1 = [(argmax(prob_3D_group1[i,:]) == imax_3D_group1) .& mask_3D_highm[i] for i in axes(prob_3D_group1, 1)]

	# use the remaining points to find the other groups
	gmm_1D_group23 = GMM(2, x_1D_gmm[.!mask_1D_group1, :], nIter=Niter, kind=:full,  nInit=Niter)
	gmm_3D_group23 = GMM(2, x_3D_gmm[.!mask_3D_group1, :], nIter=Niter, kind=:full,  nInit=Niter)

	# select goup 2 and 3 targets based on this assignment
	imax_1D_group2 = argmax(means(gmm_1D_group23)[:, 1])
	prob_1D_group23 = gmmposterior(gmm_1D_group23, x_1D_gmm)[1]
	mask_1D_group2 = [(argmax(prob_1D_group23[i,:]) == imax_1D_group2) .& (.!mask_1D_group1[i]) for i in axes(prob_1D_group1, 1)]
	mask_1D_group3 = [(argmax(prob_1D_group23[i,:]) != imax_1D_group2) .& (.!mask_1D_group1[i]) for i in axes(prob_1D_group1, 1)]

	imax_3D_group2 = argmax(means(gmm_3D_group23)[:, 1])
	prob_3D_group23 = gmmposterior(gmm_3D_group23, x_3D_gmm)[1]
	mask_3D_group2 = [(argmax(prob_3D_group23[i,:]) == imax_3D_group2) .& (.!mask_3D_group1[i]) for i in axes(prob_3D_group1, 1)]
	mask_3D_group3 = [(argmax(prob_3D_group23[i,:]) != imax_3D_group2) .& (.!mask_3D_group1[i]) for i in axes(prob_3D_group1, 1)]

	# alternatively fit 3 as it is
	gmm_1D_group123 = GMM(3, x_1D_gmm, nIter=Niter, kind=:full)
	gmm_3D_group123 = GMM(3, x_3D_gmm, nIter=Niter, kind=:full)

	# or split on C/Fe alone
	gmm_1D_calone = GMM(2, x_1D_gmm[:, 2], nIter=Niter, nInit=Niter)
	gmm_3D_calone = GMM(2, x_3D_gmm[:, 2], nIter=Niter, nInit=Niter)

	imax_1D_highC = argmax(means(gmm_1D_calone)[:, 1])
	prob_1D_calone = gmmposterior(gmm_1D_calone, x_1D_gmm[:, 2:2])[1]
	mask_1D_highC = [(argmax(prob_1D_calone[i,:]) == imax_1D_highC) for i in axes(prob_1D_calone, 1)]
	mask_1D_lowC = [(argmax(prob_1D_calone[i,:]) != imax_1D_highC) for i in axes(prob_1D_calone, 1)]

	imax_3D_highC = argmax(means(gmm_3D_calone)[:, 1])
	prob_3D_calone = gmmposterior(gmm_3D_calone, x_3D_gmm[:, 2:2])[1]
	mask_3D_highC = [(argmax(prob_3D_calone[i,:]) == imax_3D_highC) for i in axes(prob_3D_calone, 1)]
	mask_3D_lowC = [(argmax(prob_3D_calone[i,:]) != imax_3D_highC) for i in axes(prob_3D_calone, 1)]
end

# ╔═╡ d687f757-226e-4571-bb54-d6a69de69b89
get_assignment_mask(g, x, group) = begin
	prob = gmmposterior(g, x)[1]
	mask = [(argmax(prob[i, :]) == group) for i in axes(prob, 1)]
end

# ╔═╡ c20c8129-bdbe-4bb8-a7a0-36a23d8bc1ba
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(9, 5), sharex=true, sharey=true, layout="tight")

	# 1D
	ax[0].scatter(x_1D_gmm[mask_1D_group3, 1], x_1D_gmm[mask_1D_group3, 2], color="orange", alpha=0.5, s=15)
	confidence_ellipse_data!(x_1D_gmm[mask_1D_group3, 1], x_1D_gmm[mask_1D_group3, 2], ax[0], 2, facecolor="orange", ec="orange", lw=4, alpha=0.3)
	
	ax[0].scatter(x_1D_gmm[mask_1D_group2, 1], x_1D_gmm[mask_1D_group2, 2], color="green", alpha=0.5, s=15)
	confidence_ellipse_data!(x_1D_gmm[mask_1D_group2, 1], x_1D_gmm[mask_1D_group2, 2], ax[0], 2, facecolor="green", ec="green", lw=4, alpha=0.3)
	
	ax[0].scatter(x_1D_gmm[mask_1D_group1, 1], x_1D_gmm[mask_1D_group1, 2], color="blue", alpha=0.5, s=15)
	confidence_ellipse_data!(x_1D_gmm[mask_1D_group1, 1], x_1D_gmm[mask_1D_group1, 2], ax[0], 2, facecolor="blue", ec="blue", lw=4, alpha=0.3)

	#=ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.9), -4.5, 0.9, color="orange", alpha=0.3, 
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.75, 5.95), -3.5, 0.9, color="green", alpha=0.3, angle=45, 
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.52, 7.9), -2, 3, color="blue", alpha=0.3, angle=-55, 
	)
	ax[0].add_patch(ellipse1)
	ax[0].add_patch(ellipse2)
	ax[0].add_patch(ellipse3)=#
	#draw_ellipse_axes(ellipse1, ax[0], N=1, skip_horizontal=true, color="blue", ls="--", alpha=0.75)
	#draw_ellipse_axes(ellipse2, ax[0], N=1, skip_vertical=true, color="green", ls="--", alpha=0.75)
	#draw_ellipse_axes(ellipse3, ax[0], N=1, skip_vertical=true, color="orange", ls="--", alpha=0.75)
	

	# 3D
	ax[1].scatter(x_3D_gmm[mask_3D_group3, 1], x_3D_gmm[mask_3D_group3, 2], color="orange", alpha=0.5, s=15)
	confidence_ellipse_data!(x_3D_gmm[mask_3D_group3, 1], x_3D_gmm[mask_3D_group3, 2], ax[1], 2, facecolor="orange", ec="orange", lw=4, alpha=0.3)
	
	ax[1].scatter(x_3D_gmm[mask_3D_group2, 1], x_3D_gmm[mask_3D_group2, 2], color="green", alpha=0.5, s=15)
	confidence_ellipse_data!(x_3D_gmm[mask_3D_group2, 1], x_3D_gmm[mask_3D_group2, 2], ax[1], 2, facecolor="green", ec="green", lw=4, alpha=0.3)
	
	ax[1].scatter(x_3D_gmm[mask_3D_group1, 1], x_3D_gmm[mask_3D_group1, 2], color="blue", alpha=0.5, s=15)
	confidence_ellipse_data!(x_3D_gmm[mask_3D_group1, 1], x_3D_gmm[mask_3D_group1, 2], ax[1], 2, facecolor="blue", ec="blue", lw=4, alpha=0.3)
	
	#=
	ax[1].scatter(x_3D_gmm[mask_3D_highC, 1], x_3D_gmm[mask_3D_highC, 2], color="k", alpha=0.5, s=15)
	x_f = range(-6.9, 0, length=100)|>collect
	ax[1].fill_between(
		x_f, 
		mean(x_3D_gmm[mask_3D_highC, 2]) - MUST.std(x_3D_gmm[mask_3D_highC, 2]),
		mean(x_3D_gmm[mask_3D_highC, 2]) + MUST.std(x_3D_gmm[mask_3D_highC, 2]),
		color="k", alpha=0.1
	)
	ax[1].axhline(mean(x_3D_gmm[mask_3D_highC, 2]), ls="--", color="k", lw=2)
	
	
	ax[1].scatter(x_3D_gmm[mask_3D_lowC, 1], x_3D_gmm[mask_3D_lowC, 2], color="r", alpha=0.5, s=15)
	ax[1].fill_between(
		x_f, 
		mean(x_3D_gmm[mask_3D_lowC, 2]) - MUST.std(x_3D_gmm[mask_3D_lowC, 2]),
		mean(x_3D_gmm[mask_3D_lowC, 2]) + MUST.std(x_3D_gmm[mask_3D_lowC, 2]),
		color="r", alpha=0.1
	)
	ax[1].axhline(mean(x_3D_gmm[mask_3D_lowC, 2]), ls="--", color="r", lw=2)
	=#

	ax[0].text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="orange", fontsize="small")
	ax[0].text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="green", fontsize="small")
	ax[0].text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="blue", fontsize="small", ha="right")
	ax[1].text(-6.7, 7.5, L"\rm \mathbf{Group\ III}", color="orange", fontsize="small")
	ax[1].text(-3.1, 5.3, L"\rm \mathbf{Group\ II}", color="green", fontsize="small")
	ax[1].text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="blue", fontsize="small", ha="right")
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-6.9, 0, length=100)|>collect
	ax[0].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)
	ax[1].plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax[0].text(0.05, 0.87, L"\mathbf{1D}", color="tomato", transform=ax[0].transAxes, fontsize="x-large")
	ax[1].text(0.05, 0.87, L"\mathbf{3D}", color="steelblue", transform=ax[1].transAxes, fontsize="x-large")

	#ax.set_ylim(6, 8.5)
	ax[0].set_xlim(-6.9, -0.9)
	ax[0].set_ylim(4.5, 10.2)
	
	ax[0].set_ylabel(L"\rm A(C)")
	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ 2cd3aede-d981-4974-ab31-63f11a76b377
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=true, sharey=true, layout="tight")

	# 1D from Yoon directly
	@show "G1 1D"
	ellipse1_test = matplotlib.patches.Ellipse(
		means(g1_yoon)[1,:], el1_yoon[1]*1.2, el1_yoon[2]*1.2, angle= el1_yoon[3], color="steelblue", alpha=0*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse1_test)
	@show "G2 1D"
	ellipse2_test = matplotlib.patches.Ellipse(
		means(g2_yoon)[1,:], el2_yoon[1], el2_yoon[2], angle= el2_yoon[3], color="yellowgreen", alpha=0*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse2_test)
	@show "G3 1D"
	ellipse3_test = matplotlib.patches.Ellipse(
		means(g3_yoon)[1,:], el3_yoon[1], el3_yoon[2], angle= el3_yoon[3], color="tomato", alpha=0*0.5, ec="0.5"
	)
	#ax.add_patch(ellipse3_test)
	draw_ellipse_axes(ellipse1_test, ax, N=1, skip_vertical=true, color="steelblue", ls="--", alpha=1, lw=3)
	draw_ellipse_axes(ellipse2_test, ax, N=1, skip_vertical=true, color="yellowgreen", ls="--", alpha=1, lw=3)
	draw_ellipse_axes(ellipse3_test, ax, N=1, skip_vertical=true, color="tomato", ls="--", alpha=1, lw=3)
	

	# 3D
	#ax.scatter(x_3D_gmm[mask_3D_group3, 1], x_3D_gmm[mask_3D_group3, 2], color="orange", alpha=0.5, s=15)
	@show "G3 3D"
	e=confidence_ellipse_data!(x_3D_gmm[mask_3D_group3, 1], x_3D_gmm[mask_3D_group3, 2], ax, 2, facecolor="orange", ec="orange", lw=4, alpha=0.)
	draw_ellipse_axes(e, ax, color="tomato", skip_vertical=true, lw=4)
	
	#ax.scatter(x_3D_gmm[mask_3D_group2, 1], x_3D_gmm[mask_3D_group2, 2], color="green", alpha=0.5, s=15)
	@show "G2 3D"
	e=confidence_ellipse_data!(x_3D_gmm[mask_3D_group2, 1], x_3D_gmm[mask_3D_group2, 2], ax, 2, facecolor="green", ec="green", lw=4, alpha=0.)
	draw_ellipse_axes(e, ax, color="yellowgreen", skip_vertical=true, lw=4)

	@show "G1 3D"
	#ax.scatter(x_3D_gmm[mask_3D_group1, 1], x_3D_gmm[mask_3D_group1, 2], color="blue", alpha=0.5, s=15)
	e=confidence_ellipse_data!(x_3D_gmm[mask_3D_group1, 1], x_3D_gmm[mask_3D_group1, 2], ax, 2, facecolor="blue", ec="blue", lw=4, alpha=0.)
	draw_ellipse_axes(e, ax, color="steelblue", skip_vertical=true, lw=4)


	#=text_group3 = """
	$(L"\rm \mathbf{Group\ III}")
	CEMP-no
	Faint SNe
	Pop III
	"""=#
	text_group3 = L"\rm \mathbf{Group\ III}"
	ax.text(-6.7, 6.99, text_group3, color="tomato", fontsize="small", ha="left", va="bottom")
	
	
	#=text_group2 = """
	$(L"\rm \mathbf{CEMP-no}")
	Standard CCSNe
	Pop II
	"""=#
	text_group2 = L"\rm \mathbf{Group\ II}"
	ax.text(-3.1, 5.6, text_group2, color="yellowgreen", fontsize="small")
	

	#=text_group1 = """
	$(L"\rm \mathbf{CEMP-s\ +\ CEMP-r/s}")
	MRSNe, NS-NS or NS-BH merger
	AGB contribution
	"""=#
	text_group1 = L"\rm \mathbf{Group\ I}"
	ax.text(-1.4, 8.7, text_group1, color="steelblue", fontsize="small", ha="right")

	ax.text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax.text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")
	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-6.9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax.plot([],[], color="k", lw=3, ls="-", label=L"\rm 3D")
	ax.plot([],[], color="k", lw=3, ls="--", label=L"\rm 1D")

	

	#ax.set_ylim(6, 8.5)
	ax.set_xlim(-7.1, -0.9)
	ax.set_ylim(3.6, 10.4)

	ax.xaxis.set_major_locator(plt.MaxNLocator(7))
	ax.yaxis.set_major_locator(plt.MaxNLocator(7))


	ax.legend(labelspacing=0, loc="upper left")
	
	ax.set_ylabel(L"\rm A(C)")
	ax.set_xlabel(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ 02504ac9-7af2-4f1c-8cce-f6d952f26bc1
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=true, sharey=true, layout="tight")

	
	ax.scatter(x_3D_gmm[mask_3D_highC, 1], x_3D_gmm[mask_3D_highC, 2], color="steelblue", alpha=0.5, s=25)
	x_f = range(-7.1, 0, length=100)|>collect
	ax.fill_between(
		x_f, 
		mean(x_3D_gmm[mask_3D_highC, 2]) - MUST.std(x_3D_gmm[mask_3D_highC, 2]),
		mean(x_3D_gmm[mask_3D_highC, 2]) + MUST.std(x_3D_gmm[mask_3D_highC, 2]),
		color="steelblue", alpha=0.1
	)
	ax.axhline(mean(x_3D_gmm[mask_3D_highC, 2]), ls="--", color="steelblue", lw=2)
	lb = L"\rm A(C) = "*"$(round(mean(x_3D_gmm[mask_3D_highC, 2]), sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(MUST.std(x_3D_gmm[mask_3D_highC, 2]), digits=2)))"
	ax.text(-6.9, mean(x_3D_gmm[mask_3D_highC, 2])+0.1, lb, ha="left", va="bottom", color="steelblue", fontsize="small")
	
	ax.scatter(x_3D_gmm[mask_3D_lowC, 1], x_3D_gmm[mask_3D_lowC, 2], color="tomato", alpha=0.5, s=25)
	ax.fill_between(
		x_f, 
		mean(x_3D_gmm[mask_3D_lowC, 2]) - MUST.std(x_3D_gmm[mask_3D_lowC, 2]),
		mean(x_3D_gmm[mask_3D_lowC, 2]) + MUST.std(x_3D_gmm[mask_3D_lowC, 2]),
		color="tomato", alpha=0.1
	)
	ax.axhline(mean(x_3D_gmm[mask_3D_lowC, 2]), ls="--", color="tomato", lw=2)
	lb = L"\rm A(C) = "*"$(round(mean(x_3D_gmm[mask_3D_lowC, 2]), sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(MUST.std(x_3D_gmm[mask_3D_lowC, 2]), digits=2)))"
	ax.text(-6.9, mean(x_3D_gmm[mask_3D_lowC, 2])-0.18, lb, ha="left", va="top", color="tomato", fontsize="small")
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-7.1, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax.text(0.05, 0.87, L"\mathbf{3D}", color="k", transform=ax.transAxes, fontsize="x-large")

	#ax.set_ylim(6, 8.5)
	ax.set_xlim(-7.1, -0.9)
	ax.set_ylim(4.5, 9.5)
	
	ax.set_ylabel(L"\rm A(C)")
	ax.set_xlabel(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ 71dd0c9a-e623-4806-82be-73e2456b5774
let
	plt.close()
	f, ax = plt.subplots(1, 1)

	ax.hist(x_3D_gmm[mask_3D_highC, 2], histtype="step", bins=10)
	ax.hist(x_3D_gmm[mask_3D_lowC, 2], histtype="step", bins=10)


	f
end

# ╔═╡ 87e865b0-7f1c-44b5-9173-b2a65e77ba51
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5), sharex=true, sharey=true, layout="tight")

	
	ax.scatter(x_1D_gmm[mask_1D_highC, 1], x_1D_gmm[mask_1D_highC, 2], color="steelblue", alpha=0.5, s=25)
	x_f = range(-7.1, 0, length=100)|>collect
	ax.fill_between(
		x_f, 
		mean(x_1D_gmm[mask_1D_highC, 2]) - MUST.std(x_1D_gmm[mask_1D_highC, 2]),
		mean(x_1D_gmm[mask_1D_highC, 2]) + MUST.std(x_1D_gmm[mask_1D_highC, 2]),
		color="steelblue", alpha=0.1
	)
	ax.axhline(mean(x_1D_gmm[mask_1D_highC, 2]), ls="--", color="steelblue", lw=2)
	lb = L"\rm A(C) = "*"$(round(mean(x_1D_gmm[mask_1D_highC, 2]), sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(MUST.std(x_1D_gmm[mask_1D_highC, 2]), digits=2)))"
	ax.text(-6.9, mean(x_1D_gmm[mask_1D_highC, 2])+0.18, lb, ha="left", va="top", color="steelblue", fontsize="small")
	
	ax.scatter(x_1D_gmm[mask_1D_lowC, 1], x_1D_gmm[mask_1D_lowC, 2], color="tomato", alpha=0.5, s=25)
	ax.fill_between(
		x_f, 
		mean(x_1D_gmm[mask_1D_lowC, 2]) - MUST.std(x_1D_gmm[mask_1D_lowC, 2]),
		mean(x_1D_gmm[mask_1D_lowC, 2]) + MUST.std(x_1D_gmm[mask_1D_lowC, 2]),
		color="tomato", alpha=0.1
	)
	ax.axhline(mean(x_1D_gmm[mask_1D_lowC, 2]), ls="--", color="tomato", lw=2)
	lb = L"\rm A(C) = "*"$(round(mean(x_1D_gmm[mask_1D_lowC, 2]), sigdigits=3))"*""*L"\rm\ (\sigma="*"$(round(MUST.std(x_1D_gmm[mask_1D_lowC, 2]), digits=2)))"
	ax.text(-6.9, mean(x_1D_gmm[mask_1D_lowC, 2])+0.18, lb, ha="left", va="top", color="tomato", fontsize="small")
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-7.1, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls=":", lw=3)

	ax.text(0.05, 0.87, L"\mathbf{1D}", color="k", transform=ax.transAxes, fontsize="x-large")

	#ax.set_ylim(6, 8.5)
	ax.set_xlim(-7.1, -0.9)
	ax.set_ylim(4.5, 9.5)
	
	ax.set_ylabel(L"\rm A(C)")
	ax.set_xlabel(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ f3ddcd1a-d600-4008-83c4-e2b3de687b36
let
	open("Eitner2025_CEMP.txt", "w") do f
		write(f, "# object_id,reference,teff,logg,[Fe/H],[C/Fe]_1D,[C/Fe]_3D \n")

		mask = selection_mask .& (.!isnan.(parameters_all["cfe"])) .& (parameters_all["feh"] .< 0.0)
		obid = parameters_all["object_id"][mask]
		ref = replace.(parameters_all["reference"][mask], ','=>'_', ' '=>"")
		teff = parameters_all["teff"][mask]
		logg = parameters_all["logg"][mask]
		feh = parameters_all["feh"][mask]
		cfe1D = parameters_all["cfe"][mask]
		cfe3D = parameters_all["cfe"][mask] .+ corrections_saga_all[mask]
		
		for i in eachindex(obid)
			line = MUST.@sprintf(
				"%s,%s,%i,%.2f,%.2f,%.2f,%.2f\n", 
				obid[i],
				ref[i],
				teff[i],
				logg[i],
				feh[i],
				cfe1D[i],
				cfe3D[i]
			)
			write(f, line)
		end
	end
end

# ╔═╡ 83afde32-04ce-477d-844b-af40762e1557
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5))


	#CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	CEMP_mask = selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe 

	group = bins_CEMP

	c_mean =  [mean(
		parameters_all["cfe"][group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="tomato", alpha=0.1, ls="")
	ax.errorbar(
		metallicity_bin_centers, c_mean, yerr=c_sigma,
		color="tomato", lw=2.2, label=L"\rm 1D", marker="s", capsize=5,
		markersize=8, ls="--"
	)

	
	c_mean =  [mean(
		parameters_all["cfe"][group[:, i]]
		.+ corrections_saga_all[group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][group[:, i]]
		.+corrections_saga_all[group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="steelblue", alpha=0.1, ls="")
	ax.errorbar(
		metallicity_bin_centers, c_mean, yerr=c_sigma,
		color="steelblue", label=L"\rm 3D", lw=2.2, marker="s", capsize=5,
		markersize=8
	)
	
	
	cemp_line(x) = 0.7 
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=1.0, ls=":", lw=2.5)

	ax.set_xlim(-6.1, -1.5)
	#ax.set_ylim(4, 9.2)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("[C/Fe]")
	ax.legend(loc="upper right")

	f.savefig("cfe_corrected.pdf")

	f
end

# ╔═╡ 70491b18-373c-4986-be3d-34580e207769
md"## Galactic C/O ratio"

# ╔═╡ 951e2dc3-389a-4b90-a167-95ea61b2a6e6
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5))

	selection_mask = parameters_co["logg"] .>= logg_limit
	NaNmask = .!isnan.(parameters_co["ofe"]) .& .!isnan.(parameters_co["cfe"])
	CEMP_mask = (parameters_co["cfe"] .>= cfe_limit) .& selection_mask .& NaNmask
	feh = parameters_co["feh"][CEMP_mask]
	cfe = parameters_co["cfe"][CEMP_mask]
	ofe = parameters_co["ofe"][CEMP_mask]

	bins_CO_CEMP = bin_parameter(
		metallicity_bin_edges, 
		CEMP_mask,
		para=parameters_co,
		selection=trues(length(parameters_co["feh"])),
		data=saga_co_data
	)

	co_mean =  [mean(
		parameters_co["cfe"][bins_CO_CEMP[:, i]] .-
		parameters_co["ofe"][bins_CO_CEMP[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_co["cfe"][bins_CO_CEMP[:, i]] .-
		parameters_co["ofe"][bins_CO_CEMP[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]

	@info "There are $(length(feh)) CEMP stars in total with [C/Fe] and [O/Fe] measurements."
	ax.scatter(feh, cfe .- ofe, color="k", marker="s", s=50, alpha=0.5)
	
	#ax.fill_between(metallicity_bin_centers, co_mean .- c_sigma, co_mean .+ c_sigma, color="cyan", alpha=0.1, ls="")
	#ax.plot(metallicity_bin_centers, co_mean, color="steelblue", lw=3, alpha=0.7)

	ax.axhline(mean(cfe .- ofe), color="tomato", label=L"\rm \left< [C/O] \right> = "*"$(round(mean(cfe .- ofe), digits=2))"*"\n"*L"\rm \left<C/O \right> \approx"*"$(round(0.55*exp10(mean(cfe .- ofe)), digits=2))", ls="--", lw=4)

	ax.set_ylim(-2.9, 2.9)

	ax.legend(loc="lower center")
	ax.set_ylabel("[C/O]")
	ax.set_xlabel("[Fe/H]")

	f
end

# ╔═╡ Cell order:
# ╟─29b1890e-7361-40cc-a1f2-3b874faba043
# ╠═1c7ee45a-8d60-11f0-29f2-abf8c4f707d2
# ╠═f1473678-24b0-4a97-967a-4f6070de46d7
# ╟─f7e7ae1d-49dd-4c76-a9ff-cec1013c86aa
# ╟─f7fa9c00-1bc2-4246-be6a-ac30ee3d6457
# ╟─cdedf1ac-7de7-4e61-83f7-03282f747445
# ╟─a7c4ab87-e3eb-456c-a41b-0a5d72501da4
# ╟─ed0947be-f0bd-4054-86f7-afd8f36641bf
# ╟─eff9f59a-faaf-4f3c-ab1d-67c2bee0bdc0
# ╟─9b27ec05-9e7d-40e8-b004-62d45a0e8f25
# ╟─ff4b6af5-11c5-48db-95f5-48dc055e12dc
# ╟─895721c1-dc9d-47d4-b5fb-4e189b5184c0
# ╟─16f38c33-b243-45d0-8bf3-2a0b388a8d39
# ╟─550affe7-12a4-4b3f-8edb-4ff0c8c8a41b
# ╟─06b5cd98-c60f-4d61-a723-5d2b0ac6460a
# ╟─c8469b97-df26-4753-b5be-e0b2b8a02b82
# ╟─9a3a15f2-d167-47e5-85bf-327221203325
# ╟─4f163103-84ce-4e36-9b69-2b6d3e880ca1
# ╟─82ce1cc4-350c-4ca3-9c01-21321414a240
# ╟─c3f55a57-397d-4b2c-87f8-0ebc133a9ab0
# ╟─3c6ef3eb-5ef7-43b3-946f-73049a7d9515
# ╟─c508746d-879e-4e76-baed-136fa36d716c
# ╟─2fa849f5-bf68-4548-9d0e-234ba6edf104
# ╟─d5e20abd-a153-4dc6-9895-d8142fcc6acc
# ╟─ba0a854b-5fa4-468d-91c9-380523e385cd
# ╟─6d2305b8-b452-47c9-a0a2-b097029ade80
# ╟─f4e30299-95e0-4861-919e-6a767adf1194
# ╟─cfb5b8d3-4a10-44ba-9cc7-a566339a833f
# ╟─df166dde-2837-4440-a8dc-0c2a8194866e
# ╟─9834da27-c370-4795-a91e-89a382f7b81d
# ╟─53aa57c6-2b5c-4e04-94ab-c732ef7f2ec3
# ╟─64ee5874-29e7-49a3-b8b3-5f89c7156bb9
# ╟─73ebcfac-d1b0-496f-8f2d-df5ea3aa7830
# ╟─7fb3939a-1709-4a77-be21-27924e2072cb
# ╟─7188fdf5-2366-47d6-b082-b476fc30ebf3
# ╟─6d6772d6-f956-4216-8731-de2fc4b32f3e
# ╟─6174ea3c-c2e5-4c24-8299-0c7e74afa2ef
# ╟─68e0ac58-ed87-442d-a7b5-13053c9c1cdf
# ╟─52a7d86a-a6e0-4b9e-8f63-adc0cb9a0359
# ╟─e4c53afc-f220-420e-8f61-adcca2594767
# ╟─ae363cc0-a747-4e26-808b-97df38543144
# ╟─2275bda1-cd48-42c7-9462-89aefd5e3f8c
# ╟─bf6c434a-0fd9-420b-ac06-bb8e215e9fdd
# ╟─951cf3cb-8f2f-4952-96fc-0318e866724e
# ╟─d883c7c9-ef6d-49cf-b1e3-46ba8b92f2ac
# ╟─b0653e46-ae33-4725-98aa-443a2198c3a0
# ╟─e0d62ec9-04fc-4b9b-b215-a7160ba20205
# ╟─47b98cc4-0ecc-4d95-a257-f09076283b22
# ╟─75efbc2c-8be3-427d-bf62-16c9a64d2e6a
# ╟─bedda209-be73-49d4-8bc0-b1a902a1aac3
# ╟─79da8fae-65eb-4ddf-9826-93df48bb6a2f
# ╟─80b78da6-881b-4fed-9dd0-aa7d613b6d71
# ╟─7c2b57e4-1c06-4adc-9017-a176fd36fe82
# ╟─3826dca1-246c-4c61-95ee-8e497f42aa72
# ╟─3d760565-3e87-4c13-9e78-ed6394baecaf
# ╟─41af66b7-2811-4bd5-972d-5260f67fbfd7
# ╟─d6226d59-d4e3-4652-a7a9-dbe41666511f
# ╠═b606dca3-665e-4e3c-9218-eaff570551d0
# ╟─b292f8a3-9288-4851-baae-31eb00e7f571
# ╟─e2b8a4f0-e180-4f2a-8e15-793bbd320282
# ╟─387c093e-49cf-40cb-a757-b5a8d50f50fc
# ╟─01aab765-3346-4212-b28f-c8bc5b51174f
# ╟─fc4cf849-36da-4006-85c4-418e85d5e23a
# ╠═66af10ce-e791-4cf0-88ef-bceabffe0d20
# ╠═59ac0194-f85c-400a-8940-20fc1f4eef03
# ╠═5c7a35ae-86ee-4eb8-affa-7ff0d5e99f45
# ╠═85ad8731-b3a3-4960-9ebb-732cc179525e
# ╟─6aad86e1-0e90-4a1f-9ffd-fb29d9c50686
# ╟─c1b44ff5-924e-40b5-9914-0cf87825db4b
# ╟─22a6add8-76ac-4c25-8acf-c869620505a5
# ╠═9f9c5bc2-5926-4d5e-881c-1a403c931272
# ╠═4ffb9e96-89f5-458e-bb36-5abaed2ba7fa
# ╟─925ce950-46ca-45fe-8580-b2db83230a4a
# ╟─06ee2904-3e62-4d8f-a873-068384d3bfe6
# ╟─6cc7f3f0-f4f9-4f7a-aceb-c79e7a278154
# ╠═1f66327a-5581-4250-ba81-a286e0a04181
# ╠═dcf005c3-98cd-490a-874d-f922f8d49f3f
# ╟─9effae34-66bb-47a4-a1a9-c4a8e12475a3
# ╟─8be627f4-3fe6-481a-a322-8ed2d6a015aa
# ╠═addf8b42-5994-4073-a595-33df2f3b3fd4
# ╟─481caf39-bc3f-42f5-8edd-94ea4ca6f409
# ╠═54261a34-ee49-434b-8cdc-6e0e9ffe3368
# ╠═09b8ee6d-e9f5-4cf1-ab91-305952b3972d
# ╠═1a860a43-3fd2-4b5c-8c8f-bf4341cea882
# ╟─57b503a1-e7b2-464a-a49a-359576688e58
# ╠═dd18395c-e88f-4cb1-af1d-a5a347380464
# ╠═d325fc52-afd9-474d-8e4b-05455fc340d9
# ╠═5aa85e3d-8e21-47b8-a4cc-8897dc49b756
# ╠═f4c3987b-e002-4db5-b166-d80bf67ddf75
# ╟─0e30cd20-184c-4015-b780-1e94be606e00
# ╟─b377a70f-c7cf-43d8-b6c9-cd8d42fa2835
# ╟─dfe8046d-a75e-424c-bc51-a2e429b97f70
# ╟─92228544-f92c-4c50-ad50-3cb818bb4523
# ╟─a7707e78-18f9-4c96-b5ff-c51cb9311652
# ╟─dea25815-eebf-4308-a06f-09a7a49862c3
# ╠═80ba8489-3a1a-4d3b-88fe-56954397c7df
# ╠═c0a0d7cd-6b6f-47a2-a666-959e2cb68b7c
# ╠═5709600f-4817-43bd-b34b-404e173efb43
# ╟─c2cae081-827b-4606-91c6-a5610b3f72b9
# ╟─ff2843c9-a864-466c-b8fc-87b95e1da2f2
# ╟─2404c7e6-b3cb-4c92-9261-bdef297e727f
# ╟─b9b2059c-c098-48d8-94bb-587715552774
# ╟─ea48b2d1-c752-4047-8523-474aaecc82b0
# ╟─67ddc416-4cbb-4fb5-9d3c-b83ce557d141
# ╟─64b264bd-d556-4b08-8c64-b77c51849be2
# ╟─d67f968c-7405-46ff-8588-da811fc79dbc
# ╟─bd9e94c2-1b09-4cb7-82f9-df86ee21a5aa
# ╟─98316a96-6b40-4f04-8817-7bb8c9b992fc
# ╟─37c90a88-2615-461b-9e52-6a5df9e9d76a
# ╟─d6ee50f7-d83a-40ff-b6d8-d24eb55d758f
# ╟─96ca3162-7ccf-4027-978f-238613e5413e
# ╟─556f2faa-a86e-4988-b697-af42f5b9be84
# ╠═ef8e7bda-8c97-418d-864f-2e7e04c6cca9
# ╟─f7605825-8a21-4b20-87c9-cdf54f2f04de
# ╟─f9032b6e-2a5c-43d5-9ff9-1cf4c442b858
# ╟─6d3bb0c2-b17d-49e8-8141-e22131f65f33
# ╟─716d2ed3-817e-4909-92ad-a83d7b79c171
# ╟─a1b28a53-a3b5-4e87-bc2d-5ce112c1f280
# ╟─fac3be33-656b-4c3c-9fb4-eae4dc47d71b
# ╟─b6020aa9-c570-4cef-b41d-629eff19f542
# ╟─0bd75590-be29-48da-ab79-92a9be005d29
# ╟─38efffde-0ebb-4beb-97c9-8a639761a6fd
# ╟─6fe595f7-e8d0-42b9-aab4-752151128bcf
# ╟─bfb45894-1667-4582-896a-63597c6e82ab
# ╟─2065ce06-76de-40c7-bd8f-e9beaa4cc02d
# ╟─d687f757-226e-4571-bb54-d6a69de69b89
# ╟─c20c8129-bdbe-4bb8-a7a0-36a23d8bc1ba
# ╟─2cd3aede-d981-4974-ab31-63f11a76b377
# ╟─02504ac9-7af2-4f1c-8cce-f6d952f26bc1
# ╟─71dd0c9a-e623-4806-82be-73e2456b5774
# ╟─87e865b0-7f1c-44b5-9173-b2a65e77ba51
# ╟─f3ddcd1a-d600-4008-83c4-e2b3de687b36
# ╟─83afde32-04ce-477d-844b-af40762e1557
# ╟─70491b18-373c-4986-be3d-34580e207769
# ╟─951e2dc3-389a-4b90-a167-95ea61b2a6e6
