### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 3867ab6e-bcaf-4879-9eb9-f5b42dcc6709
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("..")
	using TSO
	using MUST
	using PythonPlot
	using LaTeXStrings 
	using PlutoUI
	using KernelDensity
end

# ╔═╡ dc377839-bd92-4c9c-9432-7bf08b71add6
md"# Paper v3"

# ╔═╡ 9818e49b-1a6d-44d3-9aa2-6486842b2efa
begin
	@import_dispatch "../../../../dispatch2"
	@import_m3dis "../../../../Multi3D"
end

# ╔═╡ 1f6f5fe0-3f38-4547-b9e6-0482a76ea03a
auxfuncs = MUST.pyimport("m3dis.auxfuncs")

# ╔═╡ 1150c5dd-f486-46d6-b911-5fb1d70d3792
begin
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

	pretty_from_name(name) = begin
		mi = MUST.ModelInformation(name)
		L"\rm T_{eff}="*"$(round(Int, mi.teff))"*L"\rm\ K,\ log(g)="*"$(round( mi.logg, digits=2))"*L"\rm,\ [Fe/H] ="*"$(round(mi.feh, digits=2))"
	end
	
	TableOfContents()
end

# ╔═╡ 4ea802dc-b3e3-44ab-8758-b76bef9de441
pretty_from_name_short(name) = begin
	mi = MUST.ModelInformation(name)
	"$(round(Int, mi.teff)), "*"$(round( mi.logg, digits=2)), "*"$(round(mi.feh, digits=2))"
end

# ╔═╡ 3022833b-9b2c-4e3f-ba69-0dfc10ed5e2a


# ╔═╡ ea64c5f8-21e3-4685-ba30-2fd61ab34974
md"""
# Models
We load the models automatically and get their stellar parameters from the name, as is done in cubes_web.jl
"""

# ╔═╡ 49ce048c-18ec-44d6-a782-f474ad11c83c
datadir = @in_dispatch "CEMP_models2/"

# ╔═╡ 8dfd95b6-144c-4bb0-986e-6f5cde660c3e


# ╔═╡ ccbdff2f-bc03-444b-9b0f-3ce318e10116
models = availableRuns(datadir)

# ╔═╡ 5952687e-2249-4a8e-ba83-54ab695da15a
models_info = Dict(
	m => MUST.ModelInformation(m) for m in models
)

# ╔═╡ d41e2628-ace8-483a-b3bd-8ddbb3bd9172
cemp_models = [m for m in models if MUST.ModelInformation(m).category=="CEMP"]

# ╔═╡ 2412db85-5a78-4eb0-89f3-00cad4ecdd14
scaledSolar_models = [m for m in models if MUST.ModelInformation(m).category=="ScaledSolar"]

# ╔═╡ 8c7f2268-27be-468c-9e61-32153a97e515


# ╔═╡ 9c1a9d47-0389-4e7f-ab32-fb3e1c0e761c
md"For each model, information about the spectra is stored in the models `m3dis` cubes. Each model also contains `MARCS` models that have spectra associated."

# ╔═╡ d111e9e4-ce3f-4c6b-a65b-d11f005632c9
listOfMARCS(name) = sort(MUST.glob("*.mod", joinpath(datadir, name)))

# ╔═╡ 8f903dbe-4663-42da-8002-0dce3a6a9a66
MARCSNames(name) = split(name, '/', keepempty=false) |> last

# ╔═╡ deb4aeaa-b37c-47aa-a1a6-00d6385f8318
begin
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
		marcsmodelspaths[m] = mm
	end

	marcsmodels = Dict(k=>basename.(m) for (k, m) in marcsmodelspaths)
end

# ╔═╡ 44d75cc5-2832-45ff-be3b-3a855cc1a636


# ╔═╡ c6cf7c10-a62d-4c6c-8bdc-5acd9c007f40
begin
	get_snapshot(snapid::Int, folder, tau_name="tau500") = begin
		pick_snapshot(folder, snapid, box_name="box", tau_name=tau_name)
	end
	get_snapshot(snapid::String, folder, tau_name="tau500") = begin
		MUST.marcsBox(joinpath(folder, snapid))
	end
end

# ╔═╡ 533184b3-d690-4f44-87b6-6bb49ee7e2c2
get_spectra(snapid, folder) = begin
	pick_snapshot(folder, snapid, box_name="box_m3dis") |> first
end

# ╔═╡ 5dc20e44-24eb-4e8f-bda1-fac7a1c0a5f4
same_scaled_solar(cemp_model_name) = scaledSolar_models[
	findfirst(
		MUST.same_parameters(
			MUST.ModelInformation(a), 
			MUST.ModelInformation(cemp_model_name),
		) for a in scaledSolar_models
	)
]

# ╔═╡ 4802eaa5-1c91-4813-b335-c217c5ceac87


# ╔═╡ e624d1fd-24e2-4277-8838-a9ce0026ec56
md"# Structure comparison"

# ╔═╡ 0a3ce98b-40b5-4948-87c6-57d4cf9bbba9
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

# ╔═╡ c0ef68e0-5b8a-4d6d-a529-8d6292799420


# ╔═╡ 8749760d-4857-41ed-8a0a-f7464d5fbc67
md"## 3D vs. 1D"

# ╔═╡ 6df206d1-6794-4bab-a10d-914c6f508137
md"""
3D model: $(@bind structure_3D_select Select(models, default=first(models)))
"""

# ╔═╡ e735e811-975c-49ff-b99f-7b98e63e5150
md"""
1D model: $(@bind structure_1D_select Select(marcsmodels[structure_3D_select], default=first(marcsmodels[structure_3D_select])))
"""

# ╔═╡ 6a2c237d-a0aa-4db7-9a87-b259eda25546
begin
	structure_3D1D_figsave_txt = "$(structure_3D_select)_structure"
	md"Save figure at: $(@bind structure_3D1D_figsave confirm(TextField(length(structure_3D1D_figsave_txt)+1, default=structure_3D1D_figsave_txt)))"
end

# ╔═╡ c8ad7b76-d063-457b-809e-65b818c796f4
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))
	
	# tau500 snapshot
	b, bt = get_snapshot(-1, joinpath(datadir, structure_3D_select))
	plot_profile(b, bt, :T, ax=ax)
	
	# 1D model
	bmarcs = get_snapshot(structure_1D_select, joinpath(datadir, structure_3D_select))
	ax.plot(profile(mean, bmarcs, :log10τ500, :T)..., color="tomato", lw=5)


	# inset
	left, bottom, width, height = [0.11, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")
	ax2.set_xlabel(L"\rm x\ [Mm]")

	# formation height bars
	#=ax.hlines(3100, -3.5, 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
	ax.hlines(3450, -1.2, 0.0, ls="-", color="tomato", alpha=0.8, lw=6)
	ax.text(0.1, 3100, "3D", ha="left", va="center", color="steelblue")
	ax.text(0.1, 3450, "1D", ha="left", va="center", color="tomato")=#
	
	
	ax.set_xlabel(L"\rm optical\ depth\ [\log \tau_{500}]")
	ax.set_ylabel(L"\rm temperature\ [K]")

	ax.set_xlim(-4.2, 2.3)
	ax.set_ylim(2300, 11000)

	f.savefig(structure_3D1D_figsave_txt*"_3D1D.pdf")
	f
end

# ╔═╡ cd7eb05b-0622-4e82-a141-c049ecab95e2
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))
	
	# tau500 snapshot
	b, bt = get_snapshot(-1, joinpath(datadir, structure_3D_select))
	#plot_profile(b, bt, :T, ax=ax, xvar=:log10d)
	
	_, T = profile(mean, b, :z, :T)
	_, d = profile(mean, b, :z, :log10d)
	z, tau = profile(mean, b, :z, :log10τ500)

	mask = sortperm(tau)
	tau_ip = MUST.linear_interpolation(tau[mask], d[mask], extrapolation_bc=MUST.Line())
	xlim = tau_ip.([-4.2, 2.9])
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="steelblue", lw=3, label="CEMP")


	# tau500 snapshot 2
	b, bt = get_snapshot(-1, joinpath(datadir, same_scaled_solar(structure_3D_select)))
	_, T = profile(mean, b, :z, :T)
	_, d = profile(mean, b, :z, :log10d)
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="k", lw=3.5, label="scaled-solar", ls="--")
	
	
	# 1D model
	bmarcs = get_snapshot(structure_1D_select, joinpath(datadir, structure_3D_select))
	_, T = profile(mean, bmarcs, :z, :T)
	_, d = profile(mean, bmarcs, :z, :log10d)
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="tomato", lw=3, label="1D MARCS")


	# inset
	#=left, bottom, width, height = [0.11, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")
	ax2.set_xlabel(L"\rm x\ [Mm]")=#

	# formation height bars
	#=ax.hlines(3100, -3.5, 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
	ax.hlines(3450, -1.2, 0.0, ls="-", color="tomato", alpha=0.8, lw=6)
	ax.text(0.1, 3100, "3D", ha="left", va="center", color="steelblue")
	ax.text(0.1, 3450, "1D", ha="left", va="center", color="tomato")=#
	
	
	ax.set_xlabel(L"\rm density\ [g\ cm^{-3}]")
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlim(xlim...)
	#ax.set_ylim(2300, 10000)
	ax.legend()

	#ax.set_xlim(-4.2, 2.3)
	#ax.set_ylim(2300, 11000)

	f.savefig(structure_3D1D_figsave_txt*"_Trho_3D1D.pdf")

	f
end

# ╔═╡ 1e448509-7bfe-41a5-85f1-6df811aa6548
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))
	
	# tau500 snapshot
	b, bt = get_snapshot(-1, joinpath(datadir, structure_3D_select))
	#plot_profile(b, bt, :T, ax=ax, xvar=:log10d)

	xlim = [-4.2, 2.9]
	
	_, T = profile(mean, bt, :z, :T)
	_, d = profile(mean, bt, :z, :log10τ500)
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="steelblue", lw=3, label="CEMP")


	# tau500 snapshot 2
	b, bt = get_snapshot(-1, joinpath(datadir, same_scaled_solar(structure_3D_select)))
	_, T = profile(mean, bt, :z, :T)
	_, d = profile(mean, bt, :z, :log10τ500)
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="k", lw=3.5, label="scaled-solar", ls="--")
	
	
	# 1D model
	bmarcs = get_snapshot(structure_1D_select, joinpath(datadir, structure_3D_select))
	_, T = profile(mean, bmarcs, :z, :T)
	_, d = profile(mean, bmarcs, :z, :log10τ500)
	mask = xlim[1] .< d .< xlim[2]
	ax.plot(d[mask], T[mask], color="tomato", lw=3, label="1D MARCS")


	# inset
	#=left, bottom, width, height = [0.11, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ500=0.0)[:T][:,:,1]
	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")
	ax2.set_xlabel(L"\rm x\ [Mm]")=#

	# formation height bars
	#=ax.hlines(3100, -3.5, 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
	ax.hlines(3450, -1.2, 0.0, ls="-", color="tomato", alpha=0.8, lw=6)
	ax.text(0.1, 3100, "3D", ha="left", va="center", color="steelblue")
	ax.text(0.1, 3450, "1D", ha="left", va="center", color="tomato")=#
	
	
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlim(xlim...)
	#ax.set_ylim(2300, 11000)
	ax.legend()

	#ax.set_ylim(2300, 11000)

	f.savefig(structure_3D1D_figsave_txt*"_Ttau_3D1D.pdf")

	f
end

# ╔═╡ 99e19910-17a3-4dd2-8ec0-72dd3da2726a


# ╔═╡ 01ec615e-bc60-4e50-a076-692317095591
md"## Two Model comparison"

# ╔═╡ 55b53af5-c8d3-4086-ac92-5995002d8c37
md"""
3D model (A): $(@bind structure_3D_select_A Select(models, default=first(models)))
"""

# ╔═╡ 3b07d6ba-0e3f-484c-b3f8-da3840f5d020
md"""
1D model (A): $(@bind structure_1D_select_A Select(marcsmodels[structure_3D_select_A], default=first(marcsmodels[structure_3D_select_A])))
"""

# ╔═╡ 56a77efd-66d5-47aa-934b-99e3d1aff0a8
md"""
3D model (B): $(@bind structure_3D_select_B Select(models, default=first(models)))
"""

# ╔═╡ 989f922a-9596-4045-b4fc-186a18ca5725
md"""
1D model (B): $(@bind structure_1D_select_B Select(marcsmodels[structure_3D_select_B], default=first(marcsmodels[structure_3D_select_B])))
"""

# ╔═╡ 8cebfa74-030a-4e3a-9c91-49160c1f8f64


# ╔═╡ de37b6dc-cf8b-4a1c-95b8-ac59ac04837a
begin
	structure_3D3D_figsave_txt = "$(structure_3D_select_A)_vs_$(structure_3D_select_B)"
	md"Save figure at: $(@bind structure_3D3D_figsave confirm(TextField(length(structure_3D3D_figsave_txt)+1, default=structure_3D3D_figsave_txt)))"
end

# ╔═╡ 212b2d0b-fbcc-4934-8ebb-627545fea549
let
	plt.close()
	f, ax = plt.subplots(2, 2, figsize=(6, 7), sharey=true)
	plt.subplots_adjust(wspace=0.02, hspace=0.02)

	xlim_tau = [-4.2, 0.2]

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

	panel_data(b1, b2, bmarcs, x, y; average_on=:z) = begin
		_, y1 = profile(mean, b1, average_on, y)
		_, x1 = profile(mean, b1, average_on, x)
		z, tau = profile(mean, b1, average_on, :log10τ500)
	
		mask = sortperm(tau)
		tau_ip = MUST.linear_interpolation(tau[mask], x1[mask], extrapolation_bc=MUST.Line())
		xl = tau_ip.(xlim_tau)
		
		# tau500 snapshot 2
		_, y2 = profile(mean, b2, average_on, y)
		_, x2 = profile(mean, b2, average_on, x)
		
		# 1D model
		_, y3 = profile(mean, bmarcs, average_on, y)
		_, x3 = profile(mean, bmarcs, average_on, x)

		x1, y1, x2, y2, x3, y3, xl
	end
	
	plot_panel!(ax, d, T, d2, T2, d3, T3, xlim) = begin
		if (plot_type == :carbonmarcs) | (plot_type == :carbon) | (plot_type == :all) 
			mask = xlim[1] .< d .< xlim[2]
			ax.plot(d[mask], T[mask], color="steelblue", lw=3, label="CEMP")
		end
		
		if (plot_type == :marcs) | (plot_type == :carbonmarcs) | (plot_type == :all) 
			mask = xlim[1] .< d3 .< xlim[2]
			ax.plot(d3[mask], T3[mask], color="tomato", lw=3, label="1D MARCS")
		end
		
		if (plot_type == :marcs) | (plot_type == :carbon) | (plot_type == :all) 
			mask = xlim[1] .< d2 .< xlim[2]
			ax.plot(d2[mask], T2[mask], color="k", lw=3.5, label="scaled-solar", ls="--")
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
		bmarcs = get_snapshot(
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
		bmarcs = get_snapshot(
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
		bmarcs = get_snapshot(
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
		bmarcs = get_snapshot(
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
		plot_panel!(ax[0,0], d1, T1, d2, T2, d3, T3, xlim1)
		plot_panel!(ax[0,1], d1B, T1B, d2B, T2B, d3B, T3B, xlim1B)
		plot_panel!(ax[1,0], tau1, tT1, tau2, tT2, tau3, tT3, xlim_bottom)
		plot_panel!(ax[1,1], tau1B, tT1B, tau2B, tT2B, tau3B, tT3B, xlim_bottom)

		ax[0,0].set_xlim(xlim1...)
		ax[0,1].set_xlim(xlim1B...)
		ax[1,0].set_xlim(xlim_bottom...)
		ax[1,1].set_xlim(xlim_bottom...)

		ax[0,0].set_xlabel(L"\rm density\ [g\ cm^{-3}]")
		ax[0,0].set_ylabel(L"\rm temperature\ [K]")
		ax[0,1].legend()

		ax[0,1].set_xlabel(L"\rm density\ [g\ cm^{-3}]")
		#ax[0,0].set_ylabel(L"\rm temperature\ [K]")
		#ax[0,0].legend()

		ax[1,0].set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
		ax[1,0].set_ylabel(L"\rm temperature\ [K]")

		ax[1,1].set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
		#ax[0,0].set_ylabel(L"\rm temperature\ [K]")
		#ax[0,0].legend()
	end

	ax[0,1].tick_params(labelright=false, labelleft=false, labeltop=true, labelbottom=false)
	ax[0,0].tick_params(labeltop=true, labelbottom=false)
	ax[0,1].yaxis.set_label_position("right")
	ax[0,0].xaxis.set_label_position("top")
	ax[0,1].xaxis.set_label_position("top")

	ax[1,1].tick_params(labelright=false, labelleft=false)
	ax[1,1].yaxis.set_label_position("right")
	
	ax[0,0].xaxis.set_major_locator(plt.MaxNLocator(3))
	ax[0,1].xaxis.set_major_locator(plt.MaxNLocator(3))
	
	ax[1,0].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_A),
		ha="center",va="top", transform=ax[1,0].transAxes, 
		fontsize="small", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)
	ax[1,1].text(
		0.5,0.95,pretty_from_name_short(structure_3D_select_B),
		ha="center",va="top", transform=ax[1,1].transAxes, 
		fontsize="small", bbox=Dict("facecolor"=>"none", "edgecolor"=>"k")
	)

	f.savefig(structure_3D3D_figsave_txt*fig_label)
	
	f
end

# ╔═╡ 45cd645d-055d-41a7-a7e9-0f0867d0d183


# ╔═╡ 69bf13d1-3f9a-45c6-98c0-523951ed812f
md"# Spectra"

# ╔═╡ d51a843b-168e-411f-8641-e02d96371678
md"## Comparison"

# ╔═╡ 44694a88-c977-4d56-9fba-a42da6f60682
md"""
3D CEMP model: $(@bind spectra_CEMP_select Select([nothing, cemp_models...], default=nothing))
"""

# ╔═╡ f86f0d5b-ffc1-429a-b31c-8aad4aeb396a
if !isnothing(spectra_CEMP_select)
md"""
1D model: $(@bind spectra_1D_select Select(marcsmodels[spectra_CEMP_select], default=first(marcsmodels[spectra_CEMP_select])))
"""
end

# ╔═╡ 685e9c27-dc91-4f82-8112-998383271de0
if !isnothing(spectra_CEMP_select)
	spectra_3D1D_tag_txt= "GBand_DC+0"
	md"Spectra to show: $(@bind spectra_3D1D_tag confirm(TextField(length(spectra_3D1D_tag_txt)+5, default=spectra_3D1D_tag_txt)))"
end

# ╔═╡ 136db4eb-c274-46ef-bac2-03a3bf300cac
if !isnothing(spectra_CEMP_select)
	spectra_3D1D_figsave_txt = "$(spectra_CEMP_select)_structure_3D1D.pdf"
	md"Save figure at: $(@bind spectra_3D1D_figsave confirm(TextField(length(spectra_3D1D_figsave_txt)+1, default=spectra_3D1D_figsave_txt)))"
end

# ╔═╡ 560fcc8a-a9ee-4420-bc97-eb53b9f6de94
if !isnothing(spectra_CEMP_select)
let
	f, ax = plt.subplots(1, 1, figsize=(12, 4))
	
	# 3D CEMP spectrum
	b = get_spectra(-1, joinpath(datadir, spectra_CEMP_select))
	λ, F = MUST.mean_integrated_flux(b, Symbol(spectra_3D1D_tag))
	ax.plot(λ, F, color="k", lw=1.7, label="CEMP")

	# scaled-solar
	smod = scaledSolar_models[
		findfirst(
			MUST.same_parameters(
				MUST.ModelInformation(a), 
				MUST.ModelInformation(spectra_CEMP_select),
			) for a in scaledSolar_models
		)
	]
	b = get_spectra(-1, joinpath(datadir, smod))
	λ, F = MUST.mean_integrated_flux(b, Symbol(spectra_3D1D_tag))
	ax.plot(λ, F, color="steelblue", lw=1.7, label="non-CEMP", ls="--")

	
	# 1D model
	bmarcs = get_spectra(spectra_1D_select, joinpath(datadir, spectra_CEMP_select))
	λM, FM = MUST.mean_integrated_flux(bmarcs, Symbol(spectra_3D1D_tag))
	ax.plot(λM, FM, color="tomato", lw=1.7, label="MARCS")


	ax.set_xlabel(L"\rm wavelength\ [\AA]")
	ax.set_ylabel(L"\rm flux\ [normalized]")
	ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05))

	ax.set_ylim(-0.11, 1.01)
	ax.set_xlim(4297, 4303)

	#f.savefig(spectra_3D1D_figsave)
	f
end
end

# ╔═╡ efae7789-e54d-46fa-ae43-063575e522bf


# ╔═╡ dc3f65a6-a4f5-4bae-bb83-1287384a2e10
md"## Abundance corrections"

# ╔═╡ 3389b1ca-c595-4d78-aabf-dfee9e9672a3
abundance_correction_parameters = [(5250.0, 3.0), (5750.0, 4.5)]

# ╔═╡ 8e56ef77-75ac-4774-ab13-d9849c4f6d30


# ╔═╡ 9365a10c-59e1-448f-9fc2-2dbf05bf66d6
md"We compute the abundance corrections based on the equivalent width for each of the desired models as a function of metallicity. We compute the corrections w.r.t. 1D MARCS and 3D scaled-solar reference."

# ╔═╡ cecc0bff-92f5-46d7-a523-287b13753c70
begin
	cemp_models_for_each_parameter = []
	for (teff_para, logg_para) in abundance_correction_parameters
		# find all metallicities at those parameters
		para_mask = [
			(m.teff==teff_para)&&(m.logg==logg_para)
			for m in MUST.ModelInformation.(cemp_models)
		]
		
		append!(cemp_models_for_each_parameter, [cemp_models[para_mask]])
	end

	cemp_models_for_each_parameter
end

# ╔═╡ 1c62bc8d-b0d0-4905-98f2-7c7c16cf27b7


# ╔═╡ 9be354b8-fcb0-44b5-ae8d-7ff137325a2a
md"Abundance corrections can now be computed for each of those parameters, for multiple references"

# ╔═╡ b874cdbb-7496-471d-8e2e-4e3f47a126b0
struct AbundanceCorrection
	refmodel
	model
	ab_ref
	ew_ref
	ab
	ew
	Δab
end

# ╔═╡ d967809b-c044-49a6-a392-334d40c6d89e
spectra_name_reference = "GBand_DC+0"

# ╔═╡ e6d97f9d-7069-4b09-9486-d6a6c279838c
spectra_name_models = [
	"GBand_DC+0",
	"GBand_DC+02",
	"GBand_DC+04",
	"GBand_DC+06",
	"GBand_DC+08",
	"GBand_DC-02",
	"GBand_DC-04",
	"GBand_DC-06",
	"GBand_DC-08"
]

# ╔═╡ 4ff6bac5-7f46-4f4b-a221-8b75bdb6b53f
begin
	corrections_1D = []
	for i in eachindex(cemp_models_for_each_parameter)
		metallicity_corrections = Dict()
		
		for j in eachindex(cemp_models_for_each_parameter[i])
			# current model
			model3D_name = cemp_models_for_each_parameter[i][j]

			# current metallicity
			feh = MUST.ModelInformation(model3D_name).feh

			# current model spectra
			model3D = get_spectra(-1, joinpath(datadir, model3D_name))

			metallicity_corrections[feh] = []
			for marcsmodelname in marcsmodels[model3D_name]			
				# compute corrections w.r.t all MARCS models
				marcsmodel = get_spectra(
					marcsmodelname, joinpath(datadir, model3D_name)
				)
				ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
					"[C/Fe]", 
					marcsmodel, spectra_name_reference,
					model3D, spectra_name_models,
				)
				ac = AbundanceCorrection(
					marcsmodelname,
					model3D_name,
					ab_ref,
					ew_ref,
					ab,
					ew,
					Δab,
				)
				append!(metallicity_corrections[feh], [ac])
			end
		end

		# in metallicity_corrections we have stored the abundance corrections
		# for each metallicity of a given set of parameters
		append!(corrections_1D, [metallicity_corrections])
	end
end

# ╔═╡ 258d76cb-ac6b-453e-b541-29ad51461195


# ╔═╡ f77ae3e8-374a-4e3c-b8b3-7985747bccc4
md"We can now interpolate in metallicity, and store the correction functions sorted for each reference"

# ╔═╡ a58d4a2d-45ad-4bb3-9000-6a05f43e5d71
get_info_from_name(name) = begin
	if occursin(".mod", name)
		splitname = split(name, '_')
		teff = parse(Float64, splitname[1][2:end])
		logg = parse(Float64, splitname[2][2:end])
		feh  = parse(Float64, splitname[6][2:end])
		vmic  = parse(Float64, splitname[4][2:end])

		MUST.ModelInformation(
			category="$(splitname[1][1])", 
			teff=teff, logg=logg, feh=feh, 
			extension=vmic, 
			original_name=name
		)
	else
		MUST.ModelInformation(name)
	end
end

# ╔═╡ 53fef078-5301-4772-981c-68fbb1349f8d
get_closest_model(models, reference) = begin
	models_teff = [m.teff for m in models] ./1000.
	models_logg = [m.logg for m in models]
	models_feh  = [m.feh for m in models]

	teff_select = reference.teff / 1000.
	logg_select = reference.logg 
	feh_select = reference.feh
	
	# compute the 3D distance to all MARCS models
	d = (models_teff .- teff_select) .^2 .+
		(models_logg .- logg_select) .^2 .+
		(models_feh .- feh_select) .^2

	# pick the closest 
	argmin(d)
end

# ╔═╡ 7a831dfa-f216-4d6b-9370-5f286d7f3e77
function read_z_corrections(corrs, reference; what=:Δab, vmic=nothing)
	
	metallicities = sort(keys(corrs) |> collect)
	reference_corrections = []
	reference_metallicities = []

	for z in metallicities
		# check if there is a reference with that name
		mod = if reference==:closest
			infos = [a.refmodel for a in corrs[z]] .|> get_info_from_name
			vmic_mask = if isnothing(vmic)
				trues(length(infos))
			else
				allvmics = [i.extension for i in infos]
				closest_vmic = allvmics[argmin(abs.(allvmics .- vmic))]
				[i.extension ≈ closest_vmic for i in infos]
			end
			ic = get_closest_model(
				[a.refmodel for a in corrs[z]][vmic_mask] .|> get_info_from_name, 
				corrs[z][1].model |> get_info_from_name
			) 

			corrs[z][vmic_mask][ic]
		else
			ic = findfirst([a.refmodel==reference for a in corrs[z]])
			isnothing(ic) ? nothing : corrs[z][ic]
		end
		if !isnothing(mod)
			append!(reference_corrections, [getfield(mod, what)])
			append!(reference_metallicities, [z])
		end
	end

	reference_metallicities, reference_corrections
end

# ╔═╡ 631410d5-7e7a-458f-a7fb-74b462913bd8


# ╔═╡ 1149d154-42e0-4e94-b26f-716504872afa
begin
	corrections_feh_1D = [
		read_z_corrections(
			corrections_1D[i],
			:closest, vmic=1
		) for i in eachindex(corrections_1D)
	]
	
	ab_feh_1D = [
		read_z_corrections(
			corrections_1D[i],
			:closest, vmic=1, what=:ab
		) for i in eachindex(corrections_1D)
	]
	
	ew_feh_1D = [
		read_z_corrections(
			corrections_1D[i],
			:closest, vmic=1, what=:ew
		) for i in eachindex(corrections_1D)
	]
	
	ab_ref_feh_1D = [
		read_z_corrections(
			corrections_1D[i],
			:closest, vmic=1, what=:ab_ref
		) for i in eachindex(corrections_1D)
	]
	
	ew_ref_feh_1D = [
		read_z_corrections(
			corrections_1D[i],
			:closest, vmic=1, what=:ew_ref
		) for i in eachindex(corrections_1D)
	]
	
	correctionfunctions_feh_1D = [
		MUST.linear_interpolation(
			corrections_feh_1D[i][1], corrections_feh_1D[i][2], extrapolation_bc=MUST.Line()
		) for i in eachindex(corrections_feh_1D)
	]
	@info "1D corrections computed."
end;

# ╔═╡ 87aaadf1-2d46-4d14-8723-8c843040d1d1


# ╔═╡ f3dc3826-c36e-4292-9431-a86d4325903c
let
	plt.close()

	f, ax = plt.subplots(1, 2, figsize=(8, 6), sharex=true, sharey=true)
	plt.subplots_adjust(wspace=0)

	
	ax[0].axhline(0.0, color="0.5", lw=1, ls=":")
	ax[1].axhline(0.0, color="0.5", lw=1, ls=":")
	

	color = ["tomato", "steelblue"]
	ls = ["-", "--"]
	lw = [3.5, 3.5]
	marker = ["s", "s"]
	markerfacecolor = ["tomato", "steelblue"]

	for i in eachindex(corrections_feh_1D)
		ax[0].plot(
			corrections_feh_1D[i][1], corrections_feh_1D[i][2], 
			color=color[i], marker=marker[i], ls=ls[i], lw=lw[i],
			markersize=12, markerfacecolor=markerfacecolor[i], markeredgewidth=2,
			label=L"\rm T_{eff}="*"$(round(Int,abundance_correction_parameters[i][1]))"*L"\rm \ K"*"\n"*"log(g)="*"$(round(abundance_correction_parameters[i][2], digits=2))"
			#label="$(round(Int,abundance_correction_parameters[i][1])) K, $(round(abundance_correction_parameters[i][2], digits=2))"
		)
	end

	ax[0].legend(ncol=1, loc="upper center")
	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm \Delta\ A(C)\ [dex]")
	ax[0].set_ylim(-0.79, 0.47)
	ax[0].set_xlim(-6.3, -1.7)

	ax[0].set_title(L"\rm 3D\ CEMP - 1D")
	ax[1].set_title(L"\rm 3D\ CEMP - 3D\ non\ CEMP")

	f
end

# ╔═╡ 093b3783-af82-4e34-a1e6-65a0df5bb35c


# ╔═╡ 3a59080e-f925-48c7-b01d-660c45a7f541
md"Choose parameters: $(@bind i_cog_p_select Select([i=>a for (i, a) in enumerate(abundance_correction_parameters)]))"

# ╔═╡ d8f6a261-88ab-4222-8a51-ce1b359df452
begin
	cog_p_select = abundance_correction_parameters[i_cog_p_select]
	md"Choose metallicity: $(@bind cog_z_select confirm(Slider(ab_feh_1D[i_cog_p_select][1], show_value=true)))"
end

# ╔═╡ 2dd8fd87-f543-413e-b8f3-4187ff67e52c
let
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	feh_cog_select = ab_feh_1D[i_cog_p_select][1]
	i_cog_feh_select = argmin(abs.(feh_cog_select .- cog_z_select))

	ab = ab_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ew = ew_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ab_ref = ab_ref_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ew_ref = ew_ref_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	Δab = corrections_feh_1D[i_cog_p_select][2][i_cog_feh_select]

	ew_int = MUST.linear_interpolation(ab, ew, extrapolation_bc=MUST.Line())
	da = maximum(ab) - minimum(ab)
	ab_arr = range(minimum(ab) - 0.1*da, maximum(ab) +0.1*da, length=100)
	ew_arr = ew_int.(ab_arr)
	ax.plot(
		ab_arr, ew_arr, 
		marker="", color="k", markersize=10, ls="--"
	)
	ax.plot(ab, ew, marker="s", color="k", label="3D", markersize=11, ls="", markerfacecolor="w", markeredgewidth=2)
	ax.plot([ab_ref], [ew_ref], 
		marker="X", color="tomato", ls="",
		label="1D", markersize=12
	)

	ax.set_ylabel(L"\rm EW\ [\AA]")
	ax.set_xlabel(L"\rm [C/Fe]")

	paras_as_name = L"\rm T_{eff}="*"$(round(Int,cog_p_select[1]))"*L"\rm ,\  log(g)="*"$(round(cog_p_select[2], digits=2))"
	ax.set_title(paras_as_name*L"\rm ,\ [Fe/H] ="*"$(cog_z_select)")
	ax.text(
		0.95, 0.05, 
		L"\rm \Delta [C/Fe]="*"$(round(Δab, sigdigits=2))",
		ha="right", va="bottom", transform=ax.transAxes,
		color="tomato",
	)
	ax.set_xlim(minimum(ab_arr), maximum(ab_arr))

	ax.axvline(ab_ref, color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axvline(ab_ref+ Δab, color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axhline(ew_int(ab_ref + Δab), color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axhline(ew_int(ab_ref), color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)

	ax.legend()
	ax.xaxis.set_major_locator(plt.MaxNLocator(6))

	f
end

# ╔═╡ 75f3327c-00c4-4594-86ac-e2b5237dc0bc


# ╔═╡ ef9fef62-5e50-4a1d-b1a9-4db7289663a7


# ╔═╡ 874f7823-d626-45c9-aac8-72b22fe872c6
md"# Data comparison"

# ╔═╡ 8649e73c-ddf0-42b0-b0e1-93295bcdda5f
md"## SAGA data"

# ╔═╡ 35c27a34-6c3d-4b9f-ba87-305e81eadaff
saga_all_path = "../cgisess_32aefcdc405d0a0c9fb7936f6f1324aa_data.tsv"

# ╔═╡ 0bbf88db-b471-4c9c-999f-c3ace82841be
begin
	saga_all_data = MUST.readdlm(saga_all_path, '\t', String, skipstart=1)

	teff_all, logg_all, z_all = saga_all_data[:, 5], saga_all_data[:, 6], saga_all_data[:, 7]

	cfe_all, bafe_all, eufe_all, lafe_all, feh_all = saga_all_data[:, 9], saga_all_data[:, 10], saga_all_data[:, 11], saga_all_data[:, 12], saga_all_data[:, 13]

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

	mask = (parameters_all["feh"].>0.0) 
	parameters_all["feh"][mask] .= NaN

	parameters_all
end

# ╔═╡ 8eb12511-d69c-4e4a-9eeb-052471784a32


# ╔═╡ 023e815c-5831-4445-8fb3-58981fcb09fb
md"To isolate the giant stars, one has to define a limiting log(g)."

# ╔═╡ 26a8f45f-e3fa-427a-82ab-004b598182a4
logg_limit = 2.7

# ╔═╡ bb39c096-108c-4920-8eda-0eb60e225f35


# ╔═╡ 05779478-7a59-43bb-badc-83c90a37016f
md"We can furthermore define multiple sub-limits for different CEMP groups"

# ╔═╡ 8781d794-ad28-483b-8fd1-cf8fd6d907dd
# For CEMP in general
cfe_limit = 0.7

# ╔═╡ f7bd519c-01c7-42ad-ab22-5998cd0590ff
# For CEMP-s and CEMP-r/s
ba_limit = 1.0

# ╔═╡ ea87f29c-10ba-494e-8d58-ba324c7b15e1
# For CEMP-r
eu_limit = 1.0

# ╔═╡ 96fe954a-f6e0-414e-8100-48fe0c8ae271


# ╔═╡ 94f5e285-198d-4a09-936c-5c19509eba10
md"Apply r/s cut from Yoon et al. (2016) at A(C)=7.1"

# ╔═╡ 8131eeff-aa5e-4434-a1e5-af46220b92e2
c_limit = 7.1

# ╔═╡ c5702c96-0a5d-4420-8271-fdf38d7a5a3d
parameters_all["r/s"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit

# ╔═╡ d7b304d3-79ea-42e7-817d-3f35d003653a
parameters_all["no"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit

# ╔═╡ 7ab6c379-1747-4ae6-94c1-6ef3a3c83f72


# ╔═╡ 51fa1fc7-5191-48a9-af04-648e8669ba3c
selection_mask = parameters_all["logg"] .>= logg_limit

# ╔═╡ 4edb9461-7f5b-41a0-ac30-4ae72a42d3c8


# ╔═╡ 990e7d73-3e8e-4900-a225-4e5cac6ac8c6
md"Next we can correct the data by using the abundance corrections from above"

# ╔═╡ 5e31351e-5e77-46a5-9e8c-baa53eeaee73
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

# ╔═╡ 73b3581d-a484-4ff4-8e69-445192946a3b
"""
	interpolate_corrections(teff, logg, feh; feh_interpolators=interpolators)

Interpolate the corrections in metallicity based on the closest point in the Teff, logg plane.
"""
function interpolate_corrections(teff, logg, feh; feh_interpolators=interpolators)
	logg_gr = getindex.(abundance_correction_parameters, 2)
	teff_gr = getindex.(abundance_correction_parameters, 1)

	# closest point teff-logg plane
	tefflogg = closest_point(zip(teff_gr,logg_gr)|>collect, (teff, logg))

	# interpolate the corrections in metallicity
	if !isnan(feh)
		feh_interpolators[tefflogg](feh)
	else
		feh
	end
end

# ╔═╡ 2fc5e72e-4c17-4519-915d-88479d11228f
begin
	corrections_saga_all = fill!(similar(parameters_all["cfe"]), 0.0)
	
	for i in eachindex(parameters_all["cfe"])
		logg = parameters_all["logg"][i]
		teff = parameters_all["teff"][i]
		feh = parameters_all["feh"][i]
		corrections_saga_all[i] = interpolate_corrections(
			teff, logg, feh, feh_interpolators=correctionfunctions_feh_1D
		)
	end

	@info "3D corrections applied."
end

# ╔═╡ c012c57f-ea9c-44fa-ab8d-e2794953b214


# ╔═╡ 0bb30e77-41a2-4e73-a032-cf6b408f3d7c
md"## Galactic CEMP Distributions"

# ╔═╡ 216a0e5c-12ee-4e04-9f50-cb5bd4fa510e
#metallicity_bin_edges = [-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]
metallicity_bin_edges = [-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]

# ╔═╡ 27b5837f-b21a-434b-a074-a18e9ef7f3d2
metallicity_bin_centers = (metallicity_bin_edges[2:end] .+ metallicity_bin_edges[1:end-1]) ./ 2

# ╔═╡ 952d802f-82dd-47e9-a1f9-fa98fc505d28
function bin_parameter(bins, general_mask=trues(size(saga_all_data, 1)); name="feh")
	bin_centers = (bins[2:end] .+ bins[1:end-1]) ./ 2
	bin_masks = falses(size(saga_all_data, 1), length(bin_centers))
	for i in eachindex(bin_centers)
		bin_masks[:, i] .= (bins[i] .<= parameters_all[name]) .&
		(parameters_all[name] .< bins[i+1]) .& (.!isnan.(parameters_all[name]))

		bin_masks[:, i] .= bin_masks[:, i] .& general_mask .& selection_mask
	end

	bin_masks
end

# ╔═╡ e465e7b6-91a7-4826-859b-5b60382bc96f
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

# ╔═╡ 48bf4d88-82e0-4aaa-b1e8-35d375a03910
begin
	@info "cumsum of those bins:"
	@info cumsum(count_bins_general)
	@info cumsum(count_bins_CEMP)
	@info cumsum(count_bins_CEMP_corr)
end

# ╔═╡ f07880a2-bbaa-47c9-897a-fd94d93ede4f


# ╔═╡ 6e1e3bf9-7299-4ec3-8802-bbe798089fcb
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm 1D", marker="s")
	ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 3D", marker="s")

	#ax.plot(metallicity_bin_centers, count_bins_CEMP ./ count_bins_general*100, zorder=10, color="k", lw=2.0, label=L"\rm 1D", marker="s")
	#ax.plot(metallicity_bin_centers, count_bins_CEMP_corr ./ count_bins_general*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 3D", marker="s")

	ax.set_zorder(10)
	ax.patch.set_visible(false)
	
	ax.set_xlabel(L"\rm [Fe/H]")
	ax.set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax.legend(loc="lower left")

	ax.set_xlim(-6.2, -1.8)
	ax.set_ylim(0.24*100, 1.04*100)
	
	f
end

# ╔═╡ 892a04e1-3758-435a-b8cd-e7c03b777923
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=true, sharey=true)

	plt.subplots_adjust(wspace=0)

	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, marker="", alpha=0.3, ls=":")
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, marker="", alpha=0.3, ls=":")
	
	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm CEMP-no,\ 1D", marker="s", markersize=7)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm CEMP-no,\ 3D", marker="s", markersize=7)

	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm CEMP-r/s,\ 1D", marker="o", ls="--", markeredgecolor="k", markerfacecolor="w", markersize=8)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm CEMP-r/s,\ 3D", marker="o", ls="--", markeredgecolor="tomato", markerfacecolor="w", markersize=8)


	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax[0].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.04))
	ax[1].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.04))

	ax[0].set_xlim(-6.2, -1.5)

	
	f
end

# ╔═╡ 76e3d925-6778-4d40-a02f-23e6aa6258eb
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
	ax.scatter(feh, c+corrections_saga_all[CEMP_mask], s=15, marker="o", c="None", alpha=0.6, edgecolor="tomato", rasterized=true)

	
	# add histograms
	ax_xhist.hist(feh, bins=30, color="0.5")
	ax_yhist.hist(c, orientation="horizontal", bins=30, label="1D", color="0.5")
	ax_yhist.hist(c+corrections_saga_all[CEMP_mask], orientation="horizontal", bins=30, label="3D", histtype="step", lw=3, color="tomato")
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
	ax.plot(metallicity_bin_centers, c_mean, color="k", lw=3, label=L"\rm 1D")

	
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
	ax.plot(metallicity_bin_centers, c_mean, color="tomato", label=L"\rm 3D", lw=3)
	
	
	
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


	f
end

# ╔═╡ Cell order:
# ╟─dc377839-bd92-4c9c-9432-7bf08b71add6
# ╠═3867ab6e-bcaf-4879-9eb9-f5b42dcc6709
# ╟─9818e49b-1a6d-44d3-9aa2-6486842b2efa
# ╟─1f6f5fe0-3f38-4547-b9e6-0482a76ea03a
# ╟─1150c5dd-f486-46d6-b911-5fb1d70d3792
# ╟─4ea802dc-b3e3-44ab-8758-b76bef9de441
# ╟─3022833b-9b2c-4e3f-ba69-0dfc10ed5e2a
# ╟─ea64c5f8-21e3-4685-ba30-2fd61ab34974
# ╠═49ce048c-18ec-44d6-a782-f474ad11c83c
# ╟─8dfd95b6-144c-4bb0-986e-6f5cde660c3e
# ╟─ccbdff2f-bc03-444b-9b0f-3ce318e10116
# ╟─5952687e-2249-4a8e-ba83-54ab695da15a
# ╟─d41e2628-ace8-483a-b3bd-8ddbb3bd9172
# ╟─2412db85-5a78-4eb0-89f3-00cad4ecdd14
# ╟─8c7f2268-27be-468c-9e61-32153a97e515
# ╟─9c1a9d47-0389-4e7f-ab32-fb3e1c0e761c
# ╟─d111e9e4-ce3f-4c6b-a65b-d11f005632c9
# ╟─8f903dbe-4663-42da-8002-0dce3a6a9a66
# ╟─deb4aeaa-b37c-47aa-a1a6-00d6385f8318
# ╟─44d75cc5-2832-45ff-be3b-3a855cc1a636
# ╟─c6cf7c10-a62d-4c6c-8bdc-5acd9c007f40
# ╟─533184b3-d690-4f44-87b6-6bb49ee7e2c2
# ╟─5dc20e44-24eb-4e8f-bda1-fac7a1c0a5f4
# ╟─4802eaa5-1c91-4813-b335-c217c5ceac87
# ╟─e624d1fd-24e2-4277-8838-a9ce0026ec56
# ╟─0a3ce98b-40b5-4948-87c6-57d4cf9bbba9
# ╟─c0ef68e0-5b8a-4d6d-a529-8d6292799420
# ╟─8749760d-4857-41ed-8a0a-f7464d5fbc67
# ╟─6df206d1-6794-4bab-a10d-914c6f508137
# ╟─e735e811-975c-49ff-b99f-7b98e63e5150
# ╟─6a2c237d-a0aa-4db7-9a87-b259eda25546
# ╟─c8ad7b76-d063-457b-809e-65b818c796f4
# ╟─cd7eb05b-0622-4e82-a141-c049ecab95e2
# ╟─1e448509-7bfe-41a5-85f1-6df811aa6548
# ╟─99e19910-17a3-4dd2-8ec0-72dd3da2726a
# ╟─01ec615e-bc60-4e50-a076-692317095591
# ╟─55b53af5-c8d3-4086-ac92-5995002d8c37
# ╟─3b07d6ba-0e3f-484c-b3f8-da3840f5d020
# ╟─56a77efd-66d5-47aa-934b-99e3d1aff0a8
# ╟─989f922a-9596-4045-b4fc-186a18ca5725
# ╟─8cebfa74-030a-4e3a-9c91-49160c1f8f64
# ╟─de37b6dc-cf8b-4a1c-95b8-ac59ac04837a
# ╟─212b2d0b-fbcc-4934-8ebb-627545fea549
# ╟─45cd645d-055d-41a7-a7e9-0f0867d0d183
# ╟─69bf13d1-3f9a-45c6-98c0-523951ed812f
# ╟─d51a843b-168e-411f-8641-e02d96371678
# ╟─44694a88-c977-4d56-9fba-a42da6f60682
# ╟─f86f0d5b-ffc1-429a-b31c-8aad4aeb396a
# ╟─685e9c27-dc91-4f82-8112-998383271de0
# ╟─136db4eb-c274-46ef-bac2-03a3bf300cac
# ╟─560fcc8a-a9ee-4420-bc97-eb53b9f6de94
# ╟─efae7789-e54d-46fa-ae43-063575e522bf
# ╟─dc3f65a6-a4f5-4bae-bb83-1287384a2e10
# ╠═3389b1ca-c595-4d78-aabf-dfee9e9672a3
# ╟─8e56ef77-75ac-4774-ab13-d9849c4f6d30
# ╟─9365a10c-59e1-448f-9fc2-2dbf05bf66d6
# ╟─cecc0bff-92f5-46d7-a523-287b13753c70
# ╟─1c62bc8d-b0d0-4905-98f2-7c7c16cf27b7
# ╟─9be354b8-fcb0-44b5-ae8d-7ff137325a2a
# ╟─b874cdbb-7496-471d-8e2e-4e3f47a126b0
# ╠═d967809b-c044-49a6-a392-334d40c6d89e
# ╠═e6d97f9d-7069-4b09-9486-d6a6c279838c
# ╟─4ff6bac5-7f46-4f4b-a221-8b75bdb6b53f
# ╟─258d76cb-ac6b-453e-b541-29ad51461195
# ╟─f77ae3e8-374a-4e3c-b8b3-7985747bccc4
# ╟─a58d4a2d-45ad-4bb3-9000-6a05f43e5d71
# ╟─53fef078-5301-4772-981c-68fbb1349f8d
# ╟─7a831dfa-f216-4d6b-9370-5f286d7f3e77
# ╟─631410d5-7e7a-458f-a7fb-74b462913bd8
# ╟─1149d154-42e0-4e94-b26f-716504872afa
# ╟─87aaadf1-2d46-4d14-8723-8c843040d1d1
# ╟─f3dc3826-c36e-4292-9431-a86d4325903c
# ╟─093b3783-af82-4e34-a1e6-65a0df5bb35c
# ╟─3a59080e-f925-48c7-b01d-660c45a7f541
# ╟─d8f6a261-88ab-4222-8a51-ce1b359df452
# ╟─2dd8fd87-f543-413e-b8f3-4187ff67e52c
# ╟─75f3327c-00c4-4594-86ac-e2b5237dc0bc
# ╟─ef9fef62-5e50-4a1d-b1a9-4db7289663a7
# ╟─874f7823-d626-45c9-aac8-72b22fe872c6
# ╟─8649e73c-ddf0-42b0-b0e1-93295bcdda5f
# ╠═35c27a34-6c3d-4b9f-ba87-305e81eadaff
# ╟─0bbf88db-b471-4c9c-999f-c3ace82841be
# ╟─8eb12511-d69c-4e4a-9eeb-052471784a32
# ╟─023e815c-5831-4445-8fb3-58981fcb09fb
# ╠═26a8f45f-e3fa-427a-82ab-004b598182a4
# ╟─bb39c096-108c-4920-8eda-0eb60e225f35
# ╟─05779478-7a59-43bb-badc-83c90a37016f
# ╠═8781d794-ad28-483b-8fd1-cf8fd6d907dd
# ╠═f7bd519c-01c7-42ad-ab22-5998cd0590ff
# ╠═ea87f29c-10ba-494e-8d58-ba324c7b15e1
# ╟─96fe954a-f6e0-414e-8100-48fe0c8ae271
# ╟─94f5e285-198d-4a09-936c-5c19509eba10
# ╠═8131eeff-aa5e-4434-a1e5-af46220b92e2
# ╠═c5702c96-0a5d-4420-8271-fdf38d7a5a3d
# ╠═d7b304d3-79ea-42e7-817d-3f35d003653a
# ╟─7ab6c379-1747-4ae6-94c1-6ef3a3c83f72
# ╠═51fa1fc7-5191-48a9-af04-648e8669ba3c
# ╟─4edb9461-7f5b-41a0-ac30-4ae72a42d3c8
# ╟─990e7d73-3e8e-4900-a225-4e5cac6ac8c6
# ╟─5e31351e-5e77-46a5-9e8c-baa53eeaee73
# ╟─73b3581d-a484-4ff4-8e69-445192946a3b
# ╟─2fc5e72e-4c17-4519-915d-88479d11228f
# ╟─c012c57f-ea9c-44fa-ab8d-e2794953b214
# ╟─0bb30e77-41a2-4e73-a032-cf6b408f3d7c
# ╠═216a0e5c-12ee-4e04-9f50-cb5bd4fa510e
# ╟─27b5837f-b21a-434b-a074-a18e9ef7f3d2
# ╟─952d802f-82dd-47e9-a1f9-fa98fc505d28
# ╟─e465e7b6-91a7-4826-859b-5b60382bc96f
# ╟─48bf4d88-82e0-4aaa-b1e8-35d375a03910
# ╟─f07880a2-bbaa-47c9-897a-fd94d93ede4f
# ╟─6e1e3bf9-7299-4ec3-8802-bbe798089fcb
# ╟─892a04e1-3758-435a-b8cd-e7c03b777923
# ╟─76e3d925-6778-4d40-a02f-23e6aa6258eb
