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
	"$(round(Int, mi.teff)) K, "*"$(round( mi.logg, digits=2)), "*"$(round(mi.feh, digits=2))"
end

# ╔═╡ f4a66407-d7c8-43f9-b68f-1f7f29b4eba6
pretty_from_name_no_feh(name) = begin
	mi = MUST.ModelInformation(name)
	L"\rm T_{eff}="*"$(round(Int, mi.teff))"*L"\rm\ K,\ log(g)="*"$(round( mi.logg, digits=2))"
end

# ╔═╡ 06f0b894-d945-4941-aee5-e3e699e162f0


# ╔═╡ ea64c5f8-21e3-4685-ba30-2fd61ab34974
md"""
# Models
We load the models automatically and get their stellar parameters from the name, as is done in cubes_web.jl
"""

# ╔═╡ 49ce048c-18ec-44d6-a782-f474ad11c83c
datadir = @in_dispatch "CEMP_models3/"

# ╔═╡ d30fd639-00c4-433e-b5c4-cdc59e38281e
snapshot_id = -1

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
	snapid_m3dis = MUST.snapshotid(folder |> MUST.converted_snapshots, snapid)
	pick_snapshot(folder, snapid_m3dis, box_name="box_m3dis") |> first
end

# ╔═╡ bfc09acd-8fe3-44c9-a44f-8aefe85f3d76
get_m3d_spectra(snapid, folder; extension="lam_4297-4303_contr") = begin
	i = if typeof(snapid) <: Int && (snapid<0)
		"m3dis_$(MUST.list_snapshots(MUST.converted_snapshots(folder), numbered_only=true)[end+snapid+1])"
	elseif typeof(snapid) <: Int
		"m3dis_$(snapid)"
	else
		snapid
	end
	MUST.M3DISRun(joinpath(folder, "spectra_$(i)_$(extension)"))
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

# ╔═╡ 9ccc4261-0b53-4083-b8a3-850f10994048
md"""
Formation height to overplot (3D): $(@bind formation_height_3D confirm(TextField(15)))
"""

# ╔═╡ e6f379f0-6cbf-4d74-9c5e-4fb567307eea
md"""
Formation height to overplot (1D): $(@bind formation_height_1D confirm(TextField(15)))
"""

# ╔═╡ c8ad7b76-d063-457b-809e-65b818c796f4
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))
	
	# tau500 snapshot
	b, bt = get_snapshot(snapshot_id, joinpath(datadir, structure_3D_select))
	plot_profile(b, bt, :T, ax=ax)
	#ax.plot(profile(mean, bt, :log10τ500, :T)..., color="steelblue", lw=5)
	
	# 1D model
	bmarcs = get_snapshot(structure_1D_select, joinpath(datadir, structure_3D_select))
	ax.plot(profile(mean, bmarcs, :log10τ500, :T)..., color="tomato", lw=5)


	# inset
	left, bottom, width, height = [0.15, 0.51, 0.35, 0.3]
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
	ax.set_ylim(llim, 10000)

	ax.set_title(pretty_from_name(structure_3D_select))

	@info "Model resolution (Nx, Ny, Nz): $(size(b))"
	@info "Model resolution (km): $((maximum(b.z)-minimum(b.z)) / size(b, 3) ./1e5)"

	f.savefig(structure_3D1D_figsave_txt*"_3D1D.pdf")
	f
end

# ╔═╡ 99e19910-17a3-4dd2-8ec0-72dd3da2726a


# ╔═╡ 01ec615e-bc60-4e50-a076-692317095591
md"## Two Model comparison"

# ╔═╡ 55b53af5-c8d3-4086-ac92-5995002d8c37
md"""
3D model (A): $(@bind structure_3D_select_A Select(cemp_models, default=first(cemp_models)))
"""

# ╔═╡ 3b07d6ba-0e3f-484c-b3f8-da3840f5d020
md"""
1D model (A): $(@bind structure_1D_select_A Select(marcsmodels[structure_3D_select_A], default=first(marcsmodels[structure_3D_select_A])))
"""

# ╔═╡ 56a77efd-66d5-47aa-934b-99e3d1aff0a8
md"""
3D model (B): $(@bind structure_3D_select_B Select(cemp_models, default=last(cemp_models)))
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

# ╔═╡ 687bd929-41a2-416e-ad7a-681e98096c98
begin
	test_model_marcs = MUST.readdlm("../marcs_cemp_5800_4.5_-5.txt", skipstart=2)
end;

# ╔═╡ 212b2d0b-fbcc-4934-8ebb-627545fea549
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

# ╔═╡ bb9ebb4c-4065-4b5c-82a6-05d88be2c55e


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
	md"Spectra to show (CEMP): $(@bind spectra_3D1D_tag confirm(TextField(length(spectra_3D1D_tag_txt)+5, default=spectra_3D1D_tag_txt)))"
else
	spectra_3D1D_tag = "GBand_DC+0"
	nothing
end

# ╔═╡ c7fb5cc8-5a09-47d7-bb4a-579fb5a085f9
if !isnothing(spectra_CEMP_select)
	md"Spectra to show (scaled-solar): $(@bind spectra_3D1D_ss_tag confirm(TextField(length(spectra_3D1D_tag_txt)+5, default=spectra_3D1D_tag_txt)))"
end

# ╔═╡ 136db4eb-c274-46ef-bac2-03a3bf300cac
if !isnothing(spectra_CEMP_select)
	spectra_3D1D_figsave_txt = "$(spectra_CEMP_select)_spectrum.pdf"
	md"Save figure at: $(@bind spectra_3D1D_figsave confirm(TextField(length(spectra_3D1D_figsave_txt)+1, default=spectra_3D1D_figsave_txt)))"
end

# ╔═╡ 560fcc8a-a9ee-4420-bc97-eb53b9f6de94
if !isnothing(spectra_CEMP_select)
let
	f, ax = plt.subplots(1, 1, figsize=(12, 4))
	
	# 3D CEMP spectrum
	b = get_spectra(snapshot_id, joinpath(datadir, spectra_CEMP_select))
	λ, F = MUST.mean_integrated_flux(b, Symbol(spectra_3D1D_tag))
	ax.plot(λ, F, color="k", lw=1.7, label=L"\rm 3D\ CEMP")

	# scaled-solar
	smod = scaledSolar_models[
		findfirst(
			MUST.same_parameters(
				MUST.ModelInformation(a), 
				MUST.ModelInformation(spectra_CEMP_select),
			) for a in scaledSolar_models
		)
	]
	b = get_spectra(snapshot_id, joinpath(datadir, smod))
	λ, F = MUST.mean_integrated_flux(b, Symbol(spectra_3D1D_ss_tag))
	ax.plot(λ, F, color="steelblue", lw=1.7, label=L"\rm 3D\ \text{scaled-solar}", ls="--")

	
	# 1D model
	try
		bmarcs = get_spectra(spectra_1D_select, joinpath(datadir, spectra_CEMP_select))
		λM, FM = MUST.mean_integrated_flux(bmarcs, Symbol(spectra_3D1D_tag))
		ax.plot(λM, FM, color="tomato", lw=1.7, label=L"\rm 1D")
	catch
		nothing
	end


	ax.set_xlabel(L"\rm wavelength\ [\AA]")
	ax.set_ylabel(L"\rm flux\ [normalized]")
	ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05), framealpha=0)

	ax.set_ylim(-0.14, 1.01)
	ax.set_xlim(4297, 4303)

	cs = b[MUST.spectra_key_from_tag("composition", Symbol(spectra_3D1D_tag))]
	ax.set_title(pretty_from_name(spectra_CEMP_select)*L"\rm,\ [C/Fe]="*"$(MUST.abundance_from_composition_string(cs, "[C/Fe]"))")

	f.savefig(spectra_3D1D_figsave)
	f
end
end

# ╔═╡ 3fef366c-b34b-4747-ac5c-7e3ef5377ef8


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
	#"GBand_DC+06",
	#"GBand_DC+08",
	"GBand_DC-01",
	"GBand_DC-02",
	"GBand_DC-03",
	"GBand_DC-04",
	"GBand_DC-05",
	"GBand_DC-06",
	"GBand_DC-07",
	"GBand_DC-08",
	"GBand_DC-10",
	"GBand_DC-12"
]

# ╔═╡ 4f2b6d24-9157-47ea-a11b-53a30dc34c7a
spectra_reference_interpolation = [
	"GBand_DC+0",
	"GBand_DC-01",
	"GBand_DC-02",
	"GBand_DC+02",
]

# ╔═╡ 4ff6bac5-7f46-4f4b-a221-8b75bdb6b53f
begin
	corrections_1D = []
	corrections_ss_1D = []
	corrections_3D = []
	for i in eachindex(cemp_models_for_each_parameter)
		metallicity_corrections = Dict()
		metallicity_ss_corrections = Dict()
		metallicity_corrections_3D = Dict()
		
		for j in eachindex(cemp_models_for_each_parameter[i])
			# current model
			model3D_name = cemp_models_for_each_parameter[i][j]

			# current metallicity
			feh = MUST.ModelInformation(model3D_name).feh

			metallicity_corrections[feh] = []
			metallicity_ss_corrections[feh] = []
			metallicity_corrections_3D[feh] = []

			# current model spectra
			model3D = get_spectra(snapshot_id, joinpath(datadir, model3D_name))

			# compute corrections w.r.t to the scaled solar model
			scaledsolarmodel = get_spectra(
				snapshot_id, joinpath(datadir, same_scaled_solar(model3D_name))
			)

			@show same_scaled_solar(model3D_name)
			ab_ref, ew_ref, ab, ew, Δab = try
				MUST.curve_of_growth(
					"[C/Fe]", 
					scaledsolarmodel, spectra_name_reference,
					model3D, spectra_name_models,
				)
			catch
				@warn "Problem detected. Try interpolating ther reference"
				ab_test = MUST.abundance_from_tag(model3D, spectra_name_reference, "[C/Fe]")
				MUST.curve_of_growth(
					"[C/Fe]", 
					scaledsolarmodel, spectra_reference_interpolation,
					model3D, spectra_name_models,
					interpolate_reference_to=ab_test
				)
			end
			ac = AbundanceCorrection(
				same_scaled_solar(model3D_name),
				model3D_name,
				ab_ref,
				ew_ref,
				ab,
				ew,
				Δab,
			)
			append!(metallicity_corrections_3D[feh], [ac])

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

				# same for the scaled solar model
				ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
					"[C/Fe]", 
					marcsmodel, spectra_name_reference,
					scaledsolarmodel, spectra_name_models,
				)
				ac = AbundanceCorrection(
					marcsmodelname,
					same_scaled_solar(model3D_name),
					ab_ref,
					ew_ref,
					ab,
					ew,
					Δab,
				)
				append!(metallicity_ss_corrections[feh], [ac])
			end
		end

		# in metallicity_corrections we have stored the abundance corrections
		# for each metallicity of a given set of parameters
		append!(corrections_1D, [metallicity_corrections])
		append!(corrections_ss_1D, [metallicity_ss_corrections])
		append!(corrections_3D, [metallicity_corrections_3D])
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

# ╔═╡ 8d19d381-a19e-4073-a19a-39bd2ff1be7c
begin
	corrections_feh_ss_1D = [
		read_z_corrections(
			corrections_ss_1D[i],
			:closest, vmic=1
		) for i in eachindex(corrections_ss_1D)
	]
	
	ab_feh_ss_1D = [
		read_z_corrections(
			corrections_ss_1D[i],
			:closest, vmic=1, what=:ab
		) for i in eachindex(corrections_ss_1D)
	]
	
	ew_feh_ss_1D = [
		read_z_corrections(
			corrections_ss_1D[i],
			:closest, vmic=1, what=:ew
		) for i in eachindex(corrections_ss_1D)
	]
	
	ab_ref_feh_ss_1D = [
		read_z_corrections(
			corrections_ss_1D[i],
			:closest, vmic=1, what=:ab_ref
		) for i in eachindex(corrections_ss_1D)
	]
	
	ew_ref_feh_ss_1D = [
		read_z_corrections(
			corrections_ss_1D[i],
			:closest, vmic=1, what=:ew_ref
		) for i in eachindex(corrections_ss_1D)
	]
	
	correctionfunctions_feh_ss_1D = [
		MUST.linear_interpolation(
			corrections_feh_ss_1D[i][1], corrections_feh_ss_1D[i][2], extrapolation_bc=MUST.Line()
		) for i in eachindex(corrections_feh_ss_1D)
	]
	@info "1D (scaled-solar) corrections computed."
end;

# ╔═╡ ffd0e8b7-7c9a-4e51-9a1e-1982dc59754c
begin
	corrections_feh_3D = [
		read_z_corrections(
			corrections_3D[i],
			:closest
		) for i in eachindex(corrections_3D)
	]
	
	ab_feh_3D = [
		read_z_corrections(
			corrections_3D[i],
			:closest, what=:ab
		) for i in eachindex(corrections_3D)
	]
	
	ew_feh_3D = [
		read_z_corrections(
			corrections_3D[i],
			:closest, what=:ew
		) for i in eachindex(corrections_3D)
	]
	
	ab_ref_feh_3D = [
		read_z_corrections(
			corrections_3D[i],
			:closest, what=:ab_ref
		) for i in eachindex(corrections_3D)
	]
	
	ew_ref_feh_3D = [
		read_z_corrections(
			corrections_3D[i],
			:closest, what=:ew_ref
		) for i in eachindex(corrections_3D)
	]
	
	correctionfunctions_feh_3D = [
		MUST.linear_interpolation(
			corrections_feh_3D[i][1], corrections_feh_3D[i][2], extrapolation_bc=MUST.Line()
		) for i in eachindex(corrections_feh_3D)
	]
	@info "3D corrections computed."
end;

# ╔═╡ 87aaadf1-2d46-4d14-8723-8c843040d1d1


# ╔═╡ f3dc3826-c36e-4292-9431-a86d4325903c
let
	plt.close()

	f, ax = plt.subplots(1, 3, figsize=(10, 4), sharex=true, sharey=true)
	plt.subplots_adjust(wspace=0)

	
	ax[0].axhline(0.0, color="0.5", lw=1, ls=":")
	ax[1].axhline(0.0, color="0.5", lw=1, ls=":")
	ax[2].axhline(0.0, color="0.5", lw=1, ls=":")
	

	color = ["tomato", "steelblue"]
	ls = ["-", "--"]
	lw = [2.5, 2.5]
	marker = ["s", "s"]
	markerfacecolor = ["tomato", "steelblue"]
	markersize=10

	for i in eachindex(corrections_feh_ss_1D)
		ax[2].plot(
			corrections_feh_ss_1D[i][1], corrections_feh_ss_1D[i][2], 
			color=color[i], marker=marker[i], ls=ls[i], lw=lw[i],
			markersize=markersize, markerfacecolor=markerfacecolor[i], markeredgewidth=2,
			#label=L"\rm T_{eff}="*"$(round(Int,abundance_correction_parameters[i][1]))"*L"\rm \ K"*"\n"*"log(g)="*"$(round(abundance_correction_parameters[i][2], digits=2))"
			label="$(round(Int,abundance_correction_parameters[i][1])) K, $(round(abundance_correction_parameters[i][2], digits=2))"
		)
	end

	for i in eachindex(corrections_feh_3D)
		ax[1].plot(
			corrections_feh_3D[i][1], corrections_feh_3D[i][2], 
			color=color[i], marker=marker[i], ls=ls[i], lw=lw[i],
			markersize=markersize, markerfacecolor=markerfacecolor[i], markeredgewidth=2,
			#label=L"\rm T_{eff}="*"$(round(Int,abundance_correction_parameters[i][1]))"*L"\rm \ K"*"\n"*"log(g)="*"$(round(abundance_correction_parameters[i][2], digits=2))"
			label="$(round(Int,abundance_correction_parameters[i][1])) K, $(round(abundance_correction_parameters[i][2], digits=2))"
		)
	end

	for i in eachindex(corrections_feh_1D)
		ax[0].plot(
			corrections_feh_1D[i][1], corrections_feh_1D[i][2], 
			color=color[i], marker=marker[i], ls=ls[i], lw=lw[i],
			markersize=markersize, markerfacecolor=markerfacecolor[i], markeredgewidth=2,
			#label=L"\rm T_{eff}="*"$(round(Int,abundance_correction_parameters[i][1]))"*L"\rm \ K"*"\n"*"log(g)="*"$(round(abundance_correction_parameters[i][2], digits=2))"
			label="$(round(Int,abundance_correction_parameters[i][1])) K, $(round(abundance_correction_parameters[i][2], digits=2))"
		)
	end

	# panel labels
	ax[0].text(
		0.95,0.95,"A",
		ha="right",va="top", transform=ax[0].transAxes, fontweight="bold"
	)
	ax[1].text(
		0.95,0.95,"B",
		ha="right",va="top", transform=ax[1].transAxes, fontweight="bold"
	)
	ax[2].text(
		0.95,0.95,"C",
		ha="right",va="top", transform=ax[2].transAxes, fontweight="bold"
	)

	ax[1].legend(ncol=1, loc="lower center")
	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[2].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm \Delta A(C)\ [dex]")
	ax[0].set_ylim(-0.97, 0.37)
	ax[0].set_xlim(-6.3, -1.7)

	ax[2].set_title(L"\rm 3D\ \text{scaled-solar} - 1D")
	ax[1].set_title(L"\rm 3D\ CEMP - 3D\ \text{scaled-solar}")
	ax[0].set_title(L"\rm 3D\ CEMP - 1D")

	
	#=ax[0].text(0.07, 0.5, L"\rm 3D\ \text{scaled-solar} - 1D", rotation=90, transform=ax[0].transAxes, ha="left", va="center", alpha=0.7, bbox=Dict("facecolor"=>"none", "edgecolor"=>"k"))
	ax[1].text(0.07, 0.5, L"\rm 3D\ CEMP - 3D\ \text{scaled-solar}", rotation=90, transform=ax[1].transAxes, ha="left", va="center", alpha=0.7, bbox=Dict("facecolor"=>"none", "edgecolor"=>"k"))
	ax[2].text(0.07, 0.5, L"\rm 3D\ CEMP - 1D", rotation=90, transform=ax[2].transAxes, ha="left", va="center", alpha=0.7, bbox=Dict("facecolor"=>"none", "edgecolor"=>"k"))
	=#
	f.savefig("abundance_corrections_split.pdf", bbox_inches="tight")

	f
end

# ╔═╡ 093b3783-af82-4e34-a1e6-65a0df5bb35c


# ╔═╡ c85b36d4-ff22-4ba6-a4d1-c1791fd94ee2
md"## CoG"

# ╔═╡ 3a59080e-f925-48c7-b01d-660c45a7f541
md"Choose parameters: $(@bind i_cog_p_select Select([i=>a for (i, a) in enumerate(abundance_correction_parameters)]))"

# ╔═╡ d8f6a261-88ab-4222-8a51-ce1b359df452
begin
	cog_p_select = abundance_correction_parameters[i_cog_p_select]
	md"Choose metallicity: $(@bind cog_z_select confirm(Slider(ab_feh_1D[i_cog_p_select][1], show_value=true)))"
end

# ╔═╡ 2dd8fd87-f543-413e-b8f3-4187ff67e52c
let
	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	feh_cog_select = ab_feh_1D[i_cog_p_select][1]
	i_cog_feh_select = argmin(abs.(feh_cog_select .- cog_z_select))

	ab = ab_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ew = ew_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ab_ref = ab_ref_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	ew_ref = ew_ref_feh_1D[i_cog_p_select][2][i_cog_feh_select]
	Δab = corrections_feh_1D[i_cog_p_select][2][i_cog_feh_select]

	ab_ref_3D = ab_ref_feh_3D[i_cog_p_select][2][i_cog_feh_select]
	ew_ref_3D = ew_ref_feh_3D[i_cog_p_select][2][i_cog_feh_select]
	Δab_3D = corrections_feh_3D[i_cog_p_select][2][i_cog_feh_select]

	ew_int = MUST.linear_interpolation(ab, ew, extrapolation_bc=MUST.Line())
	da = maximum(ab) - minimum(ab)
	ab_arr = range(minimum(ab) - 0.1*da, maximum(ab) +0.1*da, length=100)
	ew_arr = ew_int.(ab_arr)
	
	ax.plot(
		ab_arr, ew_arr, 
		marker="", color="k", markersize=10, ls="--"
	)
	ax.plot(ab, ew, marker="s", color="k", label="3D CEMP", markersize=11, ls="", markerfacecolor="w", markeredgewidth=2)
	
	ax.plot([ab_ref], [ew_ref], 
		marker="X", color="tomato", ls="",
		label="1D", markersize=12
	)

	ax.plot([ab_ref_3D], [ew_ref_3D], 
		marker="X", color="steelblue", ls="",
		label="3D scaled-solar", markersize=12
	)

	ax.set_ylabel(L"\rm EW\ [\AA]")
	ax.set_xlabel(L"\rm [C/Fe]")

	paras_as_name = L"\rm T_{eff}="*"$(round(Int,cog_p_select[1]))"*L"\rm ,\  log(g)="*"$(round(cog_p_select[2], digits=2))"
	ax.set_title(paras_as_name*L"\rm ,\ [Fe/H] ="*"$(cog_z_select)")
	ax.text(
		0.95, 0.05, 
		L"\rm \Delta [C/Fe]_{1D}="*"$(round(Δab, sigdigits=2))",
		ha="right", va="bottom", transform=ax.transAxes,
		color="tomato",
	)
	ax.text(
		0.95, 0.11, 
		L"\rm \Delta [C/Fe]_{CEMP}="*"$(round(Δab_3D, sigdigits=2))",
		ha="right", va="bottom", transform=ax.transAxes,
		color="steelblue",
	)
	ax.set_xlim(minimum(ab_arr), maximum(ab_arr))

	#ax.axvline(ab_ref, color="tomato", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axvline(ab_ref+ Δab, color="tomato", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axhline(ew_int(ab_ref + Δab), color="tomato", alpha=0.5, ls=":", zorder=0, lw=1.5)
	#ax.axhline(ew_int(ab_ref), color="tomato", alpha=0.5, ls=":", zorder=0, lw=1.5)

	ax.axvline(ab_ref_3D, color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axvline(ab_ref_3D+ Δab_3D, color="steelblue", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axhline(ew_int(ab_ref_3D + Δab_3D), color="steelblue", alpha=0.5, ls=":", zorder=0, lw=1.5)
	ax.axhline(ew_int(ab_ref_3D), color="k", alpha=0.5, ls=":", zorder=0, lw=1.5)


	ax.legend(loc="upper left")
	ax.xaxis.set_major_locator(plt.MaxNLocator(6))

	f
end

# ╔═╡ 75f3327c-00c4-4594-86ac-e2b5237dc0bc


# ╔═╡ 1cf58f8e-e1d0-48eb-9c00-c22fc3144688
md"## 3D - 1D evolution with metallicity"

# ╔═╡ 6c9b537b-60a8-4888-ba5d-be97621e03e8
begin
	color_index = range(0, 100) |> collect	
	cmap = plt.get_cmap("gnuplot2")
	colors = [MUST.pyconvert(Tuple, cmap((i % 10)/ 10)) for i in color_index]
	get_color!(icl) = begin
		c = colors[icl[]]
		icl[] += 1
		c
	end
	cmap
end

# ╔═╡ 537080cd-02e8-4daa-86af-2eee60d1a1ca
let
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	i_para = 1
	xlim = [-3.7, -1.2]

	get_difference(f, bt, bmarcs; same_scale=nothing) = begin
		x3D, y3D = profile(f, bt, :log10τ500, :T)
		x1D, y1D = profile(f, bmarcs, :log10τ500, :T)

		mask3D = (xlim[1] -0.1) .< x3D .< (xlim[2] +0.1)
		mask1D = (xlim[1] -0.1) .< x1D .< (xlim[2] +0.1)

		x3D, y3D = x3D[mask3D], y3D[mask3D]
		x1D, y1D = x1D[mask1D], y1D[mask1D]
		
		same_scale = if isnothing(same_scale) 
			range(
				minimum(x3D), maximum(x3D), length=300
			) |> collect
		else
			same_scale
		end

		ip3D = MUST.linear_interpolation(
			x3D[sortperm(x3D)], y3D[sortperm(x3D)], extrapolation_bc=MUST.Line()
		).(same_scale)

		ip1D = MUST.linear_interpolation(
			x1D[sortperm(x1D)], y1D[sortperm(x1D)], extrapolation_bc=MUST.Line()
		).(same_scale)

		same_scale, ip1D, ip3D
	end

	ic = Ref(1)
	FH = []
	COLS = []
	for mfp in cemp_models_for_each_parameter[i_para]
		mi = get_info_from_name(mfp)
		marcsmi = get_info_from_name.(marcsmodels[mfp])
		mclose = get_closest_model(marcsmi, mi)

		col = get_color!(ic)
		ic[] +=1
		append!(COLS, [col])
		X, Y3D, Y1D = [], [], []
		YSTD3D = []
		bmarcs = get_snapshot(marcsmodels[mfp][mclose], joinpath(datadir, mfp))
		for i in [-1, -3, -5, -7, -9]
			b, bt = get_snapshot(i, joinpath(datadir, mfp))

			same_scale, ip1D, ip3D = if i == -1
				get_difference(mean, bt, bmarcs)
			else
				get_difference(mean, bt, bmarcs, same_scale=X[1])
			end

			_, _, ipstd3D = get_difference(MUST.std, bt, bmarcs)
			

			append!(X, [same_scale])
			append!(Y3D, [ip3D])
			append!(YSTD3D, [ipstd3D])
			append!(Y1D, [ip1D])
		end

		YSTD = [MUST.std([Y3D[i][j] for i in eachindex(X)]) for j in eachindex(X[1])]
		YMEAN = [mean([Y3D[i][j] for i in eachindex(X)]) for j in eachindex(X[1])]
		YMAX = [maximum([Y3D[i][j] for i in eachindex(X)]) for j in eachindex(X[1])]
		YMIN = [minimum([Y3D[i][j] for i in eachindex(X)]) for j in eachindex(X[1])]
		#=ax.fill_between(
			X[1], 
			((Y3D[1] .- YSTD) .- Y1D[1]) ./ (Y1D[1]) .*100,  
			((Y3D[1] .+ YSTD) .- Y1D[1]) ./ (Y1D[1]) .*100
		)=#

		#=ax.fill_between(
			X[1], 
			(YMIN .- Y1D[1]) ./ (Y1D[1]) .*100,  
			(YMAX .- Y1D[1]) ./ (Y1D[1]) .*100,
			alpha=0.5
		)=#


		# we can get the contribution functions for each model like this
		s = get_m3d_spectra(-1, joinpath(datadir, mfp))
		krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
		line = s.run.atom.trans[krs[1]]
		xf, yf = line.get_cntrbf(gmean=true, df=0, norm=true)
		xf = MUST.pyconvert(Array, xf)
		yf = MUST.pyconvert(Array, yf)
		fh = xf[argmax(yf)]
		append!(FH, [fh])
		#ax.axvline(fh, color=col, ls="--", alpha=0.3, lw=5)

		ax.fill_between(
			X[1], 
			((YMEAN .- YSTD) .- Y1D[1]) ./ (Y1D[1]) .*100,  
			((YMEAN .+ YSTD) .- Y1D[1]) ./ (Y1D[1]) .*100,
			color=col,
			alpha=0.6
		)

		#=ax.fill_between(
			X[1], 
			((Y3D[1] .- YSTD3D[1]) .- Y1D[1]) ./ (Y1D[1]) .*100,  
			((Y3D[1] .+ YSTD3D[1]) .- Y1D[1]) ./ (Y1D[1]) .*100,
			color=col,
			alpha=0.2
		)=#
		ax.plot(
			[],
			[],  
			label="[Fe/H] = $(mi.feh)", color=col, 
			lw=5
		)
		
		#=ax.plot(
			X[1], 
			(Y3D[1] .- Y1D[1]) ./ (Y1D[1]) .*100,  
			label="[Fe/H] = $(mi.feh)", color=col, 
			lw=5
		)=#
	end

	ylim = MUST.pyconvert(Array, ax.get_ylim())
	ylim = [ylim[1] - 0.07*(ylim[2]-ylim[1]), max(2, ylim[2] + 0.0*(ylim[2]-ylim[1]))]
	ax.set_ylim(ylim...)
	llim = ylim[1]
	yrange = ylim[2] - ylim[1]
	for i in eachindex(FH)
		ax.hlines(
			llim +0.03*yrange + 0.035*(i-1)*yrange, FH[i], xlim[2], ls="-", color=COLS[i], alpha=0.4, lw=7
		)
		ax.text(
			FH[i]-0.1, llim +0.03*yrange + 0.035*(i-1)*yrange, 
			"$(MUST.ModelInformation(cemp_models_for_each_parameter[i_para][i]).feh)", 
			ha="right", va="center", fontsize="small",
			color=COLS[i]
		)
	end

	ax.axhline(0.0, ls=":", alpha=0.3, color="k")
	ax.set_xlim(xlim...)
	#ax.legend(labelspacing=0.1)

	ax.set_title(
		pretty_from_name_no_feh(cemp_models_for_each_parameter[i_para] |> first)
	)
	ax.set_xlabel(L"\rm optical\ depth\ [\log\tau_{500}]")
	ax.set_ylabel(L"\rm T_{3D\ CEMP}\ /\ T_{1D} - 1\ [\%]")

	tefflogg(name) = begin
		mi = MUST.ModelInformation(name)
		"$(mi.teff)_$(mi.logg)"
	end

	s = tefflogg(cemp_models_for_each_parameter[i_para] |> first)
	s *= "_difference_marcs.pdf"
	f.savefig(s)
	f
end

# ╔═╡ 87366f88-2d35-4f90-bb71-0e7c9bcb0811
molec_abund(s, name) = begin
	totn = MUST.pyconvert(
		Array, 
		s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0]
	)
	x = MUST.pyconvert(Array, s.run.ltau)
	toth = MUST.pyconvert(Array, s.run.get_toth())

	toth = if ndims(toth) < ndims(totn)
		reshape(toth, size(totn)...)
	else
		toth
	end

	x = if ndims(x) < ndims(totn)
		reshape(x, size(totn)...)
	else
		x
	end
	
	y = log10.(totn./ toth) .+ 12

	x = [MUST.mean(x[:, :, z]) for z in axes(x, 3)]
	y = [MUST.mean(y[:, :, z]) for z in axes(y, 3)]

	x, y
end

# ╔═╡ 9a3df65f-d469-42d9-b5c4-7cf66f1ed5f6
get_contr_function(s; k=1) = begin
	krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
	line = s.run.atom.trans[krs[k]]
	xf, yf = line.get_cntrbf(gmean=true, df=0, norm=true)
	xf = MUST.pyconvert(Array, xf)
	yf = MUST.pyconvert(Array, yf)

	xf, yf
end

# ╔═╡ 46d9bb97-1fd7-4d57-bcaa-8123754f195c
begin
	data_3D1DForm = []
	for i_para in eachindex(cemp_models_for_each_parameter)
		FH3D_3D1DForm = []
		FH1D_3D1DForm = []
		T3D_3D1DForm = []
		T1D_3D1DForm = []
		D3D_3D1DForm = []
		D1D_3D1DForm = []
		N3D_3D1DForm = []
		N1D_3D1DForm = []
		FEH_3D1DForm = []
		for mfp in cemp_models_for_each_parameter[i_para]
			mi = get_info_from_name(mfp)
			marcsmi = get_info_from_name.(marcsmodels[mfp])
			mclose = get_closest_model(marcsmi, mi)
			
			bmarcs = get_snapshot(marcsmodels[mfp][mclose], joinpath(datadir, mfp))
			b, bt = get_snapshot(snapshot_id, joinpath(datadir, mfp))
	
			# we can get the contribution functions for each model like this
			s = get_m3d_spectra(-1, joinpath(datadir, mfp))
			xf, yf = get_contr_function(s)
			xCH3D, NCH3D = molec_abund(s, "CH")
			fh = xf[argmax(yf)]
			append!(FH3D_3D1DForm, [fh])
	
			# for the marcs model
			s = get_m3d_spectra(marcsmodels[mfp][mclose], joinpath(datadir, mfp))
			xCH1D, NCH1D = molec_abund(s, "CH")
			xf, yf = get_contr_function(s)
			fh = xf[argmax(yf)]
			append!(FH1D_3D1DForm, [fh])
	
			# compute the mean temperature at that formation height
			b_ip = MUST.interpolate_to(b, :T, logspace=true, τ500=FH3D_3D1DForm|>last)
			append!(T3D_3D1DForm, [mean(b_ip[:T])])
			b_ip = MUST.interpolate_to(bmarcs, :T, logspace=true, τ500=FH1D_3D1DForm|>last)
			append!(T1D_3D1DForm, [mean(b_ip[:T])])

			b_ip = MUST.interpolate_to(b, :d, logspace=true, τ500=FH3D_3D1DForm|>last)
			append!(D3D_3D1DForm, [mean(b_ip[:log10d])])
			b_ip = MUST.interpolate_to(bmarcs, :d, logspace=true, τ500=FH1D_3D1DForm|>last)
			append!(D1D_3D1DForm, [mean(b_ip[:log10d])])

			# compute the CH number density
			n_ip = MUST.linear_interpolation(xCH3D[sortperm(xCH3D)], NCH3D[sortperm(xCH3D)]).(FH3D_3D1DForm|>last)
			append!(N3D_3D1DForm, [n_ip])

			n_ip = MUST.linear_interpolation(xCH1D[sortperm(xCH1D)], NCH1D[sortperm(xCH1D)]).(FH1D_3D1DForm|>last)
			append!(N1D_3D1DForm, [n_ip])
			
			append!(FEH_3D1DForm, [mi.feh])
		end

		append!(data_3D1DForm, [(
			FH3D_3D1DForm,
			FH1D_3D1DForm,
			T3D_3D1DForm,
			T1D_3D1DForm,
			D3D_3D1DForm,
			D1D_3D1DForm,
			N3D_3D1DForm,
			N1D_3D1DForm,
			FEH_3D1DForm
		)])
	end
end

# ╔═╡ b0e95f8f-0059-4575-bab6-80f086ba3e4a
let
	plt.close()
	f, ax = plt.subplots(3, 1, figsize=(5., 9), sharex=true)
	plt.subplots_adjust(wspace=0, hspace=0)

	ax[0].axhline(0.0, ls=":", alpha=0.3, color="k")
	ax[1].axhline(0.0, ls=":", alpha=0.3, color="k")
	ax[2].axhline(1.0, ls=":", alpha=0.3, color="k")

	#i_para = 1
	colors = ["tomato", "steelblue"]
	ls = ["-", "--"]
	for i_para in eachindex(cemp_models_for_each_parameter)
		FH3D_3D1DForm,FH1D_3D1DForm,T3D_3D1DForm,T1D_3D1DForm,D3D_3D1DForm,D1D_3D1DForm,N3D_3D1DForm,N1D_3D1DForm,FEH_3D1DForm = data_3D1DForm[i_para]
		#=ax.plot(
			FEH, (exp.(diss_CH ./(MUST.KBoltzmann .* T3D)) ./ exp.(diss_CH ./ (MUST.KBoltzmann .*T1D)) .* T3D.^(-3/2) ./ T1D.^(-3/2)), #./ T1D .* 100, 
			color=colors[i_para], marker="s", markersize=12, lw=3, markeredgewidth=3, ls=ls[i_para],
			label=pretty_from_name_no_feh(cemp_models_for_each_parameter[i_para] |> first)
		)=#

		mi = MUST.ModelInformation(cemp_models_for_each_parameter[i_para]|>first)
		name = "$(round(Int,mi.teff)) K, $(round(mi.logg, digits=2))"

		d_CH = 5.527509387E-12
		n_CH(T) = log10(exp(d_CH/(MUST.KBoltzmann * T)))
		ax[0].plot(
			FEH_3D1DForm, (T3D_3D1DForm .- T1D_3D1DForm) ./ T1D_3D1DForm .*100, 
			color=colors[i_para], marker="s", markersize=10, lw=2.5, markeredgewidth=3, ls=ls[i_para],
			label=name
		)
		ax[1].plot(
			FEH_3D1DForm, (D3D_3D1DForm .- D1D_3D1DForm), 
			color=colors[i_para], marker="s", markersize=10, lw=2.5, markeredgewidth=3, ls=ls[i_para],
			label=name
		)
		ax[2].plot(
			FEH_3D1DForm, exp10.(N3D_3D1DForm .- N1D_3D1DForm), 
			color=colors[i_para], marker="s", markersize=10, lw=2.5, markeredgewidth=3, ls=ls[i_para],
			label=name
		)
	end

	# panel labels
	ax[0].text(
		0.95,0.95,"A",
		ha="right",va="top", transform=ax[0].transAxes, fontweight="bold", fontsize="large"
	)
	ax[1].text(
		0.95,0.95,"B",
		ha="right",va="top", transform=ax[1].transAxes, fontweight="bold", fontsize="large"
	)
	ax[2].text(
		0.95,0.95,"C",
		ha="right",va="top", transform=ax[2].transAxes, fontweight="bold", fontsize="large"
	)
	
	ax[1].legend(labelspacing=0.02, loc="center left", bbox_to_anchor=(0.05, 0.65))

	ax[0].set_ylim(-30, 4)
	ax[1].set_ylim(-0.77, 0.13)
	ax[2].set_ylim(0.3, 4.5)

	ax[2].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm T_{3D}\ /\ T_{1D} - 1\ [\%]")
	#ax.set_ylabel(L"\rm T_{3D\ CEMP} - T_{1D}\ [K]")
	ax[1].set_ylabel(L"\rm \log_{10}(\rho_{3D}\ /\ \rho_{1D})")
	#ax[1].set_ylabel(L"\rm \log_{10}(\rho_{3D}) - \log_{10}(\rho_{1D})")
	#ax[2].set_ylabel(L"\rm A(CH)_{3D} - A(CH)_{1D}")
	#ax[2].set_ylabel(L"\rm log(n_{CH,3D}\ /\ n_{CH,1D})")
	#ax[2].set_ylabel(L"\rm log \left( \left[ \frac{n_{CH}}{n_H}\right]_{3D}\ /\ \left[ \frac{n_{CH}}{n_H}\right]_{1D} \right)")
	#ax[2].set_ylabel(L"\rm \log\left(\frac{n_{CH}}{n_{H}} \right)_{3D} - \log\left(\frac{n_{CH}}{n_{H}} \right)_{1D} ")
	ax[2].set_ylabel(L"\rm n^*_{3D}\ /\ n^*_{1D}")

	s = "differences_marcs_lte.pdf"
	f.savefig(s)
	
	f
end

# ╔═╡ ef9fef62-5e50-4a1d-b1a9-4db7289663a7


# ╔═╡ 331f1016-df68-4330-8983-a2db6219f966
md"# Molecule formation"

# ╔═╡ 924ea3f0-c822-44dc-9faa-1629c14f5232
md"## Number densities"

# ╔═╡ 494823c4-50b7-4e23-894c-caafaab48b29
"""
	plot_abund_3D(s, name; ax, cmap="Greys")

Plot abundances as a function of optical depth (500).
"""
plot_abund_3D(s, name; ax, cmap="Greys", alpha=0.8) = begin
	totn = MUST.pyconvert(
		Array, 
		s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0]
	)
	toth = MUST.pyconvert(
		Array, s.run.get_toth()
	)
	totn[totn .<= 0.0] .= minimum(totn[totn .> 0.0])
	C1 = log10.(totn ./ toth) .+ 12

	#= Plot =#
	x = MUST.pyconvert(Array, s.run.ltau)
	y = C1
	
	H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=150)

	xe = MUST.pyconvert(Array, xedges)
	ye = MUST.pyconvert(Array, yedges)
	H = MUST.pyconvert(Array, H.T)

	xx, yy = meshgrid(xe, ye)
	
	#=ax.imshow(
		H, interpolation="bicubic", origin="lower", aspect="auto",
		extent=[xe[1], xe[end], ye[1], ye[end]], cmap=cmap, norm=matplotlib.colors.LogNorm()
	)=#

	ax.imshow(
		H, 
		origin="lower",
		interpolation = "bicubic", 
		extent=[minimum(xe), maximum(xe), minimum(ye), maximum(ye)],
		cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=1), aspect="auto",
		rasterized=true,alpha=alpha
	)
	#ax.pcolormesh(xx, yy, H', cmap=cmap, norm=matplotlib.colors.LogNorm(), rasterized=true)
end

# ╔═╡ 1bd82fae-afb3-4197-999e-61fba6ccfe71
"""
	plot_abund_av_3D(s, name; ax, kwargs...)

Plot abundances (averages) as a function of optical depth (500).
"""
plot_abund_av_3D(s, name; ax, relative=true, kwargs...) = begin
	try
		totn = MUST.pyconvert(
			Array, 
			s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0]
		)
		C1 = toth = MUST.pyconvert(Array, s.run.get_toth())
		C1 = relative ? log10.(totn ./ toth) .+ 12 : log10.(totn)

		#= Plot =#
		x = MUST.pyconvert(Array, s.run.ltau)
		y = C1

		x = [MUST.mean(x[:, :, z]) for z in axes(x, 3)]
		y = [MUST.mean(y[:, :, z]) for z in axes(y, 3)]
		
		ax.plot(x, y; kwargs...)

		x, y
	catch
		totn = MUST.pyconvert(
			Array, 
			s.run.read_patch_save(
				name, lazy=false, concat=true, fdim=0, zdim=2
			)[0][0][0]
		)
		toth = MUST.pyconvert(Array, s.run.get_toth())
		C1 = relative ? log10.(totn ./ toth) .+ 12 : log10.(totn)

		#= Plot =#
		x = MUST.pyconvert(Array, s.run.ltau)
		y = C1

		ax.plot(x, y; kwargs...)

		x, y
	end
end

# ╔═╡ 731fadcd-c03a-4f7d-8f49-f7c2c29c099b


# ╔═╡ 981e31b0-7073-431c-b1a4-d82fd6e57350
md"""
3D model: $(@bind mol_3D_select Select(cemp_models, default=first(cemp_models)))
"""

# ╔═╡ 7e487942-9681-4441-8969-050c8328a19a
md"""
1D model: $(@bind mol_1D_select Select(marcsmodels[mol_3D_select], default=first(marcsmodels[mol_3D_select])))
"""

# ╔═╡ 08fba04e-f041-4c49-865e-ec254a99fc98
if !isnothing(mol_3D_select)
	mol_3D_figsave_txt = "$(mol_3D_select)_molec"
	md"Save figure at: $(@bind mol_3D_figsave confirm(TextField(length(mol_3D_figsave_txt)+5, default=mol_3D_figsave_txt)))"
end

# ╔═╡ 0e7e6f79-7225-4d29-8693-b2b93ddf4812
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 5))

	lw_CH = 4.0
	lw_CII = 10.0
	alpha_CII = 0.25
	relative=false

	# 1D MARCS model
	s = get_m3d_spectra(mol_1D_select, joinpath(datadir, mol_3D_select))
	xch, ych = plot_abund_av_3D(s, "CH", ax=ax, label=L"\rm 1D", color="tomato", ls="-", lw=lw_CH, relative=relative)
	#plot_abund_av_3D(s, "CO", ax=ax, label=L"\rm CO\ -\ 1D", color="lime", ls=":")
	#plot_abund_av_3D(s, "C_I", ax=ax, color="tomato", ls="-")
	xcii, ycii = plot_abund_av_3D(s, "C_II", ax=ax, color="tomato", ls="-", lw=lw_CII, alpha=alpha_CII, relative=relative)

	xloc_ch = -1.5
	yloc_ch = ych[argmin(abs.(xch.-xloc_ch))] - 0.1
	ax.text(xloc_ch, yloc_ch, L"\mathbf{CH}", color="tomato", ha="left", va="top", fontsize=15)

	xloc_cii = -1.5
	yloc_cii = ycii[argmin(abs.(xcii.-xloc_cii))] + 0.4
	ax.text(xloc_cii, yloc_cii, L"\mathbf{C\ II}", color="tomato", ha="left", va="bottom", fontsize=15, alpha=alpha_CII*2)


	
	
	# scaled_solar model
	s = get_m3d_spectra(-1, joinpath(datadir, same_scaled_solar(mol_3D_select)))
	#plot_abund_3D(s, "CH", ax=ax, cmap="Reds", alpha=0.6)
	xch, ych = plot_abund_av_3D(s, "CH", ax=ax, label=L"\rm 3D\ \text{scaled-solar}", color="steelblue", lw=lw_CH, relative=relative)
	#plot_abund_3D(s, "CO", ax=ax, cmap="Greys", alpha=0.6)
	#plot_abund_av_3D(s, "CO", ax=ax, label=L"\rm CO\ -\ 3D\ \text{scaled-solar}", color="lime")
	
	#plot_abund_av_3D(s, "C_I", ax=ax, color="steelblue", ls="-")
	xcii, ycii = plot_abund_av_3D(s, "C_II", ax=ax, color="steelblue", ls="-", lw=lw_CII, alpha=alpha_CII, relative=relative)

	xloc_ch = -3.9
	yloc_ch = ych[argmin(abs.(xch.-xloc_ch))] + 0.2
	ax.text(xloc_ch, yloc_ch, L"\mathbf{CH}", color="steelblue", ha="left", va="bottom", fontsize=15)

	xloc_cii = -3.9
	yloc_cii = ycii[argmin(abs.(xcii.-xloc_cii))] - 0.2
	ax.text(xloc_cii, yloc_cii, L"\mathbf{C\ II}", color="steelblue", ha="left", va="top", fontsize=15, alpha=alpha_CII*2)
	

	
	# CEMP model
	s = get_m3d_spectra(-1, joinpath(datadir, mol_3D_select))
	#plot_abund_3D(s, "CH", ax=ax, cmap="Reds", alpha=0.6)
	xch, ych = plot_abund_av_3D(s, "CH", ax=ax, label=L"\rm 3D\ CEMP", color="steelblue", ls="", marker="x", lw=lw_CH, relative=relative)
	#plot_abund_av_3D(s, "C_I", ax=ax, color="steelblue", ls="", marker="x")
	xcii, ycii = plot_abund_av_3D(s, "C_II", ax=ax, color="steelblue", ls="", lw=lw_CII, marker="x", alpha=alpha_CII, markersize=15, relative=relative)

	#=xloc_ch = -3.9
	yloc_ch = ych[argmin(abs.(xch.-xloc_ch))] - 0.1
	ax.text(xloc_ch, yloc_ch, L"\mathbf{CH}", color="steelblue", ha="left", va="top", fontsize=15)

	xloc_cii = -3.9
	yloc_cii = ycii[argmin(abs.(xcii.-xloc_cii))] + 0.1
	ax.text(xloc_cii, yloc_cii, L"\mathbf{C\ II}", color="steelblue", ha="left", va="bottom", fontsize=15, alpha=alpha_CII*2)
	=#

	ax.legend(labelspacing=0.1, loc="lower right", ncol=1)

	ax.set_xlim(-4, 0)
	if relative
		ax.set_ylabel(L"\rm log_{10}(n_X\ /\ n_H) + 12")
		ax.set_ylim(0.1, 4.9)
	else
		ax.set_ylabel(L"\rm log_{10}(n_X)\ [cm^{-3}]")
		ax.set_ylim(3.3, 9.7)
	end
	ax.set_xlabel(L"\rm optical\ depth\ [\log_{10} \tau_{500}]")
	#ax.set_ylim(3.2, 10.2)

	ax.set_title(pretty_from_name(mol_3D_select))

	f.savefig(mol_3D_figsave_txt*"3D.pdf", bbox_inches="tight")
	
	f
end

# ╔═╡ 217e0937-4b2d-4948-9a72-c236a8fc5c5a
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(7.5, 5))

	# scaled_solar model
	s = get_m3d_spectra(-1, joinpath(datadir, same_scaled_solar(mol_3D_select)))
	#plot_abund_3D(s, "CH", ax=ax, cmap="Reds", alpha=0.6)
	plot_abund_av_3D(s, "CH", ax=ax, label=L"\rm CH\ -\ 3D\ \text{scaled-solar}", color="tomato")
	#plot_abund_3D(s, "CO", ax=ax, cmap="Greys", alpha=0.6)
	#plot_abund_av_3D(s, "CO", ax=ax, label=L"\rm CO\ -\ 3D\ \text{scaled-solar}", color="lime")
	
	plot_abund_av_3D(s, "C_I", ax=ax, label=L"\rm C\ I\ -\ 3D\ \text{scaled-solar}", color="steelblue")
	plot_abund_av_3D(s, "C_II", ax=ax, label=L"\rm C\ II\ -\ 3D\ \text{scaled-solar}", color="k")

	
	# 1D MARCS model
	s = get_m3d_spectra(mol_1D_select, joinpath(datadir, mol_3D_select))
	plot_abund_av_3D(s, "CH", ax=ax, label=L"\rm CH\ -\ 1D", color="tomato", ls=":")
	#plot_abund_av_3D(s, "CO", ax=ax, label=L"\rm CO\ -\ 1D", color="lime", ls=":")
	plot_abund_av_3D(s, "C_I", ax=ax, label=L"\rm C\ I\ -\ 1D", color="steelblue", ls=":")
	plot_abund_av_3D(s, "C_II", ax=ax, label=L"\rm C\ II\ -\ 1D", color="k", ls=":")

	ax.legend(labelspacing=0.1, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2)

	ax.set_xlim(-4, 0)
	ax.set_ylabel("[X/H]")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_title(pretty_from_name(mol_3D_select))
	ax.set_ylim(0, 6.9)

	f.savefig(mol_3D_figsave_txt*"1D.pdf", bbox_inches="tight")

	f
end

# ╔═╡ b942b4f4-6979-42ee-8f18-b9e7b64c0379


# ╔═╡ a84ea668-6701-4e31-9f04-cb0feb70290d


# ╔═╡ 7dae5609-0148-438f-a207-98da5c9173d4
md"## Constribution functions"

# ╔═╡ e4232daa-7174-493d-968e-21a012204fcb
begin
	"""
		plot_contr_av(s; ax, kwargs...)
	
	Plot contribution function (averages) as a function of optical depth (500).
	"""
	plot_contr_av(s; ax, kwargs...) = begin
		krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
		line = s.run.atom.trans[krs[1]]
		x, y = line.get_cntrbf(gmean=true, df=0, norm=true)
	
		x = MUST.pyconvert(Array, x)
		y = MUST.pyconvert(Array, y)
	
		ax.plot(x, y; kwargs...)

		x, y
	end
	
	"""
		plot_contr_1D(s; ax, kwargs...)
	
	Plot contribution function (averages) as a function of optical depth (500).
	"""
	plot_contr_1D(s; ax, kwargs...) = begin
		krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
		line = s.run.atom.trans[krs[1]]
		x, y = line.get_cntrbf(df=0, norm=true)

		x = MUST.pyconvert(Array, x)
		y = MUST.pyconvert(Array, y)
	
		ax.plot(x, y; kwargs...)

		x, y
	end
end

# ╔═╡ e502770d-d765-4757-a075-3fadcdd88eed
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	# scaled_solar model
	s = get_m3d_spectra(-1, joinpath(datadir, same_scaled_solar(mol_3D_select)))
	x, y = plot_contr_av(s, ax=ax, color="steelblue", label=L"\rm 3D\ \text{scaled-solar}")
	formHeight_3DSS = x[argmax(y)]
	
	# CEMP model
	s = get_m3d_spectra(-1, joinpath(datadir, mol_3D_select))
	x, y = plot_contr_av(s, ax=ax, color="steelblue", ls="--", label=L"\rm 3D\ CEMP")
	formHeight_3DCEMP = x[argmax(y)]

	# 1D model
	s = get_m3d_spectra(mol_1D_select, joinpath(datadir, mol_3D_select))
	x, y = plot_contr_1D(s, ax=ax, color="tomato", ls="-", label=L"\rm 1D")
	formHeight_1D = x[argmax(y)]

	ax.legend(labelspacing=0.1, loc="upper left", ncol=1)

	ax.set_ylabel("contribution function")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")

	ax.set_title(pretty_from_name(mol_3D_select))

	ax.set_xlim(-6, 2)
	ax.set_ylim(-0.1, 1.1)

	
	@info "Maxima at:" formHeight_3DCEMP formHeight_3DSS formHeight_1D

	f.savefig(mol_3D_figsave_txt*"contr.pdf", bbox_inches="tight")
	
	f
end

# ╔═╡ f64c0cd9-c0ec-44e3-9e0c-440b6fdc23e8


# ╔═╡ 635d415f-aa91-4d25-8162-a550119e3f12
md"# Stellar parameter dependence"

# ╔═╡ 2caf0768-3202-4d16-8c24-4f9a6433076e
md"To investigate the effect of stellar parameters one can have a look at models with same log(g), but different temperature. Two example models are prepared for this purpose"

# ╔═╡ ac81d133-8894-46da-b85c-ab5d580ecd12
switch_dir = @in_dispatch "CEMP_models3_switch"

# ╔═╡ 01752527-37fd-4175-9982-02877ca4ae69
begin
	hot_big_CEMP_name = "CEMP_t57.50g30.00m-5.000_v1.0"
	cold_small_CEMP_name = "CEMP_t52.50g45.00m-5.000_v1.0"
	hot_small_CEMP_name = "CEMP_t57.50g45.00m-5.000_v1.0"
	cold_big_CEMP_name = "CEMP_t52.50g30.00m-5.000_v1.0"

	hot_big_ss_name = "ScaledSolar_t57.50g30.00m-5.000_v1.0"
	cold_small_ss_name = "ScaledSolar_t52.50g45.00m-5.000_v1.0"
	hot_small_ss_name = "ScaledSolar_t57.50g45.00m-5.000_v1.0"
	cold_big_ss_name = "ScaledSolar_t52.50g30.00m-5.000_v1.0"

	marcs_hot_big = "p5750_g+3.0_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
	marcs_cold_small = "p5250_g+4.5_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
	marcs_hot_small = "p5750_g+4.5_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
	marcs_cold_big = "p5250_g+3.0_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
end

# ╔═╡ 2d18e394-484d-4aef-bcb5-c45627fcf824
let
	plt.close()

	f, ax = plt.subplots(2, 1, figsize=(12, 5), sharex=true, sharey=true)
	plt.subplots_adjust(hspace=0, wspace=0)
	
	hot_big_CEMP = get_spectra(-1, joinpath(switch_dir, hot_big_CEMP_name))
	cold_small_CEMP = get_spectra(-1, joinpath(switch_dir, cold_small_CEMP_name))
	hot_small_CEMP = get_spectra(-1, joinpath(datadir, hot_small_CEMP_name))
	cold_big_CEMP = get_spectra(-1, joinpath(datadir, cold_big_CEMP_name))

	λ, F = MUST.mean_integrated_flux(hot_big_CEMP, Symbol(spectra_3D1D_tag))
	ax[0].plot(
		λ, F, 
		color="steelblue", lw=1.7, label=pretty_from_name_short(hot_big_CEMP_name)
	)
	λ, F = MUST.mean_integrated_flux(cold_big_CEMP, Symbol(spectra_3D1D_tag))
	ax[0].plot(
		λ, F, 
		color="tomato", lw=1.7, label=pretty_from_name_short(cold_big_CEMP_name)
	)

	# for legend
	ax[0].plot([],[],lw=1.7,ls="--",color="steelblue",label=pretty_from_name_short(hot_small_CEMP_name))
	ax[0].plot([],[],lw=1.7,ls="--",color="tomato",label=pretty_from_name_short(cold_small_CEMP_name))

	
	λ, F = MUST.mean_integrated_flux(hot_small_CEMP, Symbol(spectra_3D1D_tag))
	ax[1].plot(
		λ, F, 
		color="steelblue", ls ="--",lw=1.7, label=pretty_from_name_short(hot_small_CEMP_name)
	)
	λ, F = MUST.mean_integrated_flux(cold_small_CEMP, Symbol(spectra_3D1D_tag))
	ax[1].plot(
		λ, F, 
		color="tomato",ls ="--", lw=1.7, label=pretty_from_name_short(cold_small_CEMP_name)
	)

	# abundance corrections
	model3D = get_spectra(-1, joinpath(switch_dir, hot_big_CEMP_name))
	model3D_ss = get_spectra(-1, joinpath(switch_dir, hot_big_ss_name))
	marcsmodel = get_spectra(
		marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name)
	)
	ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
		"[C/Fe]", 
		model3D_ss, "GBand_DC+0",
		model3D, [
			"GBand_DC+0",
			"GBand_DC+02",
			"GBand_DC+04",
			#"GBand_DC+06",
			#"GBand_DC+08",
			"GBand_DC-01",
			"GBand_DC-02",
			"GBand_DC-03",
			"GBand_DC-04",
			"GBand_DC-05",
			"GBand_DC-06",
			"GBand_DC-07",
			"GBand_DC-08",
			"GBand_DC-10",
			"GBand_DC-12"
		],
	)
	@info "$(pretty_from_name_short(hot_big_CEMP_name)), 3D CEMP vs. 3D non-CEMP: $(Δab)"
	ab_ref, ew_ref, ab, ew, Δab = MUST.curve_of_growth(
		"[C/Fe]", 
		marcsmodel, "GBand_DC+0",
		model3D_ss, [
			"GBand_DC+0",
			"GBand_DC+02",
			"GBand_DC+04",
			#"GBand_DC+06",
			#"GBand_DC+08",
			"GBand_DC-01",
			"GBand_DC-02",
			"GBand_DC-03",
			"GBand_DC-04",
			"GBand_DC-05",
			"GBand_DC-06",
			"GBand_DC-07",
			"GBand_DC-08",
			"GBand_DC-10",
			"GBand_DC-12"
		],
	)
	@info "$(pretty_from_name_short(hot_big_CEMP_name)), 3D non-CEMP vs. 1D: $(Δab)"


	ax[0].set_ylim(-0.1, 1.1)
	ax[1].set_xlabel(L"\rm wavelength\ [\AA]")
	ax[1].set_ylabel(L"\rm flux\ [normalized]")
	ax[0].set_ylabel(L"\rm flux\ [normalized]")

	ax[0].legend(loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.9))

	f
end

# ╔═╡ 4b18b95c-e574-468a-92ce-73931621e442
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	hot_big_CEMP = get_snapshot(-1, joinpath(switch_dir, hot_big_CEMP_name)) |>last
	#cold_small_CEMP = get_snapshot(-1, joinpath(switch_dir, cold_small_CEMP_name)) |>last
	#hot_small_CEMP = get_snapshot(-1, joinpath(datadir, hot_small_CEMP_name))
	cold_big_CEMP = get_snapshot(-1, joinpath(datadir, cold_big_CEMP_name)) |> last

	hot_big_marcs = get_snapshot(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	#cold_small_marcs = get_snapshot(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	#hot_small_marcs = get_snapshot(marcs_hot_small, joinpath(datadir, hot_small_CEMP_name))
	cold_big_marcs = get_snapshot(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))

	ax.plot(
		profile(mean, hot_big_CEMP, :log10τ500, :T)..., label=pretty_from_name_short(hot_big_CEMP_name)*", 3D CEMP", color="k"
	)
	ax.plot(
		profile(mean, cold_big_CEMP, :log10τ500, :T)..., label=pretty_from_name_short(cold_big_CEMP_name)*", 3D CEMP", color="tomato"
	)
	
	ax.plot(
		profile(mean, hot_big_marcs, :log10τ500, :T)..., label=pretty_from_name_short(hot_big_CEMP_name)*", MARCS", color="k", ls="--"
	)
	ax.plot(
		profile(mean, cold_big_marcs, :log10τ500, :T)..., label=pretty_from_name_short(cold_big_CEMP_name)* ", MARCS", color="tomato", ls="--"
	)

	ax.legend()

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(2800, 9000)
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_ylabel(L"\rm temperature\ [K]")
	
	f
end

# ╔═╡ af94af20-6618-448e-b146-b3ba1b1a20f1
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	hot_big_CEMP = get_snapshot(-1, joinpath(switch_dir, hot_big_CEMP_name)) |>last
	#cold_small_CEMP = get_snapshot(-1, joinpath(switch_dir, cold_small_CEMP_name)) |>last
	hot_small_CEMP = get_snapshot(-1, joinpath(datadir, hot_small_CEMP_name)) |> last
	#cold_big_CEMP = get_snapshot(-1, joinpath(datadir, cold_big_CEMP_name)) |> last

	hot_big_marcs = get_snapshot(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	#cold_small_marcs = get_snapshot(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	hot_small_marcs = get_snapshot(marcs_hot_small, joinpath(datadir, hot_small_CEMP_name))
	#cold_big_marcs = get_snapshot(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))

	ax.plot(
		profile(mean, hot_big_CEMP, :log10τ500, :T)..., label=pretty_from_name_short(hot_big_CEMP_name)*", 3D CEMP", color="k"
	)
	ax.plot(
		profile(mean, hot_small_CEMP, :log10τ500, :T)..., label=pretty_from_name_short(hot_small_CEMP_name)*", 3D CEMP", color="tomato"
	)
	
	ax.plot(
		profile(mean, hot_big_marcs, :log10τ500, :T)..., label=pretty_from_name_short(hot_big_CEMP_name)*", MARCS", color="k", ls="--"
	)
	ax.plot(
		profile(mean, hot_small_marcs, :log10τ500, :T)..., label=pretty_from_name_short(hot_small_CEMP_name)* ", MARCS", color="tomato", ls="--"
	)

	ax.axvline(log10(2/3), color="0.5", alpha=0.5, ls=":")

	ax.legend()

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(2800, 9000)
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_ylabel(L"\rm temperature\ [K]")
	
	f
end

# ╔═╡ 22932761-ae92-4cd6-b828-07e3f86a4d0f
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	hot_big_CEMP = get_snapshot(-1, joinpath(switch_dir, hot_big_CEMP_name)) |>last
	#cold_small_CEMP = get_snapshot(-1, joinpath(switch_dir, cold_small_CEMP_name)) |>last
	#hot_small_CEMP = get_snapshot(-1, joinpath(datadir, hot_small_CEMP_name))
	cold_big_CEMP = get_snapshot(-1, joinpath(datadir, cold_big_CEMP_name)) |> last

	hot_big_marcs = get_snapshot(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	#cold_small_marcs = get_snapshot(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	#hot_small_marcs = get_snapshot(marcs_hot_small, joinpath(datadir, hot_small_CEMP_name))
	cold_big_marcs = get_snapshot(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))

	ax.plot(
		profile(mean, hot_big_CEMP, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_big_CEMP_name)*", 3D CEMP", color="k"
	)
	ax.plot(
		profile(mean, cold_big_CEMP, :log10τ500, :log10d)..., label=pretty_from_name_short(cold_big_CEMP_name)*", 3D CEMP", color="tomato"
	)
	
	ax.plot(
		profile(mean, hot_big_marcs, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_big_CEMP_name)*", MARCS", color="k", ls="--"
	)
	ax.plot(
		profile(mean, cold_big_marcs, :log10τ500, :log10d)..., label=pretty_from_name_short(cold_big_CEMP_name)* ", MARCS", color="tomato", ls="--"
	)

	ax.legend()

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(-10.2, -5.9)
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_ylabel(L"\rm density\ [g\ cm^{-3}]")
	
	f
end

# ╔═╡ a6ee120f-0dfc-433a-8c9c-0f9f828ffe29
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	hot_big_CEMP = get_snapshot(-1, joinpath(switch_dir, hot_big_CEMP_name)) |>last
	#cold_small_CEMP = get_snapshot(-1, joinpath(switch_dir, cold_small_CEMP_name)) |>last
	hot_small_CEMP = get_snapshot(-1, joinpath(datadir, hot_small_CEMP_name)) |> last
	#cold_big_CEMP = get_snapshot(-1, joinpath(datadir, cold_big_CEMP_name)) |> last

	hot_big_marcs = get_snapshot(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	#cold_small_marcs = get_snapshot(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	hot_small_marcs = get_snapshot(marcs_hot_small, joinpath(datadir, hot_small_CEMP_name))
	#cold_big_marcs = get_snapshot(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))

	ax.plot(
		profile(mean, hot_big_CEMP, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_big_CEMP_name)*", 3D CEMP", color="k"
	)
	ax.plot(
		profile(mean, hot_small_CEMP, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_small_CEMP_name)*", 3D CEMP", color="tomato"
	)
	
	ax.plot(
		profile(mean, hot_big_marcs, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_big_CEMP_name)*", MARCS", color="k", ls="--"
	)
	ax.plot(
		profile(mean, hot_small_marcs, :log10τ500, :log10d)..., label=pretty_from_name_short(hot_small_CEMP_name)* ", MARCS", color="tomato", ls="--"
	)

	ax.axvline(log10(2/3), color="0.5", alpha=0.5, ls=":")

	ax.legend()

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(-10.2, -5.9)
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")
	ax.set_ylabel(L"\rm density\ [g\ cm^{-3}]")
	
	f
end

# ╔═╡ 421da8e3-f710-46f1-a602-bdc4a3266a26
let
	plt.close()

	f, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=true)
	plt.subplots_adjust(hspace=0, wspace=0)

	ax = MUST.pyconvert(Array, ax)


	marker_alpha =0.5
	marker_lw =1.5
	line_alpha=1.0
	line_lw = 3.0
	
	mi = MUST.ModelInformation(hot_big_CEMP_name)
	mi2 = MUST.ModelInformation(cold_big_CEMP_name)
	th = "$(round(Int, mi.teff)) K"
	tc = "$(round(Int, mi2.teff)) K"
	
	s = get_m3d_spectra(-1, joinpath(switch_dir, hot_big_CEMP_name))
	x_3D_hot,y_3D_hot = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ "*th, color="steelblue", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(-1, joinpath(datadir, cold_big_CEMP_name))
	x_3D_cold,y_3D_cold = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ "*tc, color="tomato", ls="-", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	x_1D_hot,y_1D_hot = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ (MARCS)\ "*th, color="steelblue", ls="--", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))
	x_1D_cold,y_1D_cold = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ (MARCS)\ "*tc, color="tomato", ls="--", lw=line_lw, alpha=line_alpha)


	ax[1].legend(labelspacing=0.1, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2)

	ax[1].set_xlim(-4, 0)
	ax[1].set_ylabel("[X/H]")
	ax[1].set_title("log(g) = $(round(mi.logg, digits=2))")
	ax[1].set_ylim(0, 6.9)

	

	# CEMP model
	s = get_m3d_spectra(-1, joinpath(switch_dir, hot_big_CEMP_name))
	x, y = plot_contr_av(s, ax=ax[2], color="steelblue", ls="-", label=L"\rm 3D\ CEMP,\  "*pretty_from_name_short(hot_big_CEMP_name))
	formHeight_3D_hot = x[argmax(y)]

	# 1D model
	s = get_m3d_spectra(marcs_hot_big, joinpath(switch_dir, hot_big_CEMP_name))
	x, y = plot_contr_1D(s, ax=ax[2], color="steelblue", ls="--", label=L"\rm 1D,\ "*pretty_from_name_short(hot_big_CEMP_name))
	formHeight_1D_hot = x[argmax(y)]

	# CEMP model
	s = get_m3d_spectra(-1, joinpath(datadir, cold_big_CEMP_name))
	x, y = plot_contr_av(s, ax=ax[2], color="tomato", ls="-", label=L"\rm 3D\ CEMP,\ "*pretty_from_name_short(cold_big_CEMP_name))
	formHeight_3D_cold = x[argmax(y)]

	# 1D model
	s = get_m3d_spectra(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))
	x, y = plot_contr_1D(s, ax=ax[2], color="tomato", ls="--", label=L"\rm 1D,\ "*pretty_from_name_short(cold_big_CEMP_name))
	formHeight_1D_cold = x[argmax(y)]


	ax[1].axvline(formHeight_3D_hot, color="steelblue", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_3D_cold, color="tomato", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_1D_hot, color="steelblue", ls="--", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_1D_cold, color="tomato", ls="--", alpha=marker_alpha, lw=marker_lw)

	closest_temp_at_form(x, y, form) = begin
		iclose = argmin(abs.(x .- form))
		y[iclose]
	end
	ax[1].axhline(
		closest_temp_at_form(x_3D_hot, y_3D_hot, formHeight_3D_hot), 
			color="steelblue", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_3D_cold, y_3D_cold, formHeight_3D_cold), 
			color="tomato", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_1D_hot, y_1D_hot, formHeight_1D_hot), 
			color="steelblue", ls="--", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_1D_cold, y_1D_cold, formHeight_1D_cold), 
			color="tomato", ls="--", alpha=marker_alpha, lw=marker_lw
		)

	

	ax[2].legend(labelspacing=0.1, loc="upper left", ncol=2)

	ax[2].set_ylabel("contribution function")
	ax[2].set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")

	#ax.set_title(pretty_from_name(hot_big_CEMP_name))


	#ax[2].set_xlim(-6, 2)
	ax[2].set_ylim(-0.1, 1.1)
	
	f
end

# ╔═╡ 213822d8-0d1c-4ea0-905c-dd0897de9ae3
let
	plt.close()

	f, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=true)
	plt.subplots_adjust(hspace=0, wspace=0)

	ax = MUST.pyconvert(Array, ax)


	marker_alpha =0.5
	marker_lw =1.5
	line_alpha=1.0
	line_lw = 3.0
	
	mi = MUST.ModelInformation(cold_small_CEMP_name)
	mi2 = MUST.ModelInformation(cold_big_CEMP_name)
	th = "$(round(Int, mi.teff)) K"
	tc = "$(round(Int, mi2.teff)) K"
	
	s = get_m3d_spectra(-1, joinpath(switch_dir, cold_small_CEMP_name))
	x_3D_hot,y_3D_hot = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ "*th, color="steelblue", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(-1, joinpath(datadir, cold_big_CEMP_name))
	x_3D_cold,y_3D_cold = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ "*tc, color="tomato", ls="-", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	x_1D_hot,y_1D_hot = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ (MARCS)\ "*th, color="steelblue", ls="--", lw=line_lw, alpha=line_alpha)
	
	s = get_m3d_spectra(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))
	x_1D_cold,y_1D_cold = plot_abund_av_3D(s, "CH", ax=ax[1], label=L"\rm CH\ -\ (MARCS)\ "*tc, color="tomato", ls="--", lw=line_lw, alpha=line_alpha)


	ax[1].legend(labelspacing=0.1, loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2)

	ax[1].set_xlim(-4, 0)
	ax[1].set_ylabel("[X/H]")
	ax[1].set_title("log(g) = $(round(mi.logg, digits=2))")
	ax[1].set_ylim(0, 6.9)

	

	# CEMP model
	s = get_m3d_spectra(-1, joinpath(switch_dir, cold_small_CEMP_name))
	x, y = plot_contr_av(s, ax=ax[2], color="steelblue", ls="-", label=L"\rm 3D\ CEMP,\  "*pretty_from_name_short(cold_small_CEMP_name))
	formHeight_3D_hot = x[argmax(y)]

	# 1D model
	s = get_m3d_spectra(marcs_cold_small, joinpath(switch_dir, cold_small_CEMP_name))
	x, y = plot_contr_1D(s, ax=ax[2], color="steelblue", ls="--", label=L"\rm 1D,\ "*pretty_from_name_short(cold_small_CEMP_name))
	formHeight_1D_hot = x[argmax(y)]

	# CEMP model
	s = get_m3d_spectra(-1, joinpath(datadir, cold_big_CEMP_name))
	x, y = plot_contr_av(s, ax=ax[2], color="tomato", ls="-", label=L"\rm 3D\ CEMP,\ "*pretty_from_name_short(cold_big_CEMP_name))
	formHeight_3D_cold = x[argmax(y)]

	# 1D model
	s = get_m3d_spectra(marcs_cold_big, joinpath(datadir, cold_big_CEMP_name))
	x, y = plot_contr_1D(s, ax=ax[2], color="tomato", ls="--", label=L"\rm 1D,\ "*pretty_from_name_short(cold_big_CEMP_name))
	formHeight_1D_cold = x[argmax(y)]


	ax[1].axvline(formHeight_3D_hot, color="steelblue", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_3D_cold, color="tomato", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_1D_hot, color="steelblue", ls="--", alpha=marker_alpha, lw=marker_lw)
	ax[1].axvline(formHeight_1D_cold, color="tomato", ls="--", alpha=marker_alpha, lw=marker_lw)

	closest_temp_at_form(x, y, form) = begin
		iclose = argmin(abs.(x .- form))
		y[iclose]
	end
	ax[1].axhline(
		closest_temp_at_form(x_3D_hot, y_3D_hot, formHeight_3D_hot), 
			color="steelblue", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_3D_cold, y_3D_cold, formHeight_3D_cold), 
			color="tomato", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_1D_hot, y_1D_hot, formHeight_1D_hot), 
			color="steelblue", ls="--", alpha=marker_alpha, lw=marker_lw
		)
	ax[1].axhline(
		closest_temp_at_form(x_1D_cold, y_1D_cold, formHeight_1D_cold), 
			color="tomato", ls="--", alpha=marker_alpha, lw=marker_lw
		)

	

	ax[2].legend(labelspacing=0.1, loc="upper left", ncol=2)

	ax[2].set_ylabel("contribution function")
	ax[2].set_xlabel(L"\rm optical\ depth\ [\tau_{500}]")

	#ax.set_title(pretty_from_name(hot_big_CEMP_name))


	#ax[2].set_xlim(-6, 2)
	ax[2].set_ylim(-0.1, 1.1)
	
	f
end

# ╔═╡ b28a34da-219e-424e-bed3-17ac707fad74


# ╔═╡ 874f7823-d626-45c9-aac8-72b22fe872c6
md"# Data comparison"

# ╔═╡ 8649e73c-ddf0-42b0-b0e1-93295bcdda5f
md"## SAGA data"

# ╔═╡ 35c27a34-6c3d-4b9f-ba87-305e81eadaff
saga_all_path = "../cgisess_32aefcdc405d0a0c9fb7936f6f1324aa_data.tsv"

# ╔═╡ 1f1a0d0f-cc24-448b-9291-91325a9d72eb
saga_co_path = "../cgisess_co.tsv"

# ╔═╡ 0bbf88db-b471-4c9c-999f-c3ace82841be
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


# ╔═╡ 067154e6-09fa-46ba-b7dc-7e8a6757a932
let
	teff = parameters_all["teff"][selection_mask]
	logg = parameters_all["logg"][selection_mask]

	mod_teff = [a[1] for a in abundance_correction_parameters]
	mod_logg = [a[2] for a in abundance_correction_parameters]

	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 5))

	ax.scatter(teff, logg, s=15, rasterized=true, color="k", alpha=0.2, label=L"\rm SAGA")
	ax.scatter(mod_teff, mod_logg, color="red", marker="X", s=200, label=L"\rm 3D\ RHD\ models")

	ax.set_ylim(ax.get_ylim()[1], ax.get_ylim()[0])
	ax.set_xlim(ax.get_xlim()[1], ax.get_xlim()[0])

	ax.set_ylabel(L"\rm log(g)")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	#ax.legend(labelspacing=0.1, loc="lower center", ncol=2, columnspacing=0.5)

	#ax.set_ylim(5.7, 2.3)
	#ax.set_xlim(8300, 3800)

	f
end

# ╔═╡ 36cf99a4-d9cb-4282-a45c-50b7ff5833d3


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

# ╔═╡ f6d84057-c57e-4f50-90bf-f5cba97aca2a
md"Predictions from Hartwig et al. 2018, extracted from their Figure:"

# ╔═╡ dd9203c7-7c19-4e49-9458-bfa11fa089b6
begin
	hartwig18_fiducial = MUST.readdlm("../hartwig18_fiducial.csv", ',')
	hartwig18_faint20 = MUST.readdlm("../hartwig18_faint20.csv", ',')
end;

# ╔═╡ 570da95e-276c-4173-a5d1-77106cd1e8ff


# ╔═╡ 216a0e5c-12ee-4e04-9f50-cb5bd4fa510e
#metallicity_bin_edges = [-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]
metallicity_bin_edges = [-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]

# ╔═╡ 27b5837f-b21a-434b-a074-a18e9ef7f3d2
metallicity_bin_centers = (metallicity_bin_edges[2:end] .+ metallicity_bin_edges[1:end-1]) ./ 2

# ╔═╡ 952d802f-82dd-47e9-a1f9-fa98fc505d28
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

# ╔═╡ 35534a50-475a-40b9-a6f3-a352da70ab7f
get_n_below(lower_than_feh, feh, mask) = begin
	below_limit = (feh .<= lower_than_feh) .& selection_mask .& mask
	count(below_limit)
end

# ╔═╡ 7c9d2201-4078-4f42-9466-57f7bc8f4e1f
get_n_below_frac(lower_than_feh, feh, mask) = begin
	below_limit = (feh .<= lower_than_feh) .& selection_mask
	below_limit_mask = (feh .<= lower_than_feh) .& selection_mask .& mask
	count(below_limit_mask) / count(below_limit)
end

# ╔═╡ 48bf4d88-82e0-4aaa-b1e8-35d375a03910
begin
	@info "cumsum of those bins:"
	@info cumsum(count_bins_general)
	@info cumsum(count_bins_CEMP)
	@info cumsum(count_bins_CEMP_corr)
end

# ╔═╡ a625e31c-a113-45b5-b7c8-a21e31df18b1
@info "Stars with log(g) higher than limit: $(count(selection_mask .& (parameters_all["feh"] .< 0.0)))"

# ╔═╡ fd8d8aaa-7cf7-4b64-a823-cbd2f4781f05
@info "Stars with log(g) higher than limit + [C/Fe]: $(count(selection_mask .& (.!isnan.(parameters_all["cfe"])) .& (parameters_all["feh"] .< 0.0)))"

# ╔═╡ 6a49c831-5b6b-4c7c-9146-db03db278004
@info "CEMP stars: $(count(selection_mask .& (.!isnan.(parameters_all["cfe"])) .& (parameters_all["feh"] .< 0.0) .& (parameters_all["cfe"] .> 0.7)))"

# ╔═╡ f07880a2-bbaa-47c9-897a-fd94d93ede4f


# ╔═╡ 6e1e3bf9-7299-4ec3-8802-bbe798089fcb
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

# ╔═╡ 17521b74-d37b-4c84-87ad-0e4f155b3790
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

# ╔═╡ 892a04e1-3758-435a-b8cd-e7c03b777923
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

# ╔═╡ 222ffed0-a355-4e08-aa6a-fa305c078d2c
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

# ╔═╡ 52dff81a-701b-461d-aa73-791239b135af
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

# ╔═╡ 0ee1f4fc-da57-4f63-9713-c20dc2f0464d


# ╔═╡ 9c789f14-0a75-45f7-912f-9a5190a3a54a
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
	
	    # Draw major axis line
		if !skip_horizontal
	    	ax.plot([x0 - N*dx_major, x0 + N*dx_major], [y0 - N*dy_major, y0 + N*dy_major]; kwargs...)
		end

		if !skip_vertical
		    # Draw minor axis line
		    ax.plot([x0 - N*dx_minor, x0 + N*dx_minor], [y0 - N*dy_minor, y0 + N*dy_minor]; kwargs...)
		end
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

# ╔═╡ f4b9583a-9dfd-477b-beab-15ec8a5b82dd
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

# ╔═╡ e7bad683-5550-4a48-82ef-620455f9f77e


# ╔═╡ 69dc31c5-4b08-442a-ad1f-5170efa9efad
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

# ╔═╡ 1a4967b8-b19c-4e23-b6bb-6a09a30cea8e


# ╔═╡ 30bef587-80f6-4d36-a2b7-283b44f1db13
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

# ╔═╡ c08c0fd5-808f-46cd-8617-c2db2d5dc181


# ╔═╡ e897badd-4296-454d-807b-af34bfa35c9d
md"## Galactic C/O ratio"

# ╔═╡ 39de6e5e-8c11-46db-b047-c0a389596c23
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
# ╟─dc377839-bd92-4c9c-9432-7bf08b71add6
# ╠═3867ab6e-bcaf-4879-9eb9-f5b42dcc6709
# ╟─9818e49b-1a6d-44d3-9aa2-6486842b2efa
# ╟─1f6f5fe0-3f38-4547-b9e6-0482a76ea03a
# ╟─1150c5dd-f486-46d6-b911-5fb1d70d3792
# ╟─4ea802dc-b3e3-44ab-8758-b76bef9de441
# ╠═f4a66407-d7c8-43f9-b68f-1f7f29b4eba6
# ╟─06f0b894-d945-4941-aee5-e3e699e162f0
# ╟─ea64c5f8-21e3-4685-ba30-2fd61ab34974
# ╠═49ce048c-18ec-44d6-a782-f474ad11c83c
# ╠═d30fd639-00c4-433e-b5c4-cdc59e38281e
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
# ╟─bfc09acd-8fe3-44c9-a44f-8aefe85f3d76
# ╟─5dc20e44-24eb-4e8f-bda1-fac7a1c0a5f4
# ╟─4802eaa5-1c91-4813-b335-c217c5ceac87
# ╟─e624d1fd-24e2-4277-8838-a9ce0026ec56
# ╟─0a3ce98b-40b5-4948-87c6-57d4cf9bbba9
# ╟─c0ef68e0-5b8a-4d6d-a529-8d6292799420
# ╟─8749760d-4857-41ed-8a0a-f7464d5fbc67
# ╟─6df206d1-6794-4bab-a10d-914c6f508137
# ╟─e735e811-975c-49ff-b99f-7b98e63e5150
# ╟─6a2c237d-a0aa-4db7-9a87-b259eda25546
# ╟─9ccc4261-0b53-4083-b8a3-850f10994048
# ╟─e6f379f0-6cbf-4d74-9c5e-4fb567307eea
# ╟─c8ad7b76-d063-457b-809e-65b818c796f4
# ╟─99e19910-17a3-4dd2-8ec0-72dd3da2726a
# ╟─01ec615e-bc60-4e50-a076-692317095591
# ╟─55b53af5-c8d3-4086-ac92-5995002d8c37
# ╟─3b07d6ba-0e3f-484c-b3f8-da3840f5d020
# ╟─56a77efd-66d5-47aa-934b-99e3d1aff0a8
# ╟─989f922a-9596-4045-b4fc-186a18ca5725
# ╟─8cebfa74-030a-4e3a-9c91-49160c1f8f64
# ╟─de37b6dc-cf8b-4a1c-95b8-ac59ac04837a
# ╟─687bd929-41a2-416e-ad7a-681e98096c98
# ╟─212b2d0b-fbcc-4934-8ebb-627545fea549
# ╟─bb9ebb4c-4065-4b5c-82a6-05d88be2c55e
# ╟─45cd645d-055d-41a7-a7e9-0f0867d0d183
# ╟─69bf13d1-3f9a-45c6-98c0-523951ed812f
# ╟─d51a843b-168e-411f-8641-e02d96371678
# ╟─44694a88-c977-4d56-9fba-a42da6f60682
# ╟─f86f0d5b-ffc1-429a-b31c-8aad4aeb396a
# ╟─685e9c27-dc91-4f82-8112-998383271de0
# ╟─c7fb5cc8-5a09-47d7-bb4a-579fb5a085f9
# ╟─136db4eb-c274-46ef-bac2-03a3bf300cac
# ╟─560fcc8a-a9ee-4420-bc97-eb53b9f6de94
# ╟─3fef366c-b34b-4747-ac5c-7e3ef5377ef8
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
# ╠═4f2b6d24-9157-47ea-a11b-53a30dc34c7a
# ╟─4ff6bac5-7f46-4f4b-a221-8b75bdb6b53f
# ╟─258d76cb-ac6b-453e-b541-29ad51461195
# ╟─f77ae3e8-374a-4e3c-b8b3-7985747bccc4
# ╟─a58d4a2d-45ad-4bb3-9000-6a05f43e5d71
# ╟─53fef078-5301-4772-981c-68fbb1349f8d
# ╟─7a831dfa-f216-4d6b-9370-5f286d7f3e77
# ╟─631410d5-7e7a-458f-a7fb-74b462913bd8
# ╟─1149d154-42e0-4e94-b26f-716504872afa
# ╟─8d19d381-a19e-4073-a19a-39bd2ff1be7c
# ╟─ffd0e8b7-7c9a-4e51-9a1e-1982dc59754c
# ╟─87aaadf1-2d46-4d14-8723-8c843040d1d1
# ╟─f3dc3826-c36e-4292-9431-a86d4325903c
# ╟─093b3783-af82-4e34-a1e6-65a0df5bb35c
# ╟─c85b36d4-ff22-4ba6-a4d1-c1791fd94ee2
# ╟─3a59080e-f925-48c7-b01d-660c45a7f541
# ╟─d8f6a261-88ab-4222-8a51-ce1b359df452
# ╟─2dd8fd87-f543-413e-b8f3-4187ff67e52c
# ╟─75f3327c-00c4-4594-86ac-e2b5237dc0bc
# ╟─1cf58f8e-e1d0-48eb-9c00-c22fc3144688
# ╟─6c9b537b-60a8-4888-ba5d-be97621e03e8
# ╟─537080cd-02e8-4daa-86af-2eee60d1a1ca
# ╟─87366f88-2d35-4f90-bb71-0e7c9bcb0811
# ╟─9a3df65f-d469-42d9-b5c4-7cf66f1ed5f6
# ╟─46d9bb97-1fd7-4d57-bcaa-8123754f195c
# ╟─b0e95f8f-0059-4575-bab6-80f086ba3e4a
# ╟─ef9fef62-5e50-4a1d-b1a9-4db7289663a7
# ╟─331f1016-df68-4330-8983-a2db6219f966
# ╟─924ea3f0-c822-44dc-9faa-1629c14f5232
# ╟─494823c4-50b7-4e23-894c-caafaab48b29
# ╟─1bd82fae-afb3-4197-999e-61fba6ccfe71
# ╟─731fadcd-c03a-4f7d-8f49-f7c2c29c099b
# ╟─981e31b0-7073-431c-b1a4-d82fd6e57350
# ╟─7e487942-9681-4441-8969-050c8328a19a
# ╟─08fba04e-f041-4c49-865e-ec254a99fc98
# ╟─0e7e6f79-7225-4d29-8693-b2b93ddf4812
# ╟─217e0937-4b2d-4948-9a72-c236a8fc5c5a
# ╟─b942b4f4-6979-42ee-8f18-b9e7b64c0379
# ╟─a84ea668-6701-4e31-9f04-cb0feb70290d
# ╟─7dae5609-0148-438f-a207-98da5c9173d4
# ╟─e4232daa-7174-493d-968e-21a012204fcb
# ╟─e502770d-d765-4757-a075-3fadcdd88eed
# ╟─f64c0cd9-c0ec-44e3-9e0c-440b6fdc23e8
# ╟─635d415f-aa91-4d25-8162-a550119e3f12
# ╟─2caf0768-3202-4d16-8c24-4f9a6433076e
# ╠═ac81d133-8894-46da-b85c-ab5d580ecd12
# ╠═01752527-37fd-4175-9982-02877ca4ae69
# ╟─2d18e394-484d-4aef-bcb5-c45627fcf824
# ╟─4b18b95c-e574-468a-92ce-73931621e442
# ╟─af94af20-6618-448e-b146-b3ba1b1a20f1
# ╟─22932761-ae92-4cd6-b828-07e3f86a4d0f
# ╟─a6ee120f-0dfc-433a-8c9c-0f9f828ffe29
# ╟─421da8e3-f710-46f1-a602-bdc4a3266a26
# ╟─213822d8-0d1c-4ea0-905c-dd0897de9ae3
# ╟─b28a34da-219e-424e-bed3-17ac707fad74
# ╟─874f7823-d626-45c9-aac8-72b22fe872c6
# ╟─8649e73c-ddf0-42b0-b0e1-93295bcdda5f
# ╠═35c27a34-6c3d-4b9f-ba87-305e81eadaff
# ╠═1f1a0d0f-cc24-448b-9291-91325a9d72eb
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
# ╟─067154e6-09fa-46ba-b7dc-7e8a6757a932
# ╟─36cf99a4-d9cb-4282-a45c-50b7ff5833d3
# ╟─990e7d73-3e8e-4900-a225-4e5cac6ac8c6
# ╟─5e31351e-5e77-46a5-9e8c-baa53eeaee73
# ╟─73b3581d-a484-4ff4-8e69-445192946a3b
# ╟─2fc5e72e-4c17-4519-915d-88479d11228f
# ╟─c012c57f-ea9c-44fa-ab8d-e2794953b214
# ╟─0bb30e77-41a2-4e73-a032-cf6b408f3d7c
# ╟─f6d84057-c57e-4f50-90bf-f5cba97aca2a
# ╠═dd9203c7-7c19-4e49-9458-bfa11fa089b6
# ╟─570da95e-276c-4173-a5d1-77106cd1e8ff
# ╠═216a0e5c-12ee-4e04-9f50-cb5bd4fa510e
# ╟─27b5837f-b21a-434b-a074-a18e9ef7f3d2
# ╟─952d802f-82dd-47e9-a1f9-fa98fc505d28
# ╟─e465e7b6-91a7-4826-859b-5b60382bc96f
# ╠═35534a50-475a-40b9-a6f3-a352da70ab7f
# ╠═7c9d2201-4078-4f42-9466-57f7bc8f4e1f
# ╟─48bf4d88-82e0-4aaa-b1e8-35d375a03910
# ╟─a625e31c-a113-45b5-b7c8-a21e31df18b1
# ╟─fd8d8aaa-7cf7-4b64-a823-cbd2f4781f05
# ╟─6a49c831-5b6b-4c7c-9146-db03db278004
# ╟─f07880a2-bbaa-47c9-897a-fd94d93ede4f
# ╟─6e1e3bf9-7299-4ec3-8802-bbe798089fcb
# ╟─17521b74-d37b-4c84-87ad-0e4f155b3790
# ╟─892a04e1-3758-435a-b8cd-e7c03b777923
# ╟─222ffed0-a355-4e08-aa6a-fa305c078d2c
# ╟─52dff81a-701b-461d-aa73-791239b135af
# ╟─0ee1f4fc-da57-4f63-9713-c20dc2f0464d
# ╟─9c789f14-0a75-45f7-912f-9a5190a3a54a
# ╟─76e3d925-6778-4d40-a02f-23e6aa6258eb
# ╟─f4b9583a-9dfd-477b-beab-15ec8a5b82dd
# ╟─e7bad683-5550-4a48-82ef-620455f9f77e
# ╠═69dc31c5-4b08-442a-ad1f-5170efa9efad
# ╟─1a4967b8-b19c-4e23-b6bb-6a09a30cea8e
# ╟─30bef587-80f6-4d36-a2b7-283b44f1db13
# ╟─c08c0fd5-808f-46cd-8617-c2db2d5dc181
# ╟─e897badd-4296-454d-807b-af34bfa35c9d
# ╟─39de6e5e-8c11-46db-b047-c0a389596c23
