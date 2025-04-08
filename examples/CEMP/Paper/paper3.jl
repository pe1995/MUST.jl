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
		L"\rm T_{eff}="*"$(round(Int, mi.teff))"*L"\ K,\ log(g)="*"$(round( mi.logg, digits=2))"*L",\ [Fe/H] ="*"$(round(mi.feh, digits=2))"
	end
	
	TableOfContents()
end

# ╔═╡ 4ea802dc-b3e3-44ab-8758-b76bef9de441


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

# ╔═╡ 4802eaa5-1c91-4813-b335-c217c5ceac87


# ╔═╡ e624d1fd-24e2-4277-8838-a9ce0026ec56
md"# Structure comparison"

# ╔═╡ 0a3ce98b-40b5-4948-87c6-57d4cf9bbba9
"""
	plot_profile(s, q; ax, kwargs...)

Plot profile function of optical depth cube (500).
"""
plot_profile(b, bt, q; ax=nothing,cmap="GnBu", alpha=1.0, bins=200, kwargs...) = begin
	qs, logq = MUST.is_log(q)
	
	# Background distribution
	h, x, y = MUST.numpy.histogram2d(
		reshape(log10.(b[:τ500]), :), reshape(logq.(b[qs]), :), bins=bins
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
	structure_3D1D_figsave_txt = "$(structure_3D_select)_structure_3D1D.pdf"
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

	f.savefig(structure_3D1D_figsave_txt)
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

# ╔═╡ fed53f26-3d9f-4be3-957e-a7181af693d9
same_scaled_solar(cemp_model_name) = scaledSolar_models[
	findfirst(
		MUST.same_parameters(
			MUST.ModelInformation(a), 
			MUST.ModelInformation(cemp_model_name),
		) for a in scaledSolar_models
	)
]

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
					marcsmodel, "GBand_DC+0",
					model3D, [
						"GBand_DC+0",
						"GBand_DC+02",
						"GBand_DC+04",
						"GBand_DC+06",
						"GBand_DC+08",
						"GBand_DC-02",
						"GBand_DC-04",
						"GBand_DC-06",
						"GBand_DC-08"
					],
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
		marker="X", color="tomato", 
		label="1D", markersize=12
	)

	ax.set_ylabel(L"\rm EW\ [\AA]")
	ax.set_xlabel(L"\rm [C/Fe]")

	ax.set_title(L"\rm [Fe/H] ="*"$(cog_z_select)")
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

# ╔═╡ Cell order:
# ╟─dc377839-bd92-4c9c-9432-7bf08b71add6
# ╠═3867ab6e-bcaf-4879-9eb9-f5b42dcc6709
# ╟─9818e49b-1a6d-44d3-9aa2-6486842b2efa
# ╟─1f6f5fe0-3f38-4547-b9e6-0482a76ea03a
# ╟─1150c5dd-f486-46d6-b911-5fb1d70d3792
# ╟─4ea802dc-b3e3-44ab-8758-b76bef9de441
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
# ╟─4802eaa5-1c91-4813-b335-c217c5ceac87
# ╟─e624d1fd-24e2-4277-8838-a9ce0026ec56
# ╟─0a3ce98b-40b5-4948-87c6-57d4cf9bbba9
# ╟─c0ef68e0-5b8a-4d6d-a529-8d6292799420
# ╟─8749760d-4857-41ed-8a0a-f7464d5fbc67
# ╟─6df206d1-6794-4bab-a10d-914c6f508137
# ╟─e735e811-975c-49ff-b99f-7b98e63e5150
# ╟─6a2c237d-a0aa-4db7-9a87-b259eda25546
# ╟─c8ad7b76-d063-457b-809e-65b818c796f4
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
# ╟─fed53f26-3d9f-4be3-957e-a7181af693d9
# ╟─b874cdbb-7496-471d-8e2e-4e3f47a126b0
# ╠═4ff6bac5-7f46-4f4b-a221-8b75bdb6b53f
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
