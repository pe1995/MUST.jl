### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4a837d4c-724c-11ee-0b6b-21a88a56df5b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST 
	using Glob
	using PythonPlot
	using TSO 
	using LaTeXStrings
	using DelimitedFiles
	using PlutoUI
	using PlutoUI: combine
end;

# ╔═╡ bd7cdaf2-c05e-4bc9-afbb-cfb61260d62c
md"# Analysing Grid Models"

# ╔═╡ 827e54c7-d0a3-455a-8837-52255eeff202
md"## Setup"

# ╔═╡ d9b5bbbd-6229-47c3-9738-dbe9087016d8
TableOfContents()

# ╔═╡ 9e9c5903-762d-43d7-adf2-0eccee134974
begin
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	MUST.@import_dispatch "../../../dispatch2"
end;

# ╔═╡ 9ac91858-17d8-4934-8d47-fed1519efe42
availableRuns(path) = begin
	runs = glob("*/", @in_dispatch(path))
	last.(split.(runs, "/", keepempty=false))
end

# ╔═╡ 80fefa1e-de68-4448-a371-d4ae0ec70868
md"## Picking Models
You can pick a model from the output folder, that was already converted using 'convert.jl' or 'snapshot2box.jl'"

# ╔═╡ de91f663-285b-4250-9edb-41f48d5b36e1
md"Enter the data folder (relative to dispatch experiment):"

# ╔═╡ 6fea8f52-3ca6-4854-aca7-6b125d4c6542
@bind datafolder TextField(default="data/")

# ╔═╡ ac1069a6-0687-44fb-bd94-c76af8fd6b1f
md"Select one of the available models:"

# ╔═╡ 73e9bde0-c72e-4a01-aa68-486237e0c475
@bind available_runs MultiCheckBox(availableRuns(datafolder))

# ╔═╡ 73d9fd32-0dae-478a-80ff-a8e314adbd6e
function snapshots_input(simulation_names::Vector)
	snapshots_runs = Any[[] for _ in available_runs]
	names = []
	for (i, run) in enumerate(available_runs)
		o = @in_dispatch(joinpath(datafolder, run))
		snaps = MUST.list_snapshots(converted_snapshots(o))
		append!(snapshots_runs[i], snaps)
		append!(names, ["$(run)"])
	end
	
	return combine() do Child
		namesChild = [
			Child(names[i], MultiCheckBox(snapshots_runs[i]))
			for i in eachindex(names)
		]
		inputs = [
			md""" $(name): $(
				namesChild[i]
			)"""
			for (i, name) in enumerate(names)
		]
		
		md"""
		### Converted Models
		$(inputs)
		"""
	end
end

# ╔═╡ 16786647-4021-4068-9af1-2a548240758e
function snapshots_all_input(simulation_names::Vector)
	snapshots_runs = Any[[] for _ in available_runs]
	names = []
	for (i, run) in enumerate(available_runs)
		o = @in_dispatch(joinpath(datafolder, run))
		snaps = MUST.list_of_snapshots(o)
		snaps_already = MUST.list_snapshots(converted_snapshots(o))

		if !(snaps in snaps_already)
			append!(snapshots_runs[i], snaps)
			append!(names, ["$(run)"])
		end
	end
	
	return combine() do Child
		namesChild = [
			Child(names[i], Slider(snapshots_runs[i], show_value=true))
			for i in eachindex(names)
		]
		inputs = [
			md""" $(name): $(
				namesChild[i]
			)"""
			for (i, name) in enumerate(names)
		]
		
		md"""
		### Not converted snapshots
		You can convert snapshots directly in this notebook by choosing from the available snapshots and then tick the convert button.
		$(inputs)
		"""
	end
end

# ╔═╡ 68510276-b102-4a70-bd2f-a63195a4bb08


# ╔═╡ 2bf867fc-15bd-437c-b989-b1efbf2bf9d4
@bind snapshots_convert snapshots_all_input(available_runs)

# ╔═╡ fa8c564d-195b-424a-8814-0c47cd3e03a8
md"Tick box to start converting snapshots: $(@bind convert_given CheckBox(default=false))"

# ╔═╡ 2d684eb0-aa4b-40a6-9e04-ab0b42ed0516
if convert_given
	for name in keys(snapshots_convert)
		@info "Converting snapshot $(name), Number: $(snapshots_convert[name])"
		p = @in_dispatch(joinpath(datafolder, "$(name)"))
		snapshotBox(snapshots_convert[name], folder=p)
	end
end

# ╔═╡ 53706541-27d7-4093-b37f-3e384c6c8a7c


# ╔═╡ ddcbcf4d-eb85-4cde-b2ac-a4b7e93ed809
begin
	convert_given
	@bind snapshots_picks snapshots_input(available_runs)
end

# ╔═╡ 40124987-2dc6-452d-87d7-33f935854b55


# ╔═╡ 3778a669-69fc-47e6-9d0b-c514b31c8c79
md"### Loading the models"

# ╔═╡ b9f82f31-f25e-4cad-abb4-aac2d4d0e584
begin
	snapshots = Dict(name=>[] for name in keys(snapshots_picks))
	snapshots_τ = Dict(name=>[] for name in keys(snapshots_picks))
	
	for name in keys(snapshots_picks)
		for v in snapshots_picks[name]
			p = @in_dispatch(joinpath(datafolder, "$(name)"))
			snapshot, snapshot_τ = pick_snapshot(
				converted_snapshots(p), v
			)
			append!(snapshots[name], [snapshot])
			append!(snapshots_τ[name], [snapshot_τ])
		end
	end
end

# ╔═╡ cfa74be4-06c9-4911-ac29-2a17809391e6
labels = NamedTuple(
	name=>["$(name)_$(v)" for v in snapshots_picks[name]] 
	for name in keys(snapshots_picks)
)

# ╔═╡ ba8bded9-f13b-477e-a9fc-bd88e4f01e35
loggs = NamedTuple(
	name=>length(snapshots[name]) > 0 ? first(snapshots[name]).parameter.logg : -99.0
	for name in keys(snapshots_picks)
)

# ╔═╡ ac46894f-44b4-477a-81ce-687b1308c38f
md"## Picking EoS"

# ╔═╡ 574b4a88-c562-44ab-8a8e-6e5b5c5aabf0
md"Enter the folder where EoS tables can be found:"

# ╔═╡ 926aa121-ba10-4410-a91c-f36fc79a63ed
@bind eosfolder TextField(default="input_data/grd/")

# ╔═╡ 91ddaf2d-abd6-4e2d-92eb-cfa7d0ba79bc
function eos_input(simulation_names::Vector)
	all_eos = glob("*/", @in_dispatch(eosfolder))
	all_eos = convert.(String, last.(split.(all_eos, "/", keepempty=false)))
	simulation_names = convert.(String, simulation_names)

	return combine() do Child
		namesChild = [
			Child(simulation_names[i], Select(all_eos))
			for i in eachindex(simulation_names)
		]
		inputs = [
			md""" $(name): $(
				namesChild[i]
			)"""
			for (i, name) in enumerate(simulation_names)
		]
		
		md"""
		### Available EoS
		$(inputs)
		"""
	end
end

# ╔═╡ 3a94ff68-ea2f-4478-902f-4975074e8ec4
@bind eos_picks eos_input(available_runs)

# ╔═╡ ca860715-104c-4177-8249-a547863f6db4


# ╔═╡ c5f4fb05-27b6-42e8-b08a-ef60c403d9e7
eos_folders = NamedTuple(
	name=>MUST.@in_dispatch(joinpath("$(eosfolder)", "$(eos_picks[name])"))
		for name in keys(snapshots_picks)
)

# ╔═╡ 4d0c7aa4-2a13-429a-88f2-14df415aa5ec
eos = NamedTuple(
	name=>reload(SqEoS, joinpath(eos_folders[name], "eos.hdf5"))
	for name in keys(eos_folders)
)

# ╔═╡ 27d0d01a-c8c8-422b-a893-36948b908fb0
opa = NamedTuple(
	name=>reload(SqOpacity, joinpath(eos_folders[name], "binned_opacities.hdf5"))
	for name in keys(eos_picks)
)

# ╔═╡ 66f7d726-b840-484c-937f-b6a328e86b21


# ╔═╡ 1b354d9a-82f6-4540-8eaa-69218e37f87b
md"## Initial models"

# ╔═╡ f6c017c2-b296-4f82-8d2b-d97a7e36a51e
md"Pick what initial models to include in the figures"

# ╔═╡ 693c5461-e4d8-45d0-af2d-87b10f9560e3
@bind pick_initial_models MultiCheckBox(keys(snapshots_picks)|>collect)

# ╔═╡ 2bb999c6-ed87-4006-9fba-2f01fcf22f90


# ╔═╡ d55bae1d-1ca1-4537-8af1-c025a966c3b3
initial_model = NamedTuple(
	name=>@optical(
		Average3D(eos[name], joinpath(eos_folders[name], "inim.dat"), logg=loggs[name]), eos[name],
		opa[name]
	) for name in pick_initial_models
)

# ╔═╡ 271e4e9c-def2-4b3e-b383-45d95fbb309b
labels_initial = NamedTuple(
	name=>"initial: $(name)"
	for name in pick_initial_models
)

# ╔═╡ c561e5a2-dfd6-4de9-80b5-9f6ae100b217
initial_adiabats = NamedTuple(
	name=>@optical(
		TSO.adiabat(initial_model[name], eos[name]), eos[name], opa[name]
	) for name in pick_initial_models
)

# ╔═╡ 55581a87-2638-4f85-9ad7-7313c37c1971
labels_adiabat = NamedTuple(
	name=>"adiabat: $(name)"
	for name in pick_initial_models
)

# ╔═╡ e6d33b23-e77c-43d3-b401-c5620721d659
md"## Stagger models"

# ╔═╡ bd41fb98-70c4-4192-8753-6fee82f14576
function allStaggerModels(folder)
	allmodels = glob("*.mesh", folder)
	if length(allmodels) == 0
		[]
	else
		allmodels = [m for m in allmodels if !occursin("_orig", m)]
		allnames = first.(split.(basename.(allmodels), ".mesh"))
	end
end

# ╔═╡ b90c53aa-9ca6-421c-9b2f-f650a1f0832f


# ╔═╡ db7db6d6-25c3-4f43-a49b-1ef5881606ac
md"Load a Stagger model in the Multi format from folder: $(@bind staggermodelpath TextField(default=\"Stagger_grid_multi/\"))"

# ╔═╡ 80f8c89e-5c14-4df6-8a05-f5fe79a7585e
if length(allStaggerModels(staggermodelpath)) > 0
	md"Pick a Stagger model: $(@bind selectedStagger Select(allStaggerModels(staggermodelpath)))"
end

# ╔═╡ 50bc1f52-0638-427a-8cc0-f8cd49535114


# ╔═╡ 3d1d9e3e-aa7d-4ef0-9c3e-41816285a938
begin
	model_stagger = try
		MUST.multiBox(joinpath(staggermodelpath, selectedStagger))
	catch
		@warn "No stagger model could be found"
		nothing
	end
end

# ╔═╡ bb0cb1ab-7e03-4a5e-bfd8-60bf17e007d5
md"# Figure layout"

# ╔═╡ b7c0a489-b7f7-4105-bce3-ae0c1074c650
md"### Colors"

# ╔═╡ 9ea4e08e-62aa-44c3-bbf0-78b1954a9c4f
cmap = plt.get_cmap("rainbow")

# ╔═╡ c4015f47-d670-4c5d-885c-cf32b6f15829
cmap_initial = plt.get_cmap("tab20_r")

# ╔═╡ 13a647a4-3d71-4afb-ba31-867e108a8154
nmodels = length(snapshots_picks) == 0 ? 0 : sum([length(snapshots_picks[name])
	for name in keys(snapshots_picks)
])

# ╔═╡ 17ca4c8c-9563-4be4-8ab7-a7446ba28dec
begin
	colors = Dict(name=>[] for name in keys(snapshots_picks))
	colors_initial = Dict()
	ic = 1
	for name in keys(snapshots_picks)
		for j in snapshots_picks[name]
			append!(colors[name], [cmap(ic/nmodels)])
			colors_initial[name] = cmap_initial(ic/nmodels)
			ic += 1
		end
	end

	colors
end

# ╔═╡ 1be872a8-eef5-422d-b3a6-13d581a999c0


# ╔═╡ f3ffa698-bcea-4ab3-9b63-84d518c14068
md"# Average figures"

# ╔═╡ ce691c30-5025-4ffc-8185-eced3087ca13
md"Tick box to start plotting profiles: $(@bind plot_given CheckBox(default=false))"

# ╔═╡ 5e882897-d396-470c-ad46-37a39326225c
md"## Height vs. Temperature"

# ╔═╡ ef9b0e78-00f0-41d0-8429-2d36f54c0a02
plot_given && let
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for name in keys(snapshots)
		#=x = reshape(snapshots[i].z, :)
		y = reshape(snapshots[i][:T], :) 
		axA.hist2d(
			x, y, bins=200, cmin=1, cmap="Blues",
			norm=matplotlib.colors.LogNorm()
		)=#
		for (j, snap) in enumerate(snapshots[name])
			axA.plot(
				profile(mean, snap, :z, :T)..., 
				lw=2., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	for name in keys(initial_model)
		axA.plot(
			-initial_model[name].z, exp.(initial_model[name].lnT), 
			lw=1., 
			color=colors_initial[name],
			label=labels_initial[name],
			alpha=0.7,
			ls="-"
		)

		#=axA.plot(
			initial_adiabats[name].z, exp.(initial_adiabats[name].lnT), 
			lw=1., 
			color=colors_initial[name],
			label=labels_adiabat[name],
			alpha=0.7,
			ls="--"
		)=#
	end
	
	axA.legend()
	axA.set_ylabel("temperature [K]")
	axA.set_xlabel("z [cm]")

	#axA.set_xlim(-1e8, 2e8)
	#axA.set_ylim(2500, 11500)
	
	
	gcf()
end

# ╔═╡ e05efb00-233b-44b9-9204-156ae2ed0762
md"## Height vs. Density"

# ╔═╡ 36eabb43-1660-4e11-8bca-e7f09568d695
plot_given && let
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for name in keys(snapshots)
		#=x = reshape(snapshots[i].z, :)
		y = reshape(snapshots[i][:T], :) 
		axA.hist2d(
			x, y, bins=200, cmin=1, cmap="Blues",
			norm=matplotlib.colors.LogNorm()
		)=#
		for (j, snap) in enumerate(snapshots[name])
			axA.plot(
				profile(mean, snap, :z, :log10d)..., 
				lw=2., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	for name in keys(initial_model)
		axA.plot(
			-initial_model[name].z, log10.(exp.(initial_model[name].lnρ)), 
			lw=1., 
			color=colors_initial[name],
			label=labels_initial[name],
			alpha=0.7,
			ls="-"
		)
	end
	
	axA.legend()
	axA.set_ylabel(L"\rm density\ [g\ \times cm^{-3} ]")
	axA.set_xlabel("z [cm]")

	#axA.set_xlim(-1e8, 2e8)
	#axA.set_ylim(2500, 11500)
	
	gcf()
end

# ╔═╡ b34da749-03ae-4f6c-92fb-c5809abeb339
md"## Density vs. Temperature"

# ╔═╡ e99b4bed-f093-4fe6-ad1d-af0f2235470b
plot_given && let
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for name in keys(snapshots)
		for (j, snap) in enumerate(snapshots[name])
			_, T = profile(mean, snap, :z, :T)
			_, d = profile(mean, snap, :z, :log10d)
		
			axA.plot(
				d, T,
				lw=2., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	for name in keys(initial_model)
		axA.plot(
			log10.(exp.(initial_model[name].lnρ)), 
			exp.(initial_model[name].lnT),
			lw=1., 
			color=colors_initial[name],
			label=labels_initial[name],
			alpha=0.7,
			ls="-"
		)
	end
	
	axA.legend()
	axA.set_xlabel(L"\rm density\ [g\ \times\ cm^{-3}]")
	axA.set_ylabel(L"\rm temperature\ [K]")	
	
	gcf()
end

# ╔═╡ ddcf3053-c8ff-4765-9a99-2c85534485fc
md"## Height vs. velocity"

# ╔═╡ d93e8b41-6f7a-4163-82dc-80a844e57e9d
rms5(x) = sqrt.(mean(x .^2)) ./1e5

# ╔═╡ 2527e21b-e956-4d9c-956d-22db9efbaadc
plot_given && let
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for name in keys(snapshots)
		for (j, snap) in enumerate(snapshots[name])
			z, u = profile(rms5, snap, :z, :uz)
		
			axA.plot(
				z, u,
				lw=2., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	if !isnothing(model_stagger)
		z, u = profile(rms5, model_stagger, :z, :uz)
		axA.plot(
			z, u,
			lw=2., 
			color="black", 
			label="Stagger"
		)
	end

	axA.legend()
	axA.set_xlabel(L"\rm z\ [cm]")
	axA.set_ylabel(L"\rm rms(U_z)\ [km\ \times\ s^{-1}]")	
	
	gcf()
end

# ╔═╡ 14796f14-183e-496f-b5d8-41149b31d463
md"## optical depth vs. Temperature"

# ╔═╡ 6c4f5ac6-a03b-4237-aa8f-fe50f77bde6f
plot_given && let
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for name in keys(snapshots_τ)
		for (j, snap) in enumerate(snapshots_τ[name])
			axA.plot(
				profile(mean, snap, :log10τ_ross, :T)..., 
				lw=2., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	for name in keys(initial_model)
		axA.plot(
			log10.(initial_model[name].τ), exp.(initial_model[name].lnT), 
			lw=1., 
			color=colors_initial[name],
			alpha=0.7,
			ls="--"
		)
	end

	axA.legend()
	axA.set_ylabel("temperature [K]")
	axA.set_xlabel(L"\rm \log \tau_{ross}")
	
	
	gcf()
end

# ╔═╡ 4d623d6f-7366-468b-930e-71878b72aa5f
md"## RMS velocity"

# ╔═╡ 3bf11523-ad8f-49c9-84dc-5554364a8b3b
plot_given && let
	plt.close() 

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))
	
	for name in keys(snapshots_τ)
		for (j, snap) in enumerate(snapshots_τ[name])
			axA.plot(
				profile(rms5, snap, :log10τ_ross, :uz)..., 
				color=colors[name][j], 
				label=labels[name][j]
			)
		end
	end

	axA.legend()
	axA.set_ylabel(L"\rm rms(U_z)\ [km\ \times\ s^{-1}]")
	axA.set_xlabel(L"\rm \log \tau_{ross}")
	
	gcf()
end

# ╔═╡ 292e9c4b-10d9-4cf9-b5e7-2594eec4bcfd
md"# Optical surface"

# ╔═╡ 795357b7-4d73-4bb6-899a-574183fc97e1
md"## Vertical velocity"

# ╔═╡ 12446032-6cc4-4eac-94cc-6ccb8de46c5d
extent(snap) = begin
	x = MUST.axis(snap, :x) ./1e8
	y = MUST.axis(snap, :x) ./1e8

	[minimum(x), maximum(x), minimum(y), maximum(y)]
end

# ╔═╡ 488b2eec-daf8-4920-9140-7d67c6ca3de1
begin
	uz = Dict(name=>[] for name in keys(snapshots_τ))
	for name in keys(snapshots_τ)
		for (j, snap) in enumerate(snapshots_τ[name])
			uz_τ_surf = if !isnothing(snap)
				isurf = MUST.closest(log10.(MUST.axis(snap, :τ_ross, 3)), 0)
				snap[:uz][:, :, isurf]
			else
				@info "Interpolating Uz..."
				uz_τ_surf = MUST.interpolate_to(snapshots[name][j], :uz, τ_ross=1)
				MUST.axis(uz_τ_surf, :uz, 3)
			end
	
			append!(uz[name], [uz_τ_surf])
		end
	end
	uz
end

# ╔═╡ 1ffcf11f-58dd-41dd-921f-995d0a84f0d0
if (nmodels > 0) && plot_given
	fDs, axDs = [], []

	vmin = minimum([minimum(u) for name in keys(uz) for u in uz[name]]) ./1e5
	vmax = maximum([maximum(u) for name in keys(uz) for u in uz[name]]) ./1e5

	vmax = min(abs.([vmin, vmax])...)
	vmin = -vmax

	for name in keys(snapshots)
		for (j, snap) in enumerate(snapshots[name])
			plt.close()
			
			fD, axD = plt.subplots(1, 1, figsize=(5, 6))
			
			i = axD.imshow(
				uz[name][j] ./1e5,
				origin="lower",
				vmin=vmin, vmax=vmax,
				extent=extent(snap),
				cmap="coolwarm"
			)
	
			cb = fD.colorbar(i, ax=axD, fraction=0.046, pad=0.04)
			cb.set_label(L"\rm U_z\ [km\ \times\ s^{-1}]")
	
			axD.set_xlabel("x [cm]")
			axD.set_ylabel("y [cm]")	
	
			append!(fDs, [fD])
			append!(axDs, [axD])
		end
	end
	
	gcf()
end

# ╔═╡ c46c5b32-a1cd-43d3-b087-7f6af7adb88d
md"# Time evolution"

# ╔═╡ 440097e4-f518-46fe-aa1f-3d446e4cbd25
md"## Pick models"

# ╔═╡ 296cbac2-97da-401c-a038-ff2fc8770dd3
md"Pick models for time evolution:"

# ╔═╡ 146e68f1-d03f-42b5-bd43-76ce21edcfae
if nmodels>0
	@bind models confirm(Select([keys(snapshots_picks)...]))
else
	models = []
end

# ╔═╡ 38767bff-5a3b-4085-b313-ddacb941da71
md"Tick box to start time evolution: $(@bind start_timeevolution CheckBox(default=false))"

# ╔═╡ 3225b57c-be7a-438c-83d0-e9b67303b23a


# ╔═╡ 873c07cd-7f02-4144-99fd-842a0aacc6d9
md"Pick what should be plotted as a function of time"

# ╔═╡ f2f27ccc-358b-4bec-9506-b31c27d6f759
begin
	xlim = [-4.5e8, 1.7e8]
	ylim = [2500, 17000]

	xlabel = L"\rm z\ [cm]"
	ylabel = L"\rm temperature\ [K]"

	# what to plot 
	profile_to_plot(f, s, sτ) = begin
		z, T = profile(f, s, :z, :T)
		z, R = profile(f, s, :z, :log10d)

		z, T
	end

	initial_profile_to_plot(m) = begin
		log10.(exp.(m.lnρ)), exp.(m.lnT)
		-m.z, exp.(m.lnT)
	end
end;

# ╔═╡ 7d10a077-75c1-4aee-b732-35a7ff71d109
md"## Evolution from Newton to RT"

# ╔═╡ 5b856c06-da72-41be-adc7-6d4e4570ad45
cmap_time = plt.cm.get_cmap("seismic")

# ╔═╡ 34fecdab-ed33-4685-92f2-70c0e763e893
md"Time evolution of snapshots:"

# ╔═╡ 139e7d21-e102-48de-bcbd-ee1a1e78bb5d
if start_timeevolution
	model_path = joinpath(MUST.@in_dispatch(datafolder), "$(models)")
	sc = converted_snapshots(model_path)
	snaps_converted = sort(MUST.list_snapshots(sc))

	plt.close()

	fF, axF = plt.subplots(1, 1, figsize=(5, 6))

	if models in keys(initial_model)
		axF.plot(
			initial_profile_to_plot(initial_model[models])..., 
			lw=1., 
			color=colors_initial[models],
			alpha=0.7,
			ls="--"
		)
	end

	#=if models in keys(initial_adiabats)
		axF.plot(
			initial_profile_to_plot(initial_adiabats[models])..., 
			lw=1., 
			color=colors_initial[models],
			alpha=0.7,
			ls="-"
		)
	end=#

	color_sequence = [
		cmap_time(i/length(snaps_converted)) for i in eachindex(snaps_converted)
	]

	#axF.axvline(0.57e8)
	
	for (c, i) in enumerate(snaps_converted)
		try
			s, st = pick_snapshot(sc, i)
			if c == 1
				axF.plot(
					profile_to_plot(mean, s, st)..., 
					lw=1., 
					color=color_sequence[c], 
					label="$(models)"
				)
			else
				axF.plot(
					profile_to_plot(mean, s, st)..., 
					lw=1., 
					color=color_sequence[c],
					alpha=0.5
				)
			end
		catch
			@warn "snapshot $i could not be loaded."
		end
	end

	if !isnothing(xlim)
		axF.set_xlim(xlim...)
	end

	if !isnothing(ylim)
		axF.set_ylim(ylim...)
	end
	
	axF.legend()
	axF.set_ylabel(ylabel)
	axF.set_xlabel(xlabel)
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─bd7cdaf2-c05e-4bc9-afbb-cfb61260d62c
# ╟─827e54c7-d0a3-455a-8837-52255eeff202
# ╠═4a837d4c-724c-11ee-0b6b-21a88a56df5b
# ╟─d9b5bbbd-6229-47c3-9738-dbe9087016d8
# ╟─9e9c5903-762d-43d7-adf2-0eccee134974
# ╟─9ac91858-17d8-4934-8d47-fed1519efe42
# ╟─73d9fd32-0dae-478a-80ff-a8e314adbd6e
# ╟─16786647-4021-4068-9af1-2a548240758e
# ╟─91ddaf2d-abd6-4e2d-92eb-cfa7d0ba79bc
# ╟─80fefa1e-de68-4448-a371-d4ae0ec70868
# ╟─de91f663-285b-4250-9edb-41f48d5b36e1
# ╟─6fea8f52-3ca6-4854-aca7-6b125d4c6542
# ╟─ac1069a6-0687-44fb-bd94-c76af8fd6b1f
# ╟─73e9bde0-c72e-4a01-aa68-486237e0c475
# ╟─68510276-b102-4a70-bd2f-a63195a4bb08
# ╟─2bf867fc-15bd-437c-b989-b1efbf2bf9d4
# ╟─fa8c564d-195b-424a-8814-0c47cd3e03a8
# ╟─2d684eb0-aa4b-40a6-9e04-ab0b42ed0516
# ╟─53706541-27d7-4093-b37f-3e384c6c8a7c
# ╟─ddcbcf4d-eb85-4cde-b2ac-a4b7e93ed809
# ╟─40124987-2dc6-452d-87d7-33f935854b55
# ╟─3778a669-69fc-47e6-9d0b-c514b31c8c79
# ╟─b9f82f31-f25e-4cad-abb4-aac2d4d0e584
# ╟─cfa74be4-06c9-4911-ac29-2a17809391e6
# ╟─ba8bded9-f13b-477e-a9fc-bd88e4f01e35
# ╟─ac46894f-44b4-477a-81ce-687b1308c38f
# ╟─574b4a88-c562-44ab-8a8e-6e5b5c5aabf0
# ╟─926aa121-ba10-4410-a91c-f36fc79a63ed
# ╟─3a94ff68-ea2f-4478-902f-4975074e8ec4
# ╟─ca860715-104c-4177-8249-a547863f6db4
# ╟─c5f4fb05-27b6-42e8-b08a-ef60c403d9e7
# ╟─4d0c7aa4-2a13-429a-88f2-14df415aa5ec
# ╟─27d0d01a-c8c8-422b-a893-36948b908fb0
# ╟─66f7d726-b840-484c-937f-b6a328e86b21
# ╟─1b354d9a-82f6-4540-8eaa-69218e37f87b
# ╟─f6c017c2-b296-4f82-8d2b-d97a7e36a51e
# ╟─693c5461-e4d8-45d0-af2d-87b10f9560e3
# ╟─2bb999c6-ed87-4006-9fba-2f01fcf22f90
# ╟─d55bae1d-1ca1-4537-8af1-c025a966c3b3
# ╟─271e4e9c-def2-4b3e-b383-45d95fbb309b
# ╟─c561e5a2-dfd6-4de9-80b5-9f6ae100b217
# ╟─55581a87-2638-4f85-9ad7-7313c37c1971
# ╟─e6d33b23-e77c-43d3-b401-c5620721d659
# ╟─bd41fb98-70c4-4192-8753-6fee82f14576
# ╟─b90c53aa-9ca6-421c-9b2f-f650a1f0832f
# ╟─db7db6d6-25c3-4f43-a49b-1ef5881606ac
# ╟─80f8c89e-5c14-4df6-8a05-f5fe79a7585e
# ╟─50bc1f52-0638-427a-8cc0-f8cd49535114
# ╟─3d1d9e3e-aa7d-4ef0-9c3e-41816285a938
# ╟─bb0cb1ab-7e03-4a5e-bfd8-60bf17e007d5
# ╟─b7c0a489-b7f7-4105-bce3-ae0c1074c650
# ╟─9ea4e08e-62aa-44c3-bbf0-78b1954a9c4f
# ╟─c4015f47-d670-4c5d-885c-cf32b6f15829
# ╟─13a647a4-3d71-4afb-ba31-867e108a8154
# ╟─17ca4c8c-9563-4be4-8ab7-a7446ba28dec
# ╟─1be872a8-eef5-422d-b3a6-13d581a999c0
# ╟─f3ffa698-bcea-4ab3-9b63-84d518c14068
# ╟─ce691c30-5025-4ffc-8185-eced3087ca13
# ╟─5e882897-d396-470c-ad46-37a39326225c
# ╟─ef9b0e78-00f0-41d0-8429-2d36f54c0a02
# ╟─e05efb00-233b-44b9-9204-156ae2ed0762
# ╟─36eabb43-1660-4e11-8bca-e7f09568d695
# ╟─b34da749-03ae-4f6c-92fb-c5809abeb339
# ╟─e99b4bed-f093-4fe6-ad1d-af0f2235470b
# ╟─ddcf3053-c8ff-4765-9a99-2c85534485fc
# ╟─d93e8b41-6f7a-4163-82dc-80a844e57e9d
# ╟─2527e21b-e956-4d9c-956d-22db9efbaadc
# ╟─14796f14-183e-496f-b5d8-41149b31d463
# ╟─6c4f5ac6-a03b-4237-aa8f-fe50f77bde6f
# ╟─4d623d6f-7366-468b-930e-71878b72aa5f
# ╟─3bf11523-ad8f-49c9-84dc-5554364a8b3b
# ╟─292e9c4b-10d9-4cf9-b5e7-2594eec4bcfd
# ╟─795357b7-4d73-4bb6-899a-574183fc97e1
# ╟─12446032-6cc4-4eac-94cc-6ccb8de46c5d
# ╟─488b2eec-daf8-4920-9140-7d67c6ca3de1
# ╟─1ffcf11f-58dd-41dd-921f-995d0a84f0d0
# ╟─c46c5b32-a1cd-43d3-b087-7f6af7adb88d
# ╟─440097e4-f518-46fe-aa1f-3d446e4cbd25
# ╟─296cbac2-97da-401c-a038-ff2fc8770dd3
# ╟─146e68f1-d03f-42b5-bd43-76ce21edcfae
# ╟─38767bff-5a3b-4085-b313-ddacb941da71
# ╟─3225b57c-be7a-438c-83d0-e9b67303b23a
# ╟─873c07cd-7f02-4144-99fd-842a0aacc6d9
# ╠═f2f27ccc-358b-4bec-9506-b31c27d6f759
# ╟─7d10a077-75c1-4aee-b732-35a7ff71d109
# ╠═5b856c06-da72-41be-adc7-6d4e4570ad45
# ╟─34fecdab-ed33-4685-92f2-70c0e763e893
# ╠═139e7d21-e102-48de-bcbd-ee1a1e78bb5d
