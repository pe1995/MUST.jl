### A Pluto.jl notebook ###
# v0.19.37

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

# ╔═╡ 7c75aef0-dba2-11ee-2e8a-d5ef839908b5
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using PlutoUI
	using PlutoUI: combine
	using CSV
	using DataFrames
end

# ╔═╡ 625bfcd7-a29f-4512-b00b-32cfcadb11e5
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	hl = 4

	mydir = pwd()

	TableOfContents()
end

# ╔═╡ c33fd24a-6225-4d09-a37e-2f8dcb38b75c
c2m = MUST.ingredients("convert2multi.jl")

# ╔═╡ 6178d065-d16f-4f0c-ad94-c86a2441ce12
# ╠═╡ show_logs = false
MUST.@import_m3dis "../../../m3dis/experiments/Multi3D"

# ╔═╡ d6fe0fdd-0a74-4400-8821-52fa1146facd
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ 4f67ca9b-a307-47e4-8e34-1015ceaa793e
md"# The Models"

# ╔═╡ 0a1bb888-d72f-4652-aa58-286cdbc54452
begin
	md"""Select the folder where the snapshots are located:\
$(@bind datafolder confirm(TextField(default=\"./packed_models/\")))"""
end

# ╔═╡ bbd47d44-2381-4c96-a729-4833fcf379b9
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	split.(allruns, "/", keepempty=false) .|> last
end;

# ╔═╡ 0950ca67-2e63-4b84-a926-61d287f555d3
begin
	md"""Select one of the available simulations:\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""
end

# ╔═╡ 96ad0583-b4f3-4b8d-b79d-e25752fa5752
begin
	runFolder = joinpath(datafolder, selectedRun)
	runFolderAbs = abspath(joinpath(datafolder, selectedRun))
	snapsConverted = MUST.converted_snapshots(runFolder)
	listSnapsConverted = MUST.list_snapshots(snapsConverted)
	listm3dsnaps = ["m3dis_$(i)" for i in listSnapsConverted]
	maskM3D = [isfile(joinpath(runFolder, "$(i).mesh")) for i in listm3dsnaps]
end;

# ╔═╡ cbe35062-4930-4777-a8ff-3d2a3e2ed90e
md"# The Opacities
In order to compute the effective temperature we need to save the relevant EoS in a format that M3D can read. For this the EoS has to be present on a temperature grid, and not as requested in DISPATCH on internal energy. The TSO.jl opacities should always be available on the temperature grid as well. You can only pick those EoS tables for conversion here! You will be stopped otherwise."

# ╔═╡ 4abba5e5-ba99-45a4-9995-822d942b9c0f
md"__Note:__ In earlier versions there was no 'eos_T.hdf5' included when creating the EoS on the energy grid. You can uncomment the following hidden code (and set the paths to whatever you need) to artificially copy over the EoS if they are available to you. If not, please contact `eitner@mpia.de` to get the correct EoS."

# ╔═╡ 93851417-7759-4f87-b1e2-bf64cb20e1ad
#=begin
	all_eos = MUST.glob("MainSequence/*_E_*/")
	mask = [isfile(joinpath(f, "eos.hdf5")) for f in all_eos]
	all_eos = all_eos[mask]

	all_eos_T = [split.(a, "_E_", keepempty=false) for a in all_eos]
	all_eos_T = [join(a, '_') for a in all_eos_T]
	mask = [isfile(joinpath(f, "eos.hdf5")) for f in all_eos_T]

	for i in eachindex(all_eos_T)
		if mask[i]
			cp(
				joinpath(all_eos_T[i], "eos.hdf5"), 
				joinpath(all_eos[i], "eos_T.hdf5")
			)
			cp(
				joinpath(all_eos_T[i], "binned_opacities.hdf5"), 
				joinpath(all_eos[i], "binned_opacities_T.hdf5")
			)
		end
	end
end=#

# ╔═╡ f687f43d-72b3-4158-9fc2-202368da57c2
eosInPickedFolder(folder) = begin
	allFolders = MUST.glob("*/", folder)
	mask = [isfile(joinpath(f, "eos.hdf5")) for f in allFolders]
	if count(mask) > 0
		first(allFolders[mask])
	else
		""
	end
end

# ╔═╡ ea69d72f-f970-4e53-9613-6d0ac666e295
function eos_input(folder)
	path = eosInPickedFolder(folder)
	name = "Pick the EoS for the selected run"
	return combine() do Child
		namesChild = [
			Child(name, TextField(max(length(path), 20), default=path))
		]
		inputs = [
			md""" $(name): $(
				namesChild[i]
			)"""
			for i in eachindex(namesChild)
		]
		
		md"""
		$(inputs)
		"""
	end
end

# ╔═╡ 48669b88-9cd9-4b43-a019-2e853dee56d7


# ╔═╡ 7c1a244c-7b63-4b76-93b6-e45d0e577be0
md"## Pick Binned Opacities for M3D"

# ╔═╡ 78d96ea8-c8ba-48b4-b339-859934da9f97
@bind eospathSelect eos_input(runFolder)

# ╔═╡ f9b1375a-bd52-484c-ad86-653b5c792366
eospath = abspath(values(eospathSelect) |> first);

# ╔═╡ 87ae460e-99d4-4d81-98d2-8d6ea61a3100
begin
	selectedRun
	
	eos, opa = if isdir(eospath)
		# first check normal EoS
		eos = reload(SqEoS, joinpath(eospath, "eos.hdf5"))
		if is_internal_energy(@axed(eos))
			if !isfile(joinpath(eospath, "eos_T.hdf5"))
				error("The given EoS is not on the temperature grid and does not contain a 'eos_T.hdf5' alternative. If you don't have access to the correct EoS please contact eitner@mpia.de.")
			end
			e = reload(SqEoS, joinpath(eospath, "eos_T.hdf5"))
			o = reload(SqOpacity, joinpath(eospath, "binned_opacities_T.hdf5"))
	
			e, o
		else
			o = reload(SqOpacity, joinpath(eospath, "binned_opacities.hdf5"))
	
			eos, o
		end
	else
		@warn "The chosen EoS directory does not exists"
		nothing, nothing
	end
end;

# ╔═╡ 6983d7d9-4e06-45fd-abb5-29adeba5a994
if !isnothing(eos)
	selectedRun
	#newFolder = "input_multi3d/$(selectedRun)"
	#newFolder = @in_m3dis newFolder
	#!isdir(newFolder) && mkdir(newFolder)

	# save the binary eos and opacities
	f = save(eos, opa, joinpath(eospath, "eos_opa"))

	@info "Opacity table saved at $(f)"
end;

# ╔═╡ 2eff8c20-3ea5-4a72-9917-bd217e817e47
begin
	if count(.!maskM3D) > 0
		c2m.snaps2multi(
			runFolder, listSnapsConverted[.!maskM3D],
			eos=eos, label=""
		)
	end
end

# ╔═╡ 5cdd9197-bc7b-4f78-b781-c076c94b73e9
md"## Resampling Cubes for M3D"

# ╔═╡ a5633b70-aefa-485c-9fdf-3c4fa42c6301
md"Horizontal resolution: $(@bind nhorizontal TextField(default=\"30\"))"

# ╔═╡ ca55c334-5ab2-4af1-9056-aee0d43649cc
md"Upsample vertical resolution by: $(@bind xvertical TextField(default=\"2\"))"

# ╔═╡ d6ab2bc3-22a2-480d-9a31-fce57dc5682e


# ╔═╡ df7717ed-0ffc-46e9-bae6-9f12ce893d43
md"Tick to resample for M3D: $(@bind resampleM3D CheckBox(default=false))"

# ╔═╡ 1512fa8f-51c6-48ff-9ac8-5a7e21c569f0
listm3dsnapsResample = if resampleM3D
	@info "Resampling cubes..."
	snap_test = pick_snapshot(runFolderAbs, first(listSnapsConverted)) |> first
	nz = size(snap_test, 3) * parse(Int, xvertical)
	names_new = ["m3dis_$(nhorizontal)x$(nz)_$(i)" for i in listSnapsConverted]
	already_converted_m3d = isfile.(
		joinpath.(runFolderAbs, ["$(m).mesh" for m in names_new])
	)

	c2m.snaps2multi(
		runFolderAbs, 
		listSnapsConverted[.!already_converted_m3d]..., 
		n_horizontal=parse(Int, nhorizontal),
		n_vertical=nz,
		eos=eos, 
		label="$(nhorizontal)x$(nz)", 
		outfolder=runFolderAbs
	)
	
	@info "...cubes resampled."

	["m3dis_$(nhorizontal)x$(nz)_$(i)" for i in listSnapsConverted]
else
	listm3dsnaps
end;

# ╔═╡ aef4e240-d337-45f3-83ea-a260d6158b3b


# ╔═╡ 489f798b-3d3d-4398-bff3-a27129f2581a
md"# Running M3D"

# ╔═╡ b705039c-dfd8-4876-98be-598493b4b783
md"Tick to start computing: $(@bind compute CheckBox(default=false))"

# ╔═╡ 25887485-2292-4097-95e2-10bbae6190fc
begin
	if compute
		for i in eachindex(listSnapsConverted)
			MUST.heating(                                 
				listm3dsnapsResample[i], 
				joinpath(eospath, "eos_opa"),
				namelist_kwargs=(
					:model_folder=>runFolderAbs,
					:linelist=>nothing,
					:absmet=>nothing,
					:atom_params=>(:atom_file=>"", ),
					:atmos_params=>(
						:atmos_format=>"MUST", 
						:use_density=>true, 
						:use_ne=>false
					),
				), 
				slurm=false
			)
		end
	end
	if !(pwd()==mydir)
		cd(mydir)
	end
end

# ╔═╡ 55e916b7-c06f-4fb1-9a43-5e938f308e77


# ╔═╡ ea537cd6-982a-4613-ba13-a8b337e00528
md"# Effective Temperature"

# ╔═╡ 280f8e1d-4809-4063-a8a1-1df5afd3c8ed
md"## Compute"

# ╔═╡ f49c2f6c-a884-47ed-aac8-3d119efb15c2
begin
	compute
	m3dpaths = ["data/$(m)_binned" for m in listm3dsnapsResample]
end

# ╔═╡ 9de1b403-94bb-4cfb-ab36-eb19f962f141
runs = [M3DISRun(p) for p in m3dpaths]

# ╔═╡ 9dfb8acf-5615-4b43-a0bc-e762855c8d52
teff = MUST.Teff.([r.flux for r in runs])

# ╔═╡ b6149352-64ab-46a6-b839-1697ece0aa64
md"## Visualize"

# ╔═╡ 7ec96609-3ca7-46c8-8145-da360381d690
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	ax.plot(1:length(teff) |> collect, teff, marker="s", color="k", markersize=8)

	ax.set_xlabel("snapshots")
	ax.set_ylabel("effective temperature [K]")
	

	gcf()
end

# ╔═╡ 9dfd9024-3d80-4056-a943-f593177fc8ee
md"## Save"

# ╔═╡ b8bec454-0f13-472a-a473-1b9d57d64233
begin
	teff_data = Dict(
		"snaps"=>listSnapsConverted,
		"teff"=>teff
	)
	CSV.write(joinpath(runFolder, "effective_temperature.csv"), DataFrame(teff_data))
end

# ╔═╡ 371551db-dddf-4003-99d7-7f8c20e3068d


# ╔═╡ e7a96d25-27a9-4a90-932f-7b60cfef0a2c
md"Tick box to add information to data cubes: $(@bind addToCubes CheckBox(default=false))"

# ╔═╡ a1965f51-4048-4eb1-ba9e-799a3096116f
function setTeff!(name, teff)
	fid = MUST.HDF5.h5open(name, "r+")
	MUST.HDF5.delete_object(fid, "teff")
	fid["teff"] = teff
	close(fid)
end

# ╔═╡ e87f7d87-ba85-4887-b5b6-ccdeaf472c54
if addToCubes
	# path to cubes
	for (i, s) in enumerate(listSnapsConverted)
		name = joinpath(runFolderAbs, "box_sn$(s).hdf5")
		setTeff!(name, teff[i])
		
		name = joinpath(runFolderAbs, "box_tau_sn$(s).hdf5")
		setTeff!(name, teff[i])
		@info "New Teff ($(teff[i])) written to cube $(s)."
	end
end

# ╔═╡ 14a3621d-b1ce-4d27-a3b9-31a2453cac6b


# ╔═╡ Cell order:
# ╠═7c75aef0-dba2-11ee-2e8a-d5ef839908b5
# ╟─625bfcd7-a29f-4512-b00b-32cfcadb11e5
# ╟─c33fd24a-6225-4d09-a37e-2f8dcb38b75c
# ╠═6178d065-d16f-4f0c-ad94-c86a2441ce12
# ╠═d6fe0fdd-0a74-4400-8821-52fa1146facd
# ╟─4f67ca9b-a307-47e4-8e34-1015ceaa793e
# ╟─0a1bb888-d72f-4652-aa58-286cdbc54452
# ╟─bbd47d44-2381-4c96-a729-4833fcf379b9
# ╟─0950ca67-2e63-4b84-a926-61d287f555d3
# ╟─96ad0583-b4f3-4b8d-b79d-e25752fa5752
# ╟─cbe35062-4930-4777-a8ff-3d2a3e2ed90e
# ╟─4abba5e5-ba99-45a4-9995-822d942b9c0f
# ╟─93851417-7759-4f87-b1e2-bf64cb20e1ad
# ╟─f687f43d-72b3-4158-9fc2-202368da57c2
# ╟─ea69d72f-f970-4e53-9613-6d0ac666e295
# ╟─48669b88-9cd9-4b43-a019-2e853dee56d7
# ╟─7c1a244c-7b63-4b76-93b6-e45d0e577be0
# ╟─78d96ea8-c8ba-48b4-b339-859934da9f97
# ╟─f9b1375a-bd52-484c-ad86-653b5c792366
# ╟─87ae460e-99d4-4d81-98d2-8d6ea61a3100
# ╟─6983d7d9-4e06-45fd-abb5-29adeba5a994
# ╟─2eff8c20-3ea5-4a72-9917-bd217e817e47
# ╟─5cdd9197-bc7b-4f78-b781-c076c94b73e9
# ╟─a5633b70-aefa-485c-9fdf-3c4fa42c6301
# ╟─ca55c334-5ab2-4af1-9056-aee0d43649cc
# ╟─d6ab2bc3-22a2-480d-9a31-fce57dc5682e
# ╟─df7717ed-0ffc-46e9-bae6-9f12ce893d43
# ╟─1512fa8f-51c6-48ff-9ac8-5a7e21c569f0
# ╟─aef4e240-d337-45f3-83ea-a260d6158b3b
# ╟─489f798b-3d3d-4398-bff3-a27129f2581a
# ╟─b705039c-dfd8-4876-98be-598493b4b783
# ╟─25887485-2292-4097-95e2-10bbae6190fc
# ╟─55e916b7-c06f-4fb1-9a43-5e938f308e77
# ╟─ea537cd6-982a-4613-ba13-a8b337e00528
# ╟─280f8e1d-4809-4063-a8a1-1df5afd3c8ed
# ╟─f49c2f6c-a884-47ed-aac8-3d119efb15c2
# ╟─9de1b403-94bb-4cfb-ab36-eb19f962f141
# ╟─9dfb8acf-5615-4b43-a0bc-e762855c8d52
# ╟─b6149352-64ab-46a6-b839-1697ece0aa64
# ╟─7ec96609-3ca7-46c8-8145-da360381d690
# ╟─9dfd9024-3d80-4056-a943-f593177fc8ee
# ╟─b8bec454-0f13-472a-a473-1b9d57d64233
# ╟─371551db-dddf-4003-99d7-7f8c20e3068d
# ╟─e7a96d25-27a9-4a90-932f-7b60cfef0a2c
# ╟─a1965f51-4048-4eb1-ba9e-799a3096116f
# ╟─e87f7d87-ba85-4887-b5b6-ccdeaf472c54
# ╟─14a3621d-b1ce-4d27-a3b9-31a2453cac6b
