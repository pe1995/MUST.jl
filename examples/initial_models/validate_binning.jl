### A Pluto.jl notebook ###
# v0.19.38

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

# ╔═╡ 59502f3c-f1c7-11ee-1094-7f6cf8bd94f7
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot
	using LaTeXStrings
	using PlutoUI
	using Glob
end

# ╔═╡ ed8e41d3-1db2-4dbd-a0e0-d350506bde3a
md"# Binning Validation"

# ╔═╡ 51b886e1-87c7-483d-b762-9459bed3bd4e
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(TSO)), "Bergemann2023.mplstyle"))

	TableOfContents()
end

# ╔═╡ 954bff57-8ad5-4801-9f76-621bf5f723f8
md"# Selection"

# ╔═╡ 6e051fba-7269-49a2-93e2-7b699f453de9
md"## Unbinned Opacities"

# ╔═╡ f0ec0348-c6c9-4318-b76d-f927b1bc150f
md"__Select the directory where the unbinned opacity tables are stored on your system__"

# ╔═╡ 307a49fc-1592-4796-b71c-4fd900b832bd
@bind unbinnedDir confirm(
	TextField(
		length("/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/"),
		default="/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/"
	)
)

# ╔═╡ eb8b4ea1-3bb3-44ea-9768-781363e36507
findEoSAvailable(folder; kind=:unbinned) = begin
	alleos = glob("*/", folder)
	if length(alleos) == 0
		return [nothing]
	end
	alleosNames = convert.(String, last.(split.(alleos, "/", keepempty=false)))

	allfiles = [glob("*", f) for f in alleos]
	catchphrase = kind==:unbinned ? "combined_eos" : "eos"
	mask = [any(occursin.(catchphrase, allfiles[i])) for i in eachindex(allfiles)]

	mask2 = kind==:binned ? occursin.("_E_", alleosNames) : trues(length(mask))
	mask = mask .& mask2
	
	if count(mask) == 0
		return [nothing]
	end
	alleosNames[mask]
end;

# ╔═╡ 543c7cbc-a060-4ba7-a09b-685f8692864e
md"__Pick a EoS among those that are available at the given location__\
$(@bind unbinnedName Select(findEoSAvailable(unbinnedDir, kind=:unbinned)))"

# ╔═╡ e41dd9d3-5d80-4b2a-a584-69726d5e245d


# ╔═╡ 44f5edca-b77b-41e1-9a8b-964ce9e7bf9e
begin
	unbinnedEoS, unbinnedOpa = if !isnothing(unbinnedName)
		unbinnedEoSPath = joinpath(unbinnedDir, unbinnedName)
		unbinnedalleos = glob("*eos*", unbinnedEoSPath)
		unbinnedalleos = [String(split(s, '/', keepempty=false) |> last) for s in unbinnedalleos]
		unbinnedalleos = String.(unbinnedalleos)
		iselectEos = if any(occursin.("ross_", unbinnedalleos))
			findfirst(occursin.("ross_", unbinnedalleos))
		else
			1
		end
		
		unbinnedEoS = reload(
			SqEoS, 
			joinpath(unbinnedEoSPath, unbinnedalleos[iselectEos])
		)

		unbinnedalleos = glob("combined_opacities*", unbinnedEoSPath)
		unbinnedalleos = [String(split(s, '/', keepempty=false) |> last) for s in unbinnedalleos]
		unbinnedalleos = String.(unbinnedalleos)

		unbinnedOpa = reload(
			SqOpacity, 
			joinpath(unbinnedEoSPath, first(unbinnedalleos)),
			mmap=true
		)

		unbinnedEoS, unbinnedOpa
	else
		nothing, nothing
	end
end

# ╔═╡ b873758d-f946-4b6b-a3f7-e41b394fde97
md"## Binned Opacities"

# ╔═╡ 24246eac-00fa-4e08-8314-3eb7ded046a1
md"Select the binned opacities you wish to evaluate. You can pick multiple tables. For each selected table, the unbinned and binned radiative transfer will be solved for the 1D model that is given as initial condition in the selected folder. The EoS on the temperature grid is assumed to be given as `eos_T.hdf5` within the energy folder."

# ╔═╡ 801438cf-0af8-4007-a587-0c029d1c9472


# ╔═╡ 9c743d23-9316-410f-b303-07853663d681
md"__Select the directory where the binned opacity tables are stored on your system__"

# ╔═╡ 974d13f2-8a26-4621-ba3c-0fd83e2dd596
@bind binnedDir confirm(
	TextField(
		default="./"
	)
)

# ╔═╡ 13b35a40-8bee-40e1-bb8a-530ad340314c


# ╔═╡ 05f5a4ce-4f2c-4b35-82b8-a8701d7199a8
md"__Pick one (or multiple) of the available EoS__"

# ╔═╡ 0949424e-94e0-4743-87a0-27200f90037c
@bind binnedNames MultiCheckBox(sort(findEoSAvailable(binnedDir, kind=:binned)))

# ╔═╡ fa334242-46fe-4d33-b1fe-a5c255afd19d
load_eos(folder) = begin
	eosT = reload(
		SqEoS,
		joinpath(binnedDir, folder, "eos_T.hdf5")
	)

	opaT = reload(
		SqOpacity,
		joinpath(binnedDir, folder, "binned_opacities_T.hdf5"),
		mmap=false
	)

	eosT, opaT
end;

# ╔═╡ 6e0a68d5-076e-44d1-aa05-90b2c0a2c523
binnedEoS, binnedOpa = if length(binnedNames) > 0
	eosopa = load_eos.(binnedNames)

	(
		Dict(binnedNames[i]=>v for (i, v) in enumerate(first.(eosopa))), 
		Dict(binnedNames[i]=>v for (i, v) in enumerate(last.(eosopa)))
	)
else
	Dict(), Dict()
end

# ╔═╡ 96aa3d22-39f4-49b7-998b-8c97c924527a


# ╔═╡ 69fa30e6-6949-4c2c-b914-4ea9bcb8912b
md"# Radiative Transfer"

# ╔═╡ 49b18999-1ab6-4dd0-aa6a-d91b9e630d45
md"## Opacity"

# ╔═╡ 5c1c3219-4427-4368-8618-d82b9421598d
opacityUnbinned = TSO.@binned(unbinnedOpa, unbinnedEoS)

# ╔═╡ 3c50d684-0648-4512-b9dd-f81b70943628
opacityBinned = Dict(
	binnedNames[i] => TSO.@binned(binnedOpa[v]) 
	for (i, v) in enumerate(binnedNames)
)

# ╔═╡ a50d8f3d-2a94-4b5f-b58d-8e8708149209


# ╔═╡ bbd9666c-3e56-4a2c-993d-08c67c5a00a2
md"""
## 1D Models
We use the 1D given in the EoS as reference models. We compute the heating for only one model.
"""

# ╔═╡ c053cd06-3ebd-4eaf-aad7-a716ab4482f4


# ╔═╡ dd7d534e-480c-4583-82a6-2b002db20644
md"__Pick 1D model to compute the heating__\
$(@bind modelSelected Select(binnedNames))"

# ╔═╡ 6946afe6-d93f-4bd6-8adc-0bdeb3c34c84


# ╔═╡ c3e8e5f4-115c-4f8c-9373-5ea3d2b5c190
modelPath = joinpath(binnedDir, modelSelected, "inim.dat")

# ╔═╡ 7422a2bd-c8be-4475-90ec-d7cc5c688ebd
model = TSO.@optical(
	TSO.Average3D(unbinnedEoS, modelPath), unbinnedEoS, unbinnedOpa
)

# ╔═╡ 72e01a96-5117-4dc0-ab6a-d549d5754d2f


# ╔═╡ 5319cd97-9f59-4866-99b6-cbc03e7ddff2
md"## RT Solvers"

# ╔═╡ 33ab5c75-3a2e-41f8-8fd9-0f22ef407b0a
ω = weights = TSO.ω_midpoint(unbinnedOpa)

# ╔═╡ 3afd9c2f-5a3d-45a9-98f9-34129fa16518
solversBinned = Dict(
	v=>TSO.Heating(
		model, 
		eos=TSO.@axed(binnedEoS[v]), 
		opacities=opacityBinned[v]
	)
	for v in binnedNames
)

# ╔═╡ 13a71e49-797f-4321-bc9c-cac3b8169cb0
solverUnbinned = TSO.Heating(
	model, 
	eos=TSO.@axed(unbinnedEoS), 
	opacities=opacityUnbinned
)

# ╔═╡ 5c68e925-a79f-4f1d-9888-4beb63cca6a8


# ╔═╡ 2dc352a7-7fd9-4376-b1d9-1e40d3ff2900
md"## Solve RT"

# ╔═╡ 16c394c1-7d5a-4c54-976d-a85eb92ae3d6


# ╔═╡ f22615c9-4f54-4f98-b684-7f23e4793056
QBinned = Dict(
	v=>TSO.solve(solver)
	for (v, solver) in solversBinned
)

# ╔═╡ fbfc98eb-af52-4079-a997-54c362dbb7bf
begin
	unbinnedName
	QUnbinned = Dict()
end;

# ╔═╡ df6089c7-014d-4d5e-bd1e-da7e95f5b6c9
md"__Start solving monochromatic:__ $(@bind startSolving CheckBox(default=false))"

# ╔═╡ 71cae20d-f64d-414b-bc63-4cd46a182ad3
if startSolving
	if !(modelSelected in keys(QUnbinned))
		@info "Solving monochromatic RT for $(modelSelected)"
		QUnbinned[modelSelected] = TSO.solve(solverUnbinned, weights=ω)
	else
		@info "Monochromatic RT for $(modelSelected) already solved."
	end
end

# ╔═╡ f536c66c-c730-4c79-a1c7-7914d7124a1f
@info QUnbinned

# ╔═╡ bde98193-5946-4e28-91f4-7c924ef05b6b


# ╔═╡ e9d17a7f-cfe3-4ce0-a8c6-45135d532de9
md"# Heating Comparison"

# ╔═╡ 9305ee85-aa96-4182-8490-fa41c3abeef3
begin
	z, lnT, τ, lnρ = QUnbinned[modelSelected].model.z, QUnbinned[modelSelected].model.lnT, QUnbinned[modelSelected].model.τ,  QUnbinned[modelSelected].model.lnρ
end

# ╔═╡ 21050195-2a39-411d-a6a4-e96305d995d8
cmap = plt.get_cmap("rainbow")

# ╔═╡ cdfd9f56-28c5-4bb8-9aec-6873bd92ed84
begin
	colors = Dict()
	neos = length(binnedNames)
	ic = 1
	for name in binnedNames
		colors[name] = cmap(ic/neos)
		ic += 1
	end

	colors
end

# ╔═╡ f3ba0ad3-ce5a-4b87-a3a0-242490738cfb
function names_input(simulation_names::Vector)
	bins = [length(binnedOpa[s].λ) for s in simulation_names]
	return PlutoUI.combine() do Child
		namesChild = [
			Child(simulation_names[i], TextField(50, default="$(bins[i]) Bins"))
			for i in eachindex(simulation_names)
		]
		inputs = [
			md""" $(name): $(
				namesChild[i]
			)"""
			for (i, name) in enumerate(simulation_names)
		]
		
		md"""
		__Modify labels (optional)__
		$(inputs)
		"""
	end
end;

# ╔═╡ f1ee86da-5476-4e66-9ba0-3175c65798c3
@bind changedNames names_input(binnedNames)

# ╔═╡ b59bc64c-aabc-4ece-a596-ab2d17eaae99


# ╔═╡ 8e82066e-7fba-46bd-80d9-9050cebe87c0
md"## Heating throughout the atmosphere"

# ╔═╡ e9f91c4d-3094-4e84-82c5-367e439866c6
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	norm = maximum(abs.(TSO.heating(QUnbinned[modelSelected])))

	for k in binnedNames
		q = QBinned[k]
		ax.plot(
			log10.(τ), TSO.heating(q)./norm, 
			color=colors[k], label=changedNames[Symbol(k)], lw=2
	)	
	end
	
	ax.plot(
		log10.(τ), TSO.heating(QUnbinned[modelSelected])./norm, 
		marker="x", color="k", ls="",
		markerfacecolor="w", markersize=7, label="monochromatic"
	)


	ax.set_ylabel(L"\rm q_{rad}\ /\ \left|q_{rad}^{\lambda, max}\right|")	
	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_xlim(-6, 4)

	ax.legend()
	
	gcf()
end

# ╔═╡ 4deb92b6-ad19-4d86-92f5-4cf34a5ead6e
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	qu = TSO.heating(QUnbinned[modelSelected])
	norm = maximum(abs.(qu)) / 100.
	stat_label(yi) = "$(round(maximum(abs.(yi))/norm, sigdigits=3)) %"

	for k in binnedNames
		q = QBinned[k]
		y = (TSO.heating(q) .- qu)
		ax.plot(
			log10.(τ), y./norm, 
			color=colors[k], label=changedNames[Symbol(k)]*L",\ \rm max(\delta q) = "*stat_label(y), lw=2
	)	
	end
	
	ax.set_ylabel(L"\rm \left(q_{rad}^{bin} - q_{rad}^{\lambda} \right)\ /\ \left|q_{rad}^{\lambda, max}\right|\ [\%]")	
	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_xlim(-6, 4)

	ylm = max(abs.([ax.get_ylim()...])...)
	ax.set_ylim(-ylm, ylm)

	#=ax.legend(
		labelspacing=0.05, handlelength=2, loc="upper left", ncols=1,
		bbox_to_anchor=(1.0, 0.99)
	)=#
	ax.legend(
		labelspacing=0.05, handlelength=2
	)
	
	gcf()
end

# ╔═╡ adbd2b50-4976-4719-9eec-118335748672


# ╔═╡ 3528e6cf-31b8-4c08-a6d5-9d9658488c19
md"## Assignment"

# ╔═╡ 31075195-eaf0-4520-9a8c-d3283cdd6a07
bin_assignment(path) = begin
    fid = TSO.HDF5.h5open(path, "r")
    bins = TSO.HDF5.read(fid["bins"])
    λ = TSO.HDF5.read(fid["lambda"])
    close(fid)

    bins, λ
end

# ╔═╡ 58838376-2155-427e-91b7-3fde46f80cf1


# ╔═╡ 274022b1-b556-4d5e-9e70-f0a7d8f6cccf
md"__Select Opacities to show bin assignment:__ $(@bind eosBinnedPath Select(binnedNames))"

# ╔═╡ 3eb1bb0c-541e-49e7-97bc-cb9321c6ec3e


# ╔═╡ 09f5cb4a-1137-4245-850d-78508690912c
a, lam = try
	bin_assignment(joinpath(eosBinnedPath, "bin_assignment.hdf5"))
catch
	nothing, nothing
end

# ╔═╡ c8ef8c95-561d-4587-b375-5ec42fa10ddb
if !isnothing(a)
	let
		plt.close()
	
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
		im = ax.plot(
			log10.(unbinnedOpa.λ), a, ls="", color="k", lw=0.1, marker="."
		)
	
		ax.set_ylabel(L"\rm bin")
		ax.set_xlabel(L"\rm \log_{10}\lambda")
	
		gcf()
	end
else
	@info "No bin assignment found at the given location."
end

# ╔═╡ Cell order:
# ╟─ed8e41d3-1db2-4dbd-a0e0-d350506bde3a
# ╠═59502f3c-f1c7-11ee-1094-7f6cf8bd94f7
# ╟─51b886e1-87c7-483d-b762-9459bed3bd4e
# ╟─954bff57-8ad5-4801-9f76-621bf5f723f8
# ╟─6e051fba-7269-49a2-93e2-7b699f453de9
# ╟─f0ec0348-c6c9-4318-b76d-f927b1bc150f
# ╟─307a49fc-1592-4796-b71c-4fd900b832bd
# ╟─eb8b4ea1-3bb3-44ea-9768-781363e36507
# ╟─543c7cbc-a060-4ba7-a09b-685f8692864e
# ╟─e41dd9d3-5d80-4b2a-a584-69726d5e245d
# ╟─44f5edca-b77b-41e1-9a8b-964ce9e7bf9e
# ╟─b873758d-f946-4b6b-a3f7-e41b394fde97
# ╟─24246eac-00fa-4e08-8314-3eb7ded046a1
# ╟─801438cf-0af8-4007-a587-0c029d1c9472
# ╟─9c743d23-9316-410f-b303-07853663d681
# ╟─974d13f2-8a26-4621-ba3c-0fd83e2dd596
# ╟─13b35a40-8bee-40e1-bb8a-530ad340314c
# ╟─05f5a4ce-4f2c-4b35-82b8-a8701d7199a8
# ╟─0949424e-94e0-4743-87a0-27200f90037c
# ╟─fa334242-46fe-4d33-b1fe-a5c255afd19d
# ╟─6e0a68d5-076e-44d1-aa05-90b2c0a2c523
# ╟─96aa3d22-39f4-49b7-998b-8c97c924527a
# ╟─69fa30e6-6949-4c2c-b914-4ea9bcb8912b
# ╟─49b18999-1ab6-4dd0-aa6a-d91b9e630d45
# ╟─5c1c3219-4427-4368-8618-d82b9421598d
# ╟─3c50d684-0648-4512-b9dd-f81b70943628
# ╟─a50d8f3d-2a94-4b5f-b58d-8e8708149209
# ╟─bbd9666c-3e56-4a2c-993d-08c67c5a00a2
# ╟─c053cd06-3ebd-4eaf-aad7-a716ab4482f4
# ╟─dd7d534e-480c-4583-82a6-2b002db20644
# ╟─6946afe6-d93f-4bd6-8adc-0bdeb3c34c84
# ╟─c3e8e5f4-115c-4f8c-9373-5ea3d2b5c190
# ╟─7422a2bd-c8be-4475-90ec-d7cc5c688ebd
# ╟─72e01a96-5117-4dc0-ab6a-d549d5754d2f
# ╟─5319cd97-9f59-4866-99b6-cbc03e7ddff2
# ╟─33ab5c75-3a2e-41f8-8fd9-0f22ef407b0a
# ╟─3afd9c2f-5a3d-45a9-98f9-34129fa16518
# ╟─13a71e49-797f-4321-bc9c-cac3b8169cb0
# ╟─5c68e925-a79f-4f1d-9888-4beb63cca6a8
# ╟─2dc352a7-7fd9-4376-b1d9-1e40d3ff2900
# ╟─16c394c1-7d5a-4c54-976d-a85eb92ae3d6
# ╟─f22615c9-4f54-4f98-b684-7f23e4793056
# ╟─fbfc98eb-af52-4079-a997-54c362dbb7bf
# ╟─df6089c7-014d-4d5e-bd1e-da7e95f5b6c9
# ╟─71cae20d-f64d-414b-bc63-4cd46a182ad3
# ╟─f536c66c-c730-4c79-a1c7-7914d7124a1f
# ╟─bde98193-5946-4e28-91f4-7c924ef05b6b
# ╟─e9d17a7f-cfe3-4ce0-a8c6-45135d532de9
# ╠═9305ee85-aa96-4182-8490-fa41c3abeef3
# ╟─21050195-2a39-411d-a6a4-e96305d995d8
# ╟─cdfd9f56-28c5-4bb8-9aec-6873bd92ed84
# ╟─f3ba0ad3-ce5a-4b87-a3a0-242490738cfb
# ╟─f1ee86da-5476-4e66-9ba0-3175c65798c3
# ╟─b59bc64c-aabc-4ece-a596-ab2d17eaae99
# ╟─8e82066e-7fba-46bd-80d9-9050cebe87c0
# ╟─e9f91c4d-3094-4e84-82c5-367e439866c6
# ╟─4deb92b6-ad19-4d86-92f5-4cf34a5ead6e
# ╟─adbd2b50-4976-4719-9eec-118335748672
# ╟─3528e6cf-31b8-4c08-a6d5-9d9658488c19
# ╟─31075195-eaf0-4520-9a8c-d3283cdd6a07
# ╟─58838376-2155-427e-91b7-3fde46f80cf1
# ╟─274022b1-b556-4d5e-9e70-f0a7d8f6cccf
# ╟─3eb1bb0c-541e-49e7-97bc-cb9321c6ec3e
# ╟─09f5cb4a-1137-4245-850d-78508690912c
# ╟─c8ef8c95-561d-4587-b375-5ec42fa10ddb
