### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 49ce5a16-bc56-11ef-0e57-1f9a55abaaf8
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using PlutoUI
	using LaTeXStrings
	using DataFrames
	using CSV
	using StatsBase
	using FITSIO
	using DelimitedFiles
end

# ╔═╡ ad04ebd0-6bb5-47d1-8dcb-c5d9b1346d91
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))

	MUST.@import_m3dis("../../../Multi3D")
	TableOfContents()
end

# ╔═╡ 1a3653d8-74ae-4a35-a577-3994edc2adc5
md"# Models
Models can be selected from a folder containing all models as sub-folders."

# ╔═╡ 6160a040-944e-4e7d-8d4c-c12262145ad7
md"## Load Models"

# ╔═╡ 1b9d3022-c1b3-44d0-a444-a27aa902d9f0
md"Pick directory of grid models: $(@bind gridfolder confirm(TextField(default=\"./MS_models/\")))"

# ╔═╡ 4ec4a69e-0a15-4353-ac0c-bfb9e02ebdee
availableRuns(folder, inName...) = begin
	allruns = MUST.glob("*/", folder)
	n = split.(allruns, "/", keepempty=false) .|> last
	mask = falses(length(allruns))
	for inn in inName
		mask .= mask .| [occursin(inn, ni) for ni in n]
	end
	n[mask]
end

# ╔═╡ a2e3351b-fd71-4231-b010-1f55aead47a0
model_extension = ["PLATO"]#["ST20", "ST24", "P1A", "P1B", "P1C", "P1D"]

# ╔═╡ 9ed0a9ef-8f97-4aeb-9e52-c05ae6ec59c9
allRuns = availableRuns(gridfolder, model_extension...)

# ╔═╡ ac701824-7442-48ff-85bb-72315b2ad23f
begin
	grids = [
		#=MUST.Atmos1DGrid(joinpath(gridfolder, "ST20_dispatch_2024-12-13.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "ST24A_dispatch_2024-12-17.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "ST24B_dispatch_2024-12-17.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "ST24C_dispatch_2024-12-17.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "ST24D_dispatch_2024-12-17.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "P1A_dispatch_2025-02-19.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "P1B_dispatch_2025-02-19.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "P1C_dispatch_2025-02-19.mgrid")),
		MUST.Atmos1DGrid(joinpath(gridfolder, "P1D_dispatch_2025-02-19.mgrid"))=#
		MUST.Atmos1DGrid(mg) for mg in MUST.glob("*.mgrid", gridfolder)
	]

	grid = sum(grids)
end

# ╔═╡ 8c1bcdae-0d29-436f-9e81-4c7f8c763dcb


# ╔═╡ 865e1497-bc8f-403b-be99-c30d59c2f493
md"## Model Parameters
for each model we can first of all extract the parameters from the label of the folder. We could alternatively also get those values from the model cubes themselves, or from the namelist that is included."

# ╔═╡ 7d853824-7b07-4dc9-8a00-ba026efb6c1b
function parametersFromName(path; carbon=false, vmic=false)
	datadir = if occursin('/', path)
		split(path, '/', keepempty=false) .|> last
	else
		path
	end
	
	parts = split(datadir, '_', keepempty=false)
	irelevantPart = findfirst(
		x->(occursin("t", x)&occursin("g", x)&occursin("m", x)),
		parts
	)
	relevantPart = parts[irelevantPart]
	t1 = findfirst(x->x=='t', relevantPart) + 1
	t2 = findfirst(x->x=='g', relevantPart) - 1
	t3 = findfirst(x->x=='m', relevantPart) - 1

	t = parse(Float64, relevantPart[t1:t2]) .* 100
	g = parse(Float64, relevantPart[t2+2:t3]) ./10
	m = parse(Float64, relevantPart[t3+2:end])

	d = Dict("teff"=>t, "logg"=>g, "feh"=>m)
	
	# find other stuff
	if carbon
		relevantPart = parts[2]
		t1 = findfirst(x->x=='c', relevantPart) + 1
		t2 = findfirst(x->x=='v', relevantPart) - 1
		c = parse(Float64, relevantPart[t1:t2])
		d["cfe"] = c
	end

	if vmic
		relevantPart = parts[2]
		t1 = findfirst(x->x=='v', relevantPart) + 1
		t2 = length(relevantPart)
		c = parse(Float64, relevantPart[t1:t2])

		d["vmic"] = c
	end

	d
end

# ╔═╡ 48820690-b866-4971-8235-2b44c99f4be2
allParasName = parametersFromName.(allRuns)

# ╔═╡ 3d1bee96-353d-4b02-83f6-713c8b074228
teffName = [p["teff"] for p in allParasName]

# ╔═╡ 6de6198b-cc08-4fef-9c5d-61b449d96b69
loggName = [p["logg"] for p in allParasName]

# ╔═╡ 38b30580-172c-4823-b320-01c00752f735
fehName = [p["feh"] for p in allParasName]

# ╔═╡ b5846734-9647-4df0-8776-227a5fda9a80
function uniqueGrid(arrs...)
	arr1 = first(arrs)
	arrn_unique = [[] for _ in eachindex(arrs)]
	for i in eachindex(arr1)
		if !(arr1[i] in arrn_unique[1])
			for j in eachindex(arrs)
				append!(arrn_unique[j], [arrs[j][i]])
			end
		end
	end

	tuple(arrn_unique...)
end

# ╔═╡ 1e8829d9-240d-4e94-a818-d5d0dcd61b3b
function ifindTeffLogg(teff, logg, teffarr, loggarr)
	found = []
	for i in eachindex(teffarr)
		if (teff ≈ teffarr[i]) & (logg ≈ loggarr[i])
			append!(found, [i])
		end
	end

	found
end

# ╔═╡ 05783f0c-b968-47f7-9886-bcc36901a9b1
function ifindTeffLoggFeh(teff, logg, feh, teffarr, loggarr, feharr)
	found = []
	for i in eachindex(teffarr)
		if (teff ≈ teffarr[i]) & (logg ≈ loggarr[i]) & (feh ≈ feharr[i])
			append!(found, [i])
		end
	end

	found
end

# ╔═╡ 8ec1055c-7ed0-449b-b235-d0a3e8d57db9
iNameParas(allParasName, teff, logg, feh) = begin
	m = [
		(p["teff"] == teff) & 
		(p["logg"] == logg) & 
		(p["feh"] == feh) 
		for p in allParasName
	]

	findfirst(m)
end

# ╔═╡ 655c35f2-66f1-49e8-a26b-45c2bdcb44eb
teffU, loggU, fehU = uniqueGrid(teffName, loggName, fehName)

# ╔═╡ 135ebd1a-46c9-4f7f-bf8c-450a89817c81


# ╔═╡ 5d7cc279-f5db-4f55-b385-53e7a23a4808
md"Define all available models as success for the time being."

# ╔═╡ 317e9c08-7c17-4d34-adb7-88d44f410f77
sucess = [length(ifindTeffLoggFeh(t, l, f, teffName, loggName, fehName)) > 0 for (t, l, f) in zip(grid.info[!, :teff], grid.info[!, :logg], grid.info[!, :feh])]

# ╔═╡ 369fa2bc-75e0-4d44-bfad-d9ceda6ac711
grid.info[!, :completed] = sucess

# ╔═╡ f4d2ecf8-e5b0-4e3f-850c-f7873312a500


# ╔═╡ e2798255-e2a6-4dff-8af4-4a42cd5043b6
count(sucess)

# ╔═╡ 98bad81d-99cd-4a9a-b296-9952af958a4c


# ╔═╡ 735a2b14-f27b-4bf7-a96d-d0cbec8370f8
begin
	CSV.write(joinpath(gridfolder, "stellar_parameters_$(join(model_extension,"_")).csv"), select(grid.info[sucess, :], [:teff,:logg,:feh]))
end

# ╔═╡ 4f247889-b73d-4efc-a08c-ae106d9e4939


# ╔═╡ 5be1d639-2ec6-4152-9977-484597e52b77
md"# Overview"

# ╔═╡ c83842c9-7e03-49a6-b768-5e2e47101294
PLATO_targets = DataFrame(FITS("gaiadr3_tefflogg_matched.fits")[2])

# ╔═╡ 0de2ee04-d5ed-41ff-b7a8-284a59f85c3b
md"## HR-diagram"

# ╔═╡ 49d7e250-b9fc-42ff-8e31-9155f8bf5fe2
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
	ax.hist2d(
		PLATO_targets[!, :teff_gspphot], PLATO_targets[!, :logg_gspphot],
		bins=300,
		cmin=1, cmap="Greys", norm=matplotlib.colors.LogNorm(),
		rasterized=true, alpha=0.95
	)

	failed = .!grid.info[!, :completed]
	solar = grid.info[!, :feh] .== 0.0
	m1 = grid.info[!, :feh] .== -1.0
	@info "Following metallicites have been found: $(sort(unique(grid.info[!, :feh])))"
	ax.scatter(
		grid.info[.!failed .& solar, :teff], grid.info[.!failed .& solar, :logg], marker="s", color="steelblue", s=63, edgecolor="k", alpha=1, label=L"\rm [Fe/H]=0"
	)
	ax.scatter(
		grid.info[.!failed .& m1, :teff], grid.info[.!failed .& m1, :logg], marker="s", color="cyan", s=63, edgecolor="k", alpha=1, label=L"\rm [Fe/H]=-1"
	)
	ax.scatter(
		grid.info[failed, :teff], grid.info[failed, :logg], marker="s", color="white", s=63, edgecolor="k", alpha=1, label=L"\rm in\ progress"
	)

	@show count(.!failed .& solar) 
	@show count(.!failed .& m1) 
	@show count(failed)
	@show count(.!failed)
	@show count(solar)
	@show count(m1)
	@show nrow(grid.info)
	
	
	#im = ax.scatter(teffName, loggName, color="k", marker="s")
	#f.colorbar(im, ax=ax)

	ax.set_xlim(7200, 4050)
	ax.set_ylim(5.1, 2.3)
	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")

	ax.legend()

	f.savefig("model_progress.pdf")
	
	gcf()
end

# ╔═╡ a6709d15-9aa5-4717-a346-c0b6e5f022fc


# ╔═╡ 52547a03-766e-49a9-9451-87dd6ccedf94
md"# Comparison between models"

# ╔═╡ 8a44f40b-691b-4dcb-8667-10c13e55625e
paraUStrings = ["Teff $t, logg $l, [Fe/H] $f" for (t, l, f) in zip(teffU, loggU, fehU)];

# ╔═╡ 6f8bf69b-4803-4acc-959c-9c1eb9969373
md"""
	Pick model parameters

Model 1:
$(@bind pickedParameters1 Select(paraUStrings))

Model 2:
$(@bind pickedParameters2 Select(paraUStrings))

Model 3:
$(@bind pickedParameters3 Select(paraUStrings))

Model 4:
$(@bind pickedParameters4 Select(paraUStrings))
"""

# ╔═╡ 08ab61f7-734b-4791-8a58-d63566f9f49a
begin
	iselect1 = findfirst(paraUStrings .== pickedParameters1)
	iselect2 = findfirst(paraUStrings .== pickedParameters2)
	iselect3 = findfirst(paraUStrings .== pickedParameters3)
	iselect4 = findfirst(paraUStrings .== pickedParameters4)
end;

# ╔═╡ 744e44aa-289e-4b70-b416-eee6144cb40a
begin
	@info "Loading Monitoring..."
	#monitoring1 = MUST.reload!(MUST.WatchDog)
	iNameSelect1 = iNameParas(
		allParasName, teffU[iselect1], loggU[iselect1], fehU[iselect1]
	)
	iNameSelect2 = iNameParas(
		allParasName, teffU[iselect2], loggU[iselect2], fehU[iselect2]
	)
	iNameSelect3 = iNameParas(
		allParasName, teffU[iselect3], loggU[iselect3], fehU[iselect3]
	)
	iNameSelect4 = iNameParas(
		allParasName, teffU[iselect4], loggU[iselect4], fehU[iselect4]
	)
	monitoring1 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect1], folder=gridfolder), lastN=1
	)
	monitoring2 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect2], folder=gridfolder), lastN=1
	)
	monitoring3 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect3], folder=gridfolder), lastN=1
	)
	monitoring4 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect4], folder=gridfolder), lastN=1
	)

	@info "Monitoring loaded."
end

# ╔═╡ e1df03b2-93f8-417f-9bc7-f3e4f03416ec


# ╔═╡ d865e20e-43b3-41ee-9daa-98a2cd373035
name_from_index(iselect) =  L"\rm T_{eff}="*"$(round(Int, teffU[iselect])) K, "*L"\rm log(g) = "*"$(round(loggU[iselect], sigdigits=2)), "*L"\rm [Fe/H] = "*"$(fehU[iselect])"

# ╔═╡ c680a448-b739-4d17-9e03-5c5ccbe4096f
name_from_index_short(iselect) =  "$(round(Int, teffU[iselect])), $(round(loggU[iselect], sigdigits=2))"

# ╔═╡ e0efd865-da71-4947-9646-624f11a91d00


# ╔═╡ 5809a890-5a94-4b1f-9c18-405878843e04
md"## Select Spectra"

# ╔═╡ c9b2fd50-6290-4404-bf98-15216f0fa220
listSpectra(inameselect) = split.(MUST.glob("spectra*/", joinpath(gridfolder, allRuns[inameselect])), "/", keepempty=false) .|> last

# ╔═╡ 8de14495-2de8-4dda-8ebd-f9248b44cc16
md"### Plotting"

# ╔═╡ 5d73b8cf-10d4-4f58-96f3-61ba082451ab
md"""
	Pick Spectra to plot

Model 1:
$(@bind pickedSpectra1 confirm(MultiCheckBox(listSpectra(iNameSelect1))))

Model 2:
$(@bind pickedSpectra2 confirm(MultiCheckBox(listSpectra(iNameSelect2))))

Model 3:
$(@bind pickedSpectra3 confirm(MultiCheckBox(listSpectra(iNameSelect3))))

Model 4:
$(@bind pickedSpectra4 confirm(MultiCheckBox(listSpectra(iNameSelect4))))
"""

# ╔═╡ 5a404a60-c52a-4d69-90b7-95d5854987a9


# ╔═╡ 4a81eeb0-f487-456a-b206-a986b5b85139
md"### Abundance corrections"

# ╔═╡ 7ddfee1a-04f9-43e7-800d-0571c0945676
md"""
	Pick Spectra for curve of growth

Model 1:
$(@bind pickedSpectraCog1 confirm(MultiCheckBox(listSpectra(iNameSelect1))))

Model 2:
$(@bind pickedSpectraCog2 confirm(MultiCheckBox(listSpectra(iNameSelect2))))

Model 3:
$(@bind pickedSpectraCog3 confirm(MultiCheckBox(listSpectra(iNameSelect3))))

Model 4:
$(@bind pickedSpectraCog4 confirm(MultiCheckBox(listSpectra(iNameSelect4))))
"""

# ╔═╡ ccbdbc43-1178-47be-8b48-71c33d8a4238


# ╔═╡ 33b58c3f-312d-4cf1-b48d-15cc53c7d81e
loadSpectra(imodelname, name) = begin
	path = abspath(joinpath(gridfolder, allRuns[imodelname], name))
	spec =  MUST.M3DISRun(path)
end

# ╔═╡ 9941e224-8481-47f4-96d7-5e20d74c7c29
begin
	spectra1 = [loadSpectra(iNameSelect1, n) for n in  pickedSpectra1]
	spectra2 = [loadSpectra(iNameSelect2, n) for n in  pickedSpectra2]
	spectra3 = [loadSpectra(iNameSelect3, n) for n in  pickedSpectra3]
	spectra4 = [loadSpectra(iNameSelect4, n) for n in  pickedSpectra4]

	cog1 = [loadSpectra(iNameSelect1, n) for n in  pickedSpectraCog1]
	cog2 = [loadSpectra(iNameSelect2, n) for n in  pickedSpectraCog2]
	cog3 = [loadSpectra(iNameSelect3, n) for n in  pickedSpectraCog3]
	cog4 = [loadSpectra(iNameSelect4, n) for n in  pickedSpectraCog4]
end;

# ╔═╡ 5503b352-e92e-48fe-9812-a89049c67adb
md"## Curve of Growth"

# ╔═╡ f475c9c0-3e0a-4fda-a865-7c27fd79416e
let
	f, ax = plt.subplots(
		3, 1, figsize=(10, 5), sharex=true
	)
	
	plt.subplots_adjust(hspace=0.03)

	is_1d(name) = .!occursin("m3dis", name)
	wvl(line) = floor(Int, round(MUST.pyconvert(Float64, line.lam0), sigdigits=4))
	ilines = [2, 3, 4]

	if count(.!is_1d.(pickedSpectraCog1))>0
		for line in ilines
			# first model
			lte_runs = cog1[is_1d.(pickedSpectraCog1)]
			lte_lines = [run.line[line] for run in lte_runs]
		    lte_abund = MUST.getabundances.(lte_runs)
			mask = sortperm(lte_abund)
		
			nlte_run = cog1[.!is_1d.(pickedSpectraCog1)][1]
			nlte_line = nlte_run.line[line]
		    nlte_abund = MUST.getabundances(nlte_run)
		
			fLTE = MUST.equivalentwidth(lte_lines, lte_abund, LTE=true)
			fNLTE = MUST.equivalentwidth(nlte_line; LTE=false)
			aLTE = MUST.abundance(lte_lines, lte_abund, LTE=true)
			ΔNLTE = try
				nlte_abund - aLTE(fNLTE)
			catch
				@warn "problem with $(line), model 1"
				-99
			end
			
			ax[0].plot(lte_abund[mask], fLTE.(lte_abund[mask]), marker="s", label="$(wvl(lte_lines[1])) "*L"\rm \AA,\ "*L"\rm \Delta\ NLTE="*"$(round(ΔNLTE, sigdigits=2))")
			ax[0].plot([nlte_abund], [fNLTE], marker="X")
		end
	end

	if count(.!is_1d.(pickedSpectraCog2))>0
		for line in ilines
			# first model
			lte_runs = cog2[is_1d.(pickedSpectraCog2)]
			lte_lines = [run.line[line] for run in lte_runs]
		    lte_abund = MUST.getabundances.(lte_runs)
			mask = sortperm(lte_abund)
		
			nlte_run = cog2[.!is_1d.(pickedSpectraCog2)][1]
			nlte_line = nlte_run.line[line]
		    nlte_abund = MUST.getabundances(nlte_run)
		
			fLTE = MUST.equivalentwidth(lte_lines, lte_abund, LTE=true)
			fNLTE = MUST.equivalentwidth(nlte_line; LTE=false)
			aLTE = MUST.abundance(lte_lines, lte_abund, LTE=true)
			ΔNLTE = try
				nlte_abund - aLTE(fNLTE)
			catch
				@warn "problem with $(line), model 2"
				-99
			end
			
			ax[1].plot(lte_abund[mask], fLTE.(lte_abund[mask]), marker="s", label="$(wvl(lte_lines[1])) "*L"\rm \AA,\ "*L"\rm \Delta\ NLTE="*"$(round(ΔNLTE, sigdigits=2))")
			ax[1].plot([nlte_abund], [fNLTE], marker="X")
		end
	end

	if count(.!is_1d.(pickedSpectraCog3))>0
		for line in ilines
			# first model
			lte_runs = cog3[is_1d.(pickedSpectraCog3)]
			lte_lines = [run.line[line] for run in lte_runs]
		    lte_abund = MUST.getabundances.(lte_runs)
			mask = sortperm(lte_abund)
		
			nlte_run = cog3[.!is_1d.(pickedSpectraCog3)][1]
			nlte_line = nlte_run.line[line]
		    nlte_abund = MUST.getabundances(nlte_run)
		
			fLTE = MUST.equivalentwidth(lte_lines, lte_abund, LTE=true)
			fNLTE = MUST.equivalentwidth(nlte_line; LTE=false)
			aLTE = MUST.abundance(lte_lines, lte_abund, LTE=true)
			ΔNLTE = try
				nlte_abund - aLTE(fNLTE)
			catch
				@warn "problem with $(line), model 3"
				-99
			end
			
			ax[2].plot(lte_abund[mask], fLTE.(lte_abund[mask]), marker="s", label="$(wvl(lte_lines[1])) "*L"\rm \AA,\ "*L"\rm \Delta\ NLTE="*"$(round(ΔNLTE, sigdigits=2))")
			ax[2].plot([nlte_abund], [fNLTE], marker="X")
		end
	end

	
	ax[0].legend()
	ax[1].legend()
	ax[2].legend()

	f
end

# ╔═╡ 9e7219d6-3a89-4b25-bee1-531d020d895f


# ╔═╡ c64f79bb-2836-4991-aa16-f6238e63d69a
md"## Surface temperature"

# ╔═╡ 34e6cd4c-bdb6-4179-9034-b362744bf5bf
let
	f, ax = plt.subplots(
		2, 3, figsize=(10, 5), gridspec_kw=Dict("height_ratios"=> [2, 1])
	)
	
	plt.subplots_adjust(wspace=0.01, hspace=0.15)

	lw_spec=1.8

	xes = [
		monitoring1[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring2[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring3[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring4[end]["opticalSurfaces"]["x"] ./1e8
	]
	xlim = (minimum(minimum.(xes)), maximum(maximum.(xes)))
	#xlim = (-2.5, 2.5)
	
	# top left
	optSurf = monitoring1[end]["opticalSurfaces"]["Tplane"]
	x = monitoring1[end]["opticalSurfaces"]["x"] ./1e8
	y = monitoring1[end]["opticalSurfaces"]["y"] ./1e8
	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	i = ax[0,0].imshow(
		optSurf,
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		rasterized=true,
		aspect="auto"
	)
	name = name_from_index(iselect1)
	#ax[0,0].set_title(name)
	ax[0,0].text(
		0.95, 0.95, name_from_index_short(iselect1), 
		ha="right", va="top", transform=ax[0,0].transAxes,
		backgroundcolor="k", color="w"
	)
	#ax[0,0].set_xlabel(L"\rm X\ [Mm]")
	ax[0,0].set_ylabel(L"\rm Y\ [Mm]")

	# top right
	optSurf = monitoring2[end]["opticalSurfaces"]["Tplane"]
	x = monitoring2[end]["opticalSurfaces"]["x"] ./1e8
	y = monitoring2[end]["opticalSurfaces"]["y"] ./1e8
	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	i = ax[0,1].imshow(
		optSurf,
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		rasterized=true,
		aspect="auto"
	)
	name = name_from_index(iselect2)
	#ax[0,1].set_title(name)
	ax[0,1].text(
		0.95, 0.95, name_from_index_short(iselect2), 
		ha="right", va="top", transform=ax[0,1].transAxes,
		backgroundcolor="k", color="w"
	)
	#ax[0,1].set_xlabel(L"\rm X\ [Mm]")
	ax[0,1].set(yticklabels=[])  # remove the tick labels
	#ax[1].tick_params(left=false)  # remove the ticks
	

	# bottom left
	optSurf = monitoring3[end]["opticalSurfaces"]["Tplane"]
	x = monitoring3[end]["opticalSurfaces"]["x"] ./1e8
	y = monitoring3[end]["opticalSurfaces"]["y"] ./1e8
	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	i = ax[0,2].imshow(
		optSurf,
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		rasterized=true,
		aspect="auto"
	)
	name = name_from_index(iselect3)
	#ax[1,0].set_title(name)
	ax[0,2].text(
		0.95, 0.95, name_from_index_short(iselect3), 
		ha="right", va="top", transform=ax[0,2].transAxes,
		backgroundcolor="k", color="w"
	)
	#ax[0,2].set_xlabel(L"\rm X\ [Mm]")
	#ax[2].set_ylabel(L"\rm Y\ [Mm]")
	ax[0,2].set(yticklabels=[])  # remove the tick labels
	#ax[2].tick_params(left=false)  # remove the ticks



	# spectra
	for (i, spectrum) in enumerate(spectra1)
		λ, F = if occursin("m3dis", pickedSpectra1[i])
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			F_NLTE = MUST.pyconvert(Array, spectrum.run.flux)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,0].plot(λ, F_NLTE ./ c, label="3D NLTE", color="tomato", lw=lw_spec)
			λ, F_NLTE ./ c
		else
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,0].plot(λ, F_LTE ./ c, label="1D LTE", color="k", ls="--", lw=0.6*lw_spec)
			λ, F_LTE ./ c
		end
		n = join(strip.(split(name_from_index_short(iselect1), ',')), '_')
		open("$(pickedSpectra1[i])_$(n).txt", "w") do f
			writedlm(f, [λ F], '\t')
		end
	end

	for (i, spectrum) in enumerate(spectra2)
		λ, F = if occursin("m3dis", pickedSpectra2[i])
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			F_NLTE = MUST.pyconvert(Array, spectrum.run.flux)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,1].plot(λ, F_NLTE ./ c, label="3D NLTE", color="tomato", lw=lw_spec)
			λ, F_NLTE ./ c
		else
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,1].plot(λ, F_LTE ./ c, label="1D LTE", color="k", ls="--", lw=0.6*lw_spec)
			λ, F_LTE ./ c
		end
		n = join(strip.(split(name_from_index_short(iselect2), ',')), '_')
		open("$(pickedSpectra2[i])_$n.txt", "w") do f
			writedlm(f, [λ F], '\t')
		end
	end

	for (i, spectrum) in enumerate(spectra3)
		λ, F = if occursin("m3dis", pickedSpectra3[i])
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			F_NLTE = MUST.pyconvert(Array, spectrum.run.flux)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,2].plot(λ, F_NLTE ./ c, label="3D NLTE", color="tomato", lw=lw_spec)
			λ, F_NLTE ./ c			
		else
			λ = MUST.pyconvert(Array, spectrum.run.lam)
			F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
			c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
			ax[1,2].plot(λ, F_LTE ./ c, label="1D LTE", color="k", ls="--", lw=0.6*lw_spec)
			λ, F_LTE ./ c
		end
		n = join(strip.(split(name_from_index_short(iselect3), ',')), '_')
		open("$(pickedSpectra3[i])_$(n).txt", "w") do f
			writedlm(f, [λ F], '\t')
		end
	end

	# corrections Na
	#=ax[1,0].text(
		0.05, 0.93, L"\rm \Delta NLTE="*" $(0.28)", 
		ha="left", va="top", transform=ax[1,0].transAxes,
		color="k", fontsize="large"
	)
	ax[1,1].text(
		0.95, 0.05, L"\rm \Delta NLTE="*" $(-0.27)", 
		ha="right", va="bottom", transform=ax[1,1].transAxes,
		color="k", fontsize="large"
	)
	ax[1,2].text(
		0.95, 0.2, L"\rm \Delta NLTE="*" $(-0.2)", 
		ha="right", va="bottom", transform=ax[1,2].transAxes,
		color="k", fontsize="large"
	)=#

	# corrections Mg
	#=ax[1,0].text(
		0.05, 0.93, L"\rm \Delta NLTE="*" $(0.36)", 
		ha="left", va="top", transform=ax[1,0].transAxes,
		color="k", fontsize="large"
	)
	ax[1,1].text(
		0.98, 0.05, L"\rm \Delta NLTE="*" $(-0.38)", 
		ha="right", va="bottom", transform=ax[1,1].transAxes,
		color="k", fontsize="large"
	)
	ax[1,2].text(
		0.98, 0.2, L"\rm \Delta NLTE="*" $(-0.17)", 
		ha="right", va="bottom", transform=ax[1,2].transAxes,
		color="k", fontsize="large"
	)=#


	
	ax[1,2].legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.15))
	#ax[1,0].set_xlim(5888.2,5892.2)
	ax[1,0].set_xlim(5181.8,5187.2)
	#ax[1,0].set_xlim(5172.6-1.2,5172.6+1.5)
	ax[1,0].set_ylim(-0.08, 1.08)
	ax[1,1].set_xlim(ax[1,0].get_xlim()...)
	ax[1,2].set_xlim(ax[1,0].get_xlim()...)
	ax[1,1].set_ylim(ax[1,0].get_ylim()...)
	ax[1,2].set_ylim(ax[1,0].get_ylim()...)
	ax[1,1].set(yticklabels=[])  # remove the tick labels
	ax[1,2].set(yticklabels=[])  # remove the tick labels
	#ax[1,0].set(xticklabels=[])  # remove the tick labels
	#ax[1,1].set(xticklabels=[])  # remove the tick labels
	#ax[1,2].set(xticklabels=[])  # remove the tick labels
	ax[1,0].xaxis.set_major_locator(plt.MaxNLocator(4))
	ax[1,1].xaxis.set_major_locator(plt.MaxNLocator(4))
	ax[1,2].xaxis.set_major_locator(plt.MaxNLocator(4))
	ax[1,0].set_xlabel(L"\rm \lambda\ [\AA]")
	ax[1,1].set_xlabel(L"\rm \lambda\ [\AA]")
	ax[1,2].set_xlabel(L"\rm \lambda\ [\AA]")



	

	f.savefig("surface_temperature_mg_zoom2.pdf")

	f
end

# ╔═╡ a813bcff-3b68-48b3-a55e-d57964ecd2e0
let
	lw_spec=1.8

	xes = [
		monitoring1[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring2[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring3[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring4[end]["opticalSurfaces"]["x"] ./1e8
	]
	xlim = (minimum(minimum.(xes)), maximum(maximum.(xes)))
	#xlim = (-2.5, 2.5)
	
	# top left

	mons = [monitoring1, monitoring2, monitoring3]
	selects = [iselect1, iselect2, iselect3]
	f = nothing
	for imodel in 1:3
		plt.close()
		f, ax = plt.subplots(
			1, 1, figsize=(4, 4)
		)
		plt.subplots_adjust(wspace=0.01, hspace=0.15)
		optSurf = mons[imodel][end]["opticalSurfaces"]["Tplane"]
		x = mons[imodel][end]["opticalSurfaces"]["x"] ./1e8
		y = mons[imodel][end]["opticalSurfaces"]["y"] ./1e8
		extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		i = ax.imshow(
			optSurf,
			origin="lower",
			extent=extent,
			cmap="gist_heat",
			rasterized=true,
			aspect="equal"
		)
		name = name_from_index(selects[imodel])
		ax.text(
			0.95, 0.95, name_from_index_short(selects[imodel]), 
			ha="right", va="top", transform=ax.transAxes,
			backgroundcolor="k", color="w", fontsize="x-large"
		)
		ax.set_xlim(extent[1:2]...)
		ax.set_ylim(extent[1:2]...)
		ax.set_xlabel(L"\rm X\ [Mm]", fontsize="x-large")
		ax.set_ylabel(L"\rm Y\ [Mm]", fontsize="x-large")
		ax.tick_params(axis="both", labelsize=12)
		ax.xaxis.set_major_locator(plt.MaxNLocator(6))
		ax.yaxis.set_major_locator(plt.MaxNLocator(6))

		n = join(strip.(split(name_from_index_short(selects[imodel]), ',')), '_')
		f.savefig("surface_temperature_$(n).pdf")

		#=
		# spectra
		for (i, spectrum) in enumerate(spectra1)
			λ, F = if occursin("m3dis", pickedSpectra1[i])
				λ = MUST.pyconvert(Array, spectrum.run.lam)
				F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
				F_NLTE = MUST.pyconvert(Array, spectrum.run.flux)
				c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
				ax[1,0].plot(λ, F_NLTE ./ c, label="3D NLTE", color="tomato", lw=lw_spec)
				λ, F_NLTE ./ c
			else
				λ = MUST.pyconvert(Array, spectrum.run.lam)
				F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
				c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
				ax[1,0].plot(λ, F_LTE ./ c, label="1D LTE", color="k", ls="--", lw=0.6*lw_spec)
				λ, F_LTE ./ c
			end
			n = join(strip.(split(name_from_index_short(iselect1), ',')), '_')
			open("$(pickedSpectra1[i])_$(n).txt", "w") do f
				writedlm(f, [λ F], '\t')
			end
		end
	
		#ax.set_xlim(5888.2,5892.2)
		ax.set_xlim(5181.8,5187.2)
		#ax.set_xlim(5172.6-1.2,5172.6+1.5)
		ax.set_ylim(-0.08, 1.08)
		ax.xaxis.set_major_locator(plt.MaxNLocator(4))
		ax.set_xlabel(L"\rm \lambda\ [\AA]")
		=#
		f
	end

	f
end

# ╔═╡ 20ec41ae-e5d3-43ad-90ad-beb5baafb926
let
	lw_spec=2.8

	xes = [
		monitoring1[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring2[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring3[end]["opticalSurfaces"]["x"] ./1e8,
		monitoring4[end]["opticalSurfaces"]["x"] ./1e8
	]
	xlim = (minimum(minimum.(xes)), maximum(maximum.(xes)))
	#xlim = (-2.5, 2.5)
	
	# top left

	mons = [monitoring1, monitoring2, monitoring3]
	selects = [iselect1, iselect2, iselect3]
	spec = [spectra1, spectra2, spectra3]
	ps = [pickedSpectra1, pickedSpectra2, pickedSpectra3]
	f = nothing
	for imodel in 1:3
		plt.close()
		f, ax = plt.subplots(
			1, 1, figsize=(7, 3)
		)
		plt.subplots_adjust(wspace=0.01, hspace=0.15)
		
		ax.text(
			0.025, 0.05, name_from_index_short(selects[imodel]), 
			ha="left", va="bottom", transform=ax.transAxes,
			backgroundcolor="k", color="w", fontsize="x-large"
		)
		
		ax.tick_params(axis="both", labelsize=12)
		ax.xaxis.set_major_locator(plt.MaxNLocator(5))
		ax.yaxis.set_major_locator(plt.MaxNLocator(5))

		

		# spectra
		for (i, spectrum) in enumerate(spec[imodel])
			λ, F = if occursin("m3dis", ps[imodel][i])
				λ = MUST.pyconvert(Array, spectrum.run.lam)
				F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
				F_NLTE = MUST.pyconvert(Array, spectrum.run.flux)
				c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
				ax.plot(λ, F_NLTE ./ c, label="3D NLTE", color="tomato", lw=lw_spec)
				λ, F_NLTE ./ c
			else
				λ = MUST.pyconvert(Array, spectrum.run.lam)
				F_LTE = MUST.pyconvert(Array, spectrum.run.flux_lte)
				c = MUST.pyconvert(Array, spectrum.run.flux_cnt)
				ax.plot(λ, F_LTE ./ c, label="1D LTE", color="k", ls="--", lw=0.6*lw_spec)
				λ, F_LTE ./ c
			end
			
		end
	
		#ax.set_xlim(5888.2,5892.2)
		#ax.set_xlim(5880,5900)
		ax.set_xlim(5181.8,5187.2)
		#ax.set_xlim(5160.8,5200.2)
		#ax.set_xlim(5172.6-1.2,5172.6+1.5)
		ax.set_ylim(-0.08, 1.08)
		ax.set_xlabel(L"\rm \lambda\ [\AA]", fontsize="x-large")
		ax.set_ylabel(L"\rm flux\ [normalized]", fontsize="x-large")
		ax.legend(ncol=2, fontsize="large", loc="lower right", bbox_to_anchor=(1.02, -.08))
		
		n = join(strip.(split(name_from_index_short(selects[imodel]), ',')), '_')
		#f.savefig("spectrum_Mg_zoom_$(n).pdf")
		
		f
	end

	f
end

# ╔═╡ 8d19fa3a-e065-4263-9ce0-06cce2fd8afb


# ╔═╡ 65f2a089-9bdb-4b82-8537-d7b58364608c
md"## Full 3D structure"

# ╔═╡ c1c03966-366f-4edf-b57f-55d3bb0c65d3
"""
	plot_profile(s, q; ax, kwargs...)

Plot profile function of optical depth cube (500).
"""
plot_profile(s, q; ax=nothing, kwargs...) = begin
	x = MUST.pyconvert(Array, s.run.xx)
	y = MUST.pyconvert(Array, s.run.yy)
	z = MUST.pyconvert(Array, s.run.zz)
	data = Dict(
		:τ500 => exp10.(MUST.pyconvert(Array, s.run.ltau)),
		:T => MUST.pyconvert(Array, s.run.temp),
		:d => MUST.pyconvert(Array, s.run.rho),
	)

	xx, yy, zz = MUST.meshgrid(x, y, z)
	b = MUST.Box(xx, yy, zz, data, MUST.AtmosphericParameters())
	bt = MUST.height_scale_fast(b, :τ500, logspace=true)

	qs, logq = MUST.is_log(q)
	
	# Background distribution
	xt, yt = reshape(log10.(b[:τ500]), :), reshape(logq.(b[qs]), :)
	h, x, y = MUST.numpy.histogram2d(
		xt, yt, bins=150
	)
	if !isnothing(ax)
		ax.imshow(
			h.T, 
			origin="lower",
			interpolation = "bicubic", 
			extent=[minimum(x), maximum(x), minimum(y), maximum(y)],
			cmap="GnBu", norm=matplotlib.colors.LogNorm(vmin=1), aspect="auto",
			rasterized=true
		)
	end
	
	#ax.hist2d(reshape(log10.(b[:τ500]), :), reshape(b[:T], :), rasterized=true, cmap="YlGnBu", bins=300, norm=matplotlib.colors.LogNorm())

	(profile(MUST.mean, bt, :log10τ500, q)..., xt, yt)
end

# ╔═╡ 33ac4be9-57b3-48cb-8444-659a6edc0b0b
let
	specs = [spectra1, spectra2, spectra3]
	selects = [iselect1, iselect2, iselect3]
	ps = [pickedSpectra1, pickedSpectra2, pickedSpectra3]
	xlim = [-4.5, 2.2]
	ylim = [3000, 11000]
	f = nothing 
	for imodel in 1:3
		plt.close()
		f, ax = plt.subplots(1, 1, figsize=(6, 5))
		spec = specs[imodel]
		for (i, s) in enumerate(spec)
			n = join(strip.(split(name_from_index_short(selects[imodel]), ',')), '_')
			x, y = try
				x, y, xx, yy = plot_profile(s, :T, ax=ax)
				ax.plot(x, y, color="k")

				ylim = [minimum(yy[xx .> xlim[1]]), maximum(yy[xx .< xlim[2]])]

				x, y
			catch
				x = s.run.ltau
				y = s.run.temp
				ax.plot(x, y, color="tomato", lw=3.5, ls="--", label=L"\rm MARCS\ 1D")

				x, y
			end

			mn = ps[imodel][i]
			mn = split(mn[length("spectra_")+1:end], "_lam")[1]
			open("$(mn)_$(n).txt", "w") do f
				writedlm(f, [x y], '\t')
			end
		end
		

		ax.text(
			0.05, 0.95, name_from_index_short(selects[imodel]), 
			ha="left", va="top", transform=ax.transAxes,
			backgroundcolor="k", color="w", fontsize="x-large"
		)

		#ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
		#ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
		ax.set_xlim(xlim...)
		ax.set_ylim(ylim...)
		ax.tick_params(axis="both", labelsize=12)
		ax.set_ylabel(L"\rm temperature\ [K]", fontsize="x-large")
		ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]", fontsize="x-large")

		# where some data has already been plotted to ax
		handles, labels = ax.get_legend_handles_labels()
		
		# manually define a new patch 
		patch = matplotlib.patches.Patch(color="#60A6C9", label=L"\rm DISPATCH\ 3D")
		
		# handles is a list, so append manual patch
		handles.append(patch) 

		ax.legend(handles=handles, loc="center left", fontsize="large")

		n = join(strip.(split(name_from_index_short(selects[imodel]), ',')), '_')
		f.savefig("temperature3D_$(n).pdf")
		
	end

	f
end

# ╔═╡ eec1478a-d108-4adc-b67b-69a4855aa886


# ╔═╡ 42ee41d0-d723-4ce5-86bf-9ad727871082
md"## Average comparison"

# ╔═╡ c0ba2e0b-500d-4b3a-83c6-40f3fdd0be4c
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label=name_from_index(iselect1), color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label=name_from_index(iselect2), color="magenta")

	optAv = monitoring3[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label=name_from_index(iselect3), color="cyan")

	optAv = monitoring4[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label=name_from_index(iselect4), color="lime")


	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()
	f
end

# ╔═╡ 6f827de1-7306-41de-ba7c-6b310672bf6e
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label=name_from_index(iselect1), color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label=name_from_index(iselect2), color="magenta")

	optAv = monitoring3[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label=name_from_index(iselect3), color="cyan")

	optAv = monitoring4[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label=name_from_index(iselect4), color="lime")

	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm \rho\ [g\ cm^{-3}]")
	
	ax.legend()
	f
end

# ╔═╡ 666789ee-eb92-4d21-85fe-bf1a85f56a12


# ╔═╡ Cell order:
# ╠═49ce5a16-bc56-11ef-0e57-1f9a55abaaf8
# ╟─ad04ebd0-6bb5-47d1-8dcb-c5d9b1346d91
# ╟─1a3653d8-74ae-4a35-a577-3994edc2adc5
# ╟─6160a040-944e-4e7d-8d4c-c12262145ad7
# ╟─1b9d3022-c1b3-44d0-a444-a27aa902d9f0
# ╟─4ec4a69e-0a15-4353-ac0c-bfb9e02ebdee
# ╠═a2e3351b-fd71-4231-b010-1f55aead47a0
# ╠═9ed0a9ef-8f97-4aeb-9e52-c05ae6ec59c9
# ╟─ac701824-7442-48ff-85bb-72315b2ad23f
# ╟─8c1bcdae-0d29-436f-9e81-4c7f8c763dcb
# ╟─865e1497-bc8f-403b-be99-c30d59c2f493
# ╟─7d853824-7b07-4dc9-8a00-ba026efb6c1b
# ╟─48820690-b866-4971-8235-2b44c99f4be2
# ╟─3d1bee96-353d-4b02-83f6-713c8b074228
# ╟─6de6198b-cc08-4fef-9c5d-61b449d96b69
# ╟─38b30580-172c-4823-b320-01c00752f735
# ╟─b5846734-9647-4df0-8776-227a5fda9a80
# ╟─1e8829d9-240d-4e94-a818-d5d0dcd61b3b
# ╟─05783f0c-b968-47f7-9886-bcc36901a9b1
# ╟─8ec1055c-7ed0-449b-b235-d0a3e8d57db9
# ╠═655c35f2-66f1-49e8-a26b-45c2bdcb44eb
# ╟─135ebd1a-46c9-4f7f-bf8c-450a89817c81
# ╟─5d7cc279-f5db-4f55-b385-53e7a23a4808
# ╠═317e9c08-7c17-4d34-adb7-88d44f410f77
# ╠═369fa2bc-75e0-4d44-bfad-d9ceda6ac711
# ╟─f4d2ecf8-e5b0-4e3f-850c-f7873312a500
# ╠═e2798255-e2a6-4dff-8af4-4a42cd5043b6
# ╟─98bad81d-99cd-4a9a-b296-9952af958a4c
# ╠═735a2b14-f27b-4bf7-a96d-d0cbec8370f8
# ╟─4f247889-b73d-4efc-a08c-ae106d9e4939
# ╟─5be1d639-2ec6-4152-9977-484597e52b77
# ╠═c83842c9-7e03-49a6-b768-5e2e47101294
# ╟─0de2ee04-d5ed-41ff-b7a8-284a59f85c3b
# ╟─49d7e250-b9fc-42ff-8e31-9155f8bf5fe2
# ╟─a6709d15-9aa5-4717-a346-c0b6e5f022fc
# ╟─52547a03-766e-49a9-9451-87dd6ccedf94
# ╟─8a44f40b-691b-4dcb-8667-10c13e55625e
# ╟─6f8bf69b-4803-4acc-959c-9c1eb9969373
# ╟─08ab61f7-734b-4791-8a58-d63566f9f49a
# ╟─744e44aa-289e-4b70-b416-eee6144cb40a
# ╟─e1df03b2-93f8-417f-9bc7-f3e4f03416ec
# ╟─d865e20e-43b3-41ee-9daa-98a2cd373035
# ╟─c680a448-b739-4d17-9e03-5c5ccbe4096f
# ╟─e0efd865-da71-4947-9646-624f11a91d00
# ╟─5809a890-5a94-4b1f-9c18-405878843e04
# ╟─c9b2fd50-6290-4404-bf98-15216f0fa220
# ╟─8de14495-2de8-4dda-8ebd-f9248b44cc16
# ╟─5d73b8cf-10d4-4f58-96f3-61ba082451ab
# ╟─5a404a60-c52a-4d69-90b7-95d5854987a9
# ╟─4a81eeb0-f487-456a-b206-a986b5b85139
# ╟─7ddfee1a-04f9-43e7-800d-0571c0945676
# ╟─ccbdbc43-1178-47be-8b48-71c33d8a4238
# ╟─33b58c3f-312d-4cf1-b48d-15cc53c7d81e
# ╟─9941e224-8481-47f4-96d7-5e20d74c7c29
# ╟─5503b352-e92e-48fe-9812-a89049c67adb
# ╟─f475c9c0-3e0a-4fda-a865-7c27fd79416e
# ╟─9e7219d6-3a89-4b25-bee1-531d020d895f
# ╟─c64f79bb-2836-4991-aa16-f6238e63d69a
# ╟─34e6cd4c-bdb6-4179-9034-b362744bf5bf
# ╟─a813bcff-3b68-48b3-a55e-d57964ecd2e0
# ╟─20ec41ae-e5d3-43ad-90ad-beb5baafb926
# ╟─8d19fa3a-e065-4263-9ce0-06cce2fd8afb
# ╟─65f2a089-9bdb-4b82-8537-d7b58364608c
# ╟─c1c03966-366f-4edf-b57f-55d3bb0c65d3
# ╟─33ac4be9-57b3-48cb-8444-659a6edc0b0b
# ╟─eec1478a-d108-4adc-b67b-69a4855aa886
# ╟─42ee41d0-d723-4ce5-86bf-9ad727871082
# ╟─c0ba2e0b-500d-4b3a-83c6-40f3fdd0be4c
# ╟─6f827de1-7306-41de-ba7c-6b310672bf6e
# ╟─666789ee-eb92-4d21-85fe-bf1a85f56a12
