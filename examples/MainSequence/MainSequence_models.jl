### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ 30832d50-dc75-11ee-2e49-5f52349e5e50
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
end

# ╔═╡ cadc96a0-b6a8-4e46-8240-d1e72f976428
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	TableOfContents()
end

# ╔═╡ 6cb9c336-4f1d-43bd-b626-a7404a44d33b
visual = MUST.ingredients("visual.jl")

# ╔═╡ 6c5b6a85-0b7d-4afd-a068-38ce7c9ac52e
md"# Models
Models can be selected from a folder containing all models as sub-folders."

# ╔═╡ 1b3d7608-1eec-4d35-a19c-a95562828194
md"## Load Models"

# ╔═╡ e99d2039-091c-4749-a355-b5cf124c7e46
md"Pick directory of grid models: $(@bind gridfolder confirm(TextField(default=\"./packed_models/\")))"

# ╔═╡ 64d370db-ca48-44ef-9f6c-c676a45fbcd0
availableRuns(folder, inName...) = begin
	allruns = MUST.glob("*/", folder)
	n = split.(allruns, "/", keepempty=false) .|> last
	mask = falses(length(allruns))
	for inn in inName
		mask .= mask .| [occursin(inn, ni) for ni in n]
	end
	n[mask]
end

# ╔═╡ 8a760e31-7dc9-4c0e-8051-cd5d5b8f9a7e


# ╔═╡ 0f4989c9-2c2b-4733-8f35-324e74749b11
allRuns = availableRuns(gridfolder, "R5", "R6")

# ╔═╡ 6ce10bb3-1fbb-4d25-a1b9-cd25eea974de


# ╔═╡ ba036d7e-6f3f-4a4c-af22-3966e04fe74b
md"
## Stagger Grid
Load the initial Stagger grid that was used to create the models
"

# ╔═╡ e294057b-27c2-4f78-bdab-2a512ffc0d79
staggerGrid = MUST.StaggerGrid("../initial_models/stagger_grid_full_o.mgrid")

# ╔═╡ 68548a3f-5712-4ee9-99e0-6aad41624edc


# ╔═╡ 2bc84808-c748-49fb-9d27-691ddef1355f
md"## Model Parameters
for each model we can first of all extract the parameters from the label of the folder. We could alternatively also get those values from the model cubes themselves, or from the namelist that is included."

# ╔═╡ 1199b473-e011-4817-be65-0a685242cf2a
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

# ╔═╡ d749afd7-a6fa-43c1-9c6f-0a0cd2c0ec0b
allParasName = parametersFromName.(allRuns)

# ╔═╡ 0327e1f0-d924-4c7a-a1a3-296393de8c73
teffName = [p["teff"] for p in allParasName]

# ╔═╡ a043b27f-f392-4b96-b7b5-22cef4836a63
loggName = [p["logg"] for p in allParasName]

# ╔═╡ 9a8f1cf6-3577-4e6b-8bc7-354f5a13355d
fehName = [p["feh"] for p in allParasName]

# ╔═╡ 898ce0f5-7ec8-4f30-a4fe-96956ed1072f


# ╔═╡ 3b31e0f2-8acc-49dc-b69e-dd084bae4f79
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

# ╔═╡ c9cdfa8a-ee49-4e46-8e8e-4e7cd668b98a
function ifindTeffLogg(teff, logg, teffarr, loggarr)
	found = []
	for i in eachindex(teffarr)
		if (teff == teffarr[i]) & (logg == loggarr[i])
			append!(found, [i])
		end
	end

	found
end

# ╔═╡ 09505e68-b2d6-4c07-8e3f-a0e097ea3732
iNameParas(allParasName, teff, logg, feh) = begin
	m = [
		(p["teff"] == teff) & 
		(p["logg"] == logg) & 
		(p["feh"] == feh) 
		for p in allParasName
	]

	findfirst(m)
end

# ╔═╡ e28a5cc2-bd7b-4892-90a4-3ea34076d9af
teffU, loggU, fehU = uniqueGrid(teffName, loggName, fehName)

# ╔═╡ c3b2a375-9def-4e95-a10d-64c910d0c77c


# ╔═╡ 12b14982-f0c5-46fe-8c0e-379f4da4661e
md"# Overview"

# ╔═╡ c4b981a2-c759-4611-97a6-6aa731a6598b
PLATO_targets = DataFrame(FITS("gaiadr3_tefflogg_matched.fits")[2])

# ╔═╡ 81535166-89a3-4321-b620-f8e2c78d5ce4
md"## HR-diagram"

# ╔═╡ 3a0c3d77-af08-4018-825f-d055779a5911
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
	ax.hist2d(
		PLATO_targets[!, :teff_gspphot], PLATO_targets[!, :logg_gspphot],
		bins=700,
		cmin=1, cmap="Greys", norm=matplotlib.colors.LogNorm()
	)
	
	im = ax.scatter(teffName, loggName, color="k", marker="s")
	#f.colorbar(im, ax=ax)

	ax.set_xlim(7100, 4000)
	ax.set_ylim(5, 1)

	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	
	gcf()
end

# ╔═╡ 1111550d-3ec4-4974-9e72-6ce7bad32230


# ╔═╡ c1e576d3-a8f1-453e-9dd8-50c1d650e2e0
md"## Simulation Time"

# ╔═╡ 30533c62-c714-4685-bf26-c0e6d6948bfa
gettime(folder) = begin
	lastsnap = MUST.list_snapshots(MUST.converted_snapshots(folder)) |> last
	snap, _ = pick_snapshot(folder, lastsnap)
	snap.parameter.time
end

# ╔═╡ 6dc3764f-cfd5-475f-84f2-7dedb693f31e
lastsnapshot(folder) = begin
	MUST.list_snapshots(MUST.converted_snapshots(folder)) |> last
end

# ╔═╡ d0573093-6878-4d38-a2ef-0e6a54fdc878
simTime = gettime.(joinpath.(gridfolder, allRuns))

# ╔═╡ 3b2ac7f6-fe83-4465-810a-30147d785e00
lastSnap = lastsnapshot.(joinpath.(gridfolder, allRuns))

# ╔═╡ bfeafa54-beca-4ee3-88a7-48f384fe3df4


# ╔═╡ 85cd82ca-5417-4b86-80b8-e40d680b96c2
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	im = ax.scatter(
		loggName, simTime ./60 ./60, c=teffName, cmap="RdYlBu", s=70, edgecolor="k"
	)
	f.colorbar(im, ax=ax)

	ax.set_xlabel(L"\rm \log(g)\ [dex]")
	ax.set_ylabel(L"\rm stellar\ time\ [h]")
	
	gcf()
end

# ╔═╡ a392c80a-84c4-40f6-a4cd-8e090df26086
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	im = ax.scatter(
		loggName, lastSnap, c=teffName, cmap="RdYlBu", edgecolor="k", s=70
	)
	f.colorbar(im, ax=ax)

	ax.set_xlabel(L"\rm \log(g)\ [dex]")
	ax.set_ylabel(L"\rm last\ snap\ id")
	
	gcf()
end

# ╔═╡ de16e1d1-54bd-4ce5-8035-2a33389cdfe5


# ╔═╡ 4ca5a08d-2738-4e48-b84e-eb81250d615b
md"## Crashed Models"

# ╔═╡ 38273a36-ea28-4bf4-bb36-529d1a619511
crashedModels = lastSnap .< 170

# ╔═╡ 7c39ae82-d34a-4921-9b69-7f8736fcf5dc
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	im = ax.scatter(
		teffName[.!crashedModels], 
		loggName[.!crashedModels],
		c="black",
		marker="s"
	)
	im = ax.scatter(
		teffName[crashedModels], 
		loggName[crashedModels],
		c="red",
		marker="s", label="crashed"
	)
	
	ax.set_xlim(7100, 4000)
	ax.set_ylim(5, 1)
	ax.legend()

	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	
	
	gcf()
end

# ╔═╡ 6db48ee7-183e-441e-ac66-7a6156df8d9b
failedModels = DataFrame(
	Dict(
		"model"=>[r for (i, r) in enumerate(allRuns) if crashedModels[i]],
		"teff"=>[r for (i, r) in enumerate(teffName) if crashedModels[i]],
		"logg"=>[r for (i, r) in enumerate(loggName) if crashedModels[i]],
		"modelpath"=>[joinpath(gridfolder, r) for (i, r) in enumerate(allRuns) if crashedModels[i]]
	)
)

# ╔═╡ 0d77142a-d908-4ad9-b32f-a52f00b95b4b


# ╔═╡ 8950923d-e4d9-4517-83d5-95b6dbaa5b19
md"# Monitoring"

# ╔═╡ 095d5c29-d828-409c-afd9-2384d1bb3a83
begin
	# number of modes
	notCrashedRuns = allRuns[.!crashedModels]
	notCrashedteff = round.(Int, teffName[.!crashedModels])
	notCrashedlogg = round.(loggName[.!crashedModels], sigdigits=4)

	mask = sortperm((collect(zip(notCrashedteff, notCrashedlogg))))
	notCrashedRuns = notCrashedRuns[mask]
	notCrashedteff = notCrashedteff[mask]
	notCrashedlogg = notCrashedlogg[mask]

	teff = notCrashedteff
	logg = notCrashedlogg
end

# ╔═╡ 2dd0de7c-028b-43f2-960e-e0839d3fa7b7
md"## Optical Surface"

# ╔═╡ 720ff7c3-66d5-4dfc-89fb-4a385e7bfd9e
md"""
	Tick models to show optical suface
$(@bind selectForSurface MultiCheckBox(notCrashedRuns, default=sample(notCrashedRuns, 25)))
"""

# ╔═╡ c3b9f946-4fae-4aca-9722-a80dc999b95d
let
	npanles = length(selectForSurface)
	nx = ceil(Int, sqrt(npanles))
	ny = ceil(Int, npanles / nx)

	plt.close()
	f, ax = plt.subplots(nx, ny, figsize=(10, 10))

	k = 1
	for j in 0:ny-1
		for i in 0:nx-1
			if k>npanles
				ax[i][j].axis("off")
				continue
			end
			iselectForSurface = findfirst(notCrashedRuns .== selectForSurface[k])
			w = MUST.WatchDog(selectForSurface[k], folder=gridfolder)
			ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
			m = MUST.reload(w, last(ws))["opticalSurfaces"]
		
			x = m["x"] ./1e8
			y = m["y"] ./1e8
		
			extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		
			im = ax[i][j].imshow(
				m["Tplane"],
				origin="lower",
				extent=extent,
				cmap="gist_heat",
				aspect="auto"
			)

			ax[i][j].text(
				0.97, 0.97, 
				"$(teff[iselectForSurface]) K,"*" $(logg[iselectForSurface]) dex", 
				ha="right", va="top", 
				transform=ax[i][j].transAxes,
				fontsize="small",
				color="white", backgroundcolor="k"
			)
			
			#cb = f.colorbar(i, ax=ax[i][j], fraction=0.046, pad=0.04)
			k += 1
		end
	end
	#plt.set_aspect("equal")

	plt.subplots_adjust(wspace=0, hspace=0)
	f.supxlabel("X [Mm]")
	f.supylabel("Y [Mm]")

	f.savefig(joinpath(gridfolder, "optical_surface_T.pdf"), bbox_inches="tight")
	gcf()
end

# ╔═╡ 4e9f26c6-001c-4d90-b9b3-ef4f259e0ae9


# ╔═╡ 9b659650-b91d-4229-a811-6802f7d56966
which_gif = "Tplane" 

# ╔═╡ 6df1c2d5-850c-4af8-b29d-25d7a0f7666c
begin
	vminmax = []
	
	for k in 1:length(selectForSurface)
		w = MUST.WatchDog(selectForSurface[k], folder=gridfolder)
		ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
		m = MUST.reload(w, last(ws))["opticalSurfaces"]
		
		append!(
			vminmax, [(minimum(m[which_gif]), maximum(m[which_gif]))]
		)
	end

	vminmax
end

# ╔═╡ 115b47e6-5153-4940-8499-0b19de8f2e54
function optical_surface(notCrashedRuns, teff, logg, selectForSurface; which, vminmax=nothing, cmap="coolwarm", quantity="Tplane")
	npanles = length(selectForSurface)
	nx = ceil(Int, sqrt(npanles))
	ny = ceil(Int, npanles / nx)

	plt.close()
	f, ax = plt.subplots(nx, ny, figsize=(10, 10))

	k = 1
	for j in 0:ny-1
		for i in 0:nx-1
			if k>npanles
				ax[i][j].axis("off")
				continue
			end
			iselectForSurface = findfirst(notCrashedRuns .== selectForSurface[k])
			w = MUST.WatchDog(selectForSurface[k], folder=gridfolder)
			ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
			which = min(which, length(ws))
			m = MUST.reload(w, ws[which])["opticalSurfaces"]
		
			x = m["x"] ./1e8
			y = m["y"] ./1e8
		
			extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		
			im = if isnothing(vminmax) 
				ax[i][j].imshow(
					m[quantity],
					origin="lower",
					extent=extent,
					cmap=cmap,
					aspect="auto"
				)
			else
				ax[i][j].imshow(
					m[quantity],
					origin="lower",
					extent=extent,
					cmap=cmap,
					aspect="auto",
					vmin=vminmax[k][1],
					vmax=vminmax[k][2]
				)
			end

			ax[i][j].text(
				0.97, 0.97, 
				"$(teff[iselectForSurface]) K,"*" $(logg[iselectForSurface]) dex", 
				ha="right", va="top", 
				transform=ax[i][j].transAxes,
				fontsize="small",
				color="white", backgroundcolor="k"
			)
			
			#cb = f.colorbar(i, ax=ax[i][j], fraction=0.046, pad=0.04)
			k += 1
		end
	end
	#plt.set_aspect("equal")

	plt.subplots_adjust(wspace=0, hspace=0)
	f.supxlabel("X [Mm]")
	f.supylabel("Y [Mm]")

	f.savefig(
		joinpath(gridfolder, "optical_surface_$(quantity)_$(which).png"), bbox_inches="tight",
		dpi=72
	)

	joinpath(gridfolder, "optical_surface_$(quantity)_$(which).png")
end

# ╔═╡ 53892118-270c-40c7-a12b-99aef4d9e69b
let
	f = []
	for i in 1:360
		append!(f, [optical_surface(
			notCrashedRuns, notCrashedteff, notCrashedlogg, selectForSurface; 
			which=i, 
			vminmax=vminmax,
			quantity=which_gif,
			cmap="gist_heat"
		)])
	end

	visual.gifs.gifs_from_png(
		f, joinpath(gridfolder, "optical_surface_T.gif");
		duration=0.5
	)
end

# ╔═╡ 8ab32b10-f50a-4fdd-b0b6-dd0c1d154e1f


# ╔═╡ 4866dd18-e1f8-4939-83ec-f2e80407f731
md"## Mass Flux"

# ╔═╡ 73ea53c8-eb56-462c-8612-41c33d379915
begin
	runnames = split.(allRuns, "pack_", keepempty=false) .|> first
	runnames = [r*".nml" for r in runnames]
	nml = MUST.FreeNamelist.(joinpath.(gridfolder, allRuns, runnames))
	htop_scale = MUST.nmlValue.(nml, "htop_scale")
end

# ╔═╡ 44f8fb4d-0ffc-4e12-b0de-cf2d024a07f5
md"""
Click to compute mass flux: $(@bind computeMassFlux CheckBox(default=false))
"""

# ╔═╡ faba53f3-1a3b-4211-9205-928404525ed4
begin
	mass_flux = zeros(length(allRuns))
	if computeMassFlux
		for k in 1:length(allRuns)
			w = MUST.WatchDog(allRuns[k], folder=gridfolder)
			m = MUST.reload!(w, groups=["geoMassFlux", "atmosphericParameters"], lastN=30)
			mf = [mi["geoMassFlux"]["massFlux"][end] for mi in m] ./1e5
			t = [mi["atmosphericParameters"]["time"] for mi in m]
			mass_flux[k] = MUST.mean(mf) 
		end
	end
end

# ╔═╡ 97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
logg_vertical_line = 0.29

# ╔═╡ d51f056b-6ef3-430b-aacf-04c14caeb2ea
computeMassFlux && let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	vmin = minimum(mass_flux)
	vmax = maximum(mass_flux)
	#vtot = max(abs(vmin), abs(vmax))
	#vmin = -vtot
	#vmax = vtot
	
	i = ax.plot(teffName ./ exp10.(loggName), mass_flux, marker="s", ls="", color="k")

	ax.axvline(logg_vertical_line, color="red")
	ax.axhline(0.0, color="red")
	
	
	ax.set_xlabel(L"\rm T_{eff} / g")
	ax.set_ylabel(L"\rm mass\ flux\ [km/s]")
	
	f
end

# ╔═╡ 93e362e6-be03-4e9c-b11d-4002fc60211d
computeMassFlux && let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	vmin = minimum(mass_flux)
	vmax = maximum(mass_flux)
	#vtot = max(abs(vmin), abs(vmax))
	#vmin = -vtot
	#vmax = vtot

	b = 2.5
	c = 4.15
	d = 0.3
	a = 2.0 -d
	func(x) = a * exp(-b*(x-c)) + d

	b2 = 2.5
	c2 = 4.15
	d2 = 0.3
	a2 = 2.0 -d
	func2(x) = a2 * exp(-b2*(x-c2)) + d2

	x0 = 5777.0 / exp10(4.44)
	y0 = 1.1
	x1 = 0.29
	y1 = 1.5
	m = (y1 - y0) / (x1 - x0)
	func3(x) = m * x + (y1-m*x1)
	
	i = ax.plot(teffName ./ exp10.(loggName), htop_scale, marker="s", ls="", color="k", zorder=0)
	#=ax.scatter(
		[3.5, 4.15, 4.44, 5.0], 
		[7.5, 2.0, 1.1, 0.5],
		s=50, 
		color="red",
		label="usefull points"
	)=#
	ax.scatter(
		[5777.0 / exp10(4.44)], 
		[1.1],
		s=50, 
		color="red",
		label="usefull points"
	)
	ax.scatter(
		[0.29], 
		[1.52],
		s=50, 
		color="red",
		label="usefull points"
	)

	x = range(0, 1, length=100)
	#ax.plot(x, func.(x), color="magenta")
	#ax.plot(x, func2.(x), color="cyan")
	ax.plot(x, func3.(x), color="cyan")

	#ax.set_ylim(0, 3.0)
	#ax.set_xlim(3.7, 5.1)
	
	ax.legend()
	ax.axvline(logg_vertical_line, color="red")

	ax.set_xlabel(L"\rm T_{eff} / g")
	ax.set_ylabel(L"\rm top\ boundary\ scaling\ factor")
	
	f
end

# ╔═╡ 1799306b-954d-402e-8652-a66cf5f3389b


# ╔═╡ b2b24be9-6db7-46c5-a69c-3d790c4d5406
md"## Comparison between models"

# ╔═╡ ae86742d-6e97-4f32-9edb-e0feeb9f446d
paraUStrings = ["Teff $t, logg $l, [Fe/H] $f" for (t, l, f) in zip(teffU, loggU, fehU)];

# ╔═╡ 0eaa7e63-dee5-466c-8f9f-959cbab74f18
md"""
	Pick model parameters
$(@bind pickedParameters1 Select(paraUStrings))
$(@bind pickedParameters2 Select(paraUStrings))
"""

# ╔═╡ 7c27b82b-6438-4295-a001-503f2a200fbc
begin
	iselect1 = findfirst(paraUStrings .== pickedParameters1)
	iselect2 = findfirst(paraUStrings .== pickedParameters2)
end;

# ╔═╡ 835848e9-8221-40ad-8b8d-ff427660a499
begin
	@info "Loading Monitoring..."
	#monitoring1 = MUST.reload!(MUST.WatchDog)
	iNameSelect1 = iNameParas(
		allParasName, teffU[iselect1], loggU[iselect1], fehU[iselect1]
	)
	iNameSelect2 = iNameParas(
		allParasName, teffU[iselect2], loggU[iselect2], fehU[iselect2]
	)
	monitoring1 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect1], folder=gridfolder)
	)
	monitoring2 = MUST.reload!(
		MUST.WatchDog(allRuns[iNameSelect2], folder=gridfolder)
	)

	@info "Monitoring loaded."
end

# ╔═╡ 709ebdbe-e7a1-4471-a89c-d2dae3aecc36
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label="$(pickedParameters1)", color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label="$(pickedParameters2)", color="magenta")


	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()
	f
end

# ╔═╡ 995f2a3e-89ba-4afd-8845-deff05a651d2
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label="$(pickedParameters1)", color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label="$(pickedParameters2)", color="magenta")

	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm \rho\ [g\ cm^{-3}]")
	
	ax.legend()
	f
end

# ╔═╡ a901759c-addd-4485-94de-32794f4ae0f0
md"# Summary"

# ╔═╡ ab835496-0843-40c5-ba4a-37b4fd5f2e36
function add_point!(f, ax, teff, logg, segments, colors; square_logg_length=0.1, square_teff_length=0.1, alpha=0.3, add_label=false, title=nothing, zorder=1, labels=nothing)
	rect = plt.Rectangle(
		(teff - square_teff_length/2, logg - square_logg_length/2), square_teff_length, square_logg_length,
		edgecolor="black", 
		facecolor="None",
		lw=0,
		zorder=zorder
	)
	ax.add_patch(rect)
	ax.scatter(
		(teff - square_teff_length/2, logg - square_logg_length/2)..., s=0, zorder=zorder
	)

	# for each segment, fill in given color
	nseg = length(segments)
	for i in eachindex(segments)
		smin = logg - square_logg_length/2 + (i-1)*(square_logg_length/nseg)
		smax = smin + square_logg_length/nseg

		if isnothing(labels)
			ax.fill_between(
				[teff - square_teff_length/2, teff + square_teff_length/2],
				[smin, smin], [smax, smax],
				color=colors[i],
				edgecolor="black",
				alpha=alpha,
				zorder=zorder
			)
		else
			ax.fill_between(
				[teff - square_teff_length/2, teff + square_teff_length/2],
				[smin, smin], [smax, smax],
				color=colors[i],
				edgecolor="black",
				alpha=alpha,
				zorder=zorder,
				label=labels[i]
			)
		end

		if add_label
			ax.text(
				teff, smax, segments[i], va="bottom", ha="center", zorder=zorder
				
			)
		end
	end

	if !isnothing(title)
		ax.text(
			teff, logg - (square_logg_length/2), title, 
			va="bottom", ha="center", zorder=zorder
		)
	end
end

# ╔═╡ c1779c06-e40c-415c-b268-2fecb58fe087


# ╔═╡ cd1a8d0b-ea3a-43d5-b2be-c169b35e97ce
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	fehs = [0.0]

	dteff = abs(maximum(teffU)-minimum(teffU))
	dlogg = abs(maximum(loggU)-minimum(loggU))
	upperLogg = 5.1 #maximum(loggU)+0.1*dlogg
	lowerLogg = 2.8 #minimum(loggU)-0.1*dlogg
	upperTeff = 7500 #maximum(teffU)+0.1*dteff
	lowerTeff = 4000 #minimum(teffU)-0.1*dteff

	square_logg_length=0.028*(upperLogg-lowerLogg)
	square_teff_length=0.028*(upperTeff-lowerTeff)

	# PLATO targets
	ax.hist2d(
		PLATO_targets[!, :teff_gspphot], PLATO_targets[!, :logg_gspphot],
		bins=300, zorder=-1,
		cmin=1, cmap="Greys", norm=matplotlib.colors.LogNorm(),
	)

	for i in eachindex(teffU)
		iTeffLogg = ifindTeffLogg(teffU[i], loggU[i], teffName, loggName)
		segs = fehName[iTeffLogg]
		crashed = crashedModels[iTeffLogg]
		colors = [c ?  "magenta" : "lime" for c in crashed]
		zorder = any(crashed) ? 0 : 1
		colors = [
			fehs[i] in segs ? colors[findfirst(fehs[i].==segs)] : "w" 
			for i in eachindex(fehs)
		]
		#first(crashed) && continue
		add_point!(
			f, ax, teffU[i], loggU[i], fehs, colors, 
			square_logg_length=square_logg_length,
			square_teff_length=square_teff_length,
			alpha=0.8, zorder=zorder
		)
	end

	ax.set_ylim(upperLogg, lowerLogg)
	ax.set_xlim(upperTeff, lowerTeff)

	x_start = first(MUST.pyconvert(Any, ax.get_xlim()))
	y_start = first(MUST.pyconvert(Any, ax.get_ylim()))
	x_end = last(MUST.pyconvert(Any, ax.get_xlim()))
	y_end = last(MUST.pyconvert(Any, ax.get_ylim()))

	# add FeH legend
	#=ax_x_end = x_end + abs(x_end - x_start) *0.1
	ax_y_end = y_end + abs(y_end - y_start) *0.1
	ax_x_start = x_start - abs(x_end - x_start) *0.1
	ax_y_start = y_start - abs(y_end - y_start) *0.1
	add_point!(
		f, ax, ax_x_start, ax_y_start, fehs, ["w" for _ in fehs], 
		square_logg_length=square_logg_length*2.9,
		square_teff_length=square_teff_length*2.9,
		alpha=0.8, add_label=true, title="[Fe/H]"
	)=#

	# legend handles
	add_point!(
		f, ax, -100, -100, [1], ["lime"], 
		square_logg_length=square_logg_length*2.9,
		square_teff_length=square_teff_length*2.9,
		alpha=0.8, labels=["completed"]
	)
	add_point!(
		f, ax, -100, -100, [1], ["magenta"], 
		square_logg_length=square_logg_length*2.9,
		square_teff_length=square_teff_length*2.9,
		alpha=0.8, labels=["in progress"]
	)
	ax.plot([-100], [-100], color="0.5", label="PLATO targets", ls="", marker="s")


	# overplot the stagger grid
	maskStag = staggerGrid["feh"] .== 0.0
	ax.scatter(
		staggerGrid["teff"][maskStag], staggerGrid["logg"][maskStag], 
		color="k", marker="X", s=70, label="Stagger grid"
	)

	ax.legend(loc="upper left")
	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	f.savefig("DISPATCH_models_PLATO_status_overview.pdf")
	f.savefig("DISPATCH_models_PLATO_status_overview.png", dpi=300)
	
	gcf()
end

# ╔═╡ ff92f822-19b1-4c03-b29c-3227bd25de45
finalParameters = DataFrame(
	Dict(
		"model"=>[r for (i, r) in enumerate(allRuns) if !crashedModels[i]],
		"teff"=>[r for (i, r) in enumerate(teffName) if !crashedModels[i]],
		"logg"=>[r for (i, r) in enumerate(loggName) if !crashedModels[i]],
		"feh"=>[r for (i, r) in enumerate(fehName) if !crashedModels[i]]
	)
)

# ╔═╡ 7af40fdd-5e66-4802-9a11-bd5e9f50a1bf
CSV.write(joinpath(gridfolder, "stellar_parameters.csv"), finalParameters)

# ╔═╡ f85a322f-c1b7-44d1-b7f8-234c9ebce22b
CSV.write(joinpath(gridfolder, "failed_models.csv"), failedModels)

# ╔═╡ e04e4cdb-90c8-4514-a673-ef8ecb7e6507
open(joinpath(gridfolder, "failed_models.txt"), "w") do f
	write(f, join(failedModels[!, :modelpath], '\n'))
end

# ╔═╡ Cell order:
# ╠═30832d50-dc75-11ee-2e49-5f52349e5e50
# ╠═cadc96a0-b6a8-4e46-8240-d1e72f976428
# ╟─6cb9c336-4f1d-43bd-b626-a7404a44d33b
# ╟─6c5b6a85-0b7d-4afd-a068-38ce7c9ac52e
# ╟─1b3d7608-1eec-4d35-a19c-a95562828194
# ╟─e99d2039-091c-4749-a355-b5cf124c7e46
# ╟─64d370db-ca48-44ef-9f6c-c676a45fbcd0
# ╟─8a760e31-7dc9-4c0e-8051-cd5d5b8f9a7e
# ╠═0f4989c9-2c2b-4733-8f35-324e74749b11
# ╟─6ce10bb3-1fbb-4d25-a1b9-cd25eea974de
# ╟─ba036d7e-6f3f-4a4c-af22-3966e04fe74b
# ╟─e294057b-27c2-4f78-bdab-2a512ffc0d79
# ╟─68548a3f-5712-4ee9-99e0-6aad41624edc
# ╟─2bc84808-c748-49fb-9d27-691ddef1355f
# ╟─1199b473-e011-4817-be65-0a685242cf2a
# ╠═d749afd7-a6fa-43c1-9c6f-0a0cd2c0ec0b
# ╠═0327e1f0-d924-4c7a-a1a3-296393de8c73
# ╠═a043b27f-f392-4b96-b7b5-22cef4836a63
# ╠═9a8f1cf6-3577-4e6b-8bc7-354f5a13355d
# ╟─898ce0f5-7ec8-4f30-a4fe-96956ed1072f
# ╟─3b31e0f2-8acc-49dc-b69e-dd084bae4f79
# ╟─c9cdfa8a-ee49-4e46-8e8e-4e7cd668b98a
# ╟─09505e68-b2d6-4c07-8e3f-a0e097ea3732
# ╟─e28a5cc2-bd7b-4892-90a4-3ea34076d9af
# ╟─c3b2a375-9def-4e95-a10d-64c910d0c77c
# ╟─12b14982-f0c5-46fe-8c0e-379f4da4661e
# ╟─c4b981a2-c759-4611-97a6-6aa731a6598b
# ╟─81535166-89a3-4321-b620-f8e2c78d5ce4
# ╟─3a0c3d77-af08-4018-825f-d055779a5911
# ╟─1111550d-3ec4-4974-9e72-6ce7bad32230
# ╟─c1e576d3-a8f1-453e-9dd8-50c1d650e2e0
# ╟─30533c62-c714-4685-bf26-c0e6d6948bfa
# ╟─6dc3764f-cfd5-475f-84f2-7dedb693f31e
# ╟─d0573093-6878-4d38-a2ef-0e6a54fdc878
# ╟─3b2ac7f6-fe83-4465-810a-30147d785e00
# ╟─bfeafa54-beca-4ee3-88a7-48f384fe3df4
# ╟─85cd82ca-5417-4b86-80b8-e40d680b96c2
# ╟─a392c80a-84c4-40f6-a4cd-8e090df26086
# ╟─de16e1d1-54bd-4ce5-8035-2a33389cdfe5
# ╟─4ca5a08d-2738-4e48-b84e-eb81250d615b
# ╠═38273a36-ea28-4bf4-bb36-529d1a619511
# ╟─7c39ae82-d34a-4921-9b69-7f8736fcf5dc
# ╟─6db48ee7-183e-441e-ac66-7a6156df8d9b
# ╟─0d77142a-d908-4ad9-b32f-a52f00b95b4b
# ╟─8950923d-e4d9-4517-83d5-95b6dbaa5b19
# ╟─095d5c29-d828-409c-afd9-2384d1bb3a83
# ╟─2dd0de7c-028b-43f2-960e-e0839d3fa7b7
# ╟─720ff7c3-66d5-4dfc-89fb-4a385e7bfd9e
# ╟─c3b9f946-4fae-4aca-9722-a80dc999b95d
# ╟─4e9f26c6-001c-4d90-b9b3-ef4f259e0ae9
# ╠═9b659650-b91d-4229-a811-6802f7d56966
# ╠═6df1c2d5-850c-4af8-b29d-25d7a0f7666c
# ╟─115b47e6-5153-4940-8499-0b19de8f2e54
# ╠═53892118-270c-40c7-a12b-99aef4d9e69b
# ╟─8ab32b10-f50a-4fdd-b0b6-dd0c1d154e1f
# ╟─4866dd18-e1f8-4939-83ec-f2e80407f731
# ╟─73ea53c8-eb56-462c-8612-41c33d379915
# ╟─44f8fb4d-0ffc-4e12-b0de-cf2d024a07f5
# ╟─faba53f3-1a3b-4211-9205-928404525ed4
# ╠═97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
# ╟─d51f056b-6ef3-430b-aacf-04c14caeb2ea
# ╟─93e362e6-be03-4e9c-b11d-4002fc60211d
# ╟─1799306b-954d-402e-8652-a66cf5f3389b
# ╟─b2b24be9-6db7-46c5-a69c-3d790c4d5406
# ╟─ae86742d-6e97-4f32-9edb-e0feeb9f446d
# ╟─0eaa7e63-dee5-466c-8f9f-959cbab74f18
# ╟─7c27b82b-6438-4295-a001-503f2a200fbc
# ╟─835848e9-8221-40ad-8b8d-ff427660a499
# ╟─709ebdbe-e7a1-4471-a89c-d2dae3aecc36
# ╟─995f2a3e-89ba-4afd-8845-deff05a651d2
# ╟─a901759c-addd-4485-94de-32794f4ae0f0
# ╟─ab835496-0843-40c5-ba4a-37b4fd5f2e36
# ╟─c1779c06-e40c-415c-b268-2fecb58fe087
# ╟─cd1a8d0b-ea3a-43d5-b2be-c169b35e97ce
# ╟─ff92f822-19b1-4c03-b29c-3227bd25de45
# ╠═7af40fdd-5e66-4802-9a11-bd5e9f50a1bf
# ╠═f85a322f-c1b7-44d1-b7f8-234c9ebce22b
# ╠═e04e4cdb-90c8-4514-a673-ef8ecb7e6507
