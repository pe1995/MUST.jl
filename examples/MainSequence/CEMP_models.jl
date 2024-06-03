### A Pluto.jl notebook ###
# v0.19.41

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
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	split.(allruns, "/", keepempty=false) .|> last
end;

# ╔═╡ 0f4989c9-2c2b-4733-8f35-324e74749b11
allRuns = availableRuns(gridfolder)

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
allParasName = parametersFromName.(allRuns, carbon=true, vmic=true)

# ╔═╡ 0327e1f0-d924-4c7a-a1a3-296393de8c73
teffName = [p["teff"] for p in allParasName]

# ╔═╡ a043b27f-f392-4b96-b7b5-22cef4836a63
loggName = [p["logg"] for p in allParasName]

# ╔═╡ 9a8f1cf6-3577-4e6b-8bc7-354f5a13355d
fehName = [p["feh"] for p in allParasName]

# ╔═╡ d839813d-616c-4644-90b6-d8998a55767d
CFeName = [p["cfe"] for p in allParasName]

# ╔═╡ 58b73629-0d06-4194-8823-0d9857270f8f
vmicName = [p["vmic"] for p in allParasName]

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
iNameParas(allParasName, teff, logg, feh, cfe) = begin
	m = [
		(p["teff"] == teff) & 
		(p["logg"] == logg) & 
		(p["feh"] == feh) & 
		(p["cfe"] == cfe) 
		for p in allParasName
	]

	findfirst(m)
end

# ╔═╡ e28a5cc2-bd7b-4892-90a4-3ea34076d9af
teffU, loggU, fehU = uniqueGrid(teffName, loggName, fehName)

# ╔═╡ c3b2a375-9def-4e95-a10d-64c910d0c77c


# ╔═╡ 12b14982-f0c5-46fe-8c0e-379f4da4661e
md"# Overview"

# ╔═╡ 81535166-89a3-4321-b620-f8e2c78d5ce4
md"## HR-diagram"

# ╔═╡ 3a0c3d77-af08-4018-825f-d055779a5911
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

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
crashedModels = lastSnap .< 175

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
DataFrame(
	Dict(
		"model"=>[r for (i, r) in enumerate(allRuns) if crashedModels[i]],
		"teff"=>[r for (i, r) in enumerate(teffName) if crashedModels[i]],
		"logg"=>[r for (i, r) in enumerate(loggName) if crashedModels[i]]
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

# ╔═╡ c3b9f946-4fae-4aca-9722-a80dc999b95d
let
	npanles = length(notCrashedRuns)
	nx = ceil(Int, sqrt(npanles))
	ny = floor(Int, npanles / nx)

	plt.close()
	f, ax = plt.subplots(nx, ny, figsize=(10, 10))

	k = 1
	for j in 0:ny-1
		for i in 0:nx-1
			if k>npanles
				ax[i][j].axis("off")
				continue
			end
			w = MUST.WatchDog(notCrashedRuns[k], folder=gridfolder)
			ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
			m = MUST.reload(w, last(ws))["opticalSurfaces"]
		
			x = m["x"] ./1e8
			y = m["y"] ./1e8
		
			extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		
			im = ax[i][j].imshow(
				m["uzplane"] ./1e5,
				origin="lower",
				extent=extent,
				cmap="coolwarm",
				aspect="auto"
			)

			ax[i][j].text(
				0.97, 0.97, 
				"$(teff[k]) K,"*" $(logg[k]) dex", 
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

	f.savefig(joinpath(gridfolder, "optical_surface_uz.pdf"), bbox_inches="tight")
	gcf()
end

# ╔═╡ 4e9f26c6-001c-4d90-b9b3-ef4f259e0ae9


# ╔═╡ 9b659650-b91d-4229-a811-6802f7d56966
which_gif = "uzplane"

# ╔═╡ 6df1c2d5-850c-4af8-b29d-25d7a0f7666c
begin
	vminmax = []
	
	for k in 1:length(notCrashedRuns)
		w = MUST.WatchDog(notCrashedRuns[k], folder=gridfolder)
		ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
		m = MUST.reload(w, last(ws))["opticalSurfaces"]
		
		append!(
			vminmax, [(minimum(m[which_gif]), maximum(m[which_gif]))]
		)
	end

	vminmax
end

# ╔═╡ 115b47e6-5153-4940-8499-0b19de8f2e54
function optical_surface(notCrashedRuns, teff, logg; which, vminmax=nothing, cmap="coolwarm", quantity="Tplane")
	npanles = length(notCrashedRuns)
	nx = ceil(Int, sqrt(npanles))
	ny = floor(Int, npanles / nx)

	plt.close()
	f, ax = plt.subplots(nx, ny, figsize=(10, 10))

	k = 1
	for j in 0:ny-1
		for i in 0:nx-1
			if k>npanles
				ax[i][j].axis("off")
				continue
			end
			w = MUST.WatchDog(notCrashedRuns[k], folder=gridfolder)
			ws = MUST.snapshotnumber.(MUST.availableSnaps(w))
			whichRed = min(which, length(ws))
			m = MUST.reload(w, ws[whichRed])["opticalSurfaces"]
		
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
				"$(teff[k]) K,"*" $(logg[k]) dex", 
				ha="right", va="top", 
				transform=ax[i][j].transAxes,
				fontsize="small",
				color="white", backgroundcolor="k"
			)
			
			k += 1
		end
	end

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
	#=f = []
	for i in 1:360
		append!(f, [optical_surface(
			notCrashedRuns, notCrashedteff, notCrashedlogg; 
			which=i, 
			vminmax=vminmax,
			quantity=which_gif,
			cmap="coolwarm"
		)])
	end

	visual.gifs.gifs_from_png(
		f, joinpath(gridfolder, "optical_surface_uz.gif");
		duration=0.5
	)=#
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

# ╔═╡ faba53f3-1a3b-4211-9205-928404525ed4
begin
	mass_flux = zeros(length(allRuns))
	for k in 1:length(allRuns)
		w = MUST.WatchDog(allRuns[k], folder=gridfolder)
		m = MUST.reload!(w, groups=["geoMassFlux", "atmosphericParameters"], lastN=30)
		mf = [mi["geoMassFlux"]["massFlux"][end] for mi in m] ./1e5
		t = [mi["atmosphericParameters"]["time"] for mi in m]
		mass_flux[k] = MUST.mean(mf) 
	end
end

# ╔═╡ 97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
logg_vertical_line = 3.2

# ╔═╡ d51f056b-6ef3-430b-aacf-04c14caeb2ea
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	vmin = minimum(mass_flux)
	vmax = maximum(mass_flux)
	#vtot = max(abs(vmin), abs(vmax))
	#vmin = -vtot
	#vmax = vtot
	
	i = ax.plot(loggName, mass_flux, marker="s", ls="", color="k")

	ax.axvline(logg_vertical_line, color="red")
	ax.axhline(0.0, color="red")
	
	
	ax.set_xlabel(L"\rm \log(g)\ [dex]")
	ax.set_ylabel(L"\rm mass\ flux\ [km/s]")
	
	f
end

# ╔═╡ 93e362e6-be03-4e9c-b11d-4002fc60211d
let
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
	
	i = ax.plot(loggName, htop_scale, marker="s", ls="", color="k", zorder=0)
	ax.scatter(
		[3.5, 4.15, 4.44, 5.0], 
		[7.5, 2.0, 1.1, 0.5],
		s=50, 
		color="red",
		label="usefull points"
	)

	x = range(3.0, 5.5, length=100)
	ax.plot(x, func.(x), color="magenta")
	ax.plot(x, func2.(x), color="cyan")

	#ax.set_ylim(0, 3.0)
	#ax.set_xlim(3.7, 5.1)
	
	ax.legend()
	ax.axvline(logg_vertical_line, color="red")

	ax.set_xlabel(L"\rm \log(g)\ [dex]")
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
$(@bind pickedParameters Select(paraUStrings))
"""

# ╔═╡ fc8be706-2541-4e4e-a978-c657898c3c6d
iselect = findfirst(paraUStrings.==pickedParameters) ;

# ╔═╡ 1f479139-b18b-4c94-9cec-e84ea2e490d2
begin
	iTeffLoggCU = ifindTeffLogg(teffU[iselect], loggU[iselect], teffName, loggName)
	carbonsU = sort(CFeName[iTeffLoggCU])
	md"""
		Select [C/Fe] abundances to compare
	$(@bind cfeSelect1 PlutoUI.Select(carbonsU, default=first(carbonsU)))
	$(@bind cfeSelect2 PlutoUI.Select(carbonsU, default=last(carbonsU)))
	"""
end

# ╔═╡ 7c27b82b-6438-4295-a001-503f2a200fbc
begin
	iselect1 = findfirst(CFeName[iTeffLoggCU] .== cfeSelect1)
	iselect2 = findfirst(CFeName[iTeffLoggCU] .== cfeSelect2)
end;

# ╔═╡ 835848e9-8221-40ad-8b8d-ff427660a499
begin
	@info "Loading Monitoring..."
	#monitoring1 = MUST.reload!(MUST.WatchDog)
	iNameSelect1 = iNameParas(
		allParasName, teffU[iselect], loggU[iselect], fehU[iselect], carbonsU[iselect1]
	)
	iNameSelect2 = iNameParas(
		allParasName, teffU[iselect], loggU[iselect], fehU[iselect], carbonsU[iselect2]
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

	ax.set_title(pickedParameters)
	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label="[C/Fe] = $cfeSelect1", color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = optAv["T"]
	ax.plot(tau, t, label="[C/Fe] = $cfeSelect2", color="magenta")


	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm T\ [K]")
	
	
	ax.legend()
	f
end

# ╔═╡ 995f2a3e-89ba-4afd-8845-deff05a651d2
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	ax.set_title(pickedParameters)
	
	optAv = monitoring1[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label="[C/Fe] = $cfeSelect1", color="k")

	optAv = monitoring2[end]["opticalAverages"]
	tau = log10.(optAv["τ_ross"])
	t = log10.(optAv["d"])
	ax.plot(tau, t, label="[C/Fe] = $cfeSelect2", color="magenta")

	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm \rho\ [g\ cm^{-3}]")
	
	ax.legend()
	f
end

# ╔═╡ a901759c-addd-4485-94de-32794f4ae0f0
md"# Summary"

# ╔═╡ ab835496-0843-40c5-ba4a-37b4fd5f2e36
function add_point!(f, ax, teff, logg, segments, colors; square_logg_length=0.1, square_teff_length=0.1, alpha=0.3, add_label=false, title=nothing)
	rect = plt.Rectangle(
		(teff - square_teff_length/2, logg - square_logg_length/2), square_teff_length, square_logg_length,
		edgecolor="black", 
		facecolor="None",
		lw=0
	)
	ax.add_patch(rect)
	ax.scatter((teff - square_teff_length/2, logg - square_logg_length/2)..., s=0)

	# for each segment, fill in given color
	nseg = length(segments)
	for i in eachindex(segments)
		smin = logg - square_logg_length/2 + (i-1)*(square_logg_length/nseg)
		smax = smin + square_logg_length/nseg

		ax.fill_between(
			[teff - square_teff_length/2, teff + square_teff_length/2],
			[smin, smin], [smax, smax],
			color=colors[i],
			edgecolor="black",
			alpha=alpha
		)

		if add_label
			ax.text(
				teff, smax, segments[i], va="bottom", ha="center",
				
			)
		end
	end

	if !isnothing(title)
		ax.text(
			teff, logg - (square_logg_length/2), title, va="bottom", ha="center",
		)
	end
end

# ╔═╡ c1779c06-e40c-415c-b268-2fecb58fe087


# ╔═╡ cd1a8d0b-ea3a-43d5-b2be-c169b35e97ce
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	cfes = [0.0, 1.0, 2.0, 3.0]
	trad = 150
	lrad = 0.01
	damp = 0.0003
	s = 100

	dteff = abs(maximum(teffU)-minimum(teffU))
	dlogg = abs(maximum(loggU)-minimum(loggU))
	upperLogg = 4. #maximum(loggU)+0.1*dlogg
	lowerLogg = 2.8 #minimum(loggU)-0.1*dlogg
	upperTeff = 7000 #maximum(teffU)+0.1*dteff
	lowerTeff = 4000 #minimum(teffU)-0.1*dteff

	square_logg_length=0.07*(upperLogg-lowerLogg)
	square_teff_length=0.07*(upperTeff-lowerTeff)

	for i in eachindex(teffU)
		#=coords = circe_coordinates(
			teffName[i], loggName[i], cfes; 
			radius=[trad, lrad], damp=damp
		)
		#ax.scatter(teffName[i], loggName[i])

		stati = ["" for _ in cfes]
		stati[2] = if crashedModels[i]
			 "crashed"
		else
			"done"
		end

		plot_coords!(ax, coords, stati, s=s, marker="o", alpha=1, linewidth=1.1)=#

		iTeffLogg = ifindTeffLogg(teffU[i], loggU[i], teffName, loggName)
		segs = CFeName[iTeffLogg]
		crashed = crashedModels[iTeffLogg]
		colors = [c ?  "magenta" : "lime" for c in crashed]
		colors = [
			cfes[i] in segs ? colors[findfirst(cfes[i].==segs)] : "w" 
			for i in eachindex(cfes)
		]
		add_point!(
			f, ax, teffU[i], loggU[i], cfes, colors, 
			square_logg_length=square_logg_length,
			square_teff_length=square_teff_length,
			alpha=0.8
		)
	end

	ax.set_ylim(upperLogg, lowerLogg)
	ax.set_xlim(upperTeff, lowerTeff)

	x_start = first(MUST.pyconvert(Any, ax.get_xlim()))
	y_start = first(MUST.pyconvert(Any, ax.get_ylim()))
	x_end = last(MUST.pyconvert(Any, ax.get_xlim()))
	y_end = last(MUST.pyconvert(Any, ax.get_ylim()))
	
	ax_x_end = x_end + abs(x_end - x_start) *0.1
	ax_y_end = y_end + abs(y_end - y_start) *0.1
	ax_x_start = x_start - abs(x_end - x_start) *0.1
	ax_y_start = y_start - abs(y_end - y_start) *0.1

	add_point!(
		f, ax, ax_x_start, ax_y_start, cfes, ["w" for _ in cfes], 
		square_logg_length=square_logg_length*1.7,
		square_teff_length=square_teff_length*1.7,
		alpha=0.8, add_label=true, title="[C/Fe]"
	)
	
	#=coords = circe_coordinates(
		x_end, y_end, cfes; 
		radius=[trad*3, lrad*3], damp=damp
	)
	plot_coords!(ax, coords,  ["" for _ in cfes], s=s*10, marker="o")
	plot_coords_labels!(ax, coords, cfes)=#

	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	
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

# ╔═╡ Cell order:
# ╠═30832d50-dc75-11ee-2e49-5f52349e5e50
# ╠═cadc96a0-b6a8-4e46-8240-d1e72f976428
# ╟─6cb9c336-4f1d-43bd-b626-a7404a44d33b
# ╟─6c5b6a85-0b7d-4afd-a068-38ce7c9ac52e
# ╟─1b3d7608-1eec-4d35-a19c-a95562828194
# ╟─e99d2039-091c-4749-a355-b5cf124c7e46
# ╟─64d370db-ca48-44ef-9f6c-c676a45fbcd0
# ╟─0f4989c9-2c2b-4733-8f35-324e74749b11
# ╟─68548a3f-5712-4ee9-99e0-6aad41624edc
# ╟─2bc84808-c748-49fb-9d27-691ddef1355f
# ╟─1199b473-e011-4817-be65-0a685242cf2a
# ╠═d749afd7-a6fa-43c1-9c6f-0a0cd2c0ec0b
# ╠═0327e1f0-d924-4c7a-a1a3-296393de8c73
# ╠═a043b27f-f392-4b96-b7b5-22cef4836a63
# ╠═9a8f1cf6-3577-4e6b-8bc7-354f5a13355d
# ╠═d839813d-616c-4644-90b6-d8998a55767d
# ╠═58b73629-0d06-4194-8823-0d9857270f8f
# ╟─898ce0f5-7ec8-4f30-a4fe-96956ed1072f
# ╟─3b31e0f2-8acc-49dc-b69e-dd084bae4f79
# ╟─c9cdfa8a-ee49-4e46-8e8e-4e7cd668b98a
# ╟─09505e68-b2d6-4c07-8e3f-a0e097ea3732
# ╟─e28a5cc2-bd7b-4892-90a4-3ea34076d9af
# ╟─c3b2a375-9def-4e95-a10d-64c910d0c77c
# ╟─12b14982-f0c5-46fe-8c0e-379f4da4661e
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
# ╠═095d5c29-d828-409c-afd9-2384d1bb3a83
# ╟─2dd0de7c-028b-43f2-960e-e0839d3fa7b7
# ╟─c3b9f946-4fae-4aca-9722-a80dc999b95d
# ╟─4e9f26c6-001c-4d90-b9b3-ef4f259e0ae9
# ╠═9b659650-b91d-4229-a811-6802f7d56966
# ╟─6df1c2d5-850c-4af8-b29d-25d7a0f7666c
# ╟─115b47e6-5153-4940-8499-0b19de8f2e54
# ╠═53892118-270c-40c7-a12b-99aef4d9e69b
# ╟─8ab32b10-f50a-4fdd-b0b6-dd0c1d154e1f
# ╟─4866dd18-e1f8-4939-83ec-f2e80407f731
# ╟─73ea53c8-eb56-462c-8612-41c33d379915
# ╟─faba53f3-1a3b-4211-9205-928404525ed4
# ╠═97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
# ╟─d51f056b-6ef3-430b-aacf-04c14caeb2ea
# ╟─93e362e6-be03-4e9c-b11d-4002fc60211d
# ╟─1799306b-954d-402e-8652-a66cf5f3389b
# ╟─b2b24be9-6db7-46c5-a69c-3d790c4d5406
# ╟─ae86742d-6e97-4f32-9edb-e0feeb9f446d
# ╟─0eaa7e63-dee5-466c-8f9f-959cbab74f18
# ╟─fc8be706-2541-4e4e-a978-c657898c3c6d
# ╟─1f479139-b18b-4c94-9cec-e84ea2e490d2
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
