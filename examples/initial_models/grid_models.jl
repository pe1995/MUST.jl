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
function parametersFromName(path)
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

	Dict("teff"=>t, "logg"=>g, "feh"=>m)
end

# ╔═╡ d749afd7-a6fa-43c1-9c6f-0a0cd2c0ec0b
allParasName = parametersFromName.(allRuns)

# ╔═╡ 0327e1f0-d924-4c7a-a1a3-296393de8c73
teffName = [p["teff"] for p in allParasName]

# ╔═╡ a043b27f-f392-4b96-b7b5-22cef4836a63
loggName = [p["logg"] for p in allParasName]

# ╔═╡ 9a8f1cf6-3577-4e6b-8bc7-354f5a13355d
fehName = [p["feh"] for p in allParasName]

# ╔═╡ e28a5cc2-bd7b-4892-90a4-3ea34076d9af


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
mass_flux = zeros(length(allRuns))

# ╔═╡ 9a482623-0e09-419d-b2b6-29edea2dcee1
let
	mass_flux .= 0.0
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	for k in 1:length(allRuns)
		w = MUST.WatchDog(allRuns[k], folder=gridfolder)
		m = MUST.reload!(w, groups=["geoMassFlux", "atmosphericParameters"])
		mf = [mi["geoMassFlux"]["massFlux"][end] for mi in m] ./1e5
		t = [mi["atmosphericParameters"]["time"] for mi in m]
		istart = max(length(mf) - 30, 1)
		mass_flux[k] = MUST.mean(mf[istart:end]) 

		ax.plot(t, mf, label="$(teffName[k]) K,"*" $(loggName[k]) dex")
	end

	#ax.legend(bbox_to_anchor=(1.1, 1.01), loc="upper left")
	f
end

# ╔═╡ 391731ff-cef3-459b-9b48-c0f48be1b543
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	vmin = minimum(mass_flux)
	vmax = maximum(mass_flux)
	#vtot = max(abs(vmin), abs(vmax))
	#vmin = -vtot
	#vmax = vtot
	
	i = ax.scatter(teffName, loggName, c=mass_flux, vmin=vmin, vmax=vmax, cmap="RdYlBu")
	f.colorbar(i, ax=ax)

	ax.set_xlim(7100, 4000)
	ax.set_ylim(5, 3.5)

	ax.set_ylabel(L"\rm \log(g)\ [dex]")
	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	
	f
end

# ╔═╡ 97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
logg_vertical_line = 4.15

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
	#func2(x) = max(exp(-1.15*(x-5.0))-0.6, 0.1)
	
	i = ax.plot(loggName, htop_scale, marker="s", ls="", color="k", zorder=0)
	ax.scatter(
		[3.5, 4.15, 4.44, 5.0], 
		[7.5, 2.0, 1.1, 0.5],
		s=50, 
		color="red",
		label="usefull points"
	)

	x = range(3.4, 5.5, length=100)
	ax.plot(x, func.(x), color="magenta")
	#ax.plot(x, func2.(x), color="cyan")

	ax.set_ylim(0, 3.0)
	ax.set_xlim(3.7, 5.1)
	
	ax.legend()
	ax.axvline(logg_vertical_line, color="red")

	ax.set_xlabel(L"\rm \log(g)\ [dex]")
	ax.set_ylabel(L"\rm top\ boundary\ scaling\ factor")
	
	f
end

# ╔═╡ b2b24be9-6db7-46c5-a69c-3d790c4d5406


# ╔═╡ a901759c-addd-4485-94de-32794f4ae0f0
md"# Summary"

# ╔═╡ 4d74f7b2-953b-4ebd-b658-eef7ce6a4e2d
function circe_coordinates(teff, logg, fehs; radius=[60, 0.05], damp=0.0003)
	# angle of splits
	nsplits = length(fehs)
	sepangle = 2π / nsplits
	angles = [(i-1) * sepangle for i in 1:nsplits]	
	
	r = sqrt(sum(radius .^2))
	rx = [sin(a) * r for a in angles] .+ teff
	ry = [cos(a) * r * damp for a in angles] .+ logg

	zip(rx, ry) |> collect
end

# ╔═╡ 75a8cac3-a04a-4497-aa45-f303b057f1ec
function plot_coords!(ax, coords, stati; kwargs...)
	colors = []
	for i in stati
		if i == "done"
			append!(colors, ["g"])
		elseif i == "in progress"
			append!(colors, ["y"])
		elseif i == "crashed"
			append!(colors, ["r"])
		else
			append!(colors, ["w"])
		end
	end
	for j in eachindex(coords)
		ax.scatter(coords[j]...; edgecolor="k", color=colors[j], kwargs...)
	end
end

# ╔═╡ 11f4dae1-181d-4266-826f-3d6ebc3a7a8c
function plot_coords_labels!(ax, coords, fehs; kwargs...)
	for j in eachindex(coords)
		l = "$(fehs[j])"
		ax.text(coords[j]..., l; ha="center", va="center", kwargs...)
	end
end

# ╔═╡ cd1a8d0b-ea3a-43d5-b2be-c169b35e97ce
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	fehs = [-0.5, 0, 0.5]
	trad = 55
	lrad = 0.01
	damp = 0.0003
	s = 100
	
	for i in eachindex(teffName)
		coords = circe_coordinates(
			teffName[i], loggName[i], fehs; 
			radius=[trad, lrad], damp=damp
		)
		#ax.scatter(teffName[i], loggName[i])

		stati = ["" for _ in fehs]
		stati[2] = if crashedModels[i]
			 "crashed"
		else
			"done"
		end

		plot_coords!(ax, coords, stati, s=s, marker="o", alpha=1, linewidth=1.1)
	end

	ax.set_xlim(7100, 4000)
	ax.set_ylim(5, 3.9)

	x_start = first(MUST.pyconvert(Any, ax.get_xlim()))
	y_start = first(MUST.pyconvert(Any, ax.get_ylim()))
	x_end = last(MUST.pyconvert(Any, ax.get_xlim()))
	y_end = last(MUST.pyconvert(Any, ax.get_ylim()))
	
	x_end = x_end + abs(x_end - x_start) *0.15
	y_end = y_end + abs(y_end - y_start) *0.15
	
	coords = circe_coordinates(
		x_end, y_end, fehs; 
		radius=[trad*3, lrad*3], damp=damp
	)
	plot_coords!(ax, coords,  ["" for _ in fehs], s=s*10, marker="o")
	plot_coords_labels!(ax, coords, fehs)

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
# ╟─e28a5cc2-bd7b-4892-90a4-3ea34076d9af
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
# ╟─9a482623-0e09-419d-b2b6-29edea2dcee1
# ╟─391731ff-cef3-459b-9b48-c0f48be1b543
# ╠═97d8e869-ff30-4cd7-81fa-1c0fafe3ff3f
# ╟─d51f056b-6ef3-430b-aacf-04c14caeb2ea
# ╟─93e362e6-be03-4e9c-b11d-4002fc60211d
# ╟─b2b24be9-6db7-46c5-a69c-3d790c4d5406
# ╟─a901759c-addd-4485-94de-32794f4ae0f0
# ╟─4d74f7b2-953b-4ebd-b658-eef7ce6a4e2d
# ╟─75a8cac3-a04a-4497-aa45-f303b057f1ec
# ╟─11f4dae1-181d-4266-826f-3d6ebc3a7a8c
# ╟─cd1a8d0b-ea3a-43d5-b2be-c169b35e97ce
# ╟─ff92f822-19b1-4c03-b29c-3227bd25de45
# ╠═7af40fdd-5e66-4802-9a11-bd5e9f50a1bf
