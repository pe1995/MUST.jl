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

# ╔═╡ 30832d50-dc75-11ee-2e49-5f52349e5e50
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using PlutoUI
	using LaTeXStrings
end

# ╔═╡ cadc96a0-b6a8-4e46-8240-d1e72f976428
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	TableOfContents()
end

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
md"# Simulation Time"

# ╔═╡ 30533c62-c714-4685-bf26-c0e6d6948bfa
gettime(folder) = begin
	lastsnap = MUST.list_snapshots(MUST.converted_snapshots(folder)) |> last
	snap, _ = pick_snapshot(folder, lastsnap)
	snap.parameter.time
end

# ╔═╡ d0573093-6878-4d38-a2ef-0e6a54fdc878
simTime = gettime.(joinpath.(gridfolder, allRuns))

# ╔═╡ 85cd82ca-5417-4b86-80b8-e40d680b96c2
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	im = ax.scatter(loggName, simTime ./60 ./60, c=teffName)
	f.colorbar(im, ax=ax)

	ax.set_xlabel(L"\rm \log(g)\ [dex]")
	ax.set_ylabel(L"\rm stellar\ time\ [h]")
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═30832d50-dc75-11ee-2e49-5f52349e5e50
# ╟─cadc96a0-b6a8-4e46-8240-d1e72f976428
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
# ╠═d0573093-6878-4d38-a2ef-0e6a54fdc878
# ╟─85cd82ca-5417-4b86-80b8-e40d680b96c2
