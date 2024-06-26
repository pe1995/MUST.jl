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

# ╔═╡ 8aed87bc-2c8a-11ef-12b8-4733d141661b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PythonPlot
	using PlutoUI
	using ProgressLogging
end

# ╔═╡ c97b657a-cfb6-4c5c-90d9-452842a1685c
TableOfContents()

# ╔═╡ 5f35847a-2449-4c7e-b4cf-872599a2b767
@import_dispatch "../../../dispatch2"

# ╔═╡ 0a8b6595-cd8f-4345-a3b1-fb49f4b3c540
datafolder = @in_dispatch "data"

# ╔═╡ c70c90d3-5c22-4f9b-b4d3-81f89ea647cd
begin
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
end;

# ╔═╡ 4df7bf78-34b1-4758-8dfe-792d8962ec37
md"# Setup
## Select Run"

# ╔═╡ eb2799c1-c8cc-48b8-aad4-309dd5a46e57
md"Press the button to check for new monitored runs: $(@bind reload_data Button(\"Refresh\"))"

# ╔═╡ 2220e9f7-6534-48ec-8430-1b745e9ba4f6
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	split.(allruns, "/", keepempty=false) .|> last
end;

# ╔═╡ e102e1fa-4a66-4cdb-a762-8cc132f1ad0f
begin
	reload_data
	md"""
	Select one of the available monitored runs
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))
	"""
end

# ╔═╡ 315f76d5-1d3f-4854-a964-aa85af72f62f


# ╔═╡ 1ed280cd-7b36-4a29-8b84-3cedff7a72d7
wd = MUST.WatchDog(selectedRun, folder=datafolder)

# ╔═╡ 4d0c6234-caa1-4d5c-9502-be07f8ecf089
monitoring = MUST.reload!(wd)

# ╔═╡ 7f5505f8-c060-4af6-92a5-9840bbe30d23
snapshots = wd.snapshotsCompleted[1:end];

# ╔═╡ bb9a5578-0581-467c-9612-97e9f0b2ba55
md"""
	Pick the snapshot that you want to extrapolate
$(@bind timeSurface Slider(snapshots, show_value=true, default=last(snapshots)))
"""

# ╔═╡ 61fbc6ee-176e-4886-8a16-02d044550ad6
itimeSurface = findfirst(snapshots .== timeSurface);

# ╔═╡ 6bf5b73a-853a-47e1-af8e-a1111011a725
md"""
	Extrapolate by how much (%)
$(@bind extrapPercent TextField(default=\"5.0\"))
"""

# ╔═╡ eb9490a5-b22e-4912-aa00-c293af514b1c
md"# Extrapolation"

# ╔═╡ dc37f14a-bcb0-4ac5-b792-7eaeea40b27f
averageData = MUST.timeevolution(monitoring, "geometricalAverages")

# ╔═╡ 0bc765cf-4f4c-467d-9a17-dd9a5696ae07
begin
	z, T, d = averageData["z"][itimeSurface], averageData["T"][itimeSurface], averageData["d"][itimeSurface]
end

# ╔═╡ c56fe342-2a95-40ca-9aa1-e70312af7b79
begin
	fextrap = parse(Float64, extrapPercent) /100.0
	znew = range(
		minimum(z), 
		maximum(z)+abs(minimum(z)-maximum(z))*fextrap,
		length=ceil(Int, length(z)*(1.0 + fextrap))
	) |> collect
end

# ╔═╡ 0b606600-03e9-4f25-8d27-87ae2ce581e3
begin
	Tnew = MUST.linear_interpolation(z, T, extrapolation_bc=MUST.Line()).(znew)
	dnew = MUST.linear_interpolation(z, log.(d), extrapolation_bc=MUST.Line()).(znew)
end

# ╔═╡ 4c13decd-0858-4c7f-9121-543486f34529
let
	f, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=true)
	
	ax[0].plot(znew, Tnew)
	ax[1].plot(znew, dnew)

	ax[0].plot(z, T)
	ax[1].plot(z, log.(d))

	f
end

# ╔═╡ e8512de0-ee46-474f-9f0d-20aaf4e33426
open(joinpath(datafolder, selectedRun, "avextra_$(timeSurface).dat"), "w") do f
	MUST.writedlm(f, [znew Tnew dnew])
end

# ╔═╡ Cell order:
# ╠═8aed87bc-2c8a-11ef-12b8-4733d141661b
# ╟─c97b657a-cfb6-4c5c-90d9-452842a1685c
# ╠═5f35847a-2449-4c7e-b4cf-872599a2b767
# ╠═0a8b6595-cd8f-4345-a3b1-fb49f4b3c540
# ╟─c70c90d3-5c22-4f9b-b4d3-81f89ea647cd
# ╟─4df7bf78-34b1-4758-8dfe-792d8962ec37
# ╟─eb2799c1-c8cc-48b8-aad4-309dd5a46e57
# ╟─2220e9f7-6534-48ec-8430-1b745e9ba4f6
# ╟─e102e1fa-4a66-4cdb-a762-8cc132f1ad0f
# ╟─315f76d5-1d3f-4854-a964-aa85af72f62f
# ╟─1ed280cd-7b36-4a29-8b84-3cedff7a72d7
# ╟─4d0c6234-caa1-4d5c-9502-be07f8ecf089
# ╟─7f5505f8-c060-4af6-92a5-9840bbe30d23
# ╟─bb9a5578-0581-467c-9612-97e9f0b2ba55
# ╟─61fbc6ee-176e-4886-8a16-02d044550ad6
# ╟─6bf5b73a-853a-47e1-af8e-a1111011a725
# ╟─eb9490a5-b22e-4912-aa00-c293af514b1c
# ╟─dc37f14a-bcb0-4ac5-b792-7eaeea40b27f
# ╠═0bc765cf-4f4c-467d-9a17-dd9a5696ae07
# ╠═c56fe342-2a95-40ca-9aa1-e70312af7b79
# ╠═0b606600-03e9-4f25-8d27-87ae2ce581e3
# ╟─4c13decd-0858-4c7f-9121-543486f34529
# ╠═e8512de0-ee46-474f-9f0d-20aaf4e33426
