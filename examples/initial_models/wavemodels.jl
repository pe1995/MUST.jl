### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ 05f376b8-9377-11ee-1276-cb53de67918a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using PythonPlot
	using MUST
	using FFTW
	using PlutoUI
	using Glob
	using PlutoUI: combine
end

# ╔═╡ 884dbb9a-8d66-4699-b43a-f1063d3e0ff5
md"# Setup"

# ╔═╡ 8966f279-e5a1-4ce0-95a1-1e090be857ad
TableOfContents()

# ╔═╡ 303f3449-4aad-4053-9cbb-5a0e119636b0
@import_dispatch "../../../dispatch2"

# ╔═╡ a2fc92e5-63a3-4d4c-b508-10423c5eb3c1
begin
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
end

# ╔═╡ 54ebd9ae-82b8-4d08-92ed-b4345fd2cd26
availableRuns(path) = begin
	runs = glob("*/", @in_dispatch(path))
	last.(split.(runs, "/", keepempty=false))
end

# ╔═╡ 9c61e8b0-24af-4376-a60a-025e4e6a0117
md"# Data Cubes"

# ╔═╡ 508bd668-1759-430c-9558-bbe98990a361
md"## Models"

# ╔═╡ 075f2abe-e690-4635-8478-aa1fed1b1f84
datafolder = @in_dispatch("data")

# ╔═╡ ff921691-4a1b-453b-a720-58efd5c50f6e
function snapshots_all_input(simulation_names)
	o = @in_dispatch(joinpath(datafolder, simulation_names))
	snaps_already = MUST.list_snapshots(converted_snapshots(o))
end

# ╔═╡ 3ba3caf4-ed08-4134-a808-467d2ec3edb0


# ╔═╡ 2744cb7b-817b-4564-b0bb-6129eff394b0
md"Pick one of the available models: $(@bind available_runs Select(availableRuns(datafolder)))"

# ╔═╡ be4cdbef-12b3-4f80-a2f7-aae4350070cd
snapshots_converted = snapshots_all_input(available_runs)

# ╔═╡ 90aa096a-9be3-40a1-8634-8ceb800f5706
begin
	snapshots = Dict(available_runs=>[])
	snapshots_τ = Dict(available_runs=>[])
	
	for v in snapshots_converted
		p = @in_dispatch(joinpath(datafolder, "$(available_runs)"))
		snapshot, snapshot_τ = pick_snapshot(
			converted_snapshots(p), v
		)
		append!(snapshots[available_runs], [snapshot])
		append!(snapshots_τ[available_runs], [snapshot_τ])
	end
end

# ╔═╡ 5567bace-a78b-4934-b752-c8756cdf2821
md"# Fourier Transformation"

# ╔═╡ bff06662-8231-4323-b2c4-be2618f7230a
md"## FFT"

# ╔═╡ 0989c716-67b2-45a7-915f-3c6358274da3
function waves(snapshots)
	times = [snap.parameter.time for snap in snapshots]
	timemask = sortperm(times)

	velocity_field = zeros(size(snapshots[1])..., 3, length(times))
	for i in eachindex(times)
		j = timemask[i]
		velocity_field[:, :, :, 1, j] .= snapshots[j][:ux] ./1e5
		velocity_field[:, :, :, 2, j] .= snapshots[j][:uy] ./1e5
		velocity_field[:, :, :, 3, j] .= snapshots[j][:uz] ./1e5
	end
	
	# perform the Fourier transformation
	wave_field = fft(velocity_field, 4)
	frequencies = fftfreq(length(times), diff(times)[1])

	frequencies, wave_field
end

# ╔═╡ cb615a2d-1fb4-4895-9cb5-83df1674a825
f, w = waves(snapshots[available_runs])

# ╔═╡ f6341fc0-ad5e-41bc-8580-ebe428100997
md"## Components"

# ╔═╡ 291d4372-5983-4b17-8913-fc5f5d7987a2
cmap = plt.get_cmap("seismic")

# ╔═╡ 3a93111f-e883-44fb-b068-38ddc894f922
absw = abs2.(w);

# ╔═╡ d2c7819c-c913-411c-a0e4-9bee8e132c5d
imw = imag.(w);

# ╔═╡ cef22779-6abb-4afe-96b9-d80f21e1da5e
rew = real.(w);

# ╔═╡ 6edb0d17-dee4-4476-b849-58f5869880e8
let
	fig, ax = plt.subplots(1, 1)

	m = sortperm(f)
	heights = axes(w, 3)[1:10:end] |> collect
	colors = [cmap(ic/length(heights)) for ic in eachindex(heights)]

	for (i, h) in enumerate(heights)
		ax.plot(f[m], absw[10, 10, h, 3, m], color=colors[i])
	end
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─884dbb9a-8d66-4699-b43a-f1063d3e0ff5
# ╠═05f376b8-9377-11ee-1276-cb53de67918a
# ╟─8966f279-e5a1-4ce0-95a1-1e090be857ad
# ╠═303f3449-4aad-4053-9cbb-5a0e119636b0
# ╠═a2fc92e5-63a3-4d4c-b508-10423c5eb3c1
# ╟─54ebd9ae-82b8-4d08-92ed-b4345fd2cd26
# ╟─ff921691-4a1b-453b-a720-58efd5c50f6e
# ╟─9c61e8b0-24af-4376-a60a-025e4e6a0117
# ╟─508bd668-1759-430c-9558-bbe98990a361
# ╟─075f2abe-e690-4635-8478-aa1fed1b1f84
# ╟─3ba3caf4-ed08-4134-a808-467d2ec3edb0
# ╟─2744cb7b-817b-4564-b0bb-6129eff394b0
# ╟─be4cdbef-12b3-4f80-a2f7-aae4350070cd
# ╟─90aa096a-9be3-40a1-8634-8ceb800f5706
# ╟─5567bace-a78b-4934-b752-c8756cdf2821
# ╟─bff06662-8231-4323-b2c4-be2618f7230a
# ╠═0989c716-67b2-45a7-915f-3c6358274da3
# ╠═cb615a2d-1fb4-4895-9cb5-83df1674a825
# ╟─f6341fc0-ad5e-41bc-8580-ebe428100997
# ╠═291d4372-5983-4b17-8913-fc5f5d7987a2
# ╠═3a93111f-e883-44fb-b068-38ddc894f922
# ╠═d2c7819c-c913-411c-a0e4-9bee8e132c5d
# ╠═cef22779-6abb-4afe-96b9-d80f21e1da5e
# ╠═6edb0d17-dee4-4476-b849-58f5869880e8
