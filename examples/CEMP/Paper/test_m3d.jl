### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 1db442aa-a7e2-11ef-09a8-6db52d25eb54
begin
	using Pkg; Pkg.activate("..")
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ 9d4b8d47-1ff1-47d9-9e7c-245766c982bb
plt = matplotlib.pyplot

# ╔═╡ a2fadf4b-6333-47bb-8e85-6c221802f997
@import_dispatch "../../../../dispatch2"

# ╔═╡ 5e5a74a4-d8b7-4629-9601-6183d79cec97
@import_m3dis "../../../../Multi3D"

# ╔═╡ d194fe19-5409-4d35-8acc-69f29eae4113
runpath = @in_m3dis("data/test")

# ╔═╡ 1967dac1-2006-4be6-8727-9f5ba5caf095
runpath2 = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/M3D_old/m3dis/experiments/Multi3D/data/test"

# ╔═╡ 2d2a18a7-b9be-4825-9567-88b81f4d0b99
run = M3DISRun(runpath)

# ╔═╡ 34afbb37-0b1a-4ecb-bc33-ed4863720a3c
run2 = M3DISRun(runpath2)

# ╔═╡ c237f390-1c38-4af5-9638-ff166bd7bb1f
let
	plt.close()
	f, ax = plt.subplots(figsize=(10, 5))

	x, y = MUST.flux(run, norm=true)
	ax.plot(x, y, lw=2, label="new commit")

	x, y = MUST.flux(run2, norm=true)
	ax.plot(x, y, lw=2, label="old commit")

	ax.set_xlim(4297, 4303)
	ax.set_ylim(0.95, 1.001)

	ax.legend()

	f
end

# ╔═╡ Cell order:
# ╠═1db442aa-a7e2-11ef-09a8-6db52d25eb54
# ╠═9d4b8d47-1ff1-47d9-9e7c-245766c982bb
# ╠═a2fadf4b-6333-47bb-8e85-6c221802f997
# ╠═5e5a74a4-d8b7-4629-9601-6183d79cec97
# ╠═d194fe19-5409-4d35-8acc-69f29eae4113
# ╠═1967dac1-2006-4be6-8727-9f5ba5caf095
# ╠═2d2a18a7-b9be-4825-9567-88b81f4d0b99
# ╠═34afbb37-0b1a-4ecb-bc33-ed4863720a3c
# ╠═c237f390-1c38-4af5-9638-ff166bd7bb1f
