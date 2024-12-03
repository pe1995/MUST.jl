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
runpath = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/dispatch2/experiments/stellar_atmospheres/shipped_models/shipped_PS2D_E_t52.50g30.00m-2.000_v1.0/spectra_p5250_g+3.0_m0.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_3000-3200_kurucz"

# ╔═╡ 1967dac1-2006-4be6-8727-9f5ba5caf095
runpath2 = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/dispatch2/experiments/stellar_atmospheres/shipped_models/shipped_PS2D_E_t52.50g30.00m-2.000_v1.0/spectra_p5250_g+3.0_m0.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_3000-3200_bertrand"

# ╔═╡ 15a9b32f-dfc7-47c1-99da-4a1dfb268751
runpath3 = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/dispatch2/experiments/stellar_atmospheres/shipped_models/shipped_PS2D_E_t52.50g30.00m-2.000_v1.0/spectra_p5250_g+3.0_m0.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_3000-3200_HITRAN"

# ╔═╡ 2d2a18a7-b9be-4825-9567-88b81f4d0b99
run = M3DISRun(runpath)

# ╔═╡ 34afbb37-0b1a-4ecb-bc33-ed4863720a3c
run2 = M3DISRun(runpath2)

# ╔═╡ 4e6b0170-9b2f-4a08-b052-86ec0079c38f
run3 = M3DISRun(runpath3)

# ╔═╡ c237f390-1c38-4af5-9638-ff166bd7bb1f
let
	plt.close()
	f, ax = plt.subplots(figsize=(10, 5))

	x, y = MUST.flux(run, norm=true)
	ax.plot(x, y, lw=2, label="Kurucz")

	x, y = MUST.flux(run2, norm=true)
	ax.plot(x, y, lw=2, label="Bertrand")

	x, y = MUST.flux(run3, norm=true)
	ax.plot(x, y, lw=2, label="HITRAN")

	ax.legend()

	ax.set_xlim(3075, 3079)

	f
end

# ╔═╡ Cell order:
# ╠═1db442aa-a7e2-11ef-09a8-6db52d25eb54
# ╠═9d4b8d47-1ff1-47d9-9e7c-245766c982bb
# ╠═a2fadf4b-6333-47bb-8e85-6c221802f997
# ╠═5e5a74a4-d8b7-4629-9601-6183d79cec97
# ╠═d194fe19-5409-4d35-8acc-69f29eae4113
# ╠═1967dac1-2006-4be6-8727-9f5ba5caf095
# ╠═15a9b32f-dfc7-47c1-99da-4a1dfb268751
# ╠═2d2a18a7-b9be-4825-9567-88b81f4d0b99
# ╠═34afbb37-0b1a-4ecb-bc33-ed4863720a3c
# ╠═4e6b0170-9b2f-4a08-b052-86ec0079c38f
# ╠═c237f390-1c38-4af5-9638-ff166bd7bb1f
