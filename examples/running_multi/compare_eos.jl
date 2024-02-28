### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end

# ╔═╡ 789f821a-be1c-4fa3-8b75-b888f9eeed1f
eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m0_a0_v1.0/combined_eos_magg_m0_a0.hdf5")

# ╔═╡ 54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m0_a0_v1.0/combined_opacities_magg_m0_a0.hdf5", mmap=true)

# ╔═╡ bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
eos2 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_magg_m0_a0_v1.9/combined_eos_magg_m0_a0.hdf5")

# ╔═╡ 7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
opa2 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_magg_m0_a0_v1.9/combined_opacities_magg_m0_a0.hdf5", mmap=true)

# ╔═╡ dc956c31-3c2c-4ee9-8b21-11a3411cf18b
let
	plt.close()

	f, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(10, 6))
	plt.subplots_adjust(wspace=0)
	
	tt1, rr1 = TSO.meshgrid(@axed eos1)
	tt2, rr2 = TSO.meshgrid(@axed eos2)

	c1 = log10.(opa1.κ[:, :, 1000])
	c2 = log10.(opa2.κ[:, :, 1000])
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = max(maximum(c1), maximum(c2))
	
	im = ax[0].scatter(tt1, rr1, label="M3D", c=c1, vmin=vmin, vmax=vmax)
	ax[1].scatter(tt2, rr2, label="TS", c=c2, vmin=vmin, vmax=vmax)

	f.colorbar(im, ax=ax)

	ax[0].legend()
	ax[1].legend()
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═789f821a-be1c-4fa3-8b75-b888f9eeed1f
# ╠═54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
# ╠═bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
# ╠═7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
# ╠═dc956c31-3c2c-4ee9-8b21-11a3411cf18b
