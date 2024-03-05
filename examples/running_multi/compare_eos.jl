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
eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m0_a0_v1.4/combined_eos_magg_m0_a0.hdf5")

# ╔═╡ 54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m0_a0_v1.4/combined_opacities_magg_m0_a0.hdf5", mmap=true)

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

	c1 = eos1.lnEi
	c2 = eos2.lnEi
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = 40 #max(maximum(c1), maximum(c2))
	
	im = ax[0].scatter(tt1, rr1, label="M3D", c=c1, vmin=vmin, vmax=vmax)
	ax[1].scatter(tt2, rr2, label="TS", c=c2, vmin=vmin, vmax=vmax)

	f.colorbar(im, ax=ax)

	ax[0].legend()
	ax[1].legend()
	
	gcf()
end

# ╔═╡ e0f3be5e-80f3-47e6-b1a3-cdf48b412d20
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	it1 = findfirst(eos1.lnRho.>=log(1e-6))
	x1 = eos1.lnT
	y1 = eos1.lnEi[:, it1]
	
	it2 = findfirst(eos2.lnRho.>=log(1e-6))
	x2 = eos2.lnT
	y2 = eos2.lnEi[:, it2]
	
	plt.plot(x1, y1, label="M3D")
	plt.plot(x2, y2, label="TS")
	plt.axvline(log(15000))
	plt.xlim(7, 12)
	plt.ylim(28, 33)

	plt.ylabel("internal energy [ln]")
	plt.xlabel("temperature [ln]")
	plt.legend()
	
	gcf()
end

# ╔═╡ 4a29299b-0ed8-4617-b044-2a7e06afa105
ross = TSO.rosseland_opacity(eos1, opa1)

# ╔═╡ 55a6be86-5678-4d20-a8cb-3dc88e6a9a21
let
	plt.close()

	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	c1 = eos1.lnRoss
	c2 = eos2.lnRoss
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = max(maximum(c1), maximum(c2))
	
	xx, yy = TSO.meshgrid(@axed(eos1))
	ax[0].scatter(xx, yy, c=ross)

	xx, yy = TSO.meshgrid(@axed(eos2))
	im = ax[1].scatter(xx, yy, c=c2)
	plt.colorbar(im, ax=ax)
	
	gcf()
end

# ╔═╡ 382eece9-ca88-46c1-8936-871f96513358
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	it1 = findfirst(eos1.lnRho.>=log(1e-6))
	x1 = log10.(opa1.λ)
	y1 = log10.(opa1.κ[10, 10, :])
	
	it2 = findfirst(eos2.lnRho.>=log(1e-6))
	x2 = log10.(opa2.λ)
	y2 = log10.(opa2.κ[10, 10, :])
	
	plt.plot(x1, y1, label="M3D")
	plt.plot(x2, y2, label="TS")

	plt.legend()
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═789f821a-be1c-4fa3-8b75-b888f9eeed1f
# ╠═54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
# ╠═bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
# ╠═7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
# ╠═dc956c31-3c2c-4ee9-8b21-11a3411cf18b
# ╠═e0f3be5e-80f3-47e6-b1a3-cdf48b412d20
# ╠═4a29299b-0ed8-4617-b044-2a7e06afa105
# ╠═55a6be86-5678-4d20-a8cb-3dc88e6a9a21
# ╠═382eece9-ca88-46c1-8936-871f96513358
