### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ db28e099-4c6f-4786-94ad-14ab71767353
ext1 = "magg_m1_a0_c0"

# ╔═╡ 6e6d4927-afd1-450c-ad02-a62802d143c3
ext2 = "magg_m1_a0_c1"

# ╔═╡ 789f821a-be1c-4fa3-8b75-b888f9eeed1f
eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext1)_v3.0/combined_eos_$(ext1).hdf5")

# ╔═╡ 5c113664-a163-4a54-9949-21634e50df55
scat1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext1)_v3.0/combined_sopacities_$(ext1).hdf5", mmap=true)

# ╔═╡ 54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext1)_v3.0/combined_opacities_$(ext1).hdf5", mmap=true)

# ╔═╡ bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
eos2 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext2)_v3.0/combined_eos_$(ext2).hdf5")

# ╔═╡ 7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
opa2 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext2)_v3.0/combined_opacities_$(ext2).hdf5", mmap=true)

# ╔═╡ 91c07140-e22f-4522-b2c2-8e2362946ddc
scat2 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext2)_v3.0/combined_sopacities_$(ext2).hdf5", mmap=true)

# ╔═╡ dc956c31-3c2c-4ee9-8b21-11a3411cf18b
let
	plt.close()

	f, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(10, 6))
	plt.subplots_adjust(wspace=0)
	
	tt1, rr1 = TSO.meshgrid(@axed eos1)
	tt2, rr2 = TSO.meshgrid(@axed eos2)

	c1 = eos1.lnEi
	c2 = eos2.lnEi
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = maximum(eos1.lnEi)#40 #max(maximum(c1), maximum(c2))
	
	im = ax[0].scatter(tt1, rr1, label="FeH=-1, CFe=0", c=c1, vmin=vmin, vmax=vmax)
	ax[1].scatter(tt2, rr2, label="FeH=-1, CFe=1", c=c2, vmin=vmin, vmax=vmax)

	f.colorbar(im, ax=ax)

	ax[0].legend()
	ax[1].legend()
	
	gcf()
end

# ╔═╡ e54494f1-fd6c-43ba-aa75-d55323b001af
let 
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	ax.plot(eos1.lnT, eos1.lnEi[:, end], label="M3D")
	ax.plot(eos2.lnT, eos2.lnEi[:, end], label="TS")
	
	gcf()
end

# ╔═╡ 7e077307-b759-4265-8d2f-50f3479b3842
minimum(eos1.lnEi), maximum(eos1.lnEi)

# ╔═╡ 013fb284-882b-4e07-a3e1-66afb1280fa0
minimum(eos2.lnEi), maximum(eos2.lnEi)

# ╔═╡ 6588a02b-c206-4aa2-b54f-836e1778e552
exp.(eos1.lnT)

# ╔═╡ 4a29299b-0ed8-4617-b044-2a7e06afa105
#ross = TSO.rosseland_opacity(eos1, opa1)

# ╔═╡ 55a6be86-5678-4d20-a8cb-3dc88e6a9a21
let
	plt.close()

	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	c1 = eos1.lnRoss
	c2 = eos2.lnRoss
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = max(maximum(c1), maximum(c2))
	
	xx, yy = TSO.meshgrid(@axed(eos1))
	ax[0].scatter(xx, yy, c=c1)

	xx, yy = TSO.meshgrid(@axed(eos2))
	im = ax[1].scatter(xx, yy, c=c2)
	plt.colorbar(im, ax=ax)
	
	gcf()
end

# ╔═╡ 65e6937c-221a-4db5-a4f3-e2bbc45470db
begin
	ρ = 1e-6
	T = 5000

	κ1 = lookup(eos1, opa1, :κ, log.(ρ), log.(T))
	κ2 = lookup(eos2, opa2, :κ, log.(ρ), log.(T))

	s1 = lookup(eos1, scat1, :κ, log.(ρ), log.(T))
	s2 = lookup(eos2, scat2, :κ, log.(ρ), log.(T))
end

# ╔═╡ 382eece9-ca88-46c1-8936-871f96513358
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x1 = log10.(opa1.λ)
	y1 = log10.(κ1)
	x2 = log10.(opa2.λ)
	y2 = log10.(κ2)
	plt.plot(x2, y2, label="FeH=-1, CFe=1")
	plt.plot(x1, y1, label="FeH=-1, CFe=0")
	
	x1 = log10.(opa1.λ)
	y1 = log10.(s1)
	x2 = log10.(opa2.λ)
	y2 = log10.(s2)
	plt.plot(x2, y2, label="FeH=-1, CFe=1")
	plt.plot(x1, y1, label="FeH=-1, CFe=0")

	#plt.xlim(4, 4.1)

	plt.legend()
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═db28e099-4c6f-4786-94ad-14ab71767353
# ╠═6e6d4927-afd1-450c-ad02-a62802d143c3
# ╠═789f821a-be1c-4fa3-8b75-b888f9eeed1f
# ╠═5c113664-a163-4a54-9949-21634e50df55
# ╠═54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
# ╠═bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
# ╠═7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
# ╠═91c07140-e22f-4522-b2c2-8e2362946ddc
# ╠═dc956c31-3c2c-4ee9-8b21-11a3411cf18b
# ╠═e54494f1-fd6c-43ba-aa75-d55323b001af
# ╠═7e077307-b759-4265-8d2f-50f3479b3842
# ╠═013fb284-882b-4e07-a3e1-66afb1280fa0
# ╠═6588a02b-c206-4aa2-b54f-836e1778e552
# ╠═4a29299b-0ed8-4617-b044-2a7e06afa105
# ╠═55a6be86-5678-4d20-a8cb-3dc88e6a9a21
# ╠═65e6937c-221a-4db5-a4f3-e2bbc45470db
# ╠═382eece9-ca88-46c1-8936-871f96513358
