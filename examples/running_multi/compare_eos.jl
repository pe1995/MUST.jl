### A Pluto.jl notebook ###
# v0.19.41

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
ext1 = "magg_m0_a0"

# ╔═╡ 6e6d4927-afd1-450c-ad02-a62802d143c3
ext2 = "magg_m0_a0_vmic1"

# ╔═╡ 789f821a-be1c-4fa3-8b75-b888f9eeed1f
eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext1)_v1.8/combined_eos_$(ext1).hdf5")

# ╔═╡ 54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext1)_v1.8/combined_opacities_$(ext1).hdf5", mmap=true)

# ╔═╡ bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
eos2 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext2)_v5.1/combined_eos_$(ext2).hdf5")

# ╔═╡ 7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
opa2 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(ext2)_v5.1/combined_opacities_$(ext2).hdf5", mmap=true)

# ╔═╡ 5f34751e-382d-4bb3-8c00-fb6c4a4eb8bd
solar_model1 = @optical Average3D(eos1, "sun_stagger.dat") eos1 opa1

# ╔═╡ d034e8f4-9882-4e36-bb2e-7cb4fee51305
solar_model2 = @optical Average3D(eos2, "sun_stagger.dat") eos2 opa2

# ╔═╡ 9cfd2420-2351-4d45-824a-87d71a9b297d


# ╔═╡ 315cc5f0-6e5e-43cd-94b9-64f650a02792
md"# EoS"

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
	
	im = ax[0].pcolormesh(tt1, rr1, c1, vmin=vmin, vmax=vmax)
	ax[1].pcolormesh(tt2, rr2, c2, vmin=vmin, vmax=vmax)

	ax[0].set_title("MARCS")
	ax[1].set_title("M3D")
	

	f.colorbar(im, ax=ax)

	ax[0].legend()
	ax[1].legend()
	
	gcf()
end

# ╔═╡ ae45577e-01a8-4976-8c09-e33534c6d438


# ╔═╡ c3313f20-74a1-4157-9bd4-ccca1ae7ad16
md"# Opacity"

# ╔═╡ d536627e-6355-491b-8955-f7012b545b5a
begin
	optical_surface1 = TSO.optical_surface(solar_model1)
	isurf1 = TSO.interpolate_to(solar_model1, z=[optical_surface1], in_log=false)
end

# ╔═╡ 2edee1f3-d5ea-45ad-9656-ffa34b32427e
begin
	optical_surface2 = TSO.optical_surface(solar_model2)
	isurf2 = TSO.interpolate_to(solar_model2, z=[optical_surface2], in_log=false)
end

# ╔═╡ b70728c9-d046-4db2-945c-015443c00b44
κ1 = lookup(eos1, opa1, :κ, isurf1.lnρ[1], isurf1.lnT[1])

# ╔═╡ fe9ac9ee-80be-4a47-b2e4-5e7a95662ab4
κ2 = lookup(eos2, opa2, :κ, isurf2.lnρ[1], isurf2.lnT[1])

# ╔═╡ fddc345a-94ed-4bfb-b78c-b316be4390f0
let
	plt.close()
	f, ax = plt.subplots(1, 1)

	ax.plot(log10.(opa1.λ), log10.(κ1), label="MARCS")
	ax.plot(log10.(opa2.λ), log10.(κ2), label="M3D")

	ax.legend()

	gcf()
end

# ╔═╡ Cell order:
# ╠═0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═db28e099-4c6f-4786-94ad-14ab71767353
# ╠═6e6d4927-afd1-450c-ad02-a62802d143c3
# ╠═789f821a-be1c-4fa3-8b75-b888f9eeed1f
# ╠═54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
# ╠═bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
# ╠═7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
# ╠═5f34751e-382d-4bb3-8c00-fb6c4a4eb8bd
# ╠═d034e8f4-9882-4e36-bb2e-7cb4fee51305
# ╟─9cfd2420-2351-4d45-824a-87d71a9b297d
# ╟─315cc5f0-6e5e-43cd-94b9-64f650a02792
# ╟─dc956c31-3c2c-4ee9-8b21-11a3411cf18b
# ╟─ae45577e-01a8-4976-8c09-e33534c6d438
# ╟─c3313f20-74a1-4157-9bd4-ccca1ae7ad16
# ╠═d536627e-6355-491b-8955-f7012b545b5a
# ╠═2edee1f3-d5ea-45ad-9656-ffa34b32427e
# ╠═b70728c9-d046-4db2-945c-015443c00b44
# ╠═fe9ac9ee-80be-4a47-b2e4-5e7a95662ab4
# ╟─fddc345a-94ed-4bfb-b78c-b316be4390f0
