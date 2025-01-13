### A Pluto.jl notebook ###
# v0.20.0

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
#ext1 = "magg_m0_a0"

# ╔═╡ 2013feaa-c461-450f-8365-982c205fa1df
ext1 = "magg_m0_a0_vmic1"

# ╔═╡ 6e6d4927-afd1-450c-ad02-a62802d143c3
ext2 = "magg_m0_a0"

# ╔═╡ 789f821a-be1c-4fa3-8b75-b888f9eeed1f
#eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext1)_v1.8/combined_eos_$(ext1).hdf5")

# ╔═╡ 54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
#opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext1)_v1.8/combined_opacities_$(ext1).hdf5", mmap=true)

# ╔═╡ eca1793a-a4eb-4937-bf03-75d3d6472329
eos1 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/$(ext1)_v1.2/combined_eos_$(ext1).hdf5")

# ╔═╡ c67235eb-69c2-49b6-afe7-5f9d8315e480
opa1 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/$(ext1)_v1.2/combined_opacities_$(ext1).hdf5", mmap=true)

# ╔═╡ bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
eos2 = reload(SqEoS, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext2)_v1.8/ross_combined_eos_$(ext2).hdf5")

# ╔═╡ 7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
opa2 = reload(SqOpacity, "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(ext2)_v1.8/combined_opacities_$(ext2).hdf5", mmap=true)

# ╔═╡ aecb8a87-e426-40c5-bf71-1e7e858e7bd9


# ╔═╡ 5f34751e-382d-4bb3-8c00-fb6c4a4eb8bd
solar_model1 = @optical Average3D(eos1, "sun_stagger.dat") eos1 opa1

# ╔═╡ d034e8f4-9882-4e36-bb2e-7cb4fee51305
solar_model2 = @optical Average3D(eos2, "sun_stagger.dat") eos2 opa2

# ╔═╡ 9cfd2420-2351-4d45-824a-87d71a9b297d


# ╔═╡ 0d29d381-89ab-4d3d-885d-05c1febacdf4
@show size(eos1)

# ╔═╡ 1f800fdb-e62e-4685-8f0d-b57803766fb9
@show minimum(eos1.lnT) maximum(eos1.lnT)

# ╔═╡ d878ea42-2a1a-495e-9bc3-f247cfe66e11
@show length(opa1.λ)

# ╔═╡ bf86a3cb-3679-4d46-bc21-a5a62e1923fb
opa1.λ

# ╔═╡ 2c72069d-f65f-44d6-b883-0157c0f57751
eos1.lnEi

# ╔═╡ 5b987b1e-837e-4d46-b874-39d5069c62e7


# ╔═╡ 22c6f6e1-3372-4cb4-8414-c87f284dfb65
@show size(eos2)

# ╔═╡ 625c758b-025c-4e32-b154-c67200a6ea8a
@show minimum(eos2.lnT) maximum(eos2.lnT)

# ╔═╡ 51dc4c98-450c-4670-9a22-fa66c29958f9
@show length(opa2.λ)

# ╔═╡ fd5cee02-9a36-4ad1-b668-8eced9602472
opa2.λ

# ╔═╡ 315cc5f0-6e5e-43cd-94b9-64f650a02792
md"# EoS"

# ╔═╡ dc956c31-3c2c-4ee9-8b21-11a3411cf18b
let
	plt.close()

	f, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=(10, 6))
	plt.subplots_adjust(wspace=0)
	
	tt1, rr1 = TSO.meshgrid(@axed eos1)
	tt2, rr2 = TSO.meshgrid(@axed eos2)

	c1 = eos1.lnRoss
	c2 = eos2.lnRoss
	
	vmin = min(minimum(c1), minimum(c2))
	vmax = max(maximum(c1), maximum(c2))#maximum(eos1.lnEi) #40 #max(maximum(c1), maximum(c2))
	
	im = ax[0].pcolormesh(log10.(exp.(tt1)), rr1, c1, vmin=vmin, vmax=vmax)
	ax[1].pcolormesh(log10.(exp.(tt2)), rr2, c2, vmin=vmin, vmax=vmax)

	ax[0].set_title("M3D")
	ax[1].set_title("MARCS")
	

	f.colorbar(im, ax=ax)

	ax[0].legend()
	ax[1].legend()
	
	gcf()
end

# ╔═╡ ae45577e-01a8-4976-8c09-e33534c6d438


# ╔═╡ c3313f20-74a1-4157-9bd4-ccca1ae7ad16
md"# Opacity"

# ╔═╡ 253d0706-732e-4c31-9c76-21db5fa1bedb
let
	plt.close()
	xx, yy = TSO.meshgrid(@axed(eos1))
	im = plt.pcolormesh(xx, yy, log10.(opa1.κ[:,:,1593]))

	plt.colorbar(im)

	gcf()
end

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

# ╔═╡ 556ff5af-f80d-451b-b78b-9ee9079ceb7e
isurf1.lnρ[1], isurf1.lnT[1]

# ╔═╡ dd51d0ff-2c7e-4164-95fb-d10fdb473d51
 isurf2.lnρ[1], isurf2.lnT[1]

# ╔═╡ b70728c9-d046-4db2-945c-015443c00b44
κ1 = lookup(eos1, opa1, :κ, isurf1.lnρ[1], isurf1.lnT[1])

# ╔═╡ fe9ac9ee-80be-4a47-b2e4-5e7a95662ab4
κ2 = lookup(eos2, opa2, :κ, isurf1.lnρ[1], isurf1.lnT[1])

# ╔═╡ fddc345a-94ed-4bfb-b78c-b316be4390f0
let
	plt.close()
	f, ax = plt.subplots(1, 1)

	ax.plot(opa1.λ, log10.(κ1), label="neu")
	ax.plot(opa2.λ, log10.(κ2))
	#ax.plot(opa1.λ, log10.(κ2))
	ax.set_xlim(0, 25000)
	#ax.set_xlim(3500, 18000)
	#ax.set_ylabel("log10 [C/Fe]=0.4 / [C/Fe]=3.0")
	ax.set_xlabel("log10 wavelength")

	ax.legend()
	
	gcf()
end

# ╔═╡ a40184cb-edae-4d48-8ec5-7a0c1c2318f7


# ╔═╡ 4f3435a1-3983-4d12-be24-edc67ac24265
κross1 = lookup(eos1, :lnRoss, solar_model1.lnρ, solar_model1.lnT)

# ╔═╡ add01598-f142-4035-9417-d52be0313a2e
κross2 = log.(lookup(eos2, opa2, :κ_ross, solar_model1.lnρ, solar_model1.lnT))

# ╔═╡ d54e6e19-78ea-4dd9-bc37-1524d5ae406f
let
	plt.close()

	plt.plot(solar_model1.z, exp.(solar_model1.lnEi), label="neu")
	plt.plot(solar_model1.z, exp.(solar_model2.lnEi))

	#plt.xlim(-4, 4)
	#plt.ylim(3000, 30000)
	plt.legend()
	gcf()
end

# ╔═╡ Cell order:
# ╠═0dbac1c4-d566-11ee-0ad0-cf2cf182e68a
# ╠═db28e099-4c6f-4786-94ad-14ab71767353
# ╠═2013feaa-c461-450f-8365-982c205fa1df
# ╠═6e6d4927-afd1-450c-ad02-a62802d143c3
# ╠═789f821a-be1c-4fa3-8b75-b888f9eeed1f
# ╠═54cd6e9e-c4e7-4512-bd5e-59c05840dfbf
# ╠═eca1793a-a4eb-4937-bf03-75d3d6472329
# ╠═c67235eb-69c2-49b6-afe7-5f9d8315e480
# ╠═bdf0a7ba-d30d-40f2-92f4-b9dd791a9ba9
# ╠═7cd8df75-c5e5-4ce3-b3af-cb54a7ec8cce
# ╟─aecb8a87-e426-40c5-bf71-1e7e858e7bd9
# ╠═5f34751e-382d-4bb3-8c00-fb6c4a4eb8bd
# ╠═d034e8f4-9882-4e36-bb2e-7cb4fee51305
# ╟─9cfd2420-2351-4d45-824a-87d71a9b297d
# ╠═0d29d381-89ab-4d3d-885d-05c1febacdf4
# ╠═1f800fdb-e62e-4685-8f0d-b57803766fb9
# ╠═d878ea42-2a1a-495e-9bc3-f247cfe66e11
# ╠═bf86a3cb-3679-4d46-bc21-a5a62e1923fb
# ╠═2c72069d-f65f-44d6-b883-0157c0f57751
# ╟─5b987b1e-837e-4d46-b874-39d5069c62e7
# ╠═22c6f6e1-3372-4cb4-8414-c87f284dfb65
# ╠═625c758b-025c-4e32-b154-c67200a6ea8a
# ╠═51dc4c98-450c-4670-9a22-fa66c29958f9
# ╠═fd5cee02-9a36-4ad1-b668-8eced9602472
# ╟─315cc5f0-6e5e-43cd-94b9-64f650a02792
# ╠═dc956c31-3c2c-4ee9-8b21-11a3411cf18b
# ╟─ae45577e-01a8-4976-8c09-e33534c6d438
# ╟─c3313f20-74a1-4157-9bd4-ccca1ae7ad16
# ╠═253d0706-732e-4c31-9c76-21db5fa1bedb
# ╠═d536627e-6355-491b-8955-f7012b545b5a
# ╠═2edee1f3-d5ea-45ad-9656-ffa34b32427e
# ╠═556ff5af-f80d-451b-b78b-9ee9079ceb7e
# ╠═dd51d0ff-2c7e-4164-95fb-d10fdb473d51
# ╠═b70728c9-d046-4db2-945c-015443c00b44
# ╠═fe9ac9ee-80be-4a47-b2e4-5e7a95662ab4
# ╠═fddc345a-94ed-4bfb-b78c-b316be4390f0
# ╟─a40184cb-edae-4d48-8ec5-7a0c1c2318f7
# ╠═4f3435a1-3983-4d12-be24-edc67ac24265
# ╠═add01598-f142-4035-9417-d52be0313a2e
# ╠═d54e6e19-78ea-4dd9-bc37-1524d5ae406f
