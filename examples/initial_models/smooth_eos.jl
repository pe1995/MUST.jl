### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ aaead4e4-65d7-11ef-3ed8-5fe769baccf0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using ImageFiltering
end

# ╔═╡ b12df0d2-0762-4868-bf5b-266c2efa4eb8
plt = matplotlib.pyplot

# ╔═╡ c43315b6-dda2-43a2-946b-6bdb5b814306
eospath = "MainSequenceInterpolated/ST10_E_t57.77g44.40m0.000_v1.0"

# ╔═╡ 11a17019-7728-406a-a170-a2c4a64283e2
eos = reload(SqEoS, joinpath(eospath, "eos.hdf5"))

# ╔═╡ 242d254a-507e-4178-9150-5cabf4884fe7
opa = reload(SqOpacity, joinpath(eospath, "binned_opacities.hdf5"))

# ╔═╡ 535c0014-49cb-4aae-9701-8a7b7b999cd4


# ╔═╡ 4982dff0-0b1d-404b-acfa-0b26bdc7a541
model = @optical Average3D(eos, joinpath(eospath, "inim.dat")) eos

# ╔═╡ a7f5f9c8-956d-462c-9c3f-9700efbfae02


# ╔═╡ cbd08511-3443-4f11-9f50-736181d54f7c
κ = lookup(eos, opa, :κ, model.lnρ, model.lnEi)

# ╔═╡ 8ef841d8-c8d6-4f02-9feb-e697c1bc5aed
s = lookup(eos, opa, :src, model.lnρ, model.lnEi)

# ╔═╡ 5ebaf49b-93b4-44c4-9313-a24184d11f1d
let
	f, ax = plt.subplots()

	for i in 1:size(κ, 2)
		ax.plot(log10.(model.τ), log10.(κ[:, i]), label="bin $(i)")
	end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 6ae9e043-3653-4c76-9afb-28b9d164b8e1
let
	f, ax = plt.subplots()

	for i in 1:size(κ, 2)
		ax.plot(log10.(model.τ), log10.(s[:, i]), label="bin $(i)")
	end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 2d978fc9-9083-4ce6-8d05-8c776e8c0f96


# ╔═╡ 746e9e94-1779-4797-b74f-70f7212b81a0
function smooth(eos::SqEoS, radius=30)
	Δλ = diff(eos.lnEi)
	kernel = ImageFiltering.Kernel.gaussian((radius,))
	#T_conv = ImageFiltering.imfilter(eos.lnT, kernel)
	T_conv = similar(eos.lnT)
	P_conv = similar(eos.lnT)
	N_conv = similar(eos.lnT)
	R_conv = similar(eos.lnT)
	
	for i in axes(T_conv, 2)
		x=eos.lnT[:, i]
		T_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnPg[:, i]
		P_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnNe[:, i]
		N_conv[:, i] .= ImageFiltering.imfilter(x, kernel)

		x=eos.lnRoss[:, i]
		R_conv[:, i] .= ImageFiltering.imfilter(x, kernel)
	end

	SqEoS(eos.lnRho, T_conv, eos.lnEi, P_conv, R_conv, N_conv)
end

# ╔═╡ 64b85589-1974-4361-9dc5-9b9f68fe1592
function smooth(opa::SqOpacity, radius=30)
	kernel = ImageFiltering.Kernel.gaussian((radius,))
	#T_conv = ImageFiltering.imfilter(eos.lnT, kernel)
	k_conv = similar(opa.κ)
	kr_conv = similar(opa.κ_ross)
	src_conv = similar(opa.κ)
	
	for l in axes(k_conv, 3)
		for i in axes(k_conv, 2)
			x=log.(opa.κ[:, i, l])
			k_conv[:, i, l] .= ImageFiltering.imfilter(x, kernel)
			
			x=log.(opa.src[:, i, l])
			src_conv[:, i, l] .= ImageFiltering.imfilter(x, kernel)
		end
	end
	for i in axes(k_conv, 2)
		x=log.(opa.κ_ross[:, i])
		kr_conv[:, i] .= ImageFiltering.imfilter(x, kernel)
	end
	
	SqOpacity(exp.(k_conv), exp.(kr_conv), exp.(src_conv), opa.λ, opa.optical_depth)
end

# ╔═╡ f741a936-3d78-42dc-8ddd-350487f635ec
function smooth2D(opa::SqOpacity, radius1=10, radius2=10)
	kernel = ImageFiltering.Kernel.gaussian((radius1, radius2))
	#T_conv = ImageFiltering.imfilter(eos.lnT, kernel)
	k_conv = similar(opa.κ)
	kr_conv = similar(opa.κ_ross)
	src_conv = similar(opa.κ)
	
	for l in axes(k_conv, 3)
		x=log.(opa.κ[:, :, l])
		k_conv[:, :, l] .= ImageFiltering.imfilter(x, kernel)
		
		x=log.(opa.src[:, :, l])
		src_conv[:, :, l] .= ImageFiltering.imfilter(x, kernel)
	end
	
	x=log.(opa.κ_ross)
	kr_conv[:, :] .= ImageFiltering.imfilter(x, kernel)
	
	
	SqOpacity(exp.(k_conv), exp.(kr_conv), exp.(src_conv), opa.λ, opa.optical_depth)
end

# ╔═╡ 972c7faa-0849-46d3-9d14-9dc68f9bfb99


# ╔═╡ 68ab1a38-eeeb-4db0-9c78-a1e5ce2481bb
eos_filter, opa_filter = smooth(eos, 40),  smooth(opa, 40)

# ╔═╡ 44ff77bf-f4b7-4da8-9a8d-656760c07378
κ_filter = lookup(eos_filter, opa_filter, :κ, model.lnρ, model.lnEi)

# ╔═╡ 722c4ec0-90cd-49a1-a4c2-bbfbc8300081
let
	f, ax = plt.subplots()

	#for i in 1:size(κ, 2)
		i=3
		ax.plot(log10.(model.τ), log10.(κ[:, i]), label="bin $(i)")
		ax.plot(log10.(model.τ), log10.(κ_filter[:, i]), label="bin $(i) filter")
	#end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 10f3c03e-ab8f-4805-96f8-10f0798cb453


# ╔═╡ 96100d99-523b-445e-8b09-22fded071a9e
lnE, ρ = TSO.meshgrid(@axed(eos_filter));

# ╔═╡ 58b9636e-388b-4c71-9017-712dfa34c2ba
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnE, log10.(exp.(ρ)), log10.(opa_filter.κ[:, :, 3]), rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	cbar.set_label("opacity (bin $(3))")
	ax.set_xlabel("ln internal energy")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ 317e0995-e620-482f-8731-342efddd82cd
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnE, log10.(exp.(ρ)), eos_filter.lnT[:, :], rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	cbar.set_label("temperature")
	ax.set_xlabel("ln internal energy")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ 67de20cc-343f-494d-b499-176f5ecf372b


# ╔═╡ fa554409-dae1-4f63-8348-71a5908dc66e
let
	f, ax = plt.subplots()

	for i in 1:size(κ, 2)
		ax.plot(log10.(model.τ), log10.(κ_filter[:, i]), label="bin $(i)")
	end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 409542fb-b5d3-4f51-9907-e21f701ff893


# ╔═╡ 5ee29b4c-9456-469a-ae0d-b224c5abedb4
#TSO.save(eos_filter, joinpath(eospath, "eos_filter.hdf5"))

# ╔═╡ 5b1807dc-470d-4a4d-926f-18764bae76f1
#TSO.save(opa_filter, joinpath(eospath, "binned_opacities_filter.hdf5"))

# ╔═╡ 5c62f473-1fdb-4e7b-9cf7-969408ad5239
#TSO.for_dispatch(eos_filter, opa_filter, joinpath(eospath, "filter/"))

# ╔═╡ Cell order:
# ╠═aaead4e4-65d7-11ef-3ed8-5fe769baccf0
# ╠═b12df0d2-0762-4868-bf5b-266c2efa4eb8
# ╠═c43315b6-dda2-43a2-946b-6bdb5b814306
# ╠═11a17019-7728-406a-a170-a2c4a64283e2
# ╠═242d254a-507e-4178-9150-5cabf4884fe7
# ╟─535c0014-49cb-4aae-9701-8a7b7b999cd4
# ╠═4982dff0-0b1d-404b-acfa-0b26bdc7a541
# ╟─a7f5f9c8-956d-462c-9c3f-9700efbfae02
# ╠═cbd08511-3443-4f11-9f50-736181d54f7c
# ╠═8ef841d8-c8d6-4f02-9feb-e697c1bc5aed
# ╟─5ebaf49b-93b4-44c4-9313-a24184d11f1d
# ╟─6ae9e043-3653-4c76-9afb-28b9d164b8e1
# ╟─2d978fc9-9083-4ce6-8d05-8c776e8c0f96
# ╠═746e9e94-1779-4797-b74f-70f7212b81a0
# ╠═64b85589-1974-4361-9dc5-9b9f68fe1592
# ╠═f741a936-3d78-42dc-8ddd-350487f635ec
# ╟─972c7faa-0849-46d3-9d14-9dc68f9bfb99
# ╠═68ab1a38-eeeb-4db0-9c78-a1e5ce2481bb
# ╠═44ff77bf-f4b7-4da8-9a8d-656760c07378
# ╟─722c4ec0-90cd-49a1-a4c2-bbfbc8300081
# ╟─10f3c03e-ab8f-4805-96f8-10f0798cb453
# ╠═96100d99-523b-445e-8b09-22fded071a9e
# ╟─58b9636e-388b-4c71-9017-712dfa34c2ba
# ╟─317e0995-e620-482f-8731-342efddd82cd
# ╟─67de20cc-343f-494d-b499-176f5ecf372b
# ╟─fa554409-dae1-4f63-8348-71a5908dc66e
# ╟─409542fb-b5d3-4f51-9907-e21f701ff893
# ╠═5ee29b4c-9456-469a-ae0d-b224c5abedb4
# ╠═5b1807dc-470d-4a4d-926f-18764bae76f1
# ╠═5c62f473-1fdb-4e7b-9cf7-969408ad5239
