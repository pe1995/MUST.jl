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
eospath = "interpolatedModels/M14.1_E_t57.50g45.00m-7.000_v1.0"

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

# ╔═╡ 5ebaf49b-93b4-44c4-9313-a24184d11f1d
let
	f, ax = plt.subplots()

	for i in 1:size(κ, 2)
		ax.plot(log10.(model.τ), log10.(κ[:, i]), label="bin $(i)")
	end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 2d978fc9-9083-4ce6-8d05-8c776e8c0f96


# ╔═╡ 746e9e94-1779-4797-b74f-70f7212b81a0
function smooth(eos::SqEoS)
	Δλ = diff(eos.lnEi)
	kernel = ImageFiltering.Kernel.gaussian((30,))
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
function smooth(opa::SqOpacity)
	kernel = ImageFiltering.Kernel.gaussian((1,))
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

# ╔═╡ 68ab1a38-eeeb-4db0-9c78-a1e5ce2481bb
eos_filter, opa_filter = smooth(eos),  smooth(opa)

# ╔═╡ 44ff77bf-f4b7-4da8-9a8d-656760c07378
κ_filter = lookup(eos_filter, opa_filter, :κ, model.lnρ, model.lnEi)

# ╔═╡ 722c4ec0-90cd-49a1-a4c2-bbfbc8300081
let
	f, ax = plt.subplots()

	#for i in 1:size(κ, 2)
		i=4
		ax.plot(log10.(model.τ), log10.(κ[:, i]), label="bin $(i)")
		ax.plot(log10.(model.τ), log10.(κ_filter[:, i]), label="bin $(i) filter")
	#end

	ax.legend(bbox_to_anchor=(1, 1), loc="upper left", framealpha=0)

	f
end

# ╔═╡ 409542fb-b5d3-4f51-9907-e21f701ff893


# ╔═╡ 5ee29b4c-9456-469a-ae0d-b224c5abedb4
TSO.save(eos_filter, joinpath(eospath, "eos_filter.hdf5"))

# ╔═╡ 5b1807dc-470d-4a4d-926f-18764bae76f1
TSO.save(opa_filter, joinpath(eospath, "binned_opacities_filter.hdf5"))

# ╔═╡ 5c62f473-1fdb-4e7b-9cf7-969408ad5239
TSO.for_dispatch(eos_filter, opa_filter, joinpath(eospath, "filter/"))

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
# ╠═5ebaf49b-93b4-44c4-9313-a24184d11f1d
# ╟─2d978fc9-9083-4ce6-8d05-8c776e8c0f96
# ╠═746e9e94-1779-4797-b74f-70f7212b81a0
# ╠═64b85589-1974-4361-9dc5-9b9f68fe1592
# ╠═68ab1a38-eeeb-4db0-9c78-a1e5ce2481bb
# ╠═44ff77bf-f4b7-4da8-9a8d-656760c07378
# ╠═722c4ec0-90cd-49a1-a4c2-bbfbc8300081
# ╟─409542fb-b5d3-4f51-9907-e21f701ff893
# ╠═5ee29b4c-9456-469a-ae0d-b224c5abedb4
# ╠═5b1807dc-470d-4a4d-926f-18764bae76f1
# ╠═5c62f473-1fdb-4e7b-9cf7-969408ad5239
