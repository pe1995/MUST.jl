### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 50f25934-7efd-11ef-1cbb-b5dc63b61831
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end;

# ╔═╡ a8d2b460-951c-43c4-b27f-a708efe44b5c
@import_dispatch "../../../dispatch2"

# ╔═╡ 05dc1e78-efe2-4a43-8370-3b5ccc1d09b2


# ╔═╡ e4d4ef31-28a1-4c22-8d37-f850e127007a
md"""
# Problem: Model Artefacts
At the uppermost layers of DISPATCH atmospheres there are checkboard pattern artefacts, that are visible if the density drops below a certain value (sun: 1e-11 g/cm^3). Question: Are these caused by jumps in the EoS/opacity tables?
"""

# ╔═╡ 3d75b32a-38d6-4939-ab69-e31b97ad78e7


# ╔═╡ 51d78eed-96e5-4d8d-a8f8-d55d49290605
md"We first load a model with artefacts to see what we are dealing with"

# ╔═╡ 5838afbf-9aaa-4b6a-b00f-2f8d622f8764
b, bτ = pick_snapshot(@in_dispatch("data/ST1_E_t57.77g44.40m0.000_v1.0_fast2"), 3841)

# ╔═╡ fe7bfcd6-0653-4722-ab81-9078f51e72c2


# ╔═╡ 553d8c87-1250-4d9f-9186-aabbc57592cd
md"And the corresponding EoS + opacity table"

# ╔═╡ 920b14be-dc82-4a23-b93f-67efc7ba9068
begin
	eosdir = "input_data/grd/MainSequenceInterpolated/ST9_E_t57.77g44.40m0.000_v1.0"
	eospath = @in_dispatch(joinpath(eosdir, "eos.hdf5"))
	opapath = @in_dispatch(joinpath(eosdir, "binned_opacities.hdf5"))
	
	eos = reload(SqEoS, eospath)
	opa = reload(SqOpacity, opapath)
end

# ╔═╡ 73d8b920-395a-4edd-9854-eea2746b37a0
begin
	eospathT = @in_dispatch(joinpath(eosdir, "eos_T.hdf5"))
	opapathT = @in_dispatch(joinpath(eosdir, "binned_opacities_T.hdf5"))
	
	eosT = reload(SqEoS, eospathT)
	opaT = reload(SqOpacity, opapathT)
end

# ╔═╡ 5371fdf8-f55c-41da-bbb4-37c197417390


# ╔═╡ a7e25766-f9d4-420a-84ec-dccc8fa0e81e
md"To see the problem we look at the upper layers"

# ╔═╡ 8eaca024-9963-4f17-a09c-8606b2281336
let
	f, ax = plt.subplots(1, 1)

	im = ax.imshow(b[:T][:,:,end], cmap="rainbow")
	cbar = f.colorbar(im, ax=ax)
	cbar.set_label("temperature [K]")

	ax.set_xlabel("pixel")
	ax.set_ylabel("pixel")
	
	f
end

# ╔═╡ 7ad1cbde-c799-4dd5-9521-524d3ca25f03
let
	f, ax = plt.subplots(1, 1)

	im = ax.imshow(log10.(b[:d][:,:,end]), cmap="rainbow")
	cbar = f.colorbar(im, ax=ax)
	cbar.set_label("density [g/cm-3]")

	ax.set_xlabel("pixel")
	ax.set_ylabel("pixel")
	
	f
end

# ╔═╡ 528276dc-79fc-4fd5-9399-0d2549af658b


# ╔═╡ 1092219e-8215-4e35-968c-b905926b9e9b
md"Lets look at the opacity + source function an a given bin"

# ╔═╡ 8b597929-9dec-43d9-aece-7d669fd370be
bin = 7

# ╔═╡ 07b2672e-9893-4d56-a3c2-19927fdb6e4c


# ╔═╡ d7a7464a-c70f-46ea-adc2-28f2b79717a3
lnE, ρ = TSO.meshgrid(@axed(eos));

# ╔═╡ 1422e93f-ea8c-4aec-9540-c4382a20688e
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnE, ρ, log10.(opa.κ[:, :, bin]), rasterized=true)

	eb, rb = profile(MUST.mean, b, :logee, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logee, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logee, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	
	cbar = f.colorbar(im, ax=ax)

	cbar.set_label("opacity (bin $(bin))")
	ax.set_xlabel("ln internal energy")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ d287adbc-7415-42f8-9057-1a74bfb737c8


# ╔═╡ 6ef3c45a-0346-4170-b85f-c26c281e07a7
lnT, ρT = TSO.meshgrid(@axed(eosT));

# ╔═╡ dec02b9d-1072-4fa6-8e3c-4c6d1d711051
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnT, ρT, log10.(opaT.κ[:, :, bin]), rasterized=true)

	eb = log.(b[:T][:,:,end])
	rb = log.(b[:d][:,:,end])
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	
	cbar = f.colorbar(im, ax=ax)

	cbar.set_label("opacity (bin $(bin))")
	ax.set_xlabel("ln temperature")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ a5941f53-3ba0-4bec-87e1-86a4ad249147
@show minimum(eosT.lnEi)

# ╔═╡ a276d38f-a74b-4116-a251-629d9934a6b8
@show maximum(eosT.lnEi)

# ╔═╡ 217b0b1b-ab2f-470e-88fc-434ff70fe844
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnT, ρT, eosT.lnEi, rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	#eb = log.(b[:T][:,:,end])
	#rb = log.(b[:d][:,:,end])
	eb, rb = profile(MUST.mean, b, :logT, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logT, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logT, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	
	cbar.set_label("internal energy")
	ax.set_xlabel("ln temperature")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ af9a7bd0-f9a3-4987-a22b-196db027ca84
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnT, ρT, log10.(opaT.src[:,:,bin]), rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	#eb = log.(b[:T][:,:,end])
	#rb = log.(b[:d][:,:,end])
	eb, rb = profile(MUST.mean, b, :logT, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logT, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logT, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	
	cbar.set_label("source function (bin $(bin))")
	ax.set_xlabel("ln temperature")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ 2cda99b3-6679-4fc8-b2dc-982025a2c714


# ╔═╡ 20d4b759-f7ad-4bd9-8359-cbc7acc16c91


# ╔═╡ 4ca429e0-70d2-4111-bf23-edc4e4200ff0
begin
	eosT_limit = deepcopy(eosT)
	lgEE = log.(b[:ee])
	E_limit = maximum(lgEE) + (maximum(lgEE) - minimum(lgEE)) *0.2

	mask = eosT_limit.lnEi .> E_limit
	@show count(mask)/length(mask) minimum(lgEE) maximum(lgEE)
	eosT_limit.lnEi[mask] .= E_limit
end

# ╔═╡ 736d1ad8-88ac-4ee5-9879-f353c5cd681d
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnT, ρT, eosT_limit.lnEi, rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	#eb = log.(b[:T][:,:,end])
	#rb = log.(b[:d][:,:,end])
	eb, rb = profile(MUST.mean, b, :logT, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logT, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logT, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	
	cbar.set_label("internal energy")
	ax.set_xlabel("ln temperature")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ e2b9d760-1906-4cce-a8a4-208f5fcb207c


# ╔═╡ a48f4088-8372-4db3-938c-05f102e042fc
eos_test, opa_test = TSO.switch_energy(eosT_limit, opaT, upsample=2048, conservative=false)

# ╔═╡ 753a40e2-5edb-4434-92e1-530c8b1e3d22


# ╔═╡ 0a2b57f9-ece6-4a6d-b5dc-cc3dbc5f4478
lnE_test, ρ_test = TSO.meshgrid(@axed(eos_test));

# ╔═╡ 6e389367-eb59-4fee-9ccc-76277c91d714
let
	f, ax = plt.subplots(1, 1)

	#im = ax.pcolormesh(lnE_test, ρ_test, log10.(opa_test.κ[:, :, bin]), rasterized=true)
	#cbar = f.colorbar(im, ax=ax)
	#cbar.set_label("opacity (bin $(bin))")


	im2 = ax.contourf(lnE_test, ρ_test, log10.(opa_test.κ[:, :, bin]), cmap="rainbow", levels=15, alpha=1)
	cbar2 = f.colorbar(im2, ax=ax)
	cbar2.set_label("opacity (bin $(bin))")

	eb, rb = profile(MUST.mean, b, :logee, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logee, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logee, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	

	ax.set_xlabel("ln internal energy")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ 2eddd579-27c4-4f70-91ed-5605fab741a4
let
	f, ax = plt.subplots(1, 1)

	im = ax.pcolormesh(lnE_test, ρ_test, log10.(opa_test.src[:, :, bin]), rasterized=true)
	cbar = f.colorbar(im, ax=ax)

	#im2 = ax.contourf(lnE_test, ρ_test, log10.(opa_test.src[:, :, bin]), cmap="rainbow", levels=15, alpha=1)
	#cbar2 = f.colorbar(im2, ax=ax)

	cbar.set_label("src (bin $(bin))")

	eb, rb = profile(MUST.mean, b, :logee, :logd)
	ax.scatter(eb, rb, color="white", marker="x", s=10)

	eb, rb = profile(MUST.maximum, b, :logee, :logd)
	ax.scatter(eb, rb, color="cyan", marker="x", s=10)

	eb, rb = profile(MUST.minimum, b, :logee, :logd)
	ax.scatter(eb, rb, color="red", marker="x", s=10)
	

	ax.set_xlabel("ln internal energy")
	ax.set_ylabel("ln density")
	
	f
end

# ╔═╡ Cell order:
# ╠═50f25934-7efd-11ef-1cbb-b5dc63b61831
# ╠═a8d2b460-951c-43c4-b27f-a708efe44b5c
# ╟─05dc1e78-efe2-4a43-8370-3b5ccc1d09b2
# ╟─e4d4ef31-28a1-4c22-8d37-f850e127007a
# ╟─3d75b32a-38d6-4939-ab69-e31b97ad78e7
# ╟─51d78eed-96e5-4d8d-a8f8-d55d49290605
# ╠═5838afbf-9aaa-4b6a-b00f-2f8d622f8764
# ╟─fe7bfcd6-0653-4722-ab81-9078f51e72c2
# ╟─553d8c87-1250-4d9f-9186-aabbc57592cd
# ╠═920b14be-dc82-4a23-b93f-67efc7ba9068
# ╠═73d8b920-395a-4edd-9854-eea2746b37a0
# ╟─5371fdf8-f55c-41da-bbb4-37c197417390
# ╟─a7e25766-f9d4-420a-84ec-dccc8fa0e81e
# ╟─8eaca024-9963-4f17-a09c-8606b2281336
# ╟─7ad1cbde-c799-4dd5-9521-524d3ca25f03
# ╟─528276dc-79fc-4fd5-9399-0d2549af658b
# ╟─1092219e-8215-4e35-968c-b905926b9e9b
# ╠═8b597929-9dec-43d9-aece-7d669fd370be
# ╟─07b2672e-9893-4d56-a3c2-19927fdb6e4c
# ╠═d7a7464a-c70f-46ea-adc2-28f2b79717a3
# ╟─1422e93f-ea8c-4aec-9540-c4382a20688e
# ╟─d287adbc-7415-42f8-9057-1a74bfb737c8
# ╠═6ef3c45a-0346-4170-b85f-c26c281e07a7
# ╟─dec02b9d-1072-4fa6-8e3c-4c6d1d711051
# ╠═a5941f53-3ba0-4bec-87e1-86a4ad249147
# ╠═a276d38f-a74b-4116-a251-629d9934a6b8
# ╟─217b0b1b-ab2f-470e-88fc-434ff70fe844
# ╟─af9a7bd0-f9a3-4987-a22b-196db027ca84
# ╟─2cda99b3-6679-4fc8-b2dc-982025a2c714
# ╟─20d4b759-f7ad-4bd9-8359-cbc7acc16c91
# ╠═4ca429e0-70d2-4111-bf23-edc4e4200ff0
# ╟─736d1ad8-88ac-4ee5-9879-f353c5cd681d
# ╟─e2b9d760-1906-4cce-a8a4-208f5fcb207c
# ╠═a48f4088-8372-4db3-938c-05f102e042fc
# ╟─753a40e2-5edb-4434-92e1-530c8b1e3d22
# ╠═0a2b57f9-ece6-4a6d-b5dc-cc3dbc5f4478
# ╟─6e389367-eb59-4fee-9ccc-76277c91d714
# ╟─2eddd579-27c4-4f70-91ed-5605fab741a4
