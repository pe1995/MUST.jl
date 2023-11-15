### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ a5b95bce-7f4d-11ee-2c07-573697832df0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end

# ╔═╡ 369e2b95-a234-470b-bf93-ec507ac92c90
md"# Check interpolated tables"

# ╔═╡ 44c1cda6-5b01-4664-93d3-ce79d24cfbfd
matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))

# ╔═╡ b7004306-3c97-4dd6-b828-06a9af03f9e5
md"## Grid"

# ╔═╡ eb3d90d9-2324-4a9a-8485-92185dbe669b
grid = MUST.StaggerGrid("random_grid.mgrid")
#grid = MUST.StaggerGrid("dispatch_grid.mgrid")

# ╔═╡ b5a02ce1-ea93-45ab-8d71-be17308f227b
md"## Model"

# ╔═╡ e0eabd3e-769c-4c09-960a-6aba011e5795
modelname = "t58g41.79m0.000"
#modelname = "t55g40m00"

# ╔═╡ bf257ee0-d69d-409a-bee5-3185a0490dc6
idxmodel = findfirst(x->x==modelname, grid["name"])

# ╔═╡ f4af070e-9f46-4e81-a7b9-55073ad87bb7
avmodelpath = grid["av_path", idxmodel]

# ╔═╡ 7ab675cb-7dbd-43f3-b985-7558b6575e64
logg = grid["logg", idxmodel]

# ╔═╡ 15b28abc-ff82-411a-9fc0-ba4efe98e0ea
md"## EoS"

# ╔═╡ 8873e5bd-ce33-42c6-9259-2d26f19959b7
eos = reload(
	SqEoS,
	joinpath(grid["binned_tables", idxmodel], "eos.hdf5")
)

# ╔═╡ b8f5ac9c-1f18-4fd9-9a4e-15b9011115ae
opa = reload(
	SqOpacity,
	joinpath(grid["binned_tables", idxmodel], "binned_opacities.hdf5")
)

# ╔═╡ 3167d0b6-199c-485a-8d9a-f99d43483d82
eosE = reload(
	SqEoS,
	joinpath(grid["binned_E_tables", idxmodel], "eos.hdf5")
)

# ╔═╡ 9162f70c-2d36-4bd2-a827-35b3574361db
opaE = reload(
	SqOpacity,
	joinpath(grid["binned_E_tables", idxmodel], "binned_opacities.hdf5")
)

# ╔═╡ d405f353-d9e7-477e-b915-b7287d000391
md"## Average Model"

# ╔═╡ 7b317e4c-0131-48d4-b07d-f679ac54b001
avmodel = @optical TSO.flip(
	Average3D(eos, avmodelpath, logg=logg), depth=true
) eos opa

# ╔═╡ a9b0e569-a30b-41c7-9de9-4e83e6bb8038
begin
	plt.close()

	plt.plot(avmodel.z, avmodel.lnT)
	
	plt.xlabel("z")
	plt.ylabel("lnT")
	
	gcf()
end

# ╔═╡ 8bb80f78-e97e-4c09-91b4-c5a93fa22abc
begin
	plt.close()

	plt.plot(avmodel.z, exp.(avmodel.lnEi))
	
	plt.xlabel("z")
	plt.ylabel("Ei")
	
	gcf()
end

# ╔═╡ 53ef0c55-3745-46c4-a426-ae07e22e40a7
begin
	plt.close()

	plt.plot(avmodel.z, lookup(eos, :lnPg, avmodel.lnρ, avmodel.lnT))
	
	plt.xlabel("z")
	plt.ylabel("lnPg")
	
	gcf()
end

# ╔═╡ b555e3c2-cdda-4c4a-9132-c152d2305425
begin
	plt.close()

	plt.plot(log10.(avmodel.τ), avmodel.lnT)
	
	plt.xlabel("logτ")
	plt.ylabel("lnT")
	
	gcf()
end

# ╔═╡ 12550aec-f780-44ae-a3f7-2179e73b8a78
md"## Opacities"

# ╔═╡ a410367d-c065-4f75-aab5-a77bbacdab14
begin
	eE = lookup(eosE, :lnEi, avmodel.lnρ, avmodel.lnT)
	kE = lookup(eosE, opaE, :κ, avmodel.lnρ, eE)
	SE = lookup(eosE, opaE, :src, avmodel.lnρ, eE)
end;

# ╔═╡ 450232c1-0f95-4804-b6d6-af4c25616e22
begin
	k = lookup(eos, opa, :κ, avmodel.lnρ, avmodel.lnT)
	S = lookup(eos, opa, :src, avmodel.lnρ, avmodel.lnT)
end;

# ╔═╡ db5523a2-ce16-422b-98f7-ee71579f73f5
cmap = plt.cm.get_cmap("RdYlBu")

# ╔═╡ 2f000a00-673b-4ee7-b992-87af05533d0f
begin
	plt.close()

	colorsO = [MUST.pyconvert(Any, cmap(i/size(k, 2))) for i in axes(k, 2)]

	for bin in axes(k, 2)
		plt.plot(
			exp.(avmodel.lnT) ./1000.0, (k .- kE) ./ kE .*100, 
			color=colorsO[bin]
		)
	end

	plt.ylim(-5, 5)

	plt.ylabel("κ - κₑ / κₑ [%]")
	plt.xlabel("T [kK]")

	gcf()
end

# ╔═╡ d30969f0-7f69-43ab-b755-bfbdb820020b
begin
	plt.close()

	for bin in axes(k, 2)
		plt.plot(
			exp.(avmodel.lnT) ./1000.0, log10.(k), 
			color=colorsO[bin]
		)
	end

	plt.ylabel("log κ")
	plt.xlabel("T [kK]")

	gcf()
end

# ╔═╡ 69f038de-7db5-46d0-a1c0-2d6f5d04d412
begin
	plt.close()

	for bin in axes(k, 2)
		plt.plot(
			exp.(avmodel.lnT) ./1000.0, (S .- SE) ./ SE .*100, 
			color=colorsO[bin]
		)
	end

	plt.ylim(-5, 5)

	plt.ylabel("S - Sₑ / Sₑ [%]")
	plt.xlabel("T [kK]")

	gcf()
end

# ╔═╡ 13684176-15aa-4511-b5b9-bf0e8dbf05e5
begin
	plt.close()

	for bin in axes(k, 2)
		plt.plot(
			exp.(avmodel.lnT) ./1000.0, log10.(S), 
			color=colorsO[bin]
		)
	end

	plt.ylabel("log S")
	plt.xlabel("T [kK]")

	gcf()
end

# ╔═╡ Cell order:
# ╟─369e2b95-a234-470b-bf93-ec507ac92c90
# ╠═a5b95bce-7f4d-11ee-2c07-573697832df0
# ╠═44c1cda6-5b01-4664-93d3-ce79d24cfbfd
# ╟─b7004306-3c97-4dd6-b828-06a9af03f9e5
# ╠═eb3d90d9-2324-4a9a-8485-92185dbe669b
# ╟─b5a02ce1-ea93-45ab-8d71-be17308f227b
# ╠═e0eabd3e-769c-4c09-960a-6aba011e5795
# ╠═bf257ee0-d69d-409a-bee5-3185a0490dc6
# ╠═f4af070e-9f46-4e81-a7b9-55073ad87bb7
# ╠═7ab675cb-7dbd-43f3-b985-7558b6575e64
# ╟─15b28abc-ff82-411a-9fc0-ba4efe98e0ea
# ╠═8873e5bd-ce33-42c6-9259-2d26f19959b7
# ╠═b8f5ac9c-1f18-4fd9-9a4e-15b9011115ae
# ╠═3167d0b6-199c-485a-8d9a-f99d43483d82
# ╠═9162f70c-2d36-4bd2-a827-35b3574361db
# ╟─d405f353-d9e7-477e-b915-b7287d000391
# ╠═7b317e4c-0131-48d4-b07d-f679ac54b001
# ╟─a9b0e569-a30b-41c7-9de9-4e83e6bb8038
# ╟─8bb80f78-e97e-4c09-91b4-c5a93fa22abc
# ╟─53ef0c55-3745-46c4-a426-ae07e22e40a7
# ╟─b555e3c2-cdda-4c4a-9132-c152d2305425
# ╟─12550aec-f780-44ae-a3f7-2179e73b8a78
# ╠═a410367d-c065-4f75-aab5-a77bbacdab14
# ╠═450232c1-0f95-4804-b6d6-af4c25616e22
# ╠═db5523a2-ce16-422b-98f7-ee71579f73f5
# ╟─2f000a00-673b-4ee7-b992-87af05533d0f
# ╟─d30969f0-7f69-43ab-b755-bfbdb820020b
# ╟─69f038de-7db5-46d0-a1c0-2d6f5d04d412
# ╟─13684176-15aa-4511-b5b9-bf0e8dbf05e5
