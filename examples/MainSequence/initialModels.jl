### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 9983e4ec-5633-11ef-20fb-e5af68ee6ad5
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST 
	using TSO
	using PythonPlot  
end

# ╔═╡ 46147df6-e44e-4365-a1b3-cca46be052d3
begin
	MUST.@import_dispatch("../../../dispatch2/")
	plt = matplotlib.pyplot
end

# ╔═╡ 173623d0-c86f-4132-a4ad-a04dfb498a14
#stag2 = MUST.Atmos1DGrid("../initial_models/stagger_grid_full_o.mgrid") 

# ╔═╡ c0c9fa86-1f0f-42b3-b977-003bcadf254d
#stag = deepcopy(stag2) 

# ╔═╡ a48aeb43-85b6-4a4b-8395-4044fd1528bb
#MUST.absolute_path!(stag, from=stag.path) 

# ╔═╡ 90ae8ccc-760b-4210-b1b0-8030fde53ae5
#oproot = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/" 

# ╔═╡ 2cf9c2e5-a437-4d88-bc17-25315d2cb827
#eos_name(ext) = joinpath(oproot, "TSO_M3D_$(ext)_v5.0/combined_eos_$(ext).hdf5") 

# ╔═╡ 6a48ba99-242f-4ff4-a5d9-cfcdcbe89b45
#=allEoS = Dict(
	0=>eos_name("magg_m0_a0_vmic1"),
	-1=>eos_name("magg_m1_a0_vmic1"),
	-2=>eos_name("magg_m2_a0_vmic1"),
	-3=>eos_name("magg_m3_a0_vmic1"),
	-4=>eos_name("magg_m4_a0_vmic1") 
)=#

# ╔═╡ f903125c-d68c-4855-ae34-9adc3ae3c439
#=closest_eos(feh) = begin
	k = sort(keys(allEoS) |> collect)
	ik = k[argmin(abs.(k .- feh))]
	allEoS[ik] 
end=#

# ╔═╡ fe8d5cd9-1557-4db5-976b-5c92fd27ad53
#eos_stag = [closest_eos(feh) for feh in stag.info[!, "feh"]] 

# ╔═╡ 80a83e08-1049-458b-8019-9a217830458c
#stag.info[!, "matching_eos"] = eos_stag 

# ╔═╡ 6b8c8e85-933d-4bdc-98d9-4d8f8d150c81
#p = abspath(joinpath(abspath(pathof(MUST)), "../../initial_grids/Stagger/av_models"))

# ╔═╡ 25b1ae22-bc9d-4da5-b391-8ff54252910a
# copy av models to new location
#=for (i,path) in enumerate(stag.info[!, "av_path"])
	newpath = joinpath(p, basename(path))
	cp(path, newpath, force=true)
	stag.info[i, "av_path"] = newpath

	path = stag.info[i, "avo_path"]
	newpath = joinpath(p, basename(path))
	cp(path, newpath, force=true)
	stag.info[i, "avo_path"] = newpath
end=#

# ╔═╡ be9343ea-8372-44a2-a6a6-14e2d3313d8f
#MUST.relative_path!(stag) 

# ╔═╡ bb31ca87-20bb-4f80-a7f1-0524550db347
#stag.info 

# ╔═╡ 8f77314f-0ee4-41c0-a069-8448c7ac896b
#MUST.save(stag, abspath(joinpath(p, "../stagger_grid.mgrid"))) 

# ╔═╡ 600d3d9f-6a44-4912-a5c4-b96cdee2e0f8


# ╔═╡ a2eb1ff1-0d79-433a-80d6-14d58c8389eb
#stag_test = deepcopy(stag)

# ╔═╡ 9b0ce20e-7d20-4bec-bcc1-1077f7eb9e06
#MUST.absolute_path!(stag_test)

# ╔═╡ 0c474690-315a-4c5e-aec8-48ab2695f0fc
#stag_test.info

# ╔═╡ 897b9956-eef9-4472-86b2-59080a5f58e8


# ╔═╡ 0322a595-8f95-47ae-80d1-20060b4934c1


# ╔═╡ bd48a88b-48ac-458e-b9ad-7130885354c5
md"# Generating initial conditions"

# ╔═╡ 5cb85f14-fee8-4aad-a510-9ec88355ff88
iniCond = MUST.ingredients("initial_condition.jl")

# ╔═╡ 2cf223a5-073a-4e98-936a-8d33607cb773
begin
	modelFolder = "test1/S1"
	avModels = joinpath(modelFolder, "av_models")
end

# ╔═╡ 9bfce8b2-0814-4ba6-9a77-633f4fc2624e
!isdir(avModels) && mkpath(avModels)

# ╔═╡ f8de92ce-b40b-4aa2-af73-81f8301ab65a
begin
	paras = zeros(2, 3)
	paras[:, 1] = [5777.0, 6400.0]
	paras[:, 2] = [4.44, 4.0]
	paras[:, 3] = [0.0, -4.0]
end

# ╔═╡ 65cd3f65-d34c-4971-929c-c4f1b4b7d1ec
grid = iniCond.initialmodel(
	paras;
	eos=:closest, 
	savedir=avModels, 
	τbottom=2, 
	adiabatic_extrapolation=true
)

# ╔═╡ e79d89d6-87da-49a9-87d4-0924aa19f923
grid.info

# ╔═╡ cc8f4815-db8e-4f3c-abd1-dd5c6d6e0587
eos = reload.(SqEoS, grid.info[!, "matching_eos"])

# ╔═╡ f5081fff-51d1-44c0-8c51-bccb7c08f996
models = [
	@optical(Average3D(e, name), e) 
	for (name, e) in zip(grid.info[!, "av_path"], eos)
]

# ╔═╡ 0558677f-f998-4699-8b6d-00ef7a011b2f
let
	plt.close()
	f, ax = plt.subplots(1, 1)

	for i in eachindex(models)
		model = models[i]
		t = grid.info[i, "teff"]
		g = grid.info[i, "logg"]
		f = grid.info[i, "feh"]
		ax.plot(model.z, exp.(model.lnT), label="$t, $g, $f")
		#ax.plot(log10.(model.τ), exp.(model.lnT), label="$t, $g, $f")
	end

	ax.legend()
	
	gcf()
end

# ╔═╡ f1bf5692-0242-4086-9147-520055089c6d


# ╔═╡ e6b706f0-2107-4952-9b45-068ae437dfeb
MUST.save(grid, modelFolder*"_120824.mgrid")

# ╔═╡ 92ce4f62-384b-4c19-9ce2-e451abfd81ac


# ╔═╡ 834720f5-0cef-437d-b09f-939d1d763167
md"# Binning opacities"

# ╔═╡ c68896a8-6aab-4600-8a5c-52fbf9cd3651
iniCond.prepare(
	grid,
	name_extension=modelFolder
)

# ╔═╡ bfacb6e6-b445-47f7-bef1-166b67bc8eee


# ╔═╡ Cell order:
# ╠═9983e4ec-5633-11ef-20fb-e5af68ee6ad5
# ╠═46147df6-e44e-4365-a1b3-cca46be052d3
# ╠═173623d0-c86f-4132-a4ad-a04dfb498a14
# ╠═c0c9fa86-1f0f-42b3-b977-003bcadf254d
# ╠═a48aeb43-85b6-4a4b-8395-4044fd1528bb
# ╠═90ae8ccc-760b-4210-b1b0-8030fde53ae5
# ╠═2cf9c2e5-a437-4d88-bc17-25315d2cb827
# ╠═6a48ba99-242f-4ff4-a5d9-cfcdcbe89b45
# ╠═f903125c-d68c-4855-ae34-9adc3ae3c439
# ╠═fe8d5cd9-1557-4db5-976b-5c92fd27ad53
# ╠═80a83e08-1049-458b-8019-9a217830458c
# ╠═6b8c8e85-933d-4bdc-98d9-4d8f8d150c81
# ╠═25b1ae22-bc9d-4da5-b391-8ff54252910a
# ╠═be9343ea-8372-44a2-a6a6-14e2d3313d8f
# ╠═bb31ca87-20bb-4f80-a7f1-0524550db347
# ╠═8f77314f-0ee4-41c0-a069-8448c7ac896b
# ╟─600d3d9f-6a44-4912-a5c4-b96cdee2e0f8
# ╠═a2eb1ff1-0d79-433a-80d6-14d58c8389eb
# ╠═9b0ce20e-7d20-4bec-bcc1-1077f7eb9e06
# ╠═0c474690-315a-4c5e-aec8-48ab2695f0fc
# ╟─897b9956-eef9-4472-86b2-59080a5f58e8
# ╟─0322a595-8f95-47ae-80d1-20060b4934c1
# ╟─bd48a88b-48ac-458e-b9ad-7130885354c5
# ╠═5cb85f14-fee8-4aad-a510-9ec88355ff88
# ╠═2cf223a5-073a-4e98-936a-8d33607cb773
# ╠═9bfce8b2-0814-4ba6-9a77-633f4fc2624e
# ╠═f8de92ce-b40b-4aa2-af73-81f8301ab65a
# ╠═65cd3f65-d34c-4971-929c-c4f1b4b7d1ec
# ╠═e79d89d6-87da-49a9-87d4-0924aa19f923
# ╟─cc8f4815-db8e-4f3c-abd1-dd5c6d6e0587
# ╟─f5081fff-51d1-44c0-8c51-bccb7c08f996
# ╟─0558677f-f998-4699-8b6d-00ef7a011b2f
# ╟─f1bf5692-0242-4086-9147-520055089c6d
# ╠═e6b706f0-2107-4952-9b45-068ae437dfeb
# ╟─92ce4f62-384b-4c19-9ce2-e451abfd81ac
# ╟─834720f5-0cef-437d-b09f-939d1d763167
# ╠═c68896a8-6aab-4600-8a5c-52fbf9cd3651
# ╠═bfacb6e6-b445-47f7-bef1-166b67bc8eee
