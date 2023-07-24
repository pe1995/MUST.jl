### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 6f2ac3f4-d1bb-4564-b35e-20c2a98b386b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using MUST
	using TSO
	using PythonPlot	
end

# ╔═╡ f8bc39c3-c284-46ce-987f-c8fa810dc207
md"# MUST & M3D"

# ╔═╡ 44248cdb-a3e0-4a3a-862f-9384dcdc1a3d
PythonPlot.matplotlib.rcParams["font.size"] = 12

# ╔═╡ 1c1b28b5-8319-4cd1-a72a-750c7aed124a
plt = PythonPlot.matplotlib.pyplot

# ╔═╡ 24364b76-d62e-4da1-9c1f-e8b9fa53d60d
visual = MUST.ingredients("visual.jl")

# ╔═╡ 93d55469-8b09-48c1-9138-5783933a96ff
pc(x) = MUST.pyconvert(Any, x)

# ╔═╡ 78fdd0d7-2853-4f52-8a3b-39a62818779d
md"## Dispatch & M3D
We load the Python modules for handling Dispatch, M3D and loading the correct module locations."

# ╔═╡ 031de813-3a84-489b-9279-dcfb3f27a6bb
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 9ce0786f-fef0-4d1e-9a38-07ea96bd147f
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 2a19b30b-1f6c-4b50-8137-9a748d2dd028
md"## Model Atmospheres
Pick the model atmosphere. It should either be located in input_multi3d/MUST or the correct folder should be given to the namelist function."

# ╔═╡ 04f7c167-3b1d-4ad8-be9f-58fbac9a20a0
modelatmos = "m3dis_sun_magg22_10x10x230_3"

# ╔═╡ b7790d90-97bf-4f89-8c15-0d70f855471d
modelatmosfolder = "input_multi3d/magg2022_150x300/"

# ╔═╡ 6ccb7b73-2ea1-40d8-9ade-55948682cde9
md"## Line List
One can either use no line list to compute Teff, use a line list in the Turbospectrum format, or use an absmet file. The absmet file contains sampled backgournd opacities, and might be good enough for the spectrum computation in terms of Teff. To activate the absmet, one has to set the linelist to nothing, in order for the namelist to be created correctly."

# ╔═╡ 03ed4c44-da3b-45bc-946b-634ca23cf30a
linelist = "./input_multi3d/vald_2490-25540.list"

# ╔═╡ 484b7245-8e34-41b2-b9f8-2a9b3c1ba9e4
absmet = "./input_multi3d/absmet"

# ╔═╡ fe9351c2-c97b-4a66-a230-f2a0ab8067ce
md"## Opacity tables
We compute the heating from the unbinned and binned opacity table and then compare the results. Both is possible within M3D."

# ╔═╡ 5185341e-9299-40d5-9881-9e002c926bfc
eos = reload(
	SqEoS, 
	"/u/peitner/DISPATCH/binned/DIS_MARCS_v1.6.3/eos.hdf5"
)

# ╔═╡ fa7355e5-9072-4c4e-be9a-d5412f6cde5c
opa = reload(
	SqOpacity, 
	"/u/peitner/DISPATCH/binned/DIS_MARCS_v1.6.3/binned_opacities.hdf5"
)

# ╔═╡ 69ee39a1-ae7d-4ec9-bae2-b830f80836ed
md"These tables need to be combined and saved in the Multi format in order to be able to run them in M3D. We save them in common folder directly in Multi to make things easier."

# ╔═╡ fb609e98-ab06-487e-8bf1-91e232e5687b
folder_eos = MUST.@in_m3dis "input_multi3d/M3D_MARCS_v1.6.3"

# ╔═╡ 7d38779c-c2e3-4ff2-b30c-f4229f35545f
!isdir(folder_eos) && mkdir(folder_eos)

# ╔═╡ 3c76bd08-a807-4398-b277-a622e010d05e
m3d_table = save(eos, opa, joinpath(folder_eos, "eos_opa"))

# ╔═╡ 81ef72b6-eafa-49ef-9923-fb290cc4091d
md"## Running M3D
First we need to run M3D in the unbinned case, compute the heating, and the run it in the binned case for comparison. If you have computed this already, you can also read an existing output."

# ╔═╡ f9e7cd2c-7a09-4eb5-b0a1-1247ff4fbdf3
compute = true

# ╔═╡ cbe117c0-f6c7-44bb-a86e-3d57362666d7
m3d_binned = if false 
	MUST.heating(
		modelatmos, 
		m3d_table,
		namelist_kwargs=(
			:model_folder=>modelatmosfolder,
			:linelist=>nothing,
			:absmet=>nothing
		),
		slurm=false
	)
else
	MUST.M3DISRun("data/m3dis_sun_magg22_10x10x230_3_binned")
end

# ╔═╡ 4668ca71-bd41-413a-8c67-d7e7dfe9b6ba
m3d_unbinned = if compute
	MUST.heating(
		modelatmos, 
		namelist_kwargs=(
			:model_folder=>modelatmosfolder,
			:linelist=>nothing,
			:absmet=>nothing,
			:atom_params=>(:atom_file=>"", ),
			:spectrum_params=>(:daa=>1., :aa_blue=>1500, :aa_red=>9000)
		),
		slurm=true
	)
else
	MUST.M3DISRun("data/m3dis_sun_magg22_10x10x230_3_unbinned")
end

# ╔═╡ fb0a997b-f2dd-4a19-bdd6-c5293159eb2c
md"## Heating Comparison"

# ╔═╡ 4fe50a65-7729-4a13-ab39-d81c62fae60d
qrad = m3d_unbinned.read_patch_save("Qrad")

# ╔═╡ 134144b6-6f21-42d5-a816-fd5bb57b9e55
begin
	fH, axH = plt.subplots(1, 1, figsize=(7, 6))
end

# ╔═╡ Cell order:
# ╟─f8bc39c3-c284-46ce-987f-c8fa810dc207
# ╠═6f2ac3f4-d1bb-4564-b35e-20c2a98b386b
# ╠═44248cdb-a3e0-4a3a-862f-9384dcdc1a3d
# ╠═1c1b28b5-8319-4cd1-a72a-750c7aed124a
# ╠═24364b76-d62e-4da1-9c1f-e8b9fa53d60d
# ╠═93d55469-8b09-48c1-9138-5783933a96ff
# ╟─78fdd0d7-2853-4f52-8a3b-39a62818779d
# ╠═031de813-3a84-489b-9279-dcfb3f27a6bb
# ╠═9ce0786f-fef0-4d1e-9a38-07ea96bd147f
# ╟─2a19b30b-1f6c-4b50-8137-9a748d2dd028
# ╠═04f7c167-3b1d-4ad8-be9f-58fbac9a20a0
# ╠═b7790d90-97bf-4f89-8c15-0d70f855471d
# ╟─6ccb7b73-2ea1-40d8-9ade-55948682cde9
# ╠═03ed4c44-da3b-45bc-946b-634ca23cf30a
# ╠═484b7245-8e34-41b2-b9f8-2a9b3c1ba9e4
# ╟─fe9351c2-c97b-4a66-a230-f2a0ab8067ce
# ╠═5185341e-9299-40d5-9881-9e002c926bfc
# ╠═fa7355e5-9072-4c4e-be9a-d5412f6cde5c
# ╟─69ee39a1-ae7d-4ec9-bae2-b830f80836ed
# ╠═fb609e98-ab06-487e-8bf1-91e232e5687b
# ╠═7d38779c-c2e3-4ff2-b30c-f4229f35545f
# ╠═3c76bd08-a807-4398-b277-a622e010d05e
# ╟─81ef72b6-eafa-49ef-9923-fb290cc4091d
# ╠═f9e7cd2c-7a09-4eb5-b0a1-1247ff4fbdf3
# ╠═cbe117c0-f6c7-44bb-a86e-3d57362666d7
# ╠═4668ca71-bd41-413a-8c67-d7e7dfe9b6ba
# ╟─fb0a997b-f2dd-4a19-bdd6-c5293159eb2c
# ╠═4fe50a65-7729-4a13-ab39-d81c62fae60d
# ╠═134144b6-6f21-42d5-a816-fd5bb57b9e55
