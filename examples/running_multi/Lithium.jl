### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 317e15be-52df-11ee-03c7-05244dd05fc9
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using TSO
	using PythonPlot	
	using LaTeXStrings

	plt = matplotlib.pyplot
end;

# ╔═╡ f2f5e1f6-4de6-4fc3-b9b7-8eb7a168a5aa
md"# The Lithium Problem
## Setup"

# ╔═╡ 7a212f19-7007-4e2b-a2bf-bddd843a3ff8
begin
	pc(x) = MUST.pyconvert(Any, x)
	pc(T, x) = MUST.pyconvert(T, x)
end

# ╔═╡ 72d7e596-02bf-4940-ac6b-7fea8e51485d
begin
	MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
end

# ╔═╡ 183529b8-dfd5-48aa-8b34-4a97b4744b4b
extent(snap) = begin
	x = MUST.axis(snap, :x) ./1e8
	y = MUST.axis(snap, :x) ./1e8

	[minimum(x), maximum(x), minimum(y), maximum(y)]
end

# ╔═╡ 0233951f-e65e-4fb8-8597-6ba43632ee81
md"## Input Data"

# ╔═╡ 2d04602c-1a6a-4915-97c3-27b13418456a
begin
	modelatmos = "muram_m2_HDm_20x20"
	#modelatmos = "m3dis_sun_magg22_10x10x299_3"
	modelatmosfolder = "./input_multi3d/muram_m2"
	#modelatmosfolder = "./input_multi3d/magg2022_150x300"
end

# ╔═╡ 0d5f972c-f308-46a9-9b35-2abe5343329e
snapshots = [
	modelatmos
]

# ╔═╡ bfc1d715-31ac-4a1f-b1c5-9d581f8c7c10
linelist = "./input_multi3d/vald_2490-25540.list"

# ╔═╡ e0b6b3a9-dfb8-4fc4-8f87-3e82be556e27
absmet = "./input_multi3d/absmet"

# ╔═╡ d770f462-f608-4f4f-810b-e848aa580a93
md"## Model Atmospheres"

# ╔═╡ e284e0a5-3420-4cc7-b049-f66c60952d28
begin
	model_hdm = MUST.multiBox(@in_m3dis(joinpath(modelatmosfolder, modelatmos)))
end

# ╔═╡ f176583e-3008-4846-a9d6-adbafaa0c7eb
uz_surface_optical(snap) = begin
	isurf = MUST.closest(MUST.axis(model_hdm, :z), 0)
	snap[:uz][:, :, isurf] ./1e5
end

# ╔═╡ 7e1f858b-62a9-48b8-9095-104a1d3dc4e9
begin
	plt.close()

	_, d = profile(MUST.mean, model_hdm, :z, :log10d)
	_, T = profile(MUST.mean, model_hdm, :z, :T)
		
	plt.plot(d, T, marker=".")
	
	gcf()
end

# ╔═╡ 33663b6e-6f26-45bb-81fd-55d2ed0fdf42
begin
	plt.close()

	plt.imshow(uz_surface_optical(model_hdm), extent=extent(model_hdm))
	
	gcf()
end

# ╔═╡ d462e034-d131-41f7-8d1e-fb9d8364dc7a
md"## NLTE and Spectrum Synthesis"

# ╔═╡ 105e6139-7d0d-4511-b7cf-5e571f941399
compute = true

# ╔═╡ bb5109e4-ff6e-4941-8c0c-30735f7b0c30
input_parameters = (
	:model_folder=>modelatmosfolder,
	:linelist=>nothing,
	:absmet=>nothing,
	:atom_params=>(
		:atom_file=>"./input_multi3d/atoms/atom.ba06", 
		:use_atom_abnd=>false,
		:abundance=>2.1,
		:exclude_trace_cont=>false,
		:exclude_from_line_list=>false,
		:hydrogen_BPO=>false
	),
	:spectrum_params=>(
		:daa=>0.1, :aa_blue=>6700, :aa_red=>6712
	),
	:composition_params=>(
		:absdat_file=>"./input_multi3d/TS_absdat.dat", :abund_file=>"./input_multi3d/abund_asplund07"
	),
	:atmos_params=>(
		:atmos_format=>"MUST", 
		:use_density=>true, 
		:use_ne=>true,
		:FeH=>-2.0
	)
)

# ╔═╡ 99ba9ccd-98f5-4510-95e6-2e7ea0955e8f
begin
	 if compute
		m3load = MUST.spectrum(
			first(snapshots), 					# run this snapshot
			NLTE=true, 							# in NLTE
			slurm=true, 						# using Slurm (+wait)
			twostep=false, 						# first dep. => LTE + dep.
			namelist_kwargs=input_parameters    # use the input parameters from above
		)
	end
end

# ╔═╡ 8ef84170-65cc-4ac7-8e9b-e81a8f40c506
md"## Read Output"

# ╔═╡ 7ff1293e-8d31-469e-b4ab-ff5f9c3046bc
m3druns = MUST.M3DISRun.(joinpath.("data", snapshots))

# ╔═╡ 923f8164-9ace-49d0-8ab3-e274d500011d
md"## Investigate Diagnostic Li Lines"

# ╔═╡ b4d715b1-a36a-44cc-86f7-1ee10a28d3f6
hdm = m3druns[1]

# ╔═╡ 4505e740-7880-4e64-9e29-92101eaaf2c0
line_hdm = hdm.line[1]

# ╔═╡ ff26a5e4-3052-4323-8845-7952f8dc568a
begin
	plt.close()

	plt.plot(window(line_hdm, LTE=true)...)
	
	gcf()
end

# ╔═╡ 4d017237-df5f-45fd-b971-6fc990a88b03
window(line_hdm, LTE=true)

# ╔═╡ Cell order:
# ╟─f2f5e1f6-4de6-4fc3-b9b7-8eb7a168a5aa
# ╠═317e15be-52df-11ee-03c7-05244dd05fc9
# ╟─7a212f19-7007-4e2b-a2bf-bddd843a3ff8
# ╟─72d7e596-02bf-4940-ac6b-7fea8e51485d
# ╠═183529b8-dfd5-48aa-8b34-4a97b4744b4b
# ╟─0233951f-e65e-4fb8-8597-6ba43632ee81
# ╠═2d04602c-1a6a-4915-97c3-27b13418456a
# ╠═0d5f972c-f308-46a9-9b35-2abe5343329e
# ╠═bfc1d715-31ac-4a1f-b1c5-9d581f8c7c10
# ╠═e0b6b3a9-dfb8-4fc4-8f87-3e82be556e27
# ╟─d770f462-f608-4f4f-810b-e848aa580a93
# ╠═e284e0a5-3420-4cc7-b049-f66c60952d28
# ╟─f176583e-3008-4846-a9d6-adbafaa0c7eb
# ╠═7e1f858b-62a9-48b8-9095-104a1d3dc4e9
# ╟─33663b6e-6f26-45bb-81fd-55d2ed0fdf42
# ╟─d462e034-d131-41f7-8d1e-fb9d8364dc7a
# ╠═105e6139-7d0d-4511-b7cf-5e571f941399
# ╠═bb5109e4-ff6e-4941-8c0c-30735f7b0c30
# ╠═99ba9ccd-98f5-4510-95e6-2e7ea0955e8f
# ╟─8ef84170-65cc-4ac7-8e9b-e81a8f40c506
# ╠═7ff1293e-8d31-469e-b4ab-ff5f9c3046bc
# ╟─923f8164-9ace-49d0-8ab3-e274d500011d
# ╠═b4d715b1-a36a-44cc-86f7-1ee10a28d3f6
# ╠═4505e740-7880-4e64-9e29-92101eaaf2c0
# ╠═ff26a5e4-3052-4323-8845-7952f8dc568a
# ╠═4d017237-df5f-45fd-b971-6fc990a88b03
