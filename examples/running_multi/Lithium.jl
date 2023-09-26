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

# ╔═╡ 0402e3b4-cc36-4e30-952a-ea897ee6ce7c
matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))

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
	modelatmos = "muram_m2_HDm_50x50"
	modelatmosfolder = "./input_multi3d/muram_m2"
end

# ╔═╡ 0d5f972c-f308-46a9-9b35-2abe5343329e
snapshots = [
	"muram_m2_HDm_50x50",
	"muram_m2_SSDm_50x50",
	"muram_m2_HDl_50x50",
	"muram_m2_SSDl_50x50",
	"muram_m2_HDh_50x50",
	"muram_m2_SSDh_50x50"
]

# ╔═╡ f6de91e4-6de8-4502-b17d-17325464ffdf
names = [
	"HDm", 
	"SSDm"
]

# ╔═╡ bfc1d715-31ac-4a1f-b1c5-9d581f8c7c10
#linelist = "./input_multi3d/nlte_ges_linelist_jmg25jan2023_I_II"
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

	fmode, amodel = plt.subplots(1, 1, figsize=(5, 6))

	for (i, snap) in enumerate(names)
		snap = snapshots[i]
		b = MUST.multiBox(@in_m3dis(joinpath(modelatmosfolder, snap)))
		
		_, d = profile(MUST.mean, b, :z, :log10d)
		_, T = profile(MUST.mean, b, :z, :T)
			
		amodel.plot(d, T, marker="", label=names[i])
	end

	amodel.set_xlabel(L"\rm \log \rho\ [g \times cm^{-3}]")
	amodel.set_ylabel(L"\rm T\ [K]")
	amodel.set_ylim(2500, 11000)
	amodel.set_xlim(-9, -6)
	amodel.legend(handlelength=4, labelspacing=0.1, ncol=1, loc="upper left")
	
	gcf()
end

# ╔═╡ 33663b6e-6f26-45bb-81fd-55d2ed0fdf42
begin
	plt.close()

	plt.imshow(uz_surface_optical(model_hdm), extent=extent(model_hdm))

	plt.xlabel("x")
	plt.ylabel("y")
	
	gcf()
end

# ╔═╡ d462e034-d131-41f7-8d1e-fb9d8364dc7a
md"## NLTE and Spectrum Synthesis"

# ╔═╡ 105e6139-7d0d-4511-b7cf-5e571f941399
compute = false

# ╔═╡ bb5109e4-ff6e-4941-8c0c-30735f7b0c30
input_parameters = Dict(
	:model_folder=>modelatmosfolder,
	:linelist=>linelist,
	:absmet=>nothing,
	:atom_params=>(
		:atom_file=>"./input_multi3d/atoms/atom.li12_col", 
		:use_atom_abnd=>false,
		:abundance=>2.1,
		:exclude_trace_cont=>true,
		:exclude_from_line_list=>true,
		:hydrogen_BPO=>false
	),
	:spectrum_params=>(
		:daa=>0.01, :aa_blue=>6706, :aa_red=>6710
	),
	:composition_params=>(
		:absdat_file=>"./input_multi3d/TS_absdat.dat",
		:abund_file=>"./input_multi3d/abund_asplund07"
	),
	#:composition_params=>(
	#	:absdat_file=>"./input_multi3d/absdat", 
	#	:abund_file=>"./input_multi3d/abund"
	#),
	:atmos_params=>(
		:atmos_format=>"MUST", 
		:use_density=>true, 
		:use_ne=>false,
		:FeH=>-2.0,
		:dims=>16
	),
	:m3d_params=>(
		:n_nu=>1, 
		:quad_scheme=>"set_a2"
	)
)

# ╔═╡ 099cb6f0-5b6f-4c51-8b28-8cf172e064d6
md"We can run multiple instances of this setup at the same time, with e.g. different abundances"

# ╔═╡ 50a9528b-a021-4407-97de-65b5b9efcb45
begin
	ab = range(0.4, 3.5, length=10) |> collect
	append!(ab, [2.1])
	extension = "TSabsdat_noblend" 
end

# ╔═╡ 6262c95a-7be3-438e-b19f-49f36bfeadd3
begin
	params = Dict()
	output_paths = []
	abundances = []

	for abundance in ab
		abundance = round(abundance, digits=3)
		input_mod = deepcopy(input_parameters)
		
		ap_new = if haskey(input_parameters, :atom_params)
			ap = input_mod[:atom_params]
			ap_new = Tuple(a==:abundance ? a=>abundance : a=>v for (a, v) in ap)
		else
			(:abundance=>abundance, )
		end

		input_mod[:atom_params] = ap_new
		params[join(["$(extension)","$(abundance)"], "_")] = Tuple(input_mod)
		append!(abundances, [abundance])
		
		for s in snapshots
			append!(
				output_paths, [
					joinpath(
						"data", 
						join([s, "$(extension)", "$(abundance)"], "_")
					)
				]
			)
		end
	end

	abundances = sort(abundances)
end

# ╔═╡ b2731305-845b-4592-9751-70913a659787
md"Run those parameters"

# ╔═╡ 99ba9ccd-98f5-4510-95e6-2e7ea0955e8f
begin
	 if compute
		MUST.spectrum(
			snapshots, 			    # run these snapshots
			params, 				# and those those paramters
			NLTE=true, 			    # in NLTE
			slurm=true, 			# using Slurm (+ wait)
			twostep=false 			# first dep. => LTE + dep.
		)
	end
end

# ╔═╡ 8ef84170-65cc-4ac7-8e9b-e81a8f40c506
md"## Read Output"

# ╔═╡ cb8c9eac-0e1c-4cba-8e6b-754d38c7dacd
extension_results(extension, folder="data", exclude="LTE") = begin
	all_results = MUST.glob("*/", MUST.@in_m3dis(folder))
	paths = []
	
	for (i, r) in enumerate(all_results)
		if occursin(extension, r) && !occursin(exclude, r) 
			append!(paths, [r])
		end
	end

	joinpath.(folder, last.(split.(dirname.(paths), "/")))
end

# ╔═╡ 341291d8-0d37-49a5-8ecc-a62c3c24adca
getabundances(run) = MUST.pyconvert(Float64, run.atom.abnd)

# ╔═╡ 36141454-25db-4415-9d91-7ed6f626ecd8
resultpaths = extension_results(extension)

# ╔═╡ 9e93ecdb-dada-43b0-8107-bcde4702f1f5
begin
	m3druns = MUST.M3DISRun.(resultpaths) 
	abund = getabundances.(m3druns)
	abund = sort(unique(abund))
end

# ╔═╡ 923f8164-9ace-49d0-8ab3-e274d500011d
md"## Investigate Diagnostic Li Lines"

# ╔═╡ 897a718a-759f-4195-905c-afb86df6bef7
find_model(runs, kind, abundance) = begin
	models = [split(dirname(join(r.sfolder)), "/")[end-1] for r in runs]
	found = []

	for (i, model) in enumerate(models)
		
		if occursin(kind, model) 
			#=if occursin(extension, model)
				if occursin("LTE", extension)
					if !occursin("LTE", model)
						continue
					end
				else
					if occursin("LTE", model)
						continue
					end
			end=#
			abund = MUST.pyconvert(Any, runs[i].atom.abnd)
			if abund ≈ abundance
				append!(found, [runs[i]])
			end
			#end
		end
	end

	if length(found) == 0
		@show kind, abundance
	end
	
	if length(found) == 1
		first(found)
	else
		@show found
		found
	end
end

# ╔═╡ ff26a5e4-3052-4323-8845-7952f8dc568a
begin
	plt.close()

	#colors_models = ["k", "royalblue"]
	names_models = names #["HDm", "SSDm"]

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))
	
	cmap = matplotlib.cm.get_cmap("gnuplot")
	colors_models = [cmap((i-1)/length(abund)) for i in eachindex(abund)]

	# HD LTE and NLTE
	for (i, n) in enumerate(names_models)
		i==2 && break
		for (j, abundance) in enumerate(abund)
			j == 1 && continue
			
			run = find_model(m3druns, n, abundance)
			
			line1 = run.line[1]
			axA.plot(
				window(line1, LTE=false)..., 
				color=colors_models[j], 
				ls="-", 
				label="$(names_models[i]), $(abundance)"*L"\rm- NLTE", 
				lw=1
			)
		end

		run = find_model(m3druns, n, 2.1)	
		line1 = run.line[1]
		axA.plot(
			window(line1, LTE=true)..., 
			color=colors_models[i], 
			ls="--", 
			label="$(names_models[i]), $(2.1)"*L"\rm- LTE", 
			lw=1
		)
		line1 = run.line[1]
		axA.plot(
			window(line1, LTE=false)..., 
			color=colors_models[i], 
			ls="-", 
			label="$(names_models[i]), $(2.1)"*L"\rm- NLTE", 
			lw=2
		)
	end
	
		
	axA.set_xlim(6707.3, 6708.4)
	axA.set_ylabel(L"\rm flux\ [normalized]")
	axA.set_xlabel(L"\rm wavelength\ [\AA]")
	axA.legend(handlelength=4, labelspacing=0.1, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0))
	
	gcf()
end

# ╔═╡ fd8ccf89-c256-4fdf-8d7a-8388e4328bf4
EW(abund, ew) = MUST.linear_interpolation(abund, ew)

# ╔═╡ 897f9b1c-8109-47fa-866b-661f75143440
AB(abund, ew) = MUST.linear_interpolation(ew[sortperm(ew)], abund[sortperm(ew)])

# ╔═╡ 16498b25-45b0-488f-bf91-bfda9c3a59b9
begin
	plt.close()

	fB, axB = plt.subplots(1, 1, figsize=(5, 6))
	colors_modelsB = ["k", "royalblue", "tomato", "blue", "magenta", "cyan"]

	fNLTE_models = []
	fLTE_models = []
	aNLTE_models = []
	aLTE_models = []
	clabels_models = []

	box_loc = (0.03, 0.8)
	
	# HD LTE and NLTE
	for (i, n) in enumerate(names_models)
		ew_LTE = []
		ew_NLTE = []
		for (j, abundance) in enumerate(abund)
			#j == 1 && continue
			
			run = find_model(m3druns, n, abundance)
			line1 = run.line[1]

			append!(ew_LTE, [pc(line1.calc_weq(LTE=true))])
			append!(ew_NLTE, [pc(line1.calc_weq(LTE=false))])
		end

		fLTE, aLTE = EW(abund, ew_LTE), AB(abund, ew_LTE)
		fNLTE, aNLTE = EW(abund, ew_NLTE), AB(abund, ew_NLTE)

		append!(fNLTE_models, [fNLTE])
		append!(fLTE_models,  [fLTE])
		append!(aNLTE_models, [aNLTE])
		append!(aLTE_models,  [aLTE])
		
		axB.plot(
			abund, ew_LTE, 
			label="$(names_models[i])"*L"\rm- LTE",
			color=colors_modelsB[i], ls="--", lw=1.5
		)
		axB.plot(
			abund, ew_NLTE, 
			label="$(names_models[i])"*L"\rm- NLTE",
			color=colors_modelsB[i], ls="-", lw=2.5
		)
		
		ΔNLTE1 = aNLTE(2.1 |> fLTE)-2.1
		ΔNLTE2 = 2.1 - aLTE(2.1 |> fNLTE)
		#=axB.text(
			aNLTE(2.1 |> fLTE) .+ 0.03, fNLTE(2.1 |> fLTE |> aNLTE),
			L"\rm \Delta_{NLTE} = "*"$(round(ΔNLTE1, sigdigits=2))",
			ha="left", va="center", color=colors_modelsB[i]
		)
		axB.text(
			aLTE(2.1 |> fNLTE) .- 0.03, fLTE(2.1 |> fNLTE |> aLTE),
			L"\rm \Delta_{NLTE} = "*"$(round(ΔNLTE2, sigdigits=2))",
			ha="right", va="center", color=colors_modelsB[i]
		)=#
		
		axB.axvline(2.1, color="k", alpha=0.4, ls=":")
		#=axB.axvline(
			aNLTE(2.1 |> fLTE), 
			label="$(names_models[i])"*L"\rm < \Delta_{NLTE} > = "*
					"$(round(MUST.mean([ΔNLTE1, ΔNLTE2]), sigdigits=2))",
			color=colors_modelsB[i], alpha=0.4, ls=":"
		)=#
		append!(
			clabels_models, 
			["$(names_models[i]) "*L"\rm \Delta_{NLTE} = "*
					"$(round(MUST.mean([ΔNLTE1, ΔNLTE2]), sigdigits=2))"]
		)
		#axB.axvline(aLTE(2.1 |> fNLTE), color=colors_modelsB[i], alpha=0.4, ls=":")
		#axB.axhline(fNLTE(2.1 |> fLTE |> aNLTE), color="k", alpha=0.4, ls=":")
		#axB.axhline(fLTE(2.1 |> fNLTE |> aLTE), color="k", alpha=0.4, ls=":")
	end

	# Difference between the models (all compared to first)
	for i in 2:length(names_models)
		Δmodel_LTE1 = aLTE_models[i](2.1 |> fLTE_models[1])-2.1
		Δmodel_LTE2 = 2.1 - aLTE_models[1](2.1 |> fLTE_models[i])
		Δmodel_NLTE1 = aNLTE_models[i](2.1 |> fNLTE_models[1])-2.1
		Δmodel_NLTE2 = 2.1 - aNLTE_models[1](2.1 |> fNLTE_models[i])

		Δmodel_LTE = round(MUST.mean([Δmodel_LTE1, Δmodel_LTE2]), sigdigits=2)
		Δmodel_NLTE = round(MUST.mean([Δmodel_NLTE1, Δmodel_NLTE2]), sigdigits=2)

		append!(
			clabels_models, 
			[L"\rm \Delta ("*"$(names_models[i])"*L"\rm )_{LTE} = "*"$(Δmodel_LTE)"]
		)
		append!(
			clabels_models, 
			[L"\rm \Delta ("*"$(names_models[i])"*L"\rm )_{NLTE} = "*"$(Δmodel_NLTE)"]
		)
	end

	for (i, c) in enumerate(clabels_models)
		v, h = box_loc[1], box_loc[2]
		h -= 0.03*i
		axB.text(
			v, h,
			clabels_models[i],
			ha="left", va="top", color="k",
			transform=axB.transAxes
		)
	end
	
	
	axB.set_xlim(0.85, 3.2)
	axB.set_ylim(-10, 300)
	axB.set_title(L"\rm Curve\ of\ Growth\ -\ 6707 \AA")
	axB.set_ylabel(L"\rm EW\ [\AA]")
	axB.set_xlabel(L"\rm A(Li)")
	axB.legend(handlelength=4, labelspacing=0.1, ncol=1, loc="upper left")
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─f2f5e1f6-4de6-4fc3-b9b7-8eb7a168a5aa
# ╠═317e15be-52df-11ee-03c7-05244dd05fc9
# ╠═0402e3b4-cc36-4e30-952a-ea897ee6ce7c
# ╟─7a212f19-7007-4e2b-a2bf-bddd843a3ff8
# ╟─72d7e596-02bf-4940-ac6b-7fea8e51485d
# ╟─183529b8-dfd5-48aa-8b34-4a97b4744b4b
# ╟─0233951f-e65e-4fb8-8597-6ba43632ee81
# ╠═2d04602c-1a6a-4915-97c3-27b13418456a
# ╠═0d5f972c-f308-46a9-9b35-2abe5343329e
# ╠═f6de91e4-6de8-4502-b17d-17325464ffdf
# ╠═bfc1d715-31ac-4a1f-b1c5-9d581f8c7c10
# ╠═e0b6b3a9-dfb8-4fc4-8f87-3e82be556e27
# ╟─d770f462-f608-4f4f-810b-e848aa580a93
# ╟─e284e0a5-3420-4cc7-b049-f66c60952d28
# ╟─f176583e-3008-4846-a9d6-adbafaa0c7eb
# ╟─7e1f858b-62a9-48b8-9095-104a1d3dc4e9
# ╟─33663b6e-6f26-45bb-81fd-55d2ed0fdf42
# ╟─d462e034-d131-41f7-8d1e-fb9d8364dc7a
# ╠═105e6139-7d0d-4511-b7cf-5e571f941399
# ╟─bb5109e4-ff6e-4941-8c0c-30735f7b0c30
# ╟─099cb6f0-5b6f-4c51-8b28-8cf172e064d6
# ╠═50a9528b-a021-4407-97de-65b5b9efcb45
# ╟─6262c95a-7be3-438e-b19f-49f36bfeadd3
# ╟─b2731305-845b-4592-9751-70913a659787
# ╠═99ba9ccd-98f5-4510-95e6-2e7ea0955e8f
# ╟─8ef84170-65cc-4ac7-8e9b-e81a8f40c506
# ╟─cb8c9eac-0e1c-4cba-8e6b-754d38c7dacd
# ╟─341291d8-0d37-49a5-8ecc-a62c3c24adca
# ╠═36141454-25db-4415-9d91-7ed6f626ecd8
# ╟─9e93ecdb-dada-43b0-8107-bcde4702f1f5
# ╟─923f8164-9ace-49d0-8ab3-e274d500011d
# ╟─897a718a-759f-4195-905c-afb86df6bef7
# ╟─ff26a5e4-3052-4323-8845-7952f8dc568a
# ╟─fd8ccf89-c256-4fdf-8d7a-8388e4328bf4
# ╟─897f9b1c-8109-47fa-866b-661f75143440
# ╟─16498b25-45b0-488f-bf91-bfda9c3a59b9
