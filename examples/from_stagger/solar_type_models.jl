### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ dcd11d0e-312d-11ee-1d21-4dc68639a432
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using Plots
	using TSO 
	using LaTeXStrings
	using DelimitedFiles
end;

# ╔═╡ 309e7de7-5a70-41f6-8ca6-47ad4b46c214
md"# Solar-type Stars"

# ╔═╡ 41269b8e-5a49-4930-84de-aa4248531174
begin
	mean = MUST.mean
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
end

# ╔═╡ 7494384f-0a14-4ad0-95ad-e8c95fab78b7
md"Load conversion functions"

# ╔═╡ 77153ffb-00e4-449b-99b2-b377a2ce2f7a
must2multi = MUST.ingredients("convert2multi.jl")

# ╔═╡ d9a0ea6c-9197-4633-97d5-57ada5921ce4
md"## Models & Opacities"

# ╔═╡ a66ad3b4-faed-41c1-af37-9cd96bb9aa41
names = [
	"DIS_MARCS_E_t45g40m00_v0.1"
]

# ╔═╡ 6202080d-fa73-4f39-8232-914582f70d11
out_folder = [
	MUST.@in_dispatch "data/DIS_MARCS_E_t45g40m00_v0.1_2"
]

# ╔═╡ a4a4f729-f247-440d-80d2-455b30c0f63f
eos_folder = [MUST.@in_dispatch("input_data/grd/$n") for n in names]

# ╔═╡ 36b50d99-f828-42a6-a617-3ff68feb36b9
labels = [
	"magg2022_150x300_t45g40m00"
]

# ╔═╡ 46402de2-102a-41ed-ac85-2c65be8d87d4
colors = length(names) == 1 ? [:red] : palette(:roma, length(names))

# ╔═╡ ad1c4dcc-817d-4781-b23f-ed63a640f444
in_folder  = [MUST.@in_dispatch "input_data/grd/$(name)" for name in names]

# ╔═╡ b1e8f459-ed60-40ab-a8c3-c8d4c1d8fa5a
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = MUST.pick_snapshot(
			MUST.converted_snapshots(out_folder[i]), -1
		)
		
		append!(snapshots, [snapshot])
		append!(snapshots_τ, [snapshot_τ])
	end
end

# ╔═╡ 2431f1ac-80cb-47d9-a657-7c4569ea6dfd
begin
	eos = [reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5")) 
				for i in eachindex(eos_folder)]
	opa = [reload(SqOpacity, joinpath(eos_folder[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder)]
end

# ╔═╡ ab08661d-fc62-4700-b423-9b96df85d9cb
md"## Average profiles"

# ╔═╡ 42a0408c-20cd-49b0-8581-7a9d13e2361b
ls = [:solid, :dash, :solid]

# ╔═╡ daec2515-64db-45ef-8a29-ab4fc8298963
lw = [3, 3, 1]

# ╔═╡ cbf793b6-c2ba-4872-9f4c-522c3a77d5a5
begin
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)	
	for i in eachindex(names)
		pav1 = plot!(
			time_average_profile(
				mean, 
				out_folder[i], 
				:log10τ_ross, 
				:T,
				hscale=:τ
			)[1:2]..., 
			label=labels[i],
			color=:black,
			ls=ls[i], 
			lw=lw[i]
		)
	end

	#plot!(xlim=(-4,4), ylim=(3500,13300))
	plot!(ylabel="T [K]", xlabel="log τ-ross")
end

# ╔═╡ 2751d10c-3e48-4578-b330-97692eb25133
begin
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	for i in eachindex(names)
		pavr1 = plot!(
			time_average_profile(
				mean, 
				out_folder[i], 
				:log10τ_ross, 
				:log10d,
				hscale=:τ
			)[1:2]..., 
			label=labels[i],
			color=:black,
			ls=ls[i], 
			lw=lw[i]
		)
	end

	#plot!(xlim=(-4,4), ylim=(-8.7,-5.7))
	plot!(ylabel="log ρ [g × cm-3]", xlabel="log τ-ross")
end

# ╔═╡ 724239e6-aac8-4096-907e-e87deaded5c6
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ 3239e571-f59e-4648-87ff-d6003d4128b9
rms5(x) = rms(x ./ 1e5)

# ╔═╡ bc94c2a0-dbb3-4906-80d6-3bfdccd5dbb5
begin
	plot(
		framestyle=:box, 
		grid=false, 
		legendforegroundcolor=nothing,
		legendbackgroundcolor=nothing,
		legendposition=:bottomright
	)
	for i in eachindex(names)
		pavr1 = plot!(
			time_average_profile(
				rms5, 
				out_folder[i], 
				:log10τ_ross, 
				:uz,
				hscale=:τ
			)[1:2]..., 
			label=labels[i],
			color=:black,
			ls=ls[i], 
			lw=lw[i]
		)
	end

	#plot!(xlim=(-4,4), ylim=(-8.7,-5.7))
	plot!(ylabel="rms Uz [km × s-1]", xlabel="log τ-ross")
end

# ╔═╡ a8dc248c-32be-452e-94bc-ca43f6afbda7
md"## Snapshot deviations"

# ╔═╡ 42fcf382-b831-4378-bf97-f1ee6601d388
begin
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	csnaps = palette(:rainbow, length(range(-8, -1)) + 1)
	model = 1
	
	for (i, snapi) in enumerate([range(-8, -1)..., :recent])
		s = pick_snapshot(out_folder[model], snapi) |> last
		pav1 = plot!(
			profile(mean, s, :log10τ_ross, :T)..., 
			label=nothing,
			color=csnaps[i],
			ls=:solid, 
			lw=1
		)
	end

	plot!(title=labels[model])
	plot!(xlim=(-4, 5), ylim=(2500, 13500))
	plot!(ylabel="T [K]", xlabel="log τ-ross")

end

# ╔═╡ 7a945e98-c82d-48d8-978c-a5f4b392252d
md"## Surface plot"

# ╔═╡ 2b819043-697b-482c-9cf3-fbb5f1fada5f
begin
	sh = pick_snapshot(out_folder[model], :recent) |> last
	isurf = MUST.closest(log10.(MUST.axis(sh, :τ_ross, 3)), 0)	
	x, y = MUST.axis(sh, :x) ./1e8 , MUST.axis(sh, :y) ./ 1e8
	heatmap(x, y, sh[:uz][:, :, isurf] ./1e5, cmap=palette(:hot, rev=false))
end

# ╔═╡ 65ac391a-9d10-46ef-8001-a81f2a3f2dac
md"## Convert to M3DIS models"

# ╔═╡ e384f1f0-12ef-473c-b038-648bd4d43938
res = [] # reshape(collect(Iterators.product((10, 20, 80), (299))), :)

# ╔═╡ d86bda77-f24f-49df-851e-f964c1a99eab
for (x, z) in res
	for i in eachindex(names)
		must2multi.snaps2multi(
			out_folder[i], 
			[range(-8, -1)..., :recent]...,
			eos=eos[i], 
			label="magg22_$(x)x$(x)x$(z)",
			n_horizontal=x,
			n_vertical=z,
			outfolder=labels[i]
		)
	end
end

# ╔═╡ Cell order:
# ╟─309e7de7-5a70-41f6-8ca6-47ad4b46c214
# ╠═dcd11d0e-312d-11ee-1d21-4dc68639a432
# ╠═41269b8e-5a49-4930-84de-aa4248531174
# ╟─7494384f-0a14-4ad0-95ad-e8c95fab78b7
# ╠═77153ffb-00e4-449b-99b2-b377a2ce2f7a
# ╟─d9a0ea6c-9197-4633-97d5-57ada5921ce4
# ╠═a66ad3b4-faed-41c1-af37-9cd96bb9aa41
# ╠═6202080d-fa73-4f39-8232-914582f70d11
# ╠═a4a4f729-f247-440d-80d2-455b30c0f63f
# ╠═36b50d99-f828-42a6-a617-3ff68feb36b9
# ╠═46402de2-102a-41ed-ac85-2c65be8d87d4
# ╠═ad1c4dcc-817d-4781-b23f-ed63a640f444
# ╠═b1e8f459-ed60-40ab-a8c3-c8d4c1d8fa5a
# ╠═2431f1ac-80cb-47d9-a657-7c4569ea6dfd
# ╟─ab08661d-fc62-4700-b423-9b96df85d9cb
# ╠═42a0408c-20cd-49b0-8581-7a9d13e2361b
# ╠═daec2515-64db-45ef-8a29-ab4fc8298963
# ╟─cbf793b6-c2ba-4872-9f4c-522c3a77d5a5
# ╟─2751d10c-3e48-4578-b330-97692eb25133
# ╟─724239e6-aac8-4096-907e-e87deaded5c6
# ╟─3239e571-f59e-4648-87ff-d6003d4128b9
# ╟─bc94c2a0-dbb3-4906-80d6-3bfdccd5dbb5
# ╟─a8dc248c-32be-452e-94bc-ca43f6afbda7
# ╟─42fcf382-b831-4378-bf97-f1ee6601d388
# ╟─7a945e98-c82d-48d8-978c-a5f4b392252d
# ╟─2b819043-697b-482c-9cf3-fbb5f1fada5f
# ╟─65ac391a-9d10-46ef-8001-a81f2a3f2dac
# ╠═e384f1f0-12ef-473c-b038-648bd4d43938
# ╠═d86bda77-f24f-49df-851e-f964c1a99eab
