### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 1eac6457-238f-49ed-82dc-6f1c521930fc
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using Glob
	using PyPlot
	using TSO
	using LaTeXStrings
end;

# ╔═╡ e21ee951-3720-43f2-8c92-0e7907f76b3d
include_helper(name) = include(joinpath(dirname(pathof(MUST)), name))

# ╔═╡ 9456066d-36b1-4e6d-8af9-b8f134fb6e24
md"# Investigating different Models"

# ╔═╡ 468c1043-7426-457e-80d3-cd8e1b2c8e1b
begin
	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["font.size"] = 10
end

# ╔═╡ 931f7a1f-dccb-4727-9982-d04be9ffd688
begin
	mean = MUST.mean
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
end;

# ╔═╡ ee6a633b-6ccc-412f-ba74-e105aa148afc
MUST.@get_help_py gifs

# ╔═╡ fbe723d5-12de-4c20-9545-7a3e398d4546
include_helper("visual.jl")

# ╔═╡ 2e04240c-981e-443f-826b-2b40eb4b1974
md"All the models from the Stagger grid come with an average model that we can load as well. We can load the EoS from the same folder."

# ╔═╡ 486f5cac-8879-4084-9faf-6dc4ae115e43
names = [
		"DIS_MARCS_E_t5777g44m00_v0.1", 
		"DIS_MARCS_E_t5777g44m00_v0.1", 
		"DIS_MARCS_E_t5777g44m00_v0.1"
]

# ╔═╡ 9a7ce3c5-45f6-4589-a838-daaddf89e94f
out_folder = [
		MUST.@in_dispatch(
			"data_tmps/bin_tests/data/DIS_MARCS_E_t5777g44m00_v0.1_test"),
		MUST.@in_dispatch(
			"data_tmps/bin_tests/data/DIS_MARCS_E_t5777g44m00_v0.1_test2"),
		MUST.@in_dispatch(
			"data_tmps/bin_tests/data/DIS_MARCS_E_t5777g44m00_v0.1")
]

# ╔═╡ 82a51f3d-9e49-44ab-ae36-0069b6bd405c
eos_folder = [
		MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.13"),
		MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.33"),
		MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.31") 
]

# ╔═╡ fe1d7b10-88a5-46c1-a244-589bacf75970
labels = ["5 bins", "10 bins", "12 bins" ]

# ╔═╡ ee39604b-6bd0-434e-b06d-417a4ab8cb7e
colors = ["magenta", "cyan", "lime"]

# ╔═╡ 5856ad8f-b6ce-4175-a158-c415bd546a7e
in_folder  = [MUST.@in_dispatch "input_data/grd/$(name)" for name in names]

# ╔═╡ 957af12a-e77e-48a3-a2de-80b86512e5a8
function converted_snapshots(folder)
	files_converted = glob("*.hdf5", folder)
	snaps = Dict()
	for file in files_converted
		snname = basename(file)

		
		if occursin("tau", snname) 
			continue 
		end 
		
		snid   = parse(Int, snname[last(findfirst("sn", snname))+1:end-5])
		is_τ = isfile(joinpath(folder,"box_tau_sn$(snid).hdf5"))

		snname  = String(split(snname, ".hdf5") |> first)
		sntname = "box_tau_sn$(snid)"

		if is_τ 
			snaps[snid] = (snname, sntname) 
		else
			snaps[snid] = (snname, nothing) 
		end
	end

	snaps["folder"] = folder
	snaps
end

# ╔═╡ 1e0de50e-bb67-4d34-a038-e7437955ec73
list_snapshots(snapshots) = sort([k for k in keys(snapshots) if k != "folder"])

# ╔═╡ 01ab2753-515c-496f-a6fc-1c2a0e42ae25
function pick_snapshot(snapshots, i)
	i = if i == :recent 
		last(list_snapshots(snapshots))
	else 
		i
	end

	snap = snapshots[i]
	if isnothing(last(snap))
		@info "snapshot $(i) loaded."
		MUST.Box(first(snap), folder=snapshots["folder"]), nothing
	else
		@info "snapshot $(i) + τ-shot loaded."
		MUST.Box(first(snap), folder=snapshots["folder"]), 
		MUST.Box(last(snap), folder=snapshots["folder"])
	end
end

# ╔═╡ 452a144e-b725-4de2-b3e1-2f614210d62e
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = pick_snapshot(converted_snapshots(out_folder[i]),
									:recent)
		
		append!(snapshots, [snapshot])
		append!(snapshots_τ, [snapshot_τ])
	end
end

# ╔═╡ 3e747391-ba4b-47bf-b363-abcb46a9309b
begin
	eos = [reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5")) 
				for i in eachindex(eos_folder)]
	opa = [reload(SqOpacity, joinpath(eos_folder[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder)]
end

# ╔═╡ ed6c250a-84d5-4ac6-bf54-9e8fcdbe55a3
for i in eachindex(eos)
	@info "Opacity table size: $(size(opa[i].κ))"
end

# ╔═╡ d8693137-82f7-4ccb-b886-4115e3032392
# ╠═╡ show_logs = false
initial_model = [
	@optical(Average3D(eos[i], joinpath(in_folder[i], "inim.dat")), eos[i], opa[i])
		for i in eachindex(eos)
]

# ╔═╡ 4c8d357d-a41a-4274-b131-93ebb695b911
stagger_model = [
	@optical(Average3D(eos[i], MUST.@in_dispatch("input_data/solar_stagger_ext")), 
		eos[i], opa[i])
	for i in eachindex(eos)
]

# ╔═╡ 4a1730fb-1fbe-445a-b850-7539967dcbe2
begin
	folder_stagger = "/ptmp/peitner/model_grid/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
							folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t", 
							folder=folder_stagger)
end

# ╔═╡ c3ab29bc-a61b-467c-b3ab-3fa919abd83d
md"## Compare models to their initial condition"

# ╔═╡ 9fd5896c-8d4f-499d-9f97-c589d8d256c2
mesh(m::MUST.Box) = MUST.meshgrid(MUST.axis(m, :x) ./1e8, 
									MUST.axis(m, :y) ./1e8, 
									MUST.axis(m, :z) ./1e8)


# ╔═╡ 85d4b142-1529-495a-bc6d-64f5f0efa3b9
is_log(x) = begin
	sx = String(x)
	
	xnew, f = if occursin("log10", sx)
		Symbol(sx[findfirst("log10", sx)[end]+1:end]), log10
	elseif occursin("log", sx)
		Symbol(sx[findfirst("log", sx)[end]+1:end]), log
	else
		x, identity
	end

	xnew, f
end

# ╔═╡ 3d1e1cd9-5601-468d-a3f1-a3dd51177a37
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if xs == :τ_ross
		logx.(MUST.axis(model, xs, 3)), logy.(MUST.plane_statistic(f, model, ys)) 
	else
		logx.(MUST.axis(model, xs)), logy.(MUST.plane_statistic(f, model, ys))
	end
end

# ╔═╡ 8a8a40e0-5b03-4c17-93a3-4eccf12e3717
begin
	f, ax = plt.subplots(1, 1)
	basic_plot!(ax)
	
	ax.plot(profile(mean, stagger, :z, :T)..., 
				lw=1.5, color="black", label="Stagger", ls="--")

	for (i, snapshot) in enumerate(snapshots)
		ax.plot(profile(mean, snapshot, :z, :T)..., 
				lw=2., color=colors[i], label=labels[i])
	end
	

	ax.set_xlabel("z [cm]")
	ax.set_ylabel("T [K]")

	ax.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ d3ef4193-11e7-478d-aa11-ba3b8adbce55
begin
	f2, ax2 = plt.subplots(1, 1)
	basic_plot!(ax2)
	
	ax2.plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
				lw=1.5, color="black", label="Stagger", ls="--")

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		ax2.plot(profile(mean, snapshot_τ, :log10τ_ross, :T)..., 
				lw=2., color=colors[i], label=labels[i])
	end


	ax2.set_xlabel(L"\rm \tau_{ross}\ [log]")
	ax2.set_ylabel("T [K]")

	ax2.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ 05a4e3ef-58b5-4f96-94b8-51991c971451
begin
	f3, ax3 = plt.subplots(1, 1)
	basic_plot!(ax3)
	
	ax3.plot(profile(mean, stagger, :z, :d)..., 
				lw=1.5, color=:black, label="Stagger", ls="--")

	for (i, snapshot) in enumerate(snapshots)
		ax3.plot(profile(mean, snapshot, :z, :d)..., 
				lw=2., color=colors[i], label=labels[i])
	end

	ax3.set_yscale("log")
	ax3.set_xlabel("z [cm]")
	ax3.set_ylabel(L"\rm density\ [g \times cm^{-3}]")

	ax3.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ 0a726cbb-1309-4071-9df1-93cd4355fb71
begin
	f4, ax4 = plt.subplots(1, 1)
	basic_plot!(ax4)
	
	ax4.plot(profile(mean, stagger_τ, :log10τ_ross, :log10d)..., 
				lw=2., color="black", label="Stagger", ls="--")
	
	for (i, snapshot_τ) in enumerate(snapshots_τ)
		ax4.plot(profile(mean, snapshot_τ, :log10τ_ross, :log10d)..., 
				lw=1.5, color=colors[i], label=labels[i])
	end
	

	ax4.set_xlabel(L"\rm \tau_{ross}\ [log]")
	ax4.set_ylabel(L"\rm density\ [g \times cm^{-3}]")

	ax4.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ e10a8581-c0ab-4209-9e14-4d456dcf9a86
begin
	f5, ax5 = plt.subplots(1, 1)
	basic_plot!(ax5)
	
	_, dStagger = profile(mean, stagger, :z, :d)
	_, TStagger = profile(mean, stagger, :z, :T)
	ax5.plot(dStagger, TStagger, lw=1.5, color="black", label="Stagger", ls="--")
	
	for (i, snapshot) in enumerate(snapshots)
		_, d = profile(mean, snapshot, :z, :d)
		_, T = profile(mean, snapshot, :z, :T)
		
		ax5.plot(d, T, lw=2., color=colors[i], label=labels[i])
	end
	
	ax5.set_xscale("log")	
	
	ax5.set_xlabel(L"\rm density\ [g \times cm^{-3}]")
	ax5.set_ylabel("T [K]")

	ax5.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ 07815fd8-f292-4760-a950-8b56a5908acf
md"Velocity distribution"

# ╔═╡ fa8e10f7-da54-4165-887f-30e740e1f264
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ 83dc28ac-97f5-4833-8a6f-cd2d2624442a
begin
	f6, ax6 = plt.subplots(1, 1)
	basic_plot!(ax6)
	
	ax6.plot(profile(rms, stagger_τ, :log10τ_ross, :uz)..., 
				lw=1.5, color="black", label="Stagger", ls="--")

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		ax6.plot(profile(rms, snapshot_τ, :log10τ_ross, :uz)..., 
				lw=2., color=colors[i], label=labels[i])
	end


	ax6.set_xlabel("τ-ross [log]")
	ax6.set_ylabel("rms(Uz) [cm × s-1]")
	
end

# ╔═╡ e654ada6-ce50-4a2b-a51c-eae3e7aa36d3
begin
	f7, ax7 = plt.subplots(1, 1)
	basic_plot!(ax7)
	
	ax7.hist(reshape(stagger[:uz], :) ./1e5, 
				lw=2., color="black", label="Stagger", ls="-", histtype="step",
				density=true, bins=500)
	
	for (i, snapshot) in enumerate(snapshots)
		ax7.hist(reshape(snapshot[:uz], :) ./1e5, 
				lw=1.5, color=colors[i], label=labels[i], histtype="step",
				density=true, bins=500)
	end
	

	ax7.set_xlabel(L"\rm Uz\ [km \times s^{-1}]")
	ax7.set_ylabel("density function")

	ax7.set_xlim(-8, 8)

	ax7.legend(framealpha=0, loc="best")
	
	gcf()
end

# ╔═╡ 919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
md"## 3D result"

# ╔═╡ 84bd1e80-ff21-4970-a5a3-fc7452da7e6f
begin
	uz = []
	for (i,snapshot_τ) in enumerate(snapshots_τ)
		uz_τ_surf = if !isnothing(snapshot_τ)
			isurf = MUST.closest(log10.(MUST.axis(snapshot_τ, :τ_ross, 3)), 0)
			snapshot_τ[:uz][:, :, isurf]
		else
			@info "Interpolating Uz..."
			uz_τ_surf = MUST.interpolate_to(snapshots[i], :uz, τ_ross=1)
			MUST.axis(uz_τ_surf, :uz, 3)
		end

		append!(uz, [uz_τ_surf])
	end

	uz
end

# ╔═╡ 5ca0a8ee-5fbf-4d76-8bfd-91c292e08cfa
md"## Convert to M3D format"

# ╔═╡ 7a023be5-46ea-4c25-857b-3f765c044a91
labels_new = replace.(labels, " "=>"_")

# ╔═╡ b1df1501-4033-4d8a-bd67-8130c095152a
downsample=40

# ╔═╡ fdf3a692-daab-49d5-8bf6-36996617349e
output_names = ["m3dis_sun_$(label)" for (name, label) in zip(names, labels_new)]

# ╔═╡ f3e705cf-a0a5-4905-9a07-aa6535985e01
aos = [@axed(eos[i]) for i in eachindex(eos)]

# ╔═╡ 52fa7b28-d846-4d1c-8142-1cc1e6c68b92
ne = [lookup(aos[i], :lnNe, log.(model_box[:d]), log.(model_box[:ee])) 
		for (i, model_box) in enumerate(snapshots)]

# ╔═╡ b3b294ca-9ac2-4083-bcb5-789b48e9cb0b
for (i, model_box) in enumerate(snapshots)
	model_box.data[:ne] = exp.(ne[i])
	
	MUST.multiBox(model_box, output_names[i], downsample_xy=downsample)
	
	@info "New size: $(size(@view(ne[i][1:downsample:end, 1:downsample:end, :])))"
end

# ╔═╡ Cell order:
# ╟─9456066d-36b1-4e6d-8af9-b8f134fb6e24
# ╠═1eac6457-238f-49ed-82dc-6f1c521930fc
# ╠═468c1043-7426-457e-80d3-cd8e1b2c8e1b
# ╠═931f7a1f-dccb-4727-9982-d04be9ffd688
# ╠═ee6a633b-6ccc-412f-ba74-e105aa148afc
# ╠═e21ee951-3720-43f2-8c92-0e7907f76b3d
# ╠═fbe723d5-12de-4c20-9545-7a3e398d4546
# ╟─2e04240c-981e-443f-826b-2b40eb4b1974
# ╠═486f5cac-8879-4084-9faf-6dc4ae115e43
# ╠═9a7ce3c5-45f6-4589-a838-daaddf89e94f
# ╠═82a51f3d-9e49-44ab-ae36-0069b6bd405c
# ╠═fe1d7b10-88a5-46c1-a244-589bacf75970
# ╟─ee39604b-6bd0-434e-b06d-417a4ab8cb7e
# ╟─5856ad8f-b6ce-4175-a158-c415bd546a7e
# ╟─957af12a-e77e-48a3-a2de-80b86512e5a8
# ╟─1e0de50e-bb67-4d34-a038-e7437955ec73
# ╟─01ab2753-515c-496f-a6fc-1c2a0e42ae25
# ╟─452a144e-b725-4de2-b3e1-2f614210d62e
# ╟─3e747391-ba4b-47bf-b363-abcb46a9309b
# ╟─ed6c250a-84d5-4ac6-bf54-9e8fcdbe55a3
# ╟─d8693137-82f7-4ccb-b886-4115e3032392
# ╟─4c8d357d-a41a-4274-b131-93ebb695b911
# ╠═4a1730fb-1fbe-445a-b850-7539967dcbe2
# ╟─c3ab29bc-a61b-467c-b3ab-3fa919abd83d
# ╟─9fd5896c-8d4f-499d-9f97-c589d8d256c2
# ╟─85d4b142-1529-495a-bc6d-64f5f0efa3b9
# ╟─3d1e1cd9-5601-468d-a3f1-a3dd51177a37
# ╟─8a8a40e0-5b03-4c17-93a3-4eccf12e3717
# ╟─d3ef4193-11e7-478d-aa11-ba3b8adbce55
# ╟─05a4e3ef-58b5-4f96-94b8-51991c971451
# ╟─0a726cbb-1309-4071-9df1-93cd4355fb71
# ╟─e10a8581-c0ab-4209-9e14-4d456dcf9a86
# ╟─07815fd8-f292-4760-a950-8b56a5908acf
# ╠═fa8e10f7-da54-4165-887f-30e740e1f264
# ╟─83dc28ac-97f5-4833-8a6f-cd2d2624442a
# ╟─e654ada6-ce50-4a2b-a51c-eae3e7aa36d3
# ╟─919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
# ╠═84bd1e80-ff21-4970-a5a3-fc7452da7e6f
# ╟─5ca0a8ee-5fbf-4d76-8bfd-91c292e08cfa
# ╠═7a023be5-46ea-4c25-857b-3f765c044a91
# ╠═b1df1501-4033-4d8a-bd67-8130c095152a
# ╠═fdf3a692-daab-49d5-8bf6-36996617349e
# ╟─f3e705cf-a0a5-4905-9a07-aa6535985e01
# ╟─52fa7b28-d846-4d1c-8142-1cc1e6c68b92
# ╟─b3b294ca-9ac2-4083-bcb5-789b48e9cb0b
