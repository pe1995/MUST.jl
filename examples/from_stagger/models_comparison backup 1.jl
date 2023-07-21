### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1eac6457-238f-49ed-82dc-6f1c521930fc
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using Plots
	using TSO
	using LaTeXStrings
end;

# ╔═╡ 9456066d-36b1-4e6d-8af9-b8f134fb6e24
md"# Investigating different Models"

# ╔═╡ 931f7a1f-dccb-4727-9982-d04be9ffd688
begin
	mean = MUST.mean
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
end;

# ╔═╡ ee6a633b-6ccc-412f-ba74-e105aa148afc
MUST.@get_help_py gifs

# ╔═╡ 2e04240c-981e-443f-826b-2b40eb4b1974
md"All the models from the Stagger grid come with an average model that we can load as well. We can load the EoS from the same folder."

# ╔═╡ 486f5cac-8879-4084-9faf-6dc4ae115e43
names = [
		"DIS_MARCS_E_t5777g44m00_v0.1"
]

# ╔═╡ 9a7ce3c5-45f6-4589-a838-daaddf89e94f
out_folder = [
		MUST.@in_dispatch("data/sun_magg")
]

# ╔═╡ 82a51f3d-9e49-44ab-ae36-0069b6bd405c
eos_folder = [
		MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35")
]

# ╔═╡ fe1d7b10-88a5-46c1-a244-589bacf75970
labels = ["new setup (Magg)"]

# ╔═╡ ee39604b-6bd0-434e-b06d-417a4ab8cb7e
colors = ["red"] #palette(:rainbow, length(names))

# ╔═╡ 5856ad8f-b6ce-4175-a158-c415bd546a7e
in_folder  = [MUST.@in_dispatch "input_data/$(name)" for name in names]

# ╔═╡ 9ccd37a5-26ff-43b7-89a5-a214cf5995d3


# ╔═╡ 452a144e-b725-4de2-b3e1-2f614210d62e
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = MUST.pick_snapshot(
			MUST.converted_snapshots(
				out_folder[i]
			),
			:time_average
		)
		
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
	@optical(Average3D(eos[i], MUST.@in_dispatch("input_data/sun_stagger.dat")), eos[i], opa[i])
		for i in eachindex(eos)
]

# ╔═╡ 4c8d357d-a41a-4274-b131-93ebb695b911
stagger_model = [
	@optical(Average3D(eos[i], MUST.@in_dispatch("input_data/sun_stagger.dat")), 
		eos[i], opa[i])
	for i in eachindex(eos)
]

# ╔═╡ 4a1730fb-1fbe-445a-b850-7539967dcbe2
begin
	folder_stagger = "/u/peitner/DISPATCH/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
							folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t", 
							folder=folder_stagger)
end

# ╔═╡ 14dab588-7275-4a7e-bade-71375d3a16bc
for i in eachindex(snapshots)
	z = MUST.axis(snapshots[i], :z) ./1e5
	@info "vertical resolution of $(labels[i]): $(first(diff(z))) km ($(length(z)))"
end

# ╔═╡ 0d72ddb5-0096-4430-8611-e843f0aa3c61
for i in eachindex(snapshots)
	@info "simulation time of $(labels[i]): $(snapshots[i].parameter.time)"
end

# ╔═╡ c3ab29bc-a61b-467c-b3ab-3fa919abd83d
md"## Compare models to their initial condition"

# ╔═╡ fa90d92b-d1b3-491e-872b-23f57da6dace
begin
	Plots.default(fontfamily = ("Courier"), titlefont=("Courier"))
	
	function copy_ticks(sp::Plots.Subplot; minorticks=10)
		ptx = twinx(sp)
		plot!(ptx,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
		
		pty = twiny(sp)
		plot!(pty,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
	end
	
	copy_ticks(plt::Plots.Plot = current(); minorticks=10) = copy_ticks(plt[1], minorticks=minorticks)
	
	
	function basic_plot!(plot::Plots.Plot=current(); minorticks=10, 
						tickfontsize=8, 
						legendfontsize=10, 
						guidefontsize=12, 
						size=(600,400), 
						lm=2, bm=2, tm=2, rm=2, copyticks=true)
		plot!(plot, formatter=:auto)
		
		if copyticks
			copy_ticks(plot, minorticks=15)
		end
		
		plot!(plot, framestyle=:box, 
				minorticks=minorticks, 
				tickfontsize=tickfontsize, 
				legendfontsize=legendfontsize, 
				titlefontsize=legendfontsize,
				grid=false, 
				format=false,
				foreground_color_legend=nothing, 
				size=size, 
				leftmargin=lm*Plots.mm, rightmargin=rm*Plots.mm, 
				bottommargin=bm*Plots.mm, topmargin=tm*Plots.mm, 
				guidefontsize=guidefontsize)

		plot!()
	end;
end

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
	plot(-initial_model[1].z, exp.(initial_model[1].lnT), 
				lw=1.5, color=:black, label="initial condition", ls=:dash)

	for (i, snapshot) in enumerate(snapshots)
		plot!(profile(mean, snapshot, :z, :T)..., 
				lw=2., color=colors[i], label=labels[i])
	end
	

	xlabel!("z [cm]")
	ylabel!("T [K]")
	
	basic_plot!()
end

# ╔═╡ d3ef4193-11e7-478d-aa11-ba3b8adbce55
begin
	plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		plot!(profile(mean, snapshot_τ, :log10τ_ross, :T)..., 
				lw=2., color=colors[i], label=labels[i], ls=:solid)
	end


	xlabel!("τ-ross [log]")
	ylabel!("T [K]")
	
	basic_plot!()
end

# ╔═╡ 05a4e3ef-58b5-4f96-94b8-51991c971451
begin
	plot(profile(mean, stagger, :z, :d)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

	for (i, snapshot) in enumerate(snapshots)
		plot!(profile(mean, snapshot, :z, :d)..., 
				lw=2., color=colors[i], label=labels[i])
	end

	yaxis!(:log)
	xlabel!("z [cm]")
	ylabel!("density [g x cm-3]")
	
	basic_plot!()
end

# ╔═╡ 0a726cbb-1309-4071-9df1-93cd4355fb71
begin
	plot(profile(mean, stagger_τ, :log10τ_ross, :log10d)..., 
				lw=2., color=:black, label="Stagger", ls=:dash)
	
	for (i, snapshot_τ) in enumerate(snapshots_τ)
		plot!(profile(mean, snapshot_τ, :log10τ_ross, :log10d)..., 
				lw=1.5, color=colors[i], label=labels[i])
	end
	

	xlabel!("τ-ross [log]")
	ylabel!("density [g x cm-3]")
	
	basic_plot!()
end

# ╔═╡ e10a8581-c0ab-4209-9e14-4d456dcf9a86
begin
	_, dStagger = profile(mean, stagger, :z, :d)
	_, TStagger = profile(mean, stagger, :z, :T)
	plot(dStagger, TStagger, lw=1.5, color=:black, label="Stagger", ls=:dash)
	
	for (i, snapshot) in enumerate(snapshots)
		_, d = profile(mean, snapshot, :z, :d)
		_, T = profile(mean, snapshot, :z, :T)
		
		plot!(d, T, lw=2., color=colors[i], label=labels[i])
	end
	
	xaxis!(:log)	
	
	xlabel!("density [g x cm-3]")
	ylabel!("T [K]")
	
	basic_plot!()
	plot!(legend=:topleft)
end

# ╔═╡ 07815fd8-f292-4760-a950-8b56a5908acf
md"Velocity distribution"

# ╔═╡ fa8e10f7-da54-4165-887f-30e740e1f264
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ 870f696c-483b-40ef-9986-8f81ed01554c
8*35

# ╔═╡ 83dc28ac-97f5-4833-8a6f-cd2d2624442a
begin
	plot(profile(rms, stagger_τ, :log10τ_ross, :uz)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		plot!(profile(rms, snapshot_τ, :log10τ_ross, :uz)..., 
				lw=2., color=colors[i], label=labels[i], ls=:solid)
	end


	xlabel!("τ-ross [log]")
	ylabel!("rms(Uz) [cm × s-1]")
	
	basic_plot!()
end

# ╔═╡ e654ada6-ce50-4a2b-a51c-eae3e7aa36d3
begin
	stephist(reshape(stagger[:uz], :) ./1e5, 
				lw=2., color=:black, label="Stagger", 
				normalize=:probability)
	
	for (i, snapshot) in enumerate(snapshots)
		stephist!(reshape(snapshot[:uz], :) ./1e5, 
				lw=1.5, color=colors[i], label=labels[i], 
				normalize=:probability)
	end
	

	xlabel!("Uz [km × s-1]")
	ylabel!("density function")

	plot!(xlim=(-8, 8))
	
	basic_plot!()
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

# ╔═╡ e2dde825-4e66-42e2-b541-74ec0600fc81
begin
	h = []
	for (i, snapshot) in enumerate(snapshots)
		append!(h, [heatmap(MUST.axis(snapshot, :x)./1e8, 
					MUST.axis(snapshot, :y)./1e8,
					uz[i] ./1e5,
						colormap=:hot, 
						clim=(-8.0, 8.0), title="\n"*labels[i])])

		#xlabel!("x [Mm]")
		#ylabel!("y [Mm]")
	
	
		basic_plot!(copyticks=false, bm=3, lm=0, rm=0, tm=3, size=(600,500))
	end
	
	plot(h..., link=:both, size=(600, 500))
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
	
	#MUST.multiBox(model_box, output_names[i], downsample_xy=downsample)
	
	@info "New size: $(size(@view(ne[i][1:downsample:end, 1:downsample:end, :])))"
end

# ╔═╡ 5b28c9b8-a991-45d6-a2d1-15f24142d277
md"## Compute a time average"

# ╔═╡ 7837df1e-e220-4a60-944d-876b523300ad
number_of_timeslots = [5]

# ╔═╡ 09d3fbbb-c636-4f8a-8409-0eba2ea9a242
do_τ = true

# ╔═╡ Cell order:
# ╟─9456066d-36b1-4e6d-8af9-b8f134fb6e24
# ╠═1eac6457-238f-49ed-82dc-6f1c521930fc
# ╠═931f7a1f-dccb-4727-9982-d04be9ffd688
# ╠═ee6a633b-6ccc-412f-ba74-e105aa148afc
# ╟─2e04240c-981e-443f-826b-2b40eb4b1974
# ╠═486f5cac-8879-4084-9faf-6dc4ae115e43
# ╠═9a7ce3c5-45f6-4589-a838-daaddf89e94f
# ╠═82a51f3d-9e49-44ab-ae36-0069b6bd405c
# ╠═fe1d7b10-88a5-46c1-a244-589bacf75970
# ╠═ee39604b-6bd0-434e-b06d-417a4ab8cb7e
# ╟─5856ad8f-b6ce-4175-a158-c415bd546a7e
# ╟─957af12a-e77e-48a3-a2de-80b86512e5a8
# ╟─1e0de50e-bb67-4d34-a038-e7437955ec73
# ╟─01ab2753-515c-496f-a6fc-1c2a0e42ae25
# ╠═9ccd37a5-26ff-43b7-89a5-a214cf5995d3
# ╠═452a144e-b725-4de2-b3e1-2f614210d62e
# ╟─3e747391-ba4b-47bf-b363-abcb46a9309b
# ╟─ed6c250a-84d5-4ac6-bf54-9e8fcdbe55a3
# ╠═d8693137-82f7-4ccb-b886-4115e3032392
# ╠═4c8d357d-a41a-4274-b131-93ebb695b911
# ╠═4a1730fb-1fbe-445a-b850-7539967dcbe2
# ╟─14dab588-7275-4a7e-bade-71375d3a16bc
# ╟─0d72ddb5-0096-4430-8611-e843f0aa3c61
# ╟─c3ab29bc-a61b-467c-b3ab-3fa919abd83d
# ╟─fa90d92b-d1b3-491e-872b-23f57da6dace
# ╟─9fd5896c-8d4f-499d-9f97-c589d8d256c2
# ╟─85d4b142-1529-495a-bc6d-64f5f0efa3b9
# ╟─3d1e1cd9-5601-468d-a3f1-a3dd51177a37
# ╠═8a8a40e0-5b03-4c17-93a3-4eccf12e3717
# ╟─d3ef4193-11e7-478d-aa11-ba3b8adbce55
# ╟─05a4e3ef-58b5-4f96-94b8-51991c971451
# ╟─0a726cbb-1309-4071-9df1-93cd4355fb71
# ╟─e10a8581-c0ab-4209-9e14-4d456dcf9a86
# ╟─07815fd8-f292-4760-a950-8b56a5908acf
# ╠═fa8e10f7-da54-4165-887f-30e740e1f264
# ╠═870f696c-483b-40ef-9986-8f81ed01554c
# ╟─83dc28ac-97f5-4833-8a6f-cd2d2624442a
# ╟─e654ada6-ce50-4a2b-a51c-eae3e7aa36d3
# ╟─919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
# ╟─84bd1e80-ff21-4970-a5a3-fc7452da7e6f
# ╟─e2dde825-4e66-42e2-b541-74ec0600fc81
# ╟─5ca0a8ee-5fbf-4d76-8bfd-91c292e08cfa
# ╠═7a023be5-46ea-4c25-857b-3f765c044a91
# ╠═b1df1501-4033-4d8a-bd67-8130c095152a
# ╠═fdf3a692-daab-49d5-8bf6-36996617349e
# ╟─f3e705cf-a0a5-4905-9a07-aa6535985e01
# ╟─52fa7b28-d846-4d1c-8142-1cc1e6c68b92
# ╠═b3b294ca-9ac2-4083-bcb5-789b48e9cb0b
# ╟─5b28c9b8-a991-45d6-a2d1-15f24142d277
# ╠═7837df1e-e220-4a60-944d-876b523300ad
# ╠═09d3fbbb-c636-4f8a-8409-0eba2ea9a242
# ╠═7d192f93-c865-43af-a8d4-be92fb711a4e
