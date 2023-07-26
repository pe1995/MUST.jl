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
	using DelimitedFiles
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
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1"
]

# ╔═╡ 9a7ce3c5-45f6-4589-a838-daaddf89e94f
out_folder = [
	MUST.@in_dispatch("data/sun_new_magg_vres"),
	MUST.@in_dispatch("data/sun_new_magg_lres"),
	MUST.@in_dispatch("data/sun_new_magg_ires"),
	MUST.@in_dispatch("data/pretty_good_sun_new_magg5_rapidf"),	
	MUST.@in_dispatch("data/sun_new_magg_hres"),
]

# ╔═╡ 82a51f3d-9e49-44ab-ae36-0069b6bd405c
eos_folder = [
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3")
]

# ╔═╡ fe1d7b10-88a5-46c1-a244-589bacf75970
labels = [
	"magg2022_60x60",
	"magg2022_90x90",
	"magg2022_90x180",
	"magg2022_150x300",
	"magg2022_195x390"
]

# ╔═╡ ee39604b-6bd0-434e-b06d-417a4ab8cb7e
colors = length(names) == 1 ? [:red] : palette(:roma, length(names))

# ╔═╡ 5856ad8f-b6ce-4175-a158-c415bd546a7e
in_folder  = [MUST.@in_dispatch "input_data/$(name)" for name in names]

# ╔═╡ 452a144e-b725-4de2-b3e1-2f614210d62e
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = MUST.pick_snapshot(
			MUST.converted_snapshots(out_folder[i]), -2
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

# ╔═╡ c72cea6d-59bd-4721-8044-9fa28fb4039e
marcs_model = readdlm("marcs_sun.txt", skipstart=1);

# ╔═╡ 4a1730fb-1fbe-445a-b850-7539967dcbe2
begin
	folder_stagger = "/u/peitner/DISPATCH/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
							folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t", 
							folder=folder_stagger)
end

# ╔═╡ 9faf1a6b-a027-4fcc-b572-872624b3f7e8
begin
	folder_muram = "/u/peitner/DISPATCH/examples/Muram"
	muram = MUST.Box("box_MURaM_cube_small.221000_HDSun", folder=folder_muram)
	muram_τ = MUST.Box("box_MURaM_cube_small.221000_HDSun_t", folder=folder_muram)
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

# ╔═╡ 7c8ed29b-2c7a-4b4c-a60d-4de3cf082b42
teff_last = MUST.read_teff(joinpath(first(out_folder), "teff.dat"));

# ╔═╡ a0220545-f149-4862-a8ba-3a1911cb732a
#plot(teff_last[:, 1], teff_last[:, 2], xlabel="time", ylabel="Teff")

# ╔═╡ c3ab29bc-a61b-467c-b3ab-3fa919abd83d
md"## Compare models to their initial condition"

# ╔═╡ cd81869e-c954-4526-a52f-31aedb64cb5f
2.3/(7*35)

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

	plot!(marcs_model[:, 2], marcs_model[:, 5], 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS")

	xlabel!("τ-ross [log]")
	ylabel!("T [K]")
	
	basic_plot!()
end

# ╔═╡ 05a4e3ef-58b5-4f96-94b8-51991c971451
begin
	plot(profile(mean, stagger, :z, :d)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

	plot!(legendforegroundcolor=nothing, legendbackgroundcolor=nothing,
	legendposition=:bottomleft)

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

	plot!(legendforegroundcolor=nothing, legendbackgroundcolor=nothing,
		legendposition=:bottomright)
	
	for (i, snapshot_τ) in enumerate(snapshots_τ)
		plot!(profile(mean, snapshot_τ, :log10τ_ross, :log10d)..., 
				lw=1.5, color=colors[i], label=labels[i])
	end
	

	plot!(marcs_model[:, 2], log10.(marcs_model[:, 11]), 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS")

	xlabel!("τ-ross [log]")
	ylabel!("density [g x cm-3]")
	
	basic_plot!()
end

# ╔═╡ e10a8581-c0ab-4209-9e14-4d456dcf9a86
begin
	_, dStagger = profile(mean, stagger, :z, :d)
	_, TStagger = profile(mean, stagger, :z, :T)
	plot(dStagger, TStagger, lw=1.5, color=:black, label="Stagger", ls=:dash)

	plot!(legendforegroundcolor=nothing, legendbackgroundcolor=nothing)

	for (i, snapshot) in enumerate(snapshots)
		_, d = profile(mean, snapshot, :z, :d)
		_, T = profile(mean, snapshot, :z, :T)
		
		plot!(d, T, lw=2., color=colors[i], label=labels[i])
	end

	plot!(marcs_model[:, 11], marcs_model[:, 5], 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS")
	
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

# ╔═╡ 7d273e41-207c-46b9-9e7d-a3626b5717c6
md"## Differences"

# ╔═╡ 9d6e7385-c8f6-4df1-91ac-3f07f0d2567f
begin
	common_tau = range(-3.5, 1.5, length=200) |> collect
	ip(x, y) = begin
		m = sortperm(x)
		MUST.linear_interpolation(x[m], y[m]).(common_tau)
	end
end

# ╔═╡ 6e078b79-7928-4854-a07f-c5ac185caa59
begin
	xs, ys = profile(mean, stagger_τ, :log10τ_ross, :T)
	
	plot(legendforegroundcolor=nothing, legendbackgroundcolor=nothing, legendposition=:bottomleft)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		xi, yi = profile(mean, snapshot_τ, :log10τ_ross, :T)
		plot!(common_tau, ip(xi, yi) .- ip(xs, ys), 
				lw=2., color=colors[i], label="$(labels[i]) - Stagger", ls=:solid)
	end

	x, y = marcs_model[:, 2], marcs_model[:, 5]
	plot!(common_tau, ip(x, y) .- ip(xs, ys), 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS - Stagger")

	xmu, ymu = profile(mean, muram_τ, :log10τ_ross, :T)
	plot!(common_tau, ip(xmu, ymu) .- ip(xs, ys), 
		color=:black, ls=:dash, lw=2.5, label="MURaM - Stagger")

	hline!([0.0], color=:black, ls=:dot, alpha=0.5, label=nothing)

	xlabel!("τ-ross [log]")
	ylabel!("ΔT [K]")

	plot!(ylim=(-600, 600))
	
	basic_plot!()
end

# ╔═╡ c0347b03-fb0a-4660-a9f6-fc24272263ed
begin	
	plot(legendforegroundcolor=nothing, legendbackgroundcolor=nothing)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		xi, yi = profile(mean, snapshot_τ, :log10τ_ross, :T)
		plot!(common_tau, (ip(xi, yi) .- ip(xs, ys)) ./ ip(xs, ys) *100, 
				lw=2., color=colors[i], label="$(labels[i]) - Stagger", ls=:solid)
	end

	plot!(common_tau, (ip(x, y) .- ip(xs, ys)) ./ ip(xs, ys) *100, 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS - Stagger")

	plot!(common_tau, (ip(xmu, ymu) .- ip(xs, ys)) ./ ip(xs, ys) *100, 
		color=:black, ls=:dash, lw=2.5, label="MURaM - Stagger")

	hline!([0.0], color=:black, ls=:dot, alpha=0.5, label=nothing)

	xlabel!("τ-ross [log]")
	ylabel!("ΔT / T  [%]")

	plot!(ylim=(-10, 10))
	
	basic_plot!()
end

# ╔═╡ fcc3aaac-a5e2-45fa-9419-4aaa916fbe45
begin
	ipr(x, y) = 10.0 .^(ip(x, y))
	
	xs2, ys2 = profile(mean, stagger_τ, :log10τ_ross, :d)
	
	plot(legendforegroundcolor=nothing, legendbackgroundcolor=nothing)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		xi2, yi2 = profile(mean, snapshot_τ, :log10τ_ross, :d)
		plot!(common_tau, ipr(xi2, log10.(yi2)) .- ipr(xs2, log10.(ys2)), 
				lw=2., color=colors[i], label="$(labels[i]) - Stagger", ls=:solid)
	end

	x2, y2 = marcs_model[:, 2], marcs_model[:, 11]
	plot!(common_tau, ipr(x2, log10.(y2)) .- ipr(xs2, log10.(ys2)), 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS - Stagger")

	xmu2, ymu2 = profile(mean, muram_τ, :log10τ_ross, :d)
	plot!(common_tau, ipr(xmu2, log10.(ymu2)) .- ipr(xs2, log10.(ys2)), 
		color=:black, ls=:dash, lw=2.5, label="MURaM - Stagger")

	hline!([0.0], color=:black, ls=:dot, alpha=0.5, label=nothing)

	#yaxis!(:log)	

	xlabel!("τ-ross [log]")
	ylabel!("Δρ [g × cm-3]")

	#plot!(ylim=(-600, 600))
	
	basic_plot!()
end

# ╔═╡ 4c0bf79b-9b7b-4c35-a728-350a79e5eb83
begin	
	plot(legendforegroundcolor=nothing, legendbackgroundcolor=nothing)

	for (i, snapshot_τ) in enumerate(snapshots_τ)
		xi2, yi2 = profile(mean, snapshot_τ, :log10τ_ross, :d)
		plot!(common_tau, 
			(ipr(xi2, log10.(yi2)) .- ipr(xs2, log10.(ys2))) ./ ipr(xs2, log10.(ys2))*100, 
				lw=2., color=colors[i], label="$(labels[i]) - Stagger", ls=:solid)
	end

	plot!(common_tau, 
		(ipr(x2, log10.(y2)) .- ipr(xs2, log10.(ys2))) ./ ipr(xs2, log10.(ys2)) *100, 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS - Stagger")

	plot!(common_tau, 
		(ipr(xmu2, log10.(ymu2)) .- ipr(xs2, log10.(ys2))) ./ ipr(xs2, log10.(ys2)) *100, 
		color=:black, ls=:dash, lw=2.5, label="MARCS - Stagger")

	hline!([0.0], color=:black, ls=:dot, alpha=0.5, label=nothing)

	xlabel!("τ-ross [log]")
	ylabel!("Δρ / ρ [%]")

	#yaxis!(:log)	

	#plot!(ylim=(-600, 600))
	
	basic_plot!()
end

# ╔═╡ a10e407c-70f8-4809-a0e4-eb334447d623
begin	
	plot(legendforegroundcolor=nothing, legendbackgroundcolor=nothing, legendposition=:bottomleft)

	snapsis = [MUST.pick_snapshot(
			MUST.converted_snapshots(out_folder[4]), j
	)[2]
		for j in [-8, -7, -6, -5, -4, -3, -2, -1, :recent]
	]
	
	for (i, snapshot_τ) in enumerate(snapsis)
		xi, yi = profile(mean, snapshot_τ, :log10τ_ross, :T)
		plot!(common_tau, ip(xi, yi) .- ip(xs, ys), 
				lw=2., ls=:solid)
	end

	
	plot!(common_tau, ip(x, y) .- ip(xs, ys), 
		color=:royalblue, ls=:dot, lw=2.5, label="MARCS - Stagger")

	plot!(common_tau, ip(xmu, ymu) .- ip(xs, ys), 
		color=:black, ls=:dash, lw=2.5, label="MURaM - Stagger")

	hline!([0.0], color=:black, ls=:dot, alpha=0.5, label=nothing)

	xlabel!("τ-ross [log]")
	ylabel!("ΔT [K]")

	plot!(ylim=(-600, 600), legendposition=:top, legendcolumns=2)
	
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
	
	plot(h..., link=:both, size=(2100, 1400))
end

# ╔═╡ 3615d990-9c9a-4e69-8f96-0e3c2ac6896b
md"## Upsample for M3D
With the `g-` methods one can upsample the cubes pretty easily"

# ╔═╡ a33fbc44-eade-4ef2-9a6d-87047b5cdb1f
#snapshotsresample = gresample.(snapshots, nz=280, method=:pchip)

# ╔═╡ 87e8b92e-dd35-4024-b795-467a96660e7b
#snapshotsresample[1]

# ╔═╡ 5ca0a8ee-5fbf-4d76-8bfd-91c292e08cfa
md"## Convert to M3D format"

# ╔═╡ fc232599-e4e9-42d0-b108-316897c363ce
begin
	function snaps2multi(folder, snaps...; 
						eos, label, 
						n_horizontal=10, n_vertical=280, outfolder="", method=:linear)
		i = 0
		for snap in snaps
			snapshot, snapshot_τ = try
				pick_snapshot(folder, snap)
			catch
				continue
			end

			aos = @axed eos
			ne = lookup(
				aos, :lnNe, 
				log.(snapshot[:d]), 
				log.(snapshot[:ee])
			) 
			snapshot.data[:ne] = exp.(ne)
		
			snapshotsresample = gresample(
				snapshot, 
				nz=n_vertical, 
				nx=n_horizontal, 
				ny=n_horizontal,
				method=method
			)
	
			if length(outfolder) > 1
				(!isdir(outfolder)) && mkdir(outfolder)
			end
	
			i += 1
			output_name = joinpath(outfolder, "m3dis_sun_$(label)_$(i)")
		
			MUST.multiBox(
				snapshotsresample, 
				output_name
			)
	
			#@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
		end
	end
	
	function snaps2multi(snaps::MUST.Box...; 
						eos, label, n_horizontal=10, n_vertical=280, outfolder="")
		i = 0
		for snap in snaps
			snapshot = snap
			snapshotsresample = gresample(
				snapshot, 
				nz=n_vertical, 
				nx=n_horizontal, 
				ny=n_horizontal
			)
	
			if length(outfolder) > 1
				(!isdir(outfolder)) && mkdir(outfolder)
			end
	
			i += 1
			output_name = joinpath(outfolder, "m3dis_sun_$(label)_$(i)")
			aos = @axed eos
	
			ne = lookup(
				aos, :lnNe, 
				log.(snapshotsresample[:d]), log.(snapshotsresample[:ee])
			) 
			
			snapshotsresample.data[:ne] = exp.(ne)
		
			MUST.multiBox(
				snapshotsresample, 
				output_name
			)
	
			#@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
		end
	end
end

# ╔═╡ 7a023be5-46ea-4c25-857b-3f765c044a91
labels_new = replace.(labels, " "=>"_")

# ╔═╡ b1df1501-4033-4d8a-bd67-8130c095152a
downsample=25

# ╔═╡ fdf3a692-daab-49d5-8bf6-36996617349e
output_names = ["m3dis_sun_$(label)" for (name, label) in zip(names, labels_new)]

# ╔═╡ f3e705cf-a0a5-4905-9a07-aa6535985e01
aos = [@axed(eos[i]) for i in eachindex(eos)]

# ╔═╡ f2be5e98-6545-4f37-9a20-026f4914c9f7
md"One can alternatively also convert multiple snapshots"

# ╔═╡ 2aab4558-dc94-475d-83d2-696d513c2aec
labels

# ╔═╡ 8c3db52c-75f5-4051-ae57-07dabebaf021
begin
	target_x = 10
	target_z = 280
end

# ╔═╡ 389e61c2-0dd4-458a-98a5-5b787fa0e957
label = "magg22_$(target_x)x$(target_x)x$(target_z)"

# ╔═╡ 055f1b98-53d7-4368-bc75-d0073985d47d
for i in eachindex(labels)	
	@info labels[i]
	
	snaps2multi(
		out_folder[i], -8, -7, -6, -5, -4, -3, -2, -1, :recent,
		eos=eos[i], 
		label=label,
		n_horizontal=target_x, 
		n_vertical=target_z, 
		outfolder=labels[i], 
		method=:linear
	)
end

# ╔═╡ 313bb855-4d09-4df7-ba9d-f261cff27794
#label_stagger = "stagger_10x10x230"

# ╔═╡ 69ec4993-ddfd-496f-87d9-eddf9ab97714
#=eos_stagger = reload(
	SqEoS, 
	joinpath(MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35"), "eos.hdf5")
)=#

# ╔═╡ 14275f0e-62e2-4a6f-98e6-dcbed3a1e158
#=snaps2multi(
	stagger,
	eos=eos_stagger, 
	label=label_stagger,
	n_horizontal=10, 
	n_vertical=230, outfolder="stagger_sun"
)=#

# ╔═╡ 5b28c9b8-a991-45d6-a2d1-15f24142d277
md"## Compute a time average"

# ╔═╡ 7837df1e-e220-4a60-944d-876b523300ad
number_of_timeslots = [6]

# ╔═╡ 09d3fbbb-c636-4f8a-8409-0eba2ea9a242
do_time = false

# ╔═╡ b564b2bf-61d9-44e4-a785-fa66db6b10e5
if do_time
	for i in eachindex(names)
		eos_i = reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5"))
		#opa = reload(SqOpacity, joinpath(eos_folder[i], "binned_opacity.hdf5"))
		
		snaps = MUST.converted_snapshots(out_folder[i])
		isnap = MUST.list_snapshots(snaps)
	
		time_isnaps = isnap[end-number_of_timeslots[i]:end-1]
		
		time_snaps = [pick_snapshot(snaps, j) |> first for j in time_isnaps]
		time_av = MUST.time_statistic(mean, time_snaps)
		MUST.save(time_av, folder=out_folder[i], name="box_tav")
	
	
		time_snaps = [pick_snapshot(snaps, j) |> last for j in time_isnaps]
		time_av = MUST.time_statistic(mean, time_snaps)
		MUST.save(time_av, folder=out_folder[i], name="box_tau_tav")
	end
end

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
# ╠═452a144e-b725-4de2-b3e1-2f614210d62e
# ╟─3e747391-ba4b-47bf-b363-abcb46a9309b
# ╟─ed6c250a-84d5-4ac6-bf54-9e8fcdbe55a3
# ╠═d8693137-82f7-4ccb-b886-4115e3032392
# ╠═4c8d357d-a41a-4274-b131-93ebb695b911
# ╠═c72cea6d-59bd-4721-8044-9fa28fb4039e
# ╠═4a1730fb-1fbe-445a-b850-7539967dcbe2
# ╠═9faf1a6b-a027-4fcc-b572-872624b3f7e8
# ╟─14dab588-7275-4a7e-bade-71375d3a16bc
# ╟─0d72ddb5-0096-4430-8611-e843f0aa3c61
# ╠═7c8ed29b-2c7a-4b4c-a60d-4de3cf082b42
# ╠═a0220545-f149-4862-a8ba-3a1911cb732a
# ╟─c3ab29bc-a61b-467c-b3ab-3fa919abd83d
# ╠═cd81869e-c954-4526-a52f-31aedb64cb5f
# ╟─fa90d92b-d1b3-491e-872b-23f57da6dace
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
# ╟─7d273e41-207c-46b9-9e7d-a3626b5717c6
# ╠═9d6e7385-c8f6-4df1-91ac-3f07f0d2567f
# ╟─6e078b79-7928-4854-a07f-c5ac185caa59
# ╟─c0347b03-fb0a-4660-a9f6-fc24272263ed
# ╟─fcc3aaac-a5e2-45fa-9419-4aaa916fbe45
# ╟─4c0bf79b-9b7b-4c35-a728-350a79e5eb83
# ╠═a10e407c-70f8-4809-a0e4-eb334447d623
# ╟─919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
# ╟─84bd1e80-ff21-4970-a5a3-fc7452da7e6f
# ╟─e2dde825-4e66-42e2-b541-74ec0600fc81
# ╟─3615d990-9c9a-4e69-8f96-0e3c2ac6896b
# ╠═a33fbc44-eade-4ef2-9a6d-87047b5cdb1f
# ╠═87e8b92e-dd35-4024-b795-467a96660e7b
# ╟─5ca0a8ee-5fbf-4d76-8bfd-91c292e08cfa
# ╟─fc232599-e4e9-42d0-b108-316897c363ce
# ╠═7a023be5-46ea-4c25-857b-3f765c044a91
# ╟─b1df1501-4033-4d8a-bd67-8130c095152a
# ╟─fdf3a692-daab-49d5-8bf6-36996617349e
# ╟─f3e705cf-a0a5-4905-9a07-aa6535985e01
# ╟─f2be5e98-6545-4f37-9a20-026f4914c9f7
# ╠═2aab4558-dc94-475d-83d2-696d513c2aec
# ╠═8c3db52c-75f5-4051-ae57-07dabebaf021
# ╠═389e61c2-0dd4-458a-98a5-5b787fa0e957
# ╠═055f1b98-53d7-4368-bc75-d0073985d47d
# ╠═313bb855-4d09-4df7-ba9d-f261cff27794
# ╠═69ec4993-ddfd-496f-87d9-eddf9ab97714
# ╠═14275f0e-62e2-4a6f-98e6-dcbed3a1e158
# ╟─5b28c9b8-a991-45d6-a2d1-15f24142d277
# ╠═7837df1e-e220-4a60-944d-876b523300ad
# ╠═09d3fbbb-c636-4f8a-8409-0eba2ea9a242
# ╠═b564b2bf-61d9-44e4-a785-fa66db6b10e5
