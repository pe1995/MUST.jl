### A Pluto.jl notebook ###
# v0.19.25

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

	pythonplot()
end;

# ╔═╡ 9456066d-36b1-4e6d-8af9-b8f134fb6e24
md"# Investigating different Models"

# ╔═╡ 931f7a1f-dccb-4727-9982-d04be9ffd688
begin
	mean = MUST.mean
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
end;

# ╔═╡ 2e04240c-981e-443f-826b-2b40eb4b1974
md"All the models from the Stagger grid come with an average model that we can load as well. We can load the EoS from the same folder."

# ╔═╡ 87b753e7-78e5-437a-9e5e-51693641bee7
begin
	name = "DIS_MARCS_E_t5777g44m00_v0.1"
	in_folder = MUST.@in_dispatch "input_data/grd/$(name)"
	out_folder = MUST.@in_dispatch "data/$(name)_test2"
	eos_folder = MUST.@in_dispatch "input_data/DIS_MARCS_E_v1.4.33"
end;

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
snapshot, snapshot_τ = pick_snapshot(converted_snapshots(out_folder), :recent)

# ╔═╡ 3e747391-ba4b-47bf-b363-abcb46a9309b
begin
	eos = reload(SqEoS, joinpath(eos_folder, "eos.hdf5"))
	opa = reload(SqOpacity, joinpath(eos_folder, "binned_opacities.hdf5"))
end

# ╔═╡ 1c406c40-8619-4016-b701-adb5b12a22d4
begin
	title_s = split(name[last(findfirst("DIS_MARCS_E_", name))+1:end], "_") |> first
	title = String(title_s)*", $(length(opa.λ)) Bins"
end

# ╔═╡ d8693137-82f7-4ccb-b886-4115e3032392
# ╠═╡ show_logs = false
initial_model = @optical Average3D(eos, joinpath(in_folder, "inim.dat")) eos opa

# ╔═╡ 4c8d357d-a41a-4274-b131-93ebb695b911
stagger_model = @optical Average3D(eos, 
					MUST.@in_dispatch("input_data/solar_stagger_ext")) eos opa

# ╔═╡ 4a1730fb-1fbe-445a-b850-7539967dcbe2
begin
	folder_stagger = "/ptmp/peitner/model_grid/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
							folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t", 
							folder=folder_stagger)
end

# ╔═╡ 6184a62b-4562-4c1c-bf9c-7f2cc57afe52
md"We can compute the optical depth, if it is not present already"

# ╔═╡ 53ced3c5-71a1-424b-ac93-9c5c47ca9ef8
if !(:τ_ross in keys(snapshot.data))
	@info "Looking up Rosseland opacity and computing optical depth"
	
	# add the rosseland opacity
	MUST.add!(snapshot, lookup(eos, opa, :κ_ross, snapshot[:d], snapshot[:ee]), :kr)

	# compute optical depth and add it
	MUST.add!(snapshot, MUST.optical_depth(b_s, opacity=:kr, density=:d), :τ_ross)
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
	plot(-initial_model.z, exp.(initial_model.lnT), 
			lw=3., color=:black, label="initial condition", ls=:solid)

	plot!(profile(mean, snapshot, :z, :T)..., 
			lw=3., color=:red, label="DISPATCH")
	

	xlabel!("z [cm]")
	ylabel!("T [K]")

	plot!(title="\n"*title)
	
	basic_plot!()
end

# ╔═╡ d3ef4193-11e7-478d-aa11-ba3b8adbce55
begin
	#plot(log10.(initial_model.τ), exp.(initial_model.lnT), 
	#		lw=3., color=:black, label="initial condition", ls=:solid)

	plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
			lw=2., color=:black, label="Stagger", ls=:dot)

	plot!(profile(mean, snapshot_τ, :log10τ_ross, :T)..., 
			lw=3., color=:red, label="DISPATCH", ls=:solid)


	xlabel!("τ-ross [log]")
	ylabel!("T [K]")

	plot!(title="\n"*title)
	
	basic_plot!()
end

# ╔═╡ 05a4e3ef-58b5-4f96-94b8-51991c971451
begin
	plot(-initial_model.z, exp.(initial_model.lnρ), 
			lw=3., color=:black, label="initial condition", ls=:solid)

	plot!(profile(mean, snapshot, :z, :d)..., 
			lw=3., color=:red, label="DISPATCH")
	

	yaxis!(:log)
	xlabel!("z [cm]")
	ylabel!("density [g x cm-3]")

	plot!(title="\n"*title)
	
	basic_plot!()
end

# ╔═╡ 0a726cbb-1309-4071-9df1-93cd4355fb71
begin
	#plot(log10.(initial_model.τ), initial_model.lnρ, 
	#		lw=3., color=:black, label="initial condition", ls=:solid)
	
	plot(profile(mean, stagger_τ, :log10τ_ross, :log10d)..., 
			lw=2., color=:black, label="Stagger", ls=:dot)
	
	plot!(profile(mean, snapshot_τ, :log10τ_ross, :log10d)..., 
			lw=3., color=:red, label="DISPATCH")
	

	xlabel!("τ-ross [log]")
	ylabel!("density [g x cm-3]")

	plot!(title="\n"*title)
	
	basic_plot!()
end

# ╔═╡ e10a8581-c0ab-4209-9e14-4d456dcf9a86
begin
	plot(exp.(initial_model.lnρ), exp.(initial_model.lnT), 
			lw=3., color=:black, label="initial condition", ls=:solid)
	

	_, d = profile(mean, snapshot, :z, :d)
	_, T = profile(mean, snapshot, :z, :T)
	
	plot!(d, T, 
			lw=3., color=:red, label="DISPATCH")
	
	xaxis!(:log)	
	
	xlabel!("density [g x cm-3]")
	ylabel!("T [K]")

	plot!(title="\n"*title)
	
	basic_plot!()
	plot!(legend=:topleft)
end

# ╔═╡ 919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
md"## 3D result"

# ╔═╡ 84bd1e80-ff21-4970-a5a3-fc7452da7e6f
uz_τ_surf = if !isnothing(snapshot_τ)
	isurf = MUST.closest(log10.(MUST.axis(snapshot_τ, :τ_ross, 3)), 0)
	snapshot_τ[:uz][:, :, isurf]
else
	@info "Interpolating Uz..."
	uz_τ_surf = MUST.interpolate_to(snapshot, :uz, τ_ross=1)
	MUST.axis(uz_τ_surf, :uz, 3)
end

# ╔═╡ e2dde825-4e66-42e2-b541-74ec0600fc81
begin
	heatmap(MUST.axis(snapshot, :x)./1e8, 
			MUST.axis(snapshot, :y)./1e8,
			uz_τ_surf ./1e5,
				colormap=:hot,
				colorbar_title= "\nvertical velocity [km/s]")
	
	xlabel!("x [Mm]")
	ylabel!("y [Mm]")

	plot!(title="\n"*title)

	basic_plot!(copyticks=false, bm=2, lm=2, rm=6, size=(600,500))
end

# ╔═╡ Cell order:
# ╟─9456066d-36b1-4e6d-8af9-b8f134fb6e24
# ╠═1eac6457-238f-49ed-82dc-6f1c521930fc
# ╠═931f7a1f-dccb-4727-9982-d04be9ffd688
# ╟─2e04240c-981e-443f-826b-2b40eb4b1974
# ╠═87b753e7-78e5-437a-9e5e-51693641bee7
# ╠═957af12a-e77e-48a3-a2de-80b86512e5a8
# ╠═1e0de50e-bb67-4d34-a038-e7437955ec73
# ╠═01ab2753-515c-496f-a6fc-1c2a0e42ae25
# ╠═452a144e-b725-4de2-b3e1-2f614210d62e
# ╠═3e747391-ba4b-47bf-b363-abcb46a9309b
# ╟─1c406c40-8619-4016-b701-adb5b12a22d4
# ╟─d8693137-82f7-4ccb-b886-4115e3032392
# ╟─4c8d357d-a41a-4274-b131-93ebb695b911
# ╠═4a1730fb-1fbe-445a-b850-7539967dcbe2
# ╟─6184a62b-4562-4c1c-bf9c-7f2cc57afe52
# ╠═53ced3c5-71a1-424b-ac93-9c5c47ca9ef8
# ╟─c3ab29bc-a61b-467c-b3ab-3fa919abd83d
# ╟─fa90d92b-d1b3-491e-872b-23f57da6dace
# ╟─9fd5896c-8d4f-499d-9f97-c589d8d256c2
# ╟─85d4b142-1529-495a-bc6d-64f5f0efa3b9
# ╟─3d1e1cd9-5601-468d-a3f1-a3dd51177a37
# ╟─8a8a40e0-5b03-4c17-93a3-4eccf12e3717
# ╟─d3ef4193-11e7-478d-aa11-ba3b8adbce55
# ╟─05a4e3ef-58b5-4f96-94b8-51991c971451
# ╟─0a726cbb-1309-4071-9df1-93cd4355fb71
# ╠═e10a8581-c0ab-4209-9e14-4d456dcf9a86
# ╟─919c0ab0-23ad-4acd-ab2a-1c90cf0aa8e1
# ╠═84bd1e80-ff21-4970-a5a3-fc7452da7e6f
# ╟─e2dde825-4e66-42e2-b541-74ec0600fc81
