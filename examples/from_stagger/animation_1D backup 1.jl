### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ d8e58bc8-efe1-4fe0-90c2-68e04c3cdf5b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using Glob
	using Plots
	using TSO
	using LaTeXStrings
	using Printf
end

# ╔═╡ 6331d0c9-04e5-4f73-96d3-63e973d8ea36
include_helper(name) = include(joinpath(dirname(pathof(MUST)), name))

# ╔═╡ e138b8b9-ecd5-4023-b594-eb57b350e354
md"# Paper I: Validation & the Sun"

# ╔═╡ 7e5757f5-63e1-4cad-90b4-fc21b9a8ce79
md"## Code Setup"

# ╔═╡ 338a10a4-5283-4b5b-9989-080d28ccb60c
md"### Plotting defaults"

# ╔═╡ 1781d460-5fce-4a88-81c5-d9ba09971c7c
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

# ╔═╡ 1c220eed-89ad-4298-bdca-5570dbfcad66
MUST.@get_help_py gifs

# ╔═╡ 761025f4-b29a-4685-b2ff-c750331450fb
md"### Dispatch module"

# ╔═╡ 41e28c24-b0d3-4f04-9403-4a68ddaf2a07
mean = MUST.mean

# ╔═╡ 6e7d0b12-93db-4767-8ca5-04c95060ade3
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ db5ddb63-de7e-4a20-871e-324fd963c185
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 7630938a-0aac-4841-a5fc-cfb2eb836b2d
mesh(m::MUST.Box) = MUST.meshgrid(MUST.axis(m, :x) ./1e8, 
									MUST.axis(m, :y) ./1e8, 
									MUST.axis(m, :z) ./1e8)

# ╔═╡ 90686e9d-99fc-453c-afdb-ddff63a573b5
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

# ╔═╡ add9584d-b0f4-42a0-9414-5a9a9eb5b206
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if xs == :τ_ross
		logx.(MUST.axis(model, xs, 3)), logy.(MUST.plane_statistic(f, model, ys)) 
	else
		logx.(MUST.axis(model, xs)), logy.(MUST.plane_statistic(f, model, ys))
	end
end

# ╔═╡ ae145d37-f463-49aa-914c-a3cd8fd1a650
names_res = "DIS_MARCS_E_t5777g44m00_v0.1"

# ╔═╡ d97ab00a-ba3f-4457-a97b-2ae22a90f21f
out_folder_res = MUST.@in_dispatch "data/sun_magg"

# ╔═╡ 9a6c9eb2-f0c0-41aa-851d-a2038896bc97
in_folder_res = MUST.@in_dispatch "input_data/"

# ╔═╡ 35de1998-70fc-4685-ad37-91c19a79c714
eos_folder_res = MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35")

# ╔═╡ 3d439950-7377-4d01-94bf-ee37dff19579
snapshot, snapshot_τ = pick_snapshot(out_folder_res, :recent)

# ╔═╡ d9f80869-a8ae-4ace-b92d-61a54419e717
begin
	eos_res = reload(SqEoS, joinpath(eos_folder_res, "eos.hdf5"))
	opa_res = reload(SqOpacity, joinpath(eos_folder_res, "binned_opacities.hdf5"))
end

# ╔═╡ 870e623b-9c32-43d3-b3af-f8bb75aa8caf
@info "Opacity table size: $(size(opa_res.κ))"

# ╔═╡ f38476e4-d2eb-46fb-8a07-c7e4e2bd4e46
# ╠═╡ show_logs = false
initial_model = @optical(
	Average3D(eos_res, joinpath(in_folder_res, "sun_stagger.dat")),
				eos_res, opa_res)

# ╔═╡ df40139a-6c67-412a-8fa8-8e45ff4dba23
begin
	folder_stagger = "/u/peitner/DISPATCH/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
							folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t", 
							folder=folder_stagger)
end

# ╔═╡ 5fd204d9-ccd3-4683-aa0c-03dbf958e435
resolution(snap) = @sprintf "%.1f km" first(diff(MUST.axis(snap, :z) ./1e5))

# ╔═╡ 9e5d5453-9a6e-457d-a50f-c5558b10562d
@info "Resolution $(out_folder_res): $(resolution(snapshot))"

# ╔═╡ 3f3da68e-37ae-4044-a215-d6ab44789d35
labels_res = "M3DIS - $(resolution(snapshot))"

# ╔═╡ b3aa6e35-24dc-4370-9b42-7992bb742277
md"## Animation"

# ╔═╡ f37aca05-14c5-4803-ac6f-526bcc3453c0
begin
	snaps = MUST.converted_snapshots(out_folder_res)
	isnaps = MUST.list_snapshots(snaps)
	f = []
	
	for i in isnaps			
		s, st = pick_snapshot(snaps, i)

		isnothing(st) && continue 
		
		plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

		
		plot!(profile(mean, st, :log10τ_ross, :T)..., 
				lw=2., color="red", label="t: $(s.parameter.time)")
		
	
		xlabel!("log τ")
		ylabel!("T [K]")

		plot!(ylim=(4000,22000))
		basic_plot!()
		
		savefig("profile_$(i).png")
		append!(f, ["profile_$(i).png"])
	end
	
	gifs.gifs_from_png(f, "average_taut_rtres.gif", duration=2.0)
end

# ╔═╡ 8e362634-766e-4688-9995-054ecd5764f7
begin
	f2 = []
	
	for i in isnaps			
		s, st = pick_snapshot(snaps, i)

		isnothing(st) && continue 
		
		plot(profile(rms, stagger_τ, :log10τ_ross, :uz)..., 
				lw=1.5, color=:black, label="Stagger", ls=:dash)

		
		plot!(profile(rms, st, :log10τ_ross, :uz)..., 
				lw=2., color="red", label="t: $(s.parameter.time)")
		
	
		xlabel!("log τ")
		ylabel!("rms(Uz) [cm × s-1]")

		plot!(ylim=(0,3.5e5))
		basic_plot!()
		
		savefig("profile_$(i).png")
		append!(f2, ["profile_$(i).png"])
	end
	
	gifs.gifs_from_png(f, "average_tauuz_rtres.gif", duration=2.0)
end

# ╔═╡ c9369187-a6c1-4b86-83f3-371f7d46e0f8
begin
	f3 = []
	
	for i in isnaps			
		s, st = pick_snapshot(snaps, i)

		isnothing(st) && continue 
		
		plot(
			log10.(MUST.plane_statistic(mean, stagger_τ, :d)),
			MUST.plane_statistic(mean, stagger_τ, :T),
			lw=1.5, color=:black, label="Stagger", ls=:dash
		)
		
		plot!(
			log10.(MUST.plane_statistic(mean, st, :d)),
			MUST.plane_statistic(mean, st, :T),
			lw=2., color="red", label="t: $(s.parameter.time)"
		)
		
	
		ylabel!("T")
		xlabel!("log ρ")

		basic_plot!()
		
		savefig("profile_$(i).png")
		append!(f3, ["profile_$(i).png"])
	end
	
	gifs.gifs_from_png(f, "average_trho_rtres.gif", duration=2.0)
end

# ╔═╡ Cell order:
# ╠═d8e58bc8-efe1-4fe0-90c2-68e04c3cdf5b
# ╠═6331d0c9-04e5-4f73-96d3-63e973d8ea36
# ╟─e138b8b9-ecd5-4023-b594-eb57b350e354
# ╟─7e5757f5-63e1-4cad-90b4-fc21b9a8ce79
# ╟─338a10a4-5283-4b5b-9989-080d28ccb60c
# ╠═e7b79c35-7df1-49c0-96b4-dd25262b52f8
# ╟─1781d460-5fce-4a88-81c5-d9ba09971c7c
# ╠═1c220eed-89ad-4298-bdca-5570dbfcad66
# ╟─761025f4-b29a-4685-b2ff-c750331450fb
# ╟─41e28c24-b0d3-4f04-9403-4a68ddaf2a07
# ╟─6e7d0b12-93db-4767-8ca5-04c95060ade3
# ╠═db5ddb63-de7e-4a20-871e-324fd963c185
# ╟─7630938a-0aac-4841-a5fc-cfb2eb836b2d
# ╟─90686e9d-99fc-453c-afdb-ddff63a573b5
# ╟─add9584d-b0f4-42a0-9414-5a9a9eb5b206
# ╠═ae145d37-f463-49aa-914c-a3cd8fd1a650
# ╠═d97ab00a-ba3f-4457-a97b-2ae22a90f21f
# ╠═9a6c9eb2-f0c0-41aa-851d-a2038896bc97
# ╠═35de1998-70fc-4685-ad37-91c19a79c714
# ╠═3d439950-7377-4d01-94bf-ee37dff19579
# ╠═d9f80869-a8ae-4ace-b92d-61a54419e717
# ╟─870e623b-9c32-43d3-b3af-f8bb75aa8caf
# ╠═f38476e4-d2eb-46fb-8a07-c7e4e2bd4e46
# ╠═df40139a-6c67-412a-8fa8-8e45ff4dba23
# ╠═5fd204d9-ccd3-4683-aa0c-03dbf958e435
# ╟─9e5d5453-9a6e-457d-a50f-c5558b10562d
# ╠═3f3da68e-37ae-4044-a215-d6ab44789d35
# ╟─b3aa6e35-24dc-4370-9b42-7992bb742277
# ╠═f37aca05-14c5-4803-ac6f-526bcc3453c0
# ╠═8e362634-766e-4688-9995-054ecd5764f7
# ╠═c9369187-a6c1-4b86-83f3-371f7d46e0f8
