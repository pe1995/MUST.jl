### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 1fd2045b-9560-4351-beb9-d857b3b4de20
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using PyPlot
	using TSO
	using LaTeXStrings
	using Printf
	using KernelDensity
end

# ╔═╡ 8af8db87-db68-4262-b9a0-4169b18b27f9
include_helper(name) = include(joinpath(dirname(pathof(MUST)), name))

# ╔═╡ 19dacb2e-1d6b-4e7f-9c29-54c44051b5cb
md"# Paper I: Validation & the Sun"

# ╔═╡ 8f754b21-bfdb-43b6-874f-bb14759fae73
md"## Code Setup"

# ╔═╡ be161252-813b-42ad-94ba-aa51b236b5dd
md"### Plotting defaults"

# ╔═╡ 08f225ec-7e5c-471c-b482-3c4c813b9af9
begin
	plt.style.use("dark_background")
	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["font.size"] = 10
end

# ╔═╡ 07dc6754-ea9f-4e5b-8e96-5e30a9e3a9f2
include_helper("visual.jl")

# ╔═╡ e3a96c0b-3ca1-4d82-bb6a-ec74520ee142
MUST.@get_help_py gifs

# ╔═╡ c6e8f4d6-e7ca-468f-b846-30599ec400a3
md"### Dispatch module"

# ╔═╡ 53e3c235-e619-4983-88f8-b1871ce6e3f4
mean = MUST.mean

# ╔═╡ 53603ac2-bc15-4357-b40f-b1f148aaebb8
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ f33c9020-fb45-48ac-9575-174e71c808bd
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 49609c42-692f-4aae-9301-400ac08a42e3
mesh(m::MUST.Box) = MUST.meshgrid(MUST.axis(m, :x) ./1e8, 
									MUST.axis(m, :y) ./1e8, 
									MUST.axis(m, :z) ./1e8)

# ╔═╡ 8388262a-9cbb-4a3c-bc50-788dd27b3d18
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

# ╔═╡ 2bc9b4c1-953c-4f9b-89c7-b8cf987da505
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if xs == :τ_ross
		logx.(MUST.axis(model, xs, 3)), logy.(MUST.plane_statistic(f, model, ys)) 
	else
		logx.(MUST.axis(model, xs)), logy.(MUST.plane_statistic(f, model, ys))
	end
end

# ╔═╡ d2e8e1c7-6012-496c-8020-baa77a424b7d
names_res = "DIS_MARCS_E_t5777g44m00_v0.1"

# ╔═╡ 9796598b-65d4-48b2-a54f-c3abe861a1e8
out_folder_res = MUST.@in_dispatch "data/t5777g44m00_v0.1_test12"

# ╔═╡ 1944aa8a-9605-4783-8d22-5fa0208d0196
in_folder_res = MUST.@in_dispatch "input_data/grd/$(names_res)"

# ╔═╡ 67656849-b6d9-4579-a504-b573f4df297b
eos_folder_res = MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35")

# ╔═╡ c0c80317-8ca1-4b67-8786-5a9a9e43cb42
snapshot, snapshot_τ = pick_snapshot(out_folder_res, :recent)

# ╔═╡ 061a9f26-439d-4416-82a9-1d6487d97e66
begin
	eos_res = reload(SqEoS, joinpath(eos_folder_res, "eos.hdf5"))
	opa_res = reload(SqOpacity, joinpath(eos_folder_res, "binned_opacities.hdf5"))
end

# ╔═╡ 11ff8209-70a1-41f6-bae2-5b379b5d9c2b
@info "Opacity table size: $(size(opa_res.κ))"

# ╔═╡ d8028fe5-13d3-4a91-b037-1e6a7a6f2862
initial_model = @optical(Average3D(eos_res, joinpath(in_folder_res, "inim.dat")),
							eos_res, opa_res)

# ╔═╡ 70bec586-ff55-480e-bb2e-6efcb8ed3e7c
resolution(snap) = @sprintf "%.1f km" first(diff(MUST.axis(snap, :z) ./1e5))

# ╔═╡ b010fe01-2c36-40a6-a4f5-6cc0f0212583
@info "Resolution $(out_folder_res): $(resolution(snapshot))"

# ╔═╡ acf0be2d-138e-41a9-9c16-73bc29a3e391
labels_res = "M3DIS - $(resolution(snapshot))"

# ╔═╡ 51c29d9e-94ea-41c4-9eee-595f6dcf6b9b
md"## Animation"

# ╔═╡ f9d5c234-15cd-4eaa-8a17-2968b08f8f36
begin
	snaps = MUST.converted_snapshots(out_folder_res)
	isnaps = MUST.list_snapshots(snaps)
	f = []
	
	for i in isnaps	
		s, _ = pick_snapshot(snaps, i)

		close()
		
		ft, axt = cube_with_velocities(s, 
			vmin_3d=3500, 
			vmax_3d=15500, 
			s_3d=12, 
			arrow_length_ratio=0.2,
			skipv=3,
			xoff=5,
			yoff=5,
			zoff=5,
			len_vec=0.65, 
			cmap="RdYlBu",
			show_time=true
		)
		
		ft.savefig("cube_$(i).png", bbox_inches="tight")
		append!(f, ["cube_$(i).png"])
	end
	
	gifs.gifs_from_png(f, "cube_velocity_animation_lres.gif", duration=0.2)
end

# ╔═╡ Cell order:
# ╟─19dacb2e-1d6b-4e7f-9c29-54c44051b5cb
# ╟─8f754b21-bfdb-43b6-874f-bb14759fae73
# ╠═1fd2045b-9560-4351-beb9-d857b3b4de20
# ╟─be161252-813b-42ad-94ba-aa51b236b5dd
# ╠═08f225ec-7e5c-471c-b482-3c4c813b9af9
# ╠═8af8db87-db68-4262-b9a0-4169b18b27f9
# ╠═07dc6754-ea9f-4e5b-8e96-5e30a9e3a9f2
# ╠═e3a96c0b-3ca1-4d82-bb6a-ec74520ee142
# ╟─c6e8f4d6-e7ca-468f-b846-30599ec400a3
# ╠═53e3c235-e619-4983-88f8-b1871ce6e3f4
# ╠═53603ac2-bc15-4357-b40f-b1f148aaebb8
# ╠═f33c9020-fb45-48ac-9575-174e71c808bd
# ╠═49609c42-692f-4aae-9301-400ac08a42e3
# ╠═8388262a-9cbb-4a3c-bc50-788dd27b3d18
# ╠═2bc9b4c1-953c-4f9b-89c7-b8cf987da505
# ╟─d2e8e1c7-6012-496c-8020-baa77a424b7d
# ╠═9796598b-65d4-48b2-a54f-c3abe861a1e8
# ╠═1944aa8a-9605-4783-8d22-5fa0208d0196
# ╠═67656849-b6d9-4579-a504-b573f4df297b
# ╠═c0c80317-8ca1-4b67-8786-5a9a9e43cb42
# ╠═061a9f26-439d-4416-82a9-1d6487d97e66
# ╠═11ff8209-70a1-41f6-bae2-5b379b5d9c2b
# ╠═d8028fe5-13d3-4a91-b037-1e6a7a6f2862
# ╠═70bec586-ff55-480e-bb2e-6efcb8ed3e7c
# ╠═b010fe01-2c36-40a6-a4f5-6cc0f0212583
# ╠═acf0be2d-138e-41a9-9c16-73bc29a3e391
# ╟─51c29d9e-94ea-41c4-9eee-595f6dcf6b9b
# ╠═f9d5c234-15cd-4eaa-8a17-2968b08f8f36
