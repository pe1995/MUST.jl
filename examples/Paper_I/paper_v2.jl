### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 7be90c6e-260a-11ee-12d4-93fd873d1015
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using PythonPlot
	using TSO
	using LaTeXStrings
	using Printf
	using KernelDensity
	using DelimitedFiles
	using PlutoUI

	plt = matplotlib.pyplot
end;

# ╔═╡ 485a693d-4872-47d4-975d-e91a7bbc0be7
md"# Paper I: Validation & The Sun (v2)"

# ╔═╡ a61a48b3-cfe9-4d0b-a0f7-7008add4e1d9
TableOfContents()

# ╔═╡ feefb5e4-c4aa-428f-87ea-82b32072fb26
include_helper(name) = joinpath(dirname(pathof(MUST)), name)

# ╔═╡ b20211ed-1d69-49f0-9af3-002affd58ff2
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	path = include_helper(path)
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ c4846fea-8b7a-4338-bfa9-df9b9c0aff6f
begin
	#PythonPlot.matplotlib.rcParams["font.size"] = 12
	np = MUST.pyimport("numpy")
	ndimage = MUST.pyimport("scipy.ndimage")
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	hl = 4
end;

# ╔═╡ 5e79333c-94c1-47cd-86d5-a6ee22b4a62e
visual = ingredients("visual.jl")

# ╔═╡ 55bbe1f3-3ced-438b-9b57-777eb2fbc41d
md"## Dispatch Module"

# ╔═╡ 0c1ca595-6c09-48f8-ac15-9a499b9072c4
mean = MUST.mean

# ╔═╡ 22bacb57-9f7f-446b-93c6-4cf0fca5f731
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ 0b02d772-a51b-4aa5-9d73-23c85bac5a6a
rms5(x) = rms(x) ./ 1e5

# ╔═╡ 911bed42-e814-4a42-9be4-148334315fe2
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ eb51ab2b-d6e5-45a8-9131-5383545d02f0
mesh(m::MUST.Box) = MUST.meshgrid(
	MUST.axis(m, :x) ./1e8, 		
	MUST.axis(m, :y) ./1e8, 								
	MUST.axis(m, :z) ./1e8
)

# ╔═╡ 1e9b3054-c949-4a94-9db5-cba24602a470
"""
	time_average_profiles(model_folder)

Time average the ```plane_statistics``` of all snapshots located in that folder.
"""
function time_average_profiles(f, model_folder, args...; which=first, kwargs...)
	x, y = [], []
	for i in MUST.list_snapshots(MUST.converted_snapshots(model_folder))
		snap = which(MUST.pick_snapshot(model_folder, i))
		xi, yi = profile(f, snap, args...; kwargs...)

		append!(x, [xi])
		append!(y, [yi])
	end

	# interpolate to common x axis
	x_common = range(
		maximum(minimum.(x)), 
		minimum(maximum.(x)), 
		length=length(x |> first)
	)
	y_common = []
	for i in eachindex(x)
		m = sortperm(x[i])
		append!(y_common, [MUST.linear_interpolation(x[i][m], y[i][m]).(x_common)])
	end
	
	# compute the mean and the std
	y_mean = zeros(length(x_common))
	y_std  = zeros(length(x_common))
	yi     = zeros(length(x))

	for i in eachindex(x_common)
		for j in eachindex(yi)
			yi[j] = y_common[j][i]
		end
		y_mean[i] = MUST.mean(yi)
		y_std[i]  = MUST.std(yi)
	end

	x_common, y_mean, y_std
end

# ╔═╡ 152c58ee-6568-4fd3-ae3d-2f8e397bbd04
md"## General Models"

# ╔═╡ ec5e9d64-e461-44ec-9753-c306cc8d29dd
md"Solar stagger model (converted to optical depth with current EoS."

# ╔═╡ 33690f11-20d9-48ee-82b9-bc0df0f9e560
md"### Stagger"

# ╔═╡ 8e521645-e67e-463d-ad8b-9a9a8dd9056e
begin
	folder_stagger = "/u/peitner/DISPATCH/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
		folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t",
		folder=folder_stagger)	
end

# ╔═╡ b8bd4aa7-8e9e-40de-8188-0ddb8b233b5e


# ╔═╡ f6edb03c-5177-4777-b3f2-fc69513a7249
md"### Muram"

# ╔═╡ f0327365-b8ff-44ac-9aa1-d42ecd5a0a3f
begin
	folder_muram = "/u/peitner/DISPATCH/examples/Muram"
	muram = MUST.Box("box_MURaM_cube_small.221000_HDSun", folder=folder_muram)
	muram_τ = MUST.Box("box_MURaM_cube_small.221000_HDSun_t", folder=folder_muram)
end

# ╔═╡ b67ca362-66e0-4e6b-ac1e-094515498379
keys(muram.data)

# ╔═╡ 9a504c51-4e79-4ef4-8a4f-570f549422c0
#muram_τ = MUST.height_scale(muram, :τ_ross)

# ╔═╡ 0a2286e2-29e4-454b-b3a7-6a33a6828760
#MUST.save(muram_τ, folder=folder_muram, name="box_MURaM_cube_small.221000_HDSun_t")

# ╔═╡ 2e44bee0-da4c-47cf-8c7e-e8c46e696bd2


# ╔═╡ 9465b4d5-08e2-453a-b758-6d4d6744c20b
md"### MARCS"

# ╔═╡ d87a8635-7527-4584-81f6-eef0da1101d6
marcs_model = readdlm(
	"/u/peitner/DISPATCH/examples/from_stagger/marcs_sun.txt", 
	skipstart=1
)

# ╔═╡ a22bef7c-6409-4eb5-a17e-66ea1a42ab73


# ╔═╡ 0016626d-bb90-47eb-83f8-ec1f9d905a03
md"We can also read other marcs models. The one above I modified to make it better readable as csv, but I now added a reader for the true MARCS models."

# ╔═╡ dfb2ddf1-10ec-426d-8907-eec2b694fe5a
marcs_models = MUST.MARCSModel.(
	joinpath.(
	"models/marcs/", [
	"p5500_g+4.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
	"p6000_g+4.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
	"p6500_g+4.5_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
	"p4500_g+4.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
	"p5777_g+4.4_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"
	])
)

# ╔═╡ 4eecac61-40f9-45d4-b789-d017aef6f149


# ╔═╡ 2270321b-6581-4d45-a00f-0e728022bc56
md"MARCS models with different alpha values"

# ╔═╡ 5b42ec94-a9fc-4ada-aa7e-5fb6bac656b1
begin
	marcs_alpha = readdlm.(
		[
			"marcs_alpha/soleil-alpha0.5.txt", 
			"marcs_alpha/soleil-alpha1.5-nopturb.txt", 
			"marcs_alpha/soleil-alpha1.5.txt"
		],
		skipstart=1
	)
	marcs_alpha = Dict(
		"alpha0.5"=>marcs_alpha[1],
		"alpha1.5n"=>marcs_alpha[2],
		"alpha1.5"=>marcs_alpha[3]
	)
end

# ╔═╡ 0f310354-ca88-468f-a760-8194322855f3


# ╔═╡ c9a52157-81c4-4b8c-8dcd-b728a0890b00
md"### Co5bold"

# ╔═╡ 515e5910-4bf1-4769-84b9-cff2a08751a0
co5boldpath = glob("*", "co5bold_tauR_av-on_tauR")

# ╔═╡ fa715a70-036e-420a-a0fa-06cefc3ef3c0
co5bolds_τ = readdlm.(co5boldpath, skipstart=1)

# ╔═╡ 87cc7333-4559-467a-be0f-0f2ec1dc0611
begin
	co5boldmean_τ = zeros(size(co5bolds_τ[1], 1), 4)
	co5boldstd_τ = zeros(size(co5bolds_τ[1], 1), 4)
	co5boldrms_τ = zeros(size(co5bolds_τ[1], 1), 4)
	
	vals = zeros(length(co5bolds_τ), 4)
	for i in axes(co5boldmean_τ, 1)
		for j in eachindex(co5bolds_τ)
			vals[j, 1] = co5bolds_τ[j][i, 1]
			vals[j, 2] = co5bolds_τ[j][i, 2]
			vals[j, 3] = log10.(co5bolds_τ[j][i, 3])
			vals[j, 4] = co5bolds_τ[j][i, 6]
		end

		mv = mean(vals, dims=1)
		co5boldmean_τ[i, 1] = mv[1]
		co5boldmean_τ[i, 2] = mv[2]
		co5boldmean_τ[i, 3] = mv[3]
		co5boldmean_τ[i, 4] = mv[4]./1e5
		

		mv = MUST.std(vals, dims=1)
		co5boldstd_τ[i, 1] = mv[1]
		co5boldstd_τ[i, 2] = mv[2]
		co5boldstd_τ[i, 3] = mv[3]
		co5boldstd_τ[i, 4] = mv[4]./1e5
		

		co5boldrms_τ[i, 1] = rms(vals[:, 1])
		co5boldrms_τ[i, 2] = rms(vals[:, 2])
		co5boldrms_τ[i, 3] = rms(vals[:, 3])
		co5boldrms_τ[i, 4] = rms5(vals[:, 4])
	end
end

# ╔═╡ 3babe5e7-82ee-4d5b-a871-060b6050a53f
md"## Dispatch models"

# ╔═╡ c3b8db0b-c9ee-4f4d-a1a4-75e3551d7675
names = [
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t55g45m00_v0.1",
	"DIS_MARCS_E_t60g45m00_v0.1",
	"DIS_MARCS_E_t65g45m00_v0.1",
	"DIS_MARCS_E_t45g40m00_v0.1"
]

# ╔═╡ 6c3227a7-993b-45bc-809e-be6a0907f384
out_folder = [
	"models/sun_magg_150x300",
	"models/sun_magg_90x90",
	"models/sun_magg_90x180",
	"models/sun_magg_195x390",
	"models/t55g45m00_magg_150x300",
	"models/t60g45m00_magg_150x300",
	"models/t65g45m00_magg_150x300",
	"models/t45g40m00_magg_150x300"
]

# ╔═╡ 239ce3e3-5522-4d03-96d7-607c69418781
in_folder = [
	MUST.@in_dispatch "input_data/grd/$(name)" for name in names
]

# ╔═╡ f78e1d7f-6756-47d9-bb6c-5b6c623dc899
eos_folder = [
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_t55g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_t60g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_t65g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_t45g40m00_v0.1")
]

# ╔═╡ 8b046999-9214-468b-8b6f-90df998ecebd
eos_folder_T = [
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_v1.6.3"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_t55g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_t60g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_t65g45m00_v0.1"),
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_t45g40m00_v0.1")
]

# ╔═╡ f6c1c537-9a8d-4746-abcd-17fd0a712353
colors = [
	"tomato",
	"cyan",
	"magenta",
	"lime",
	"black",
	"black",
	"black",
	"black"
]

# ╔═╡ cff24989-efe0-4c50-8163-7260869895d2
ls = [
	"-",
	"-",
	"-",
	"-",
	"-",
	"--",
	"-",
	":"
]

# ╔═╡ 4fb7c54f-7f7a-48b7-a6a6-f80cb9c31360
labels = [
	"M3DIS",
	"M3DIS - lres",
	"M3DIS - ires",
	"M3DIS - hres",
	"M3DIS - t55g45m00",
	"M3DIS - t60g45m00",
	"M3DIS - t65g45m00",
	"M3DIS - t45g40m00"
]

# ╔═╡ 1005536b-f9c3-4526-be24-17ad9b063524
rt_resolution = [
	300,
	90,
	180,
	390,
	300,
	300,
	300, 
	300
]

# ╔═╡ 68e9c5c9-34df-41f6-a4e6-7f5deb935d59
md"Specific snapshots (generic)"

# ╔═╡ dd5a0ce9-2daa-4d99-9dbb-443c113e3af5
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = pick_snapshot(out_folder[i], :recent)
		
		append!(snapshots, [snapshot])
		append!(snapshots_τ, [snapshot_τ])
	end
end

# ╔═╡ 089e4d19-22ad-4da4-8886-6980d70b5f31
md"Handy shortcuts to pick certain models when needed."

# ╔═╡ d9272f17-8f18-457e-80d3-aa7dd61831a9
models = Dict(
	"best"      => 1,
	"lres"      => 2,
	"ires"      => 3,
	"hres"      => 4,
	"t55g45m00" => 5,
	"t60g45m00" => 6,
	"t65g45m00" => 7,
	"t45g40m00" => 8
)

# ╔═╡ 418c34ff-1c08-4ae6-82bf-7ac2fced7244
begin
	eos = [reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5")) 
				for i in eachindex(eos_folder)]
	opa = [reload(SqOpacity, joinpath(eos_folder[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder)]
end

# ╔═╡ fd9d2f80-3b17-4c3a-b8e8-37de365c0b2a
begin
	eos_T = [reload(SqEoS, joinpath(eos_folder_T[i], "eos.hdf5")) 
				for i in eachindex(eos_folder_T)]
	opa_T = [reload(SqOpacity, joinpath(eos_folder_T[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder_T)]
end

# ╔═╡ 23728e99-134a-445f-802b-f01467a03e10
for i in eachindex(eos)
	@info "Opacity table size: $(size(opa[i].κ))"
end

# ╔═╡ 114947c5-43c4-4d82-9031-8f27903d0415
resolution(snap) = @sprintf "%.1f km" first(diff(MUST.axis(snap, :z) ./1e5))

# ╔═╡ 26ad98c3-8c6a-42f9-944b-035691fa9ad0
resolution_ross(snap, rt_res, eos_snap) = begin
	snap_re = gresample(snap; nz=rt_res)
	κ = exp.(lookup(eos_snap, :lnRoss, log.(snap_re[:d]), log.(snap_re[:ee])))
	add!(snap_re, κ, :kross)
	o = MUST.optical_depth(snap_re, opacity=:kross, density=:d)
	
	res = []
	for j in axes(snap_re.y, 2) 
		for i in axes(snap_re.x, 1) 
			column = log10.(o[i, j, :])
			#column = log10.(snap[:τ_ross][i, j, :])
			i0 = argmin(abs.(column))
			dt = (abs(column[i0-1] - column[i0]) + abs(column[i0] - column[i0+1]))/2
			append!(res, [dt])
		end
	end
			
	@sprintf "%.3f Δlgτ(lgτ=0)" mean(res)
end

# ╔═╡ c5e1dc32-de92-4deb-92a9-06818d898307
for (i,snap) in enumerate(snapshots)
	@info "Resolution $(out_folder[i]): $(resolution(snap))"
end

# ╔═╡ 574fe0c0-fdcc-493f-8aad-bbfdc5b55c6a
for (i, snap) in enumerate(snapshots)
	@info "Resolution (optical) $(out_folder[i]): $(resolution_ross(snap, rt_resolution[i], eos[i]))"
end

# ╔═╡ d1163018-162d-4489-aa43-c423bcfd64f7
for (i,snap) in enumerate(snapshots)
	@info "Resolution $(out_folder[i]): $(snap.parameter.time)"
end

# ╔═╡ 01f1b55d-ba4b-4d74-8eba-12ede107cfbf
md"## Conversions"

# ╔═╡ 3fb509e5-c611-45e2-b475-f43bd28ac53a
must2multi = MUST.ingredients("convert2multi.jl")

# ╔═╡ 84fc6b08-0825-41ab-b214-5fef7e2a2dc0
res = [] #reshape(collect(Iterators.product((150, 300), (299))), :)

# ╔═╡ 35d2d2c5-6b70-4bb9-acdc-216557859282
# ╠═╡ show_logs = false
begin
	modelsConversion = [
		models["best"],
		models["lres"],
		models["ires"],
		models["hres"],
		models["t55g45m00"], 
		models["t60g45m00"], 
		models["t65g45m00"]
	]
	for (x, z) in res
		for i in modelsConversion
			zres = rt_resolution[i]-1 # z
			must2multi.snaps2multi(
				out_folder[i], 
				[range(-8, -1)..., :recent]...,
				eos=eos[i], 
				label="magg22_$(x)x$(x)x$(zres)",
				n_horizontal=x,
				n_vertical=zres, 
				outfolder= replace(labels[i], " "=>"")
			)
		end
	end
end

# ╔═╡ ae765c38-1c13-4914-85a1-41c41bff4f3f
md"For computing the effective temperature with Multi3D, we furhtermore save the corresponding opacities in Multi-readable format in the same folder."

# ╔═╡ 0ef33323-1d6a-484c-9292-3485b4871b38
for i in modelsConversion
	if isdir(replace(labels[i], " "=>""))
		save(eos_T[i], opa_T[i], joinpath(replace(labels[i], " "=>""), "eos_opa"))
	end
end

# ╔═╡ 381be17e-b257-4e51-aa79-941db153f398
md"## Figures
### (A) Internal Comparison"

# ╔═╡ e1517fd6-523f-4083-bb74-ebf666813002
begin
	modelA = models["best"]
	plt.close()
	
	fA, axA = plt.subplots(3, 1, sharex=true, figsize=(5, 9))
	plt.subplots_adjust(wspace=0, hspace=0)
	#visual.basic_plot!.(axA)

	m3disA = pick_snapshot(out_folder[modelA], :recent) |> last
	m3disAz = pick_snapshot(out_folder[modelA], :recent) |> first
	m3disA_x, m3disA_y, m3disA_yerr = time_average_profiles(
		mean, 
		out_folder[modelA], 
		:log10τ_ross, 
		:T, 
		which=last
	)
	m3disA_xr, m3disA_yr, m3disA_yerrr = time_average_profiles(
		mean, 
		out_folder[modelA], 
		:log10τ_ross, 
		:log10d, 
		which=last
	)
	m3disA_xu, m3disA_yu, m3disA_yerru = time_average_profiles(
		rms5, 
		out_folder[modelA], 
		:log10τ_ross, 
		:uz, 
		which=last
	)

	
	# logτ vs. T
	axA[0].plot(
		m3disA_x, m3disA_y,
		color="cornflowerblue",
		ls="-",
		label=L"\rm M3DIS\ (D)",
		lw=4
	)
	axA[0].plot(
		profile(mean, stagger_τ, :log10τ_ross, :T)...,
		color="k",
		ls="-",
		label=L"\rm Stagger"
	)
	axA[0].plot(
		co5boldmean_τ[:, 1], co5boldmean_τ[:, 2],
		color="k",
		ls="--",
		label=L"\rm Co5bold"
	)
	axA[0].plot(
		marcs_models[5].structure["lgTauR"],marcs_models[5].structure["T"],
		color="k",
		ls=":",
		label=L"\rm MARCS"
	)
	

	# logτ vs. logρ
	axA[1].plot(
		m3disA_xr, m3disA_yr,
		color="cornflowerblue",
		ls="-",
		label=L"\rm M3DIS\ (D)",
		lw=4
	)
	axA[1].plot(
		profile(mean, stagger_τ, :log10τ_ross, :log10d)...,
		color="k",
		ls="-",
		label=L"\rm Stagger"
	)
	axA[1].plot(
		co5boldmean_τ[:, 1], co5boldmean_τ[:, 3],
		color="k",
		ls="--",
		label=L"\rm Co5bold"
	)
	axA[1].plot(
		marcs_models[5].structure["lgTauR"], log10.(marcs_models[5].structure["Density"]),
		color="k",
		ls=":",
		label=L"\rm MARCS"
	)
	

	# logτ vs. uz
	axA[2].plot(
		m3disA_xu, m3disA_yu,
		color="cornflowerblue",
		ls="-",
		label=labels[modelA],
		lw=4
	)
	axA[2].plot(
		profile(rms5, stagger_τ, :log10τ_ross, :uz)...,
		color="k",
		ls="-",
		label="Stagger"
	)
	#=axA[2].plot(
		co5boldmean_τ[:, 1], co5boldmean_τ[:, 4],
		color="k",
		ls="--",
		label="Co5bold"
	)=#
	

	axA[0].legend(framealpha=0, labelspacing=0.01, handlelength=hl)
	
	axA[0].set_ylabel(L"\rm T\ [K]")
	axA[1].set_ylabel(L"\rm \log \rho\ [g \ cm^{-3}]")
	axA[2].set_ylabel(L"\rm v_{z}\ [km \ s^{-1}]")
	
	axA[2].set_xlabel(L"\rm \log \tau_{ross}")
	#axA[2].set_xlim(-4, 2)
	#axA[0].set_ylim(4000, 10500)
	#axA[1].set_ylim(-8.75, -6.25)
	
	
	gcf()
end

# ╔═╡ e73c5a0e-ae3b-4bbd-b098-d8e6ff215c2e
begin
	ftest, axtest = plt.subplots(1, 1)
	m3disAtest = pick_snapshot(out_folder[models["t45g40m00"]], :recent) |> first

	axtest.plot(profile(mean, m3disAtest, :z, :T)...)

	axtest.plot(
		-marcs_models[4].structure["Depth"], marcs_models[4].structure["T"],
		color="k",
		ls=":",
		label=L"\rm MARCS"
	)

	gcf()
end

# ╔═╡ f6f45939-55f8-42c3-b0d2-79c1cefb645c
begin
	_, m3disA_yTr = profile(
		mean,
		m3disA,
		:log10τ_ross,
		:T
	)
	_, m3disA_xTr = profile(
		mean,
		m3disA,
		:log10τ_ross,
		:log10d
	)
	
	_, staggerA_yTr = profile(
		mean,
		stagger_τ,
		:log10τ_ross,
		:T
	)
	_, staggerA_xTr = profile(
		mean,
		stagger_τ,
		:log10τ_ross,
		:log10d
	)
end

# ╔═╡ d77216a7-038d-491f-b7c2-ba74d627e2a2
begin
	mergewithfolder(path) = joinpath(
		MUST.@in_dispatch("input_data/grd/DIS_MARCS_v1.6.3"), path
	)
	eosT = reload(SqEoS, mergewithfolder("eos.hdf5"))
	opaT = reload(SqOpacity, mergewithfolder("binned_opacities.hdf5"))

	
	MUST.add!(
		m3disA, 
		exp.(lookup(eosT, :lnPg, log.(m3disA[:d]),  log.(m3disA[:T]))),
		:pg
	)
	MUST.add!(
		m3disA, 
		exp.(lookup(eosT, :lnNe, log.(m3disA[:d]),  log.(m3disA[:T]))),
		:ne
	)
	MUST.add!(
		m3disAz, 
		exp.(lookup(eosT, :lnPg, log.(m3disAz[:d]),  log.(m3disAz[:T]))),
		:pg
	)
	MUST.add!(
		m3disAz, 
		exp.(lookup(eosT, :lnNe, log.(m3disAz[:d]),  log.(m3disAz[:T]))),
		:ne
	)

	MUST.add!(
		stagger_τ, 
		exp.(lookup(eosT, :lnPg, log.(stagger_τ[:d]),  log.(stagger_τ[:T]))),
		:pg
	)
	
	_, m3disA_yTp = profile(
		mean,
		m3disA,
		:log10τ_ross,
		:T
	)
	_, m3disA_xTp = profile(
		mean,
		m3disA,
		:log10τ_ross,
		:log10pg
	)

	_, staggerA_yTp = profile(
		mean,
		stagger_τ,
		:log10τ_ross,
		:T
	)
	_, staggerA_xTp = profile(
		mean,
		stagger_τ,
		:log10τ_ross,
		:log10pg
	)

	co5boldA_xTp = exp.(lookup(
		eosT, :lnPg, log.(exp10.(co5boldmean_τ[:, 3])),  log.(co5boldmean_τ[:, 2])
	))
end

# ╔═╡ 518657e8-1865-4ecf-bb57-73cb502b5b19
begin
	plt.close()
	
	fA2, axA2 = plt.subplots(1, 1, sharex=true, figsize=(5, 6))
	plt.subplots_adjust(wspace=0, hspace=0)
	#visual.basic_plot!.(axA)
	
	# logτ vs. T
	axA2.plot(
		m3disA_yr, m3disA_y,
		color="cornflowerblue",
		ls="-",
		label=L"\rm M3DIS\ (D)",
		lw=3
	)
	axA2.plot(
		staggerA_xTr, staggerA_yTr,
		color="k",
		ls="-",
		label=L"\rm Stagger"
	)
	axA2.plot(
		co5boldmean_τ[:, 3], co5boldmean_τ[:, 2],
		color="k",
		ls="--",
		label=L"\rm Co5bold"
	)
	axA2.plot(
		log10.(marcs_models[5].structure["Density"]), marcs_models[5].structure["T"],
		color="k",
		ls=":",
		label=L"\rm MARCS"
	)

	axA2.legend(framealpha=0, labelspacing=0.01, handlelength=hl)

	axA2.set_ylabel(L"\rm T\ [K]")
	axA2.set_xlabel(L"\rm \log \rho\ [g \times cm^{-3}]")
	
	axA2.set_ylim(3600, 14500)
	axA2.set_xlim(-8.75, -5.9)
	
	
	gcf()
end

# ╔═╡ 65a010db-d74f-4326-b909-8586c4d033a7
begin
	plt.close()
	
	fA3, axA3 = plt.subplots(1, 1, sharex=true, figsize=(5, 6))
	plt.subplots_adjust(wspace=0, hspace=0)
	#visual.basic_plot!.(axA)
	
	# logτ vs. T
	axA3.plot(
		m3disA_xTp, m3disA_yTp,
		color="cornflowerblue",
		ls="-",
		label=L"\rm M3DIS\ (D)",
		lw=4.5
	)
	axA3.plot(
		staggerA_xTp, staggerA_yTp,
		color="k",
		ls="-",
		label=L"\rm Stagger"
	)
	axA3.plot(
		log10.(co5boldA_xTp), co5boldmean_τ[:, 2],
		color="k",
		ls="--",
		label=L"\rm Co5bold"
	)
	#=axA3.plot(
		log10.(marcs_models[5].structure["Pg"]), marcs_models[5].structure["T"],
		color="k",
		ls="-.",
		label=L"\rm MARCS"
	)=#

	# add MARCS models with different alpha
	axA3.plot(
		log10.(marcs_alpha["alpha1.5n"][:, 7]), marcs_alpha["alpha1.5n"][:, 5],
		color="k",
		ls="-.", lw=2.5,
		label=L"\rm MARCS\ \alpha=1.5\ no\ p_{turb}", zorder=0
	)
	axA3.plot(
		log10.(marcs_alpha["alpha1.5"][:, 7]), marcs_alpha["alpha1.5"][:, 5],
		color="cyan",
		ls="-.", lw=2.5,
		label=L"\rm MARCS\ \alpha=1.5 + p_{turb}", zorder=0
	)
	axA3.plot(
		log10.(marcs_alpha["alpha0.5"][:, 7]), marcs_alpha["alpha0.5"][:, 5],
		color="magenta",
		ls="-.", lw=2.5,
		label=L"\rm MARCS\ \alpha=0.5 + p_{turb}", zorder=0
	)

	axA3.legend(framealpha=0, labelspacing=0.01, handlelength=hl)

	axA3.set_ylabel(L"\rm T\ [K]")
	axA3.set_xlabel(L"\rm \log P_g\ [g \times cm \times s^{-2}]")
	
	#axA3.set_ylim(3600, 14500)
	#xA3.set_xlim(3.2, 6.1)

	axA3.set_ylim(5000, 11100)
	axA3.set_xlim(4.8, 5.3)
	

	fA3.savefig("T-Pg_3D_MARCS.pdf", bbox_inches="tight")
	fA3.savefig("T-Pg_3D_MARCS.png", bbox_inches="tight", dpi=300)
	
	gcf()
end

# ╔═╡ 13b57468-6d36-4413-95c6-0743e2633924
let
	open("average_sun_log.txt", "w") do f
		z, t = profile(mean, m3disA, :log10τ_ross, :T)
		z, r = profile(mean, m3disA, :log10τ_ross, :log10d)
		z, p = profile(mean, m3disA, :log10τ_ross, :log10pg)
		z, ne = profile(mean, m3disA, :log10τ_ross, :log10ne)
		z, uz = profile(rms5, m3disA, :log10τ_ross, :uz)
		
		MUST.writedlm(f, [z t exp10.(r) exp10.(p) exp10.(ne) uz], ',')
	end
	open("average_sun.txt", "w") do f
		z, t = profile(mean, m3disA, :log10τ_ross, :T)
		z, r = profile(mean, m3disA, :log10τ_ross, :d)
		z, p = profile(mean, m3disA, :log10τ_ross, :pg)
		z, ne = profile(mean, m3disA, :log10τ_ross, :ne)
		z, uz = profile(rms5, m3disA, :log10τ_ross, :uz)

		MUST.writedlm(f, [z t r p ne uz], ',')
	end
	open("average_sun_log_z.txt", "w") do f
		z, t = profile(mean, m3disAz, :z, :T)
		z, r = profile(mean, m3disAz, :z, :log10d)
		z, p = profile(mean, m3disAz, :z, :log10pg)
		z, ne = profile(mean, m3disAz, :z, :log10ne)
		z, uz = profile(rms5, m3disAz, :z, :uz)
		
		MUST.writedlm(f, [z t exp10.(r) exp10.(p) exp10.(ne) uz], ',')
	end
	open("average_sun_z.txt", "w") do f
		z, t = profile(mean, m3disAz, :z, :T)
		z, r = profile(mean, m3disAz, :z, :d)
		z, p = profile(mean, m3disAz, :z, :pg)
		z, ne = profile(mean, m3disAz, :z, :ne)
		z, uz = profile(rms5, m3disAz, :z, :uz)

		MUST.writedlm(f, [z t r p ne uz], ',')
	end
end

# ╔═╡ 58bb3285-544f-4e8a-a27e-9b7e90804fe8
md"### (B) Comparison with other Models"

# ╔═╡ dd667e0e-e554-4dfe-9493-aa38a505210f
begin
	common_tauB = range(-4, 2.5, length=100) |> collect
	
	ip(x, y) = begin
		m = sortperm(x)
		TSO.linear_interpolation(
			x[m], 
			y[m], 
			extrapolation_bc=MUST.Line()
		).(common_tauB)
	end
	
	absolute_difference(snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(mean, snapA, x, y)
		xB, yB = MUST.profile(mean, snapB, x, y)
		common_tauB, ip(xB, yB) .- ip(xA, yA)
	end
	absolute_difference(f, snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(f, snapA, x, y)
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, ip(xB, yB) .- ip(xA, yA)
	end
	absolute_difference(f, xA, yA, snapB, x, y) = begin
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, ip(xB, yB) .- ip(xA, yA)
	end
	absolute_difference(f, xA, yA, xB, yB, x, y) = begin
		common_tauB, ip(xB, yB) .- ip(xA, yA)
	end
end

# ╔═╡ 6de3b5df-9766-41ca-8d51-f18023a419a3
begin
	relative_difference(snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(mean, snapA, x, y)
		xB, yB = MUST.profile(mean, snapB, x, y)
		common_tauB, (ip(xB, yB) .- ip(xA, yA)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(f, snapA, x, y)
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, (ip(xB, yB) .- ip(xA, yA)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, xA, yA, snapB, x, y) = begin
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, (ip(xB, yB) .- ip(xA, yA)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, xA, yA, xB, yB, x, y) = begin
		common_tauB, (ip(xB, yB) .- ip(xA, yA)) ./ ip(xA, yA) * 100.0
	end
end

# ╔═╡ 8743067b-00ff-43c6-a6fd-800063b1d9de
begin
	modelB = models["best"]
	
	# Figure
	plt.close()
	
	fB, axB = plt.subplots(3, 1, sharex=true, figsize=(6, 10))
	plt.subplots_adjust(wspace=0.25, hspace=0.0)
	#visual.basic_plot!.(axB)
	

	m3disB = pick_snapshot(out_folder[modelB], :recent) |> last
	m3disB_x, m3disB_y, m3disB_yerr = time_average_profiles(
		mean, 
		out_folder[modelB], 
		:log10τ_ross, 
		:T, 
		which=last
	)
	m3disB_xr, m3disB_yr, m3disB_yerrr = time_average_profiles(
		mean, 
		out_folder[modelB], 
		:log10τ_ross, 
		:log10d, 
		which=last
	)
	m3disB_xu, m3disB_yu, m3disB_yerru = time_average_profiles(
		rms5, 
		out_folder[modelB], 
		:log10τ_ross, 
		:uz, 
		which=last
	)

	

	lwD = 2.0
	# relative comparison
	begin
		## logτ vs. T
		axB[0].axhline(0.0, color="0.5", alpha=0.2)
		axB[0].plot(
			relative_difference(
				mean,
				m3disB_x, 
				m3disB_y,
				stagger_τ,
				:log10τ_ross, :T
			)...,
			color="k",
			label=L"\rm Stagger - M3DIS\ (D)",
			lw=lwD
		)
		axB[0].fill_between( 
			relative_difference(
				mean,
				m3disB_x, 
				m3disB_y.-m3disB_yerr,
				stagger_τ,
				:log10τ_ross, :T
			)...,
			relative_difference(
				mean,
				m3disB_x, 
				m3disB_y.+m3disB_yerr,
				stagger_τ,
				:log10τ_ross, :T
			)[2],
			color="0.5", alpha=0.3
		)

		
		#=axB[0, 1].plot(
			absolute_difference(
				m3disB,
				muram_τ,
				:log10τ_ross, :T
			)...,
			color="k",
			label="M3DIS - MURaM",
			ls="--",
			lw=lwD
		)=#

		
		axB[0].plot(
			common_tauB, 
			(ip(marcs_models[5].structure["lgTauR"], marcs_models[5].structure["T"]) .- ip(m3disB_x, m3disB_y)) ./ ip(m3disB_x, m3disB_y) *100.0,
			color="k",
			label=L"\rm MARCS - M3DIS\ (D)",
			ls=":",
			lw=lwD
		)
		axB[0].fill_between( 
			common_tauB, 
			(ip(marcs_models[5].structure["lgTauR"], marcs_models[5].structure["T"]) .- ip(m3disB_x, m3disB_y.-m3disB_yerr)) ./ ip(m3disB_x, m3disB_y.-m3disB_yerr) *100.0,
			(ip(marcs_models[5].structure["lgTauR"], marcs_models[5].structure["T"]) .- ip(m3disB_x, m3disB_y.+m3disB_yerr)) ./ ip(m3disB_x, m3disB_y.+m3disB_yerr) *100.0,
			color="0.5", alpha=0.3
		)

		
		axB[0].plot(
			common_tauB, 
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 2]) .- ip(m3disB_x, m3disB_y)) ./ ip(m3disB_x, m3disB_y) *100.0,
			color="k",
			label=L"\rm Co5bold - M3DIS\ (D)",
			ls="--",
			lw=lwD
		)
		#stdB = √(m3disB_yerr .^2 .+ co5boldstd_τ[:, 2] .^2)
		axB[0].fill_between( 
			common_tauB, 
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 2].+co5boldstd_τ[:, 2]) .- ip(m3disB_x, m3disB_y)) ./ ip(m3disB_x, m3disB_y) *100.0,
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 2].-co5boldstd_τ[:, 2]) .- ip(m3disB_x, m3disB_y)) ./ ip(m3disB_x, m3disB_y) *100.0,
			color="red", alpha=0.3
		)
		axB[0].fill_between( 
			common_tauB, 
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 2]) .- ip(m3disB_x, m3disB_y.-m3disB_yerr)) ./ ip(m3disB_x, m3disB_y.-m3disB_yerr) *100.0,
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 2]) .- ip(m3disB_x, m3disB_y.+m3disB_yerr)) ./ ip(m3disB_x, m3disB_y.+m3disB_yerr) *100.0,
			color="0.5", alpha=0.3
		)

		

		
		## logτ vs. logρ
		axB[1].axhline(0.0, color="0.5", alpha=0.2)
		axB[1].plot(
			absolute_difference(
				mean,
				m3disB_xr, m3disB_yr,
				stagger_τ,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label=L"\rm Stagger - M3DIS\ (D)",
			lw=lwD
		)
		axB[1].fill_between( 
			absolute_difference(
				mean,
				m3disB_xr, 
				m3disB_yr.-m3disB_yerrr,
				stagger_τ,
				:log10τ_ross, :log10d
			)...,
			absolute_difference(
				mean,
				m3disB_xr, 
				m3disB_yr.+m3disB_yerrr,
				stagger_τ,
				:log10τ_ross, :log10d
			)[2],
			color="0.5", alpha=0.3
		)
		
		#=axB[1, 1].plot(
			absolute_difference(
				m3disB,
				muram_τ,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label="M3DIS - MURaM",
			ls="--",
			lw=lwD
		)=#
		axB[1].plot(
			common_tauB, 
			ip(marcs_models[5].structure["lgTauR"], log10.(marcs_models[5].structure["Density"])) .- ip(m3disB_xr, m3disB_yr),
			color="k",
			label=L"\rm MARCS - M3DIS\ (D)",
			ls=":",
			lw=lwD
		)
		axB[1].fill_between( 
			common_tauB, 
			(ip(marcs_models[5].structure["lgTauR"], log10.(marcs_models[5].structure["Density"])) .- ip(m3disB_xr, m3disB_yr.-m3disB_yerrr)),
			(ip(marcs_models[5].structure["lgTauR"], log10.(marcs_models[5].structure["Density"])) .- ip(m3disB_xr, m3disB_yr.+m3disB_yerrr)),
			color="0.5", alpha=0.3
		)


		axB[1].plot(
			common_tauB, 
			ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 3]) .- ip(m3disB_xr, m3disB_yr),
			color="k",
			label=L"\rm CO5BOLD - M3DIS\ (D)",
			ls="--",
			lw=lwD
		)
		axB[1].fill_between( 
			common_tauB, 
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 3].+co5boldstd_τ[:, 3]) .- ip(m3disB_xr, m3disB_yr)),
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 3].-co5boldstd_τ[:, 3]) .- ip(m3disB_xr, m3disB_yr)),
			color="red", alpha=0.3
		)
		axB[1].fill_between( 
			common_tauB, 
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 3]) .- ip(m3disB_xr, m3disB_yr.-m3disB_yerrr)),
			(ip(co5boldmean_τ[:, 1], co5boldmean_τ[:, 3]) .- ip(m3disB_xr, m3disB_yr.+m3disB_yerrr)),
			color="0.5", alpha=0.3
		)
		

		
		## logτ vs. rms
		axB[2].axhline(0.0, color="0.5", alpha=0.2)
		axB[2].plot(
			absolute_difference(
				rms5,
				m3disB_xu, m3disB_yu,
				stagger_τ,
				:log10τ_ross, :uz
			)...,
			color="k",
			label="M3DIS - Stagger",
			lw=lwD
		)  
		axB[2].fill_between( 
			absolute_difference(
				rms5,
				m3disB_xu, 
				m3disB_yu.-m3disB_yerru,
				stagger_τ,
				:log10τ_ross, :uz
			)...,
			absolute_difference(
				rms5,
				m3disB_xu, 
				m3disB_yu.+m3disB_yerru,
				stagger_τ,
				:log10τ_ross, :uz
			)[2],
			color="0.5", alpha=0.3
		)
		
		#=axB[2].plot(
			common_tauB, 
			ip(m3disB_xu, m3disB_yu) .- 
				ip(co5boldmean_τ[:, 1], co5boldrms_τ[:, 4]),
			color="k",
			label="M3DIS - Co5bold",
			ls="--",
			lw=lwD
		)
		axB[2].fill_between( 
			common_tauB, 
			(ip(m3disB_xu, m3disB_yu.-m3disB_yerru) .- 
				ip(co5boldmean_τ[:, 1], co5boldrms_τ[:, 4])),
			(ip(m3disB_xu, m3disB_yu.+m3disB_yerru) .- 
				ip(co5boldmean_τ[:, 1], co5boldrms_τ[:, 4])),
			color="0.5", alpha=0.3
		)=#
		
		#=axB[2, 1].plot(
			absolute_difference(
				rms5,
				m3disB,
				muram_τ,
				:log10τ_ross, :uz
			)...,
			color="k",
			label="M3DIS - MURaM",
			ls="--",
			lw=lwD
		)=#
	

		axB[0].set_xlim(-4, 2)
		axB[0].set_ylim(-6.25, 6.25)
		axB[1].set_ylim(-0.17, 0.17)
		axB[2].set_ylim(-0.75, 0.75)
	
		
		axB[1].legend(framealpha=0, ncol=1, labelspacing=0.01, handlelength=hl, loc="lower center")
		axB[0].set_ylabel(L"\rm \Delta\ \left(<T> \right)\ /\ <T>\ [\%]")
		axB[1].set_ylabel(
			L"\rm \Delta\ \left( <\log \rho>\right)"
		)
		axB[2].set_ylabel(
			L"\rm \Delta\ \left( rms\ v_{z}\right)\ [km\ s^{-1}]"
		)
	end

	axB[2].set_xlabel(L"\rm \log \tau_{ross}")
	axB[2].set_xlabel(L"\rm \log \tau_{ross}")
		

	fB.savefig("comparison_other_models.pdf", bbox_inches="tight")
	fB.savefig("comparison_other_models.png", bbox_inches="tight", dpi=600)

	gcf()
end

# ╔═╡ 6b63bde1-9a9a-4178-ac06-630ae811964b
md"### (C) Resolution Comparison"

# ╔═╡ 792cd263-3322-4a1f-9c1c-babcdcd88fe1
begin
	lresmodelC = models["lres"]
	iresmodelC = models["ires"]
	bestmodelC = models["best"]
	
	lresC_geo, lresC = pick_snapshot(out_folder[lresmodelC], :recent) 
	iresC_geo, iresC = pick_snapshot(out_folder[iresmodelC], :recent) 
	bestC_geo, bestC = pick_snapshot(out_folder[bestmodelC], :recent) 
	
	# lres
	lresC_x, lresC_y, lresC_z = time_average_profiles(
		mean, 
		out_folder[lresmodelC], 
		:log10τ_ross, 
		:T, 
		which=last
	)
	lresC_xr, lresC_yr, lresC_zr = time_average_profiles(
		mean, 
		out_folder[lresmodelC], 
		:log10τ_ross, 
		:log10d, 
		which=last
	)
	lresC_xu, lresC_yu, lresC_zu = time_average_profiles(
		rms5, 
		out_folder[lresmodelC], 
		:log10τ_ross, 
		:uz, 
		which=last
	)
	
	
	
	# ires
	iresC_x, iresC_y, iresC_z = time_average_profiles(
		mean, 
		out_folder[iresmodelC], 
		:log10τ_ross, 
		:T, 
		which=last
	)
	iresC_xr, iresC_yr, iresC_zr = time_average_profiles(
		mean, 
		out_folder[iresmodelC], 
		:log10τ_ross, 
		:log10d, 
		which=last
	)
	iresC_xu, iresC_yu, iresC_zu = time_average_profiles(
		rms5, 
		out_folder[iresmodelC], 
		:log10τ_ross, 
		:uz, 
		which=last
	)
	
	
	
	
	# best
	bestC_x, bestC_y, bestC_z = time_average_profiles(
		mean, 
		out_folder[bestmodelC], 
		:log10τ_ross, 
		:T, 
		which=last
	)
	bestC_xr, bestC_yr, bestC_zr = time_average_profiles(
		mean, 
		out_folder[bestmodelC], 
		:log10τ_ross, 
		:log10d, 
		which=last
	)
	bestC_xu, bestC_yu, bestC_zu = time_average_profiles(
		rms5, 
		out_folder[bestmodelC], 
		:log10τ_ross, 
		:uz, 
		which=last
	)

	brC = gresample(bestC_geo, nz=299)
	bestC_xH, bestC_yH = profile(MUST.mean, brC, :z, :log10τ_ross)

	brC = gresample(iresC_geo, nz=180)
	iresC_xH, iresC_yH = profile(MUST.mean, brC, :z, :log10τ_ross)

	brC = lresC_geo
	lresC_xH, lresC_yH = profile(MUST.mean, brC, :z, :log10τ_ross)
end

# ╔═╡ e02162d4-6345-43e2-bb3f-2a453cc1b9eb
begin
	# Figure
	plt.close()
	
	fC, axC = plt.subplots(3, 1, sharex=true, figsize=(6, 10))
	plt.subplots_adjust(wspace=0.25, hspace=0.0)
	#visual.basic_plot!.(axC)
	

	lwC = 2


	# relative comparison
	begin
		## logτ vs. T
		axC[0].axhline(0.0, color="0.5", alpha=0.2)
		#=axC[0].plot(
			relative_difference(
				mean,
				iresC_x, iresC_y,
				lresC_x, lresC_y,
				:log10τ_ross, :T
			)...,
			color="r",
			label=L"\rm\ low\ (B)\ -\ interm.\ (D)",#@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls="-"
		)=#
		axC[0].plot(
			relative_difference(
				mean,
				bestC_x, bestC_y,
				iresC_x, iresC_y,
				:log10τ_ross, :T
			)...,
			color="k",
			label=L"\rm\ interm.\ (B)\ -\ high\ (D)",#@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[0].fill_between( 
			relative_difference(
				mean,
				bestC_x, bestC_y.-bestC_z,
				iresC_x, iresC_y,#.+iresC_z,
				:log10τ_ross, :T
			)...,
			relative_difference(
				mean,
				bestC_x, bestC_y.+bestC_z,
				iresC_x, iresC_y,#.-iresC_z,
				:log10τ_ross, :T
			)[2],
			color="0.5", alpha=0.2
		)
		axC[0].plot(
			relative_difference(
				mean,
				bestC_x, bestC_y,
				lresC_x, lresC_y,
				:log10τ_ross, :T
			)...,
			color="k",
			label=L"\rm low\ (A)\ -\ high\ (D)",#@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)
		axC[0].fill_between( 
			relative_difference(
				mean,
				bestC_x, bestC_y.-bestC_z,
				lresC_x, lresC_y,#.+lresC_z,
				:log10τ_ross, :T
			)...,
			relative_difference(
				mean,
				bestC_x, bestC_y.+bestC_z,
				lresC_x, lresC_y,#.-lresC_z,
				:log10τ_ross, :T
			)[2],
			color="0.5", alpha=0.2
		)
		axC[0].axhline(0.0, color="0.5", alpha=0.2)
		

		
		
		## logτ vs. logρ
		axC[1].axhline(0.0, color="0.5", alpha=0.2)
		#=axC[1].plot(
			absolute_difference(
				mean,
				iresC_xr, iresC_yr,
				lresC_xr, lresC_yr,
				:log10τ_ross, :log10d
			)...,
			color="r",
			label=label=L"\rm\ low\ (B)\ -\ interm.\ (D)",
			lw=lwC,
			ls="-"
		)=#
		axC[1].plot(
			absolute_difference(
				mean,
				bestC_xr, bestC_yr,
				iresC_xr, iresC_yr,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label=label=L"\rm\ interm.\ (B)\ -\ high\ (D)",
			lw=lwC,
			ls=":"
		)
		axC[1].fill_between( 
			absolute_difference(
				mean,
				bestC_xr, bestC_yr.-bestC_zr,
				iresC_xr, iresC_yr,#.+iresC_zr,
				:log10τ_ross, :log10d
			)...,
			absolute_difference(
				mean,
				bestC_xr, bestC_yr.+bestC_zr,
				iresC_xr, iresC_yr,#.-iresC_zr,
				:log10τ_ross, :log10d
			)[2],
			color="0.5", alpha=0.2
		)
		axC[1].plot(
			absolute_difference(
				mean,
				bestC_xr, bestC_yr,
				lresC_xr, lresC_yr,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label=L"\rm low\ (A)\ -\ high\ (D)",
			lw=lwC,
			ls="--"
		)
		axC[1].fill_between( 
			absolute_difference(
				mean,
				bestC_xr, bestC_yr.-bestC_zr,
				lresC_xr, lresC_yr,#.+lresC_zr,
				:log10τ_ross, :log10d
			)...,
			absolute_difference(
				mean,
				bestC_xr, bestC_yr.+bestC_zr,
				lresC_xr, lresC_yr,#.-lresC_zr,
				:log10τ_ross, :log10d
			)[2],
			color="0.5", alpha=0.2
		)




		
		## logτ vs. rms
		axC[2].axhline(0.0, color="0.5", alpha=0.2)
		#=axC[2].plot(
			absolute_difference(
				rms5,
				iresC_xu, iresC_yu,
				lresC_xu, lresC_yu,
				:log10τ_ross, :uz
			)...,
			color="r",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls="-"
		)=#
		axC[2].plot(
			absolute_difference(
				rms5,
				bestC_xu, bestC_yu,
				iresC_xu, iresC_yu,
				:log10τ_ross, :uz
			)...,
			color="k",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[2].fill_between( 
			absolute_difference(
				mean,
				bestC_xu, bestC_yu.-bestC_zu,
				iresC_xu, iresC_yu,#.+iresC_zu,
				:log10τ_ross, :log10d
			)...,
			absolute_difference(
				mean,
				bestC_xu, bestC_yu.+bestC_zu,
				iresC_xu, iresC_yu,#.-iresC_zu,
				:log10τ_ross, :log10d
			)[2],
			color="0.5", alpha=0.2
		)
		axC[2].axhline(0.0, color="0.5", alpha=0.2)
		axC[2].plot(
			absolute_difference(
				rms5,
				bestC_xu, bestC_yu,
				lresC_xu, lresC_yu,
				:log10τ_ross, :uz
			)...,
			color="k",
			label=@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)
		axC[2].fill_between( 
			absolute_difference(
				mean,
				bestC_xu, bestC_yu.-bestC_zu,
				lresC_xu, lresC_yu,#.+lresC_zu,
				:log10τ_ross, :log10d
			)...,
			absolute_difference(
				mean,
				bestC_xu, bestC_yu.+bestC_zu,
				lresC_xu, lresC_yu,#.-lresC_zu,
				:log10τ_ross, :log10d
			)[2],
			color="0.5", alpha=0.2
		)

		
		

		axC[0].set_xlim(-4, 2)
		axC[0].set_ylim(-6.25, 6.25)
		axC[1].set_ylim(-0.075, 0.075)
		axC[2].set_ylim(-0.47, 0.47)


		# Test: add geometrical scale on top
		whereC = axC[0]
		refz = 0.7
		zdiffC = MUST.pyconvert(Vector, whereC.get_ylim())

		whereC.text(
			0.43, 1-refz, 
			"average geometrical scale", 
			ha="right", va="center", transform=whereC.transAxes,
			color="k", alpha=0.5
		)
		whereC.axhline(
			zdiffC[2] - (zdiffC[2] - zdiffC[1]) * refz,
			0.45, 1.0, color="k", alpha=0.1, lw=0.85
		)
		whereC.axhline(
			zdiffC[2] - (zdiffC[2] - zdiffC[1]) * (refz+0.2),
			0.0, 1.0, color="k", alpha=0.1, lw=0.85
		)
		
		constant_bestC_yH = fill!(
			similar(bestC_yH), zdiffC[2] - (zdiffC[2] - zdiffC[1]) * (refz+0.05)
		)
		whereC.plot(
			bestC_yH, constant_bestC_yH, 
			color="k", 
			ls="", 
			marker="|", 
			markersize=3,
			alpha=0.5
		)

		constant_iresC_yH = fill!(
			similar(iresC_yH), zdiffC[2] - (zdiffC[2] - zdiffC[1]) * (refz+0.1)
		)
		whereC.plot(
			iresC_yH, constant_iresC_yH, 
			color="k", 
			ls="", 
			marker="|", 
			markersize=3,
			alpha=0.5
		)

		constant_lresC_yH = fill!(
			similar(lresC_yH), zdiffC[2] - (zdiffC[2] - zdiffC[1]) * (refz+0.15)
		)
		whereC.plot(
			lresC_yH, constant_lresC_yH, 
			color="k", 
			ls="", 
			marker="|", 
			markersize=3,
			alpha=0.5
		)
		
		
		axC[1].legend(framealpha=0, loc="upper center", ncol=1, labelspacing=0.01, handlelength=hl)
		axC[0].set_ylabel(L"\rm \Delta\ \left(<T> \right)\ /\ <T>\ [\%]")
		#axC[1].set_ylabel(
		#	L"\rm \Delta\ \left( <\log \rho>\right)\ [g \times cm^{-3}]"
		#)
		axC[1].set_ylabel(
			L"\rm \Delta\ \left( <\log \rho>\right)"
		)
		axC[2].set_ylabel(
			L"\rm \Delta\ \left( rms\ v_{z}\right)\ [km\ s^{-1}]"
		)
	end

	axC[2].set_xlabel(L"\rm \log \tau_{ross}")
	axC[2].set_xlabel(L"\rm \log \tau_{ross}")
		

	fC.savefig("comparison_resolution.pdf", bbox_inches="tight")
	fC.savefig("comparison_resolution.png", bbox_inches="tight", dpi=600)

	gcf()
end

# ╔═╡ b9890085-fcca-472d-b192-da7abf349b45
md"### (D) Opacity Table"

# ╔═╡ b14da0d1-63da-44aa-b1fd-4340ec64528e
begin
	modelD = models["best"]
	eosD = eos[modelD]
	opaD = opa[modelD]

	# the formation opacities are stored in the unbinned table
	fopaD = reload(
		SqOpacity, 
		joinpath(
			"/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6/",
			"combined_formation_opacities.hdf5"
		),
		mmap=true
	)
	

	# the bin assignment should be saved for this opacity table
	fid = TSO.HDF5.h5open(
		joinpath(eos_folder[modelD], "bin_assignment.hdf5"), 
		"r"
	)
	
	binassignment = TSO.HDF5.read(fid["bins"])
	λ = TSO.HDF5.read(fid["lambda"])
	close(fid)


	
	# Figure
	cD = plt.cm.gnuplot2
	normD = matplotlib.colors.BoundaryNorm(
		range(0.5, length(opaD.λ)+0.5, step=1) |> collect, 
		cD.N - 50
	)
	
	plt.close()
	fD, axD = plt.subplots(1, 1, figsize=(5,6))
	#visual.basic_plot!(axD)
	
	im = axD.scatter(
		log10.(λ), -log10.(fopaD.κ_ross), 
		c=binassignment,
		s=1,
		cmap=cD,
		norm=normD,
		rasterized=true
	)
	cbarD = fD.colorbar(
		im, 
		ax=axD, 
		ticks=range(1, length(opaD.λ), step=1) |> collect
	)


	axD.set_ylabel(
		L"\rm - \log \tau_{ross}\ (\log \tau_{\lambda}=1)"
	)
	axD.set_xlabel(
		L"\rm \log \lambda\ [\AA]"
	)
	cbarD.set_label(
		L"\rm opacity\ bins"
	)
	
	fD.savefig("formation_opacities.pdf", bbox_inches="tight", dpi=600)
	fD.savefig("formation_opacities.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ 516fe726-9d67-44dd-bbdb-f0ec17841185
md"### (E) Effective Temperature"

# ╔═╡ dc93e1ef-c947-465c-ab37-b1fd251b695d
md"### (F) Radiative Heating (Binning)"

# ╔═╡ 64213279-f47a-4b46-b1ac-730945d5b8f1
md"### (G) 2D surface"

# ╔═╡ da45e33e-9013-443e-be60-e0b53d23d305
surface_optical_ip(snap, what) = begin
	#MUST.closest(log10.(MUST.axis(snap, :τ_ross, 3)), 0)
	MUST.interpolate_to(snap, what, τ_ross=0.0, logspace=true)
end

# ╔═╡ d59e82db-89f3-46fd-8af4-6bf96fab14cf
surface_optical(snap, what) = begin
	isurf = surface_optical_ip(snap, what)
	isurf[what][:, :, 1] 
end

# ╔═╡ c9004d63-d6e5-415a-87ef-3db717a3f788
uz_surface_optical(snap) = surface_optical(snap, :uz)

# ╔═╡ b214c18c-58ee-474f-809d-d07582ab5647
rho_surface_optical(snap) = surface_optical(snap, :d)

# ╔═╡ 4935697c-974b-4de8-8182-b5c33f3628f4
temp_surface_optical(snap) = surface_optical(snap, :T)

# ╔═╡ acf85ab7-3d8b-4e83-8a73-1ed00598882f
extent(snap) = begin
	x = MUST.axis(snap, :x) ./1e8
	y = MUST.axis(snap, :x) ./1e8

	[minimum(x), maximum(x), minimum(y), maximum(y)]
end

# ╔═╡ 4ae6b55e-952a-4b70-937f-c81a2f790e83
begin
	modelG = models["best"]
	m3disG = pick_snapshot(out_folder[modelG], -1) |> last
	
	fG, axG = plt.subplots(1, 1, figsize=(5, 6))

	# limits for color bar
	uzG  = uz_surface_optical(m3disG) ./1e5
	#uzG = surface_optical(m3disG, :d)
	vmin = minimum([minimum(u) for u in uzG])
	vmax = maximum([maximum(u) for u in uzG])

	vmax = min(abs.([vmin, vmax])...)
	vmin = -vmax

	# plot the surface in uz
	imG = axG.imshow(
		uzG,
		origin="lower",
		vmin=vmin, vmax=vmax,
		cmap="coolwarm",
		extent=extent(m3disG)
	)

	# setup colorbar
	cG = fG.colorbar(imG, ax=axG, fraction=0.046, pad=0.04)
	cG.set_label(L"\rm v_z\ [km \ s^{-1}]")

	# Other labels
	axG.set_xlabel("x [Mm]")
	axG.set_ylabel("y [Mm]")

	#fG.savefig("vertical_velocity_surface.pdf", bbox_inches="tight")
	fG.savefig("vertical_velocity_surface.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ 8324d6c9-f46a-4724-be76-5eb92f753e21
let
	modelG = models["lres"]
	m3disG = pick_snapshot(out_folder[modelG], -1) |> first
	
	fG, axG = plt.subplots(1, 1, figsize=(5, 6))

	# plot the surface in uz
	imG = axG.imshow(
		m3disG[:T][:,:,end],
		origin="lower",
		cmap="coolwarm",
		extent=extent(m3disG)
	)

	gcf()
end

# ╔═╡ faccedff-2a3a-40b1-ae88-c42a69112d16
md"### (H) 3D cube with velocities"

# ╔═╡ 1634883c-2a93-4b31-bc3a-662a894733c4
begin
	plt.close()

	modelH = models["best"]
	m3disH = pick_snapshot(out_folder[modelH], -1) |> first
	
	fH, axH = visual.cube_with_velocities(
		m3disH, vmax_3d=17000, cmap="RdYlBu_r"
	)


	#fG.savefig("temperature_cube.pdf", bbox_inches="tight")
	fH.savefig("temperature_cube.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ f409f3e8-ef97-4fb9-a8e6-7a8f1e2b2d22
md"### (I) M3D resampling effect"

# ╔═╡ 759b0406-1a86-42d9-b8d7-8784a1574c02
begin
	modelI = models["best"]
	m3disI = pick_snapshot(out_folder[modelI], -1) |> first
	m3disI_τ = pick_snapshot(out_folder[modelI], -1) |> last

	snap_m3disI = gresample(
		m3disI, 
		nz=299
	)
	snap_m3disI_τ = MUST.height_scale_fast(snap_m3disI, :τ_ross)
	

	# We look at the surface and average profile of this snapshot for different 
	# resamplings
	samplings = reshape(collect(Iterators.product((5, 10, 80), (299,))), :)

	fI, axI = plt.subplots(3, 2, figsize=(10,15))
	for j in 0:1
		for i in 0:2
			#visual.basic_plot!(axI[i, j])
		end
	end
	plt.subplots_adjust(wspace=0.3, hspace=0.1)

	vminI = 5500
	vmaxI = 6900
	i0I_full = argmin(abs.(log.(MUST.axis(m3disI_τ, :τ_ross, 3))))
	TI_full = m3disI[:T][:, :, i0I_full]
	
	for (i, resolution) in enumerate(samplings)
		snap = gresample(
			m3disI, 
			nx=resolution[1], 
			ny=resolution[1],
			nz=resolution[2],
			method=:linear
		)

		snap_r = MUST.height_scale_fast(snap, :τ_ross)
		
		i0I_i = argmin(abs.(log.(MUST.axis(snap_r, :τ_ross, 3)) .- 0.0))
		TI_i  = snap_r[:T][:, :, i0I_i]
		
		imI_i = axI[i-1, 0].imshow(
			TI_i,
			origin="lower",
			vmin=vminI, vmax=vmaxI,
			cmap="hot",
			extent=extent(snap_r), 
			aspect="auto"
		)
		@info minimum(TI_i) maximum(TI_i)

		#=imI_i = axI[i-1, 1].imshow(
			TI_full,
			origin="lower",
			vmin=vminI, vmax=vmaxI,
			cmap="hot",
			extent=extent(m3disI), 
			aspect="auto"
		)=#

		#cI = fI.colorbar(imI_i, ax=axI[i-1, 0], fraction=0.046, pad=0.04)
		#cI.set_label(L"\rm T\ [K]", fontsize="large")

		xI_i , yI_i = profile(MUST.mean, snap, :z, :T)		
		axI[i-1, 1].plot(
			xI_i, yI_i,
			color="k", 
			marker=".", 
			markersize=5,
			ls="",
			label="$(resolution[1])x$(resolution[1])x$(resolution[2])"
		)
		
		xI_i , yI_i = profile(MUST.mean, snap_m3disI, :z, :T)		
		axI[i-1, 1].plot(
			xI_i, yI_i,
			color="tomato", 
			marker="",
			ls="-",
			lw=2,
			label="original"
		)

		#axI[i-1, 1].set_ylim(4000, 10000)
		#axI[i-1, 1].set_xlim(-4.0, 1.75)
		axI[i-1, 1].legend(framealpha=0, labelspacing=0.15, handlelength=hl)
		axI[i-1, 1].set_ylabel(L"\rm T\ [K]")
	end

	for (i, resolution) in enumerate(samplings)
		axI[i-1, 0].set_ylabel(L"\rm Y\ [Mm]")	
	end
	axI[2, 0].set_xlabel(L"\rm X\ [Mm]")
	axI[2, 1].set_xlabel(L"\rm \log \tau_{ross}")	


	fI.savefig("resampling_effect.png", bbox_inches="tight", dpi=600)
	fI.savefig("resampling_effect.pdf", bbox_inches="tight")
	
	gcf()
end

# ╔═╡ 82b1bd73-d81a-404c-aefd-643b7008d2b7
md"### (J) Average profiles of other stars"

# ╔═╡ 306dd6bc-c27d-4af4-9188-6e96981541fb
begin
	plt.close()

	modelsJ = [
		models["t55g45m00"], 
		models["t60g45m00"], 
		models["t65g45m00"],
		#models["t45g40m00"]
	]
	lsJ     = ["-", "--", "-", ":"]
	lwJ     = [1.5, 1.5, 3, 2.0]
	labelsJ = [
		L"\rm 5500\ K, 4.5\ dex", 
		L"\rm 6000\ K, 4.5\ dex", 
		L"\rm 6500\ K, 4.5\ dex",
		#L"\rm 4500\ K, 4.0\ dex"
	]
	
	
	fJ, axJ = plt.subplots(3, 1, sharex=true, figsize=(5, 9))
	plt.subplots_adjust(wspace=0, hspace=0)
	#visual.basic_plot!.(axJ)

	for (i, model) in enumerate(modelsJ)
		m3disJ = pick_snapshot(out_folder[model], :recent) |> last
		m3disJ_x, m3disJ_y, m3disJ_yerr = time_average_profile(
			mean, 
			out_folder[model], 
			:log10τ_ross, 
			:T, 
			hscale=:τ
		)
		m3disJ_xr, m3disJ_yr, m3disJ_yerrr = time_average_profile(
			mean, 
			out_folder[model], 
			:log10τ_ross, 
			:log10d, 
			hscale=:τ
		)
		m3disJ_xu, m3disJ_yu, m3disJ_yerru = time_average_profile(
			rms5, 
			out_folder[model], 
			:log10τ_ross, 
			:uz, 
			hscale=:τ
		)

		# logτ vs. T
		axJ[0].plot(
			m3disJ_x, m3disJ_y,
			color="k",
			ls=lsJ[i],
			label=labelsJ[i],
			lw=lwJ[i]
		)

		# logτ vs. rho
		axJ[1].plot(
			m3disJ_xr, m3disJ_yr,
			color="k",
			ls=lsJ[i],
			label=labelsJ[i],
			lw=lwJ[i]
		)

		# logτ vs. u
		axJ[2].plot(
			m3disJ_xu, m3disJ_yu,
			color="k",
			ls=lsJ[i],
			label=labelsJ[i],
			lw=lwJ[i]
		)
		
	end

	#axJ[0].plot(
	#	marcs_models[4].structure["lgTauR"], marcs_models[4].structure["T"]
	#)

	axJ[0].legend(framealpha=0, labelspacing=0.01, handlelength=hl)
	
	axJ[0].set_ylabel(L"\rm T\ [K]")
	axJ[1].set_ylabel(L"\rm \log \rho\ [g \ cm^{-3}]")
	axJ[2].set_ylabel(L"\rm rms\ v_{z}\ [km \ s^{-1}]")
	
	axJ[2].set_xlabel(L"\rm \log \tau_{ross}")
	axJ[2].set_xlim(-3.75, 3)
	axJ[0].set_ylim(2600, 13500)
	axJ[1].set_ylim(-8.75, -5.85)

	fJ.savefig("average_other_models.pdf", bbox_inches="tight")
	fJ.savefig("average_other_models.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ 521f69d2-c668-4e95-b422-bc95a858286c
md"### (K) 2D surfaces"

# ╔═╡ bec4d310-281f-4a74-bd61-cebbfd1872f5
begin
	modelK = models["best"]

	modelsK = [
		modelK, 
		models["t55g45m00"], 
		models["t60g45m00"], 
		models["t65g45m00"], 
		models["t45g40m00"]
	]
	teffK   = [5777, 5500, 6000, 6500, 4500]
	loggK   = [4.44, 4.5, 4.5, 4.5, 4.0]
	fsK = 12

	for (i, model) in enumerate(modelsK)
		plt.close()
		
		fK, axK = plt.subplots(1, 1, figsize=(5, 6))
		#visual.basic_plot!(axK)
		
		m3disK = pick_snapshot(out_folder[model], -2) |> last

		# limits for color bar
		uzK  = uz_surface_optical(m3disK) ./1e5
		vmin = minimum(uzK)
		vmax = maximum(uzK)

		vmax = min(abs.([vmin, vmax])...)
		vmin = -vmax
	
		# plot the surface in uz
		imK = axK.imshow(
			uzK,
			origin="lower",
			vmin=vmin, vmax=vmax,
			cmap="seismic_r",
			extent=extent(m3disK)
		)
		
		#axK[i-1].set_title(labels[model])
		
		# setup colorbar
		cK = fK.colorbar(imK, ax=axK, fraction=0.046, pad=0.04)
		cK.set_label(L"\rm v_z\ [km\ s^{-1}]", fontsize=fsK)
		axK.tick_params(axis="both", labelsize=fsK-2)
		cK.ax.tick_params(axis="both", labelsize=fsK-2)
		axK.set_ylabel("y [Mm]", fontsize=fsK)	
		axK.set_xlabel("x [Mm]", fontsize=fsK)
		

		#=axK.text(
			0.97, 0.97, 
			L"\rm T_{eff}"*" $(teffK[i]) K, "*"log(g)"*" $(loggK[i]) dex", 
			ha="right", va="top", 
			transform=axK.transAxes,
			color="white", fontsize="large", backgroundcolor="k"
		)=#
		axK.set_title(
			L"\rm T_{eff}"*" $(teffK[i]) K, "*"log(g)"*" $(loggK[i]) dex",
			fontsize=fsK
		)

		fK.savefig("vertical_velocity_surface_$(labels[model]).pdf",) #bbox_inches="tight")
		fK.savefig("vertical_velocity_surface_$(labels[model]).png",)
		#bbox_inches="tight", dpi=600)
	end
	
	gcf()
end

# ╔═╡ 3c3d3319-1e82-4027-8ac4-a460acb15518
let
	for (i, model) in enumerate(modelsK)
		plt.close()
		
		fK, axK = plt.subplots(1, 1, figsize=(5, 6))
		#visual.basic_plot!(axK)
		
		m3disK = pick_snapshot(out_folder[model], -2) |> last

		# limits for color bar
		uzK  = temp_surface_optical(m3disK)
		vmin = minimum(uzK)
		vmax = maximum(uzK)
	
		# plot the surface in uz
		imK = axK.imshow(
			uzK,
			origin="lower",
			vmin=vmin, vmax=vmax,
			cmap="gist_heat",
			extent=extent(m3disK)
		)
		
		#axK[i-1].set_title(labels[model])
		
		# setup colorbar
		cK = fK.colorbar(imK, ax=axK, fraction=0.046, pad=0.04)
		cK.set_label(L"\rm T\ [K]", fontsize=fsK)
		axK.tick_params(axis="both", labelsize=fsK-2)
		cK.ax.tick_params(axis="both", labelsize=fsK-2)
		axK.set_ylabel("y [Mm]", fontsize=fsK)	
		axK.set_xlabel("x [Mm]", fontsize=fsK)
		#=axK.text(
			0.97, 0.97, 
			L"\rm T_{eff}"*" $(teffK[i]) K, "*"log(g)"*" $(loggK[i]) dex", 
			ha="right", va="top", 
			transform=axK.transAxes,
			color="white", fontsize="large", backgroundcolor="k"
		)=#
		axK.set_title(
			L"\rm T_{eff}"*" $(teffK[i]) K, "*"log(g)"*" $(loggK[i]) dex", 
			fontsize=fsK
		)

		fK.savefig("temperature_surface_$(labels[model]).pdf")
		fK.savefig("temperature_surface_$(labels[model]).png")
		#bbox_inches="tight", dpi=600)
	end
	
	gcf()
end

# ╔═╡ 04543c5d-d519-4ed2-9d98-050fbd1361a1
md"### (L) Oszillations"

# ╔═╡ 6dc706bc-0fbd-4277-948f-df7b4306bf4d
scipyFFT = MUST.pyimport("scipy.fft")

# ╔═╡ 6a2d0807-a437-40db-9ed3-934b4a1018b9
fft(args...; kwargs...) = MUST.pyconvert(Array, scipyFFT.fftn(args...; kwargs...))

# ╔═╡ 848900a2-3a56-48e2-b771-45ec8dc3b392
fftfreq(axis) = begin
	s = length(axis)
	d = abs(diff(axis) |> first)
	MUST.pyconvert(Array, scipyFFT.fftfreq(s, d=d))
end

# ╔═╡ ba8f433c-400a-4829-b534-1b88e35808f9
begin
	modelL = models["best"]
	m3disL = pick_snapshot(out_folder[modelL], :recent) |> first

	zL, velocityL = profile(MUST.mean, m3disL, :z, :uz)
	velocityL ./= 1e5
	velocityFFTL = fft(velocityL)
	velocityFFTabsL = sqrt.(real(velocityFFTL).^2 .+ imag(velocityFFTL).^2)
end

# ╔═╡ 1977a610-41f0-4288-9c5e-58a140acdab3
begin
	plt.close()
	fL, axL = plt.subplots(1, 1, figsize=(5, 6))
	#visual.basic_plot!(axL)
	
	x_fft = 1 ./ fftfreq(MUST.axis(m3disL, :z) ./1e8)
	y_fft = velocityFFTabsL
	axL.plot(x_fft, y_fft, lw=1, label=L"\rm M3DIS\ (D)", color="k")

	axL.set_xlim(-0.2, 0.2)
	axL.set_ylim(0.05, 1.1)
	axL.legend(framealpha=0, handlelength=hl)
	gcf()
end

# ╔═╡ 1f458a9f-a0dd-459b-8a60-aded8d8e4617
md"### (M) Vertical Slices"

# ╔═╡ d602d302-5964-4525-b62e-d7eb7aec5182
begin
	modelM = models["best"]
	m3disM = pick_snapshot(out_folder[modelM], -1) |> first
	
	plt.close()
	fM, axM = plt.subplots(1, 1, figsize=(10, 6))
	#visual.basic_plot!(axM)
	pickeveryM = 3

	# We plot the velocity in x and z direction in a plane of contant y in the center
	y0 = 200 #argmin(abs.(MUST.axis(m3disM, :y)))
	y_y0 = MUST.axis(m3disM, :y)[y0] ./1e8
	ux = m3disM[:ux][1:pickeveryM:end, y0, 1:pickeveryM:end] ./1e5
	uz = m3disM[:uz][1:pickeveryM:end, y0, 1:pickeveryM:end] ./1e5

	xxM = m3disM.x[1:pickeveryM:end, y0, 1:pickeveryM:end] ./1e8
	zzM = m3disM.z[1:pickeveryM:end, y0, 1:pickeveryM:end] ./1e8

	TM = m3disM[:T][1:pickeveryM:end, y0, 1:pickeveryM:end]
	
	csM = axM.contour(xxM, zzM, TM, cmap="seismic", levels=10, linewidths=4, alpha=0.5)
	plt.clabel(csM, inline=1)

	norm = matplotlib.colors.Normalize(vmin=csM.cvalues.min(), vmax=csM.cvalues.max())
	sm = plt.cm.ScalarMappable(norm=norm, cmap = csM.cmap)
	#fM.colorbar(sm, ax=axM, fraction=visual.cbar_fraction, pad=visual.cbar_pad)
	
	axM.quiver(xxM, zzM, ux, uz, color="k", scale=290, zorder=100)
	axM.axhline(0.0, color="k", ls="-", lw=2)

	axM.set_ylabel(L"\rm Z\ [Mm]", fontsize="large")
	axM.set_xlabel(L"\rm X\ [Mm]", fontsize="large")
	
	
	axM.set_xlim(-2.3, 2.3)
	axM.set_ylim(-1.6, 0.7)

	fM.savefig("vertical_slice_velocity_sun.png", dpi=600, bbox_inches="tight")
	fM.savefig("vertical_slice_velocity_sun.pdf", bbox_inches="tight")
	
	gcf()
end

# ╔═╡ 78ce2d34-9898-4140-a7fd-8fe2e807d972
scaleheight(p, ρ, logg) = begin
	p / (ρ * exp10(logg)) ./ 1e5
end

# ╔═╡ 94c55529-523a-4a44-9ed8-1caaf34e4da5
convective_turnover(hp, v) = abs(hp ./ v)

# ╔═╡ 3466e7cb-17ab-4a1a-ac4f-7b3f1264584d
md"We can measure the average vertical velocity and average scale height within one granular slice, as seen in the previous figure."

# ╔═╡ cb14997c-ed1f-4be0-80ce-bba1d94765d4
function granularstatistics(model; 
			x_limits=[-2.0, -1], 
			y_limits=[y_y0-0.5, y_y0+0.5], 
			z_limits=[-1.0, 0.1])
	
	m3disM = pick_snapshot(out_folder[model], -1) |> first
	pgM = exp.(lookup(eos[model], :lnPg, log.(m3disM[:d]), log.(m3disM[:ee])))
	loggM = 4.44

	xgranM = MUST.axis(m3disM, :x) ./1e8
	ygranM = MUST.axis(m3disM, :y) ./1e8
	zgranM = MUST.axis(m3disM, :z) ./1e8
	
	granularfilter(arr, low, high) = low .< arr .< high
	mask = [
		granularfilter(xgranM, x_limits...),
		granularfilter(ygranM, y_limits...),
		granularfilter(zgranM, z_limits...)
	]

	pgranM = pgM[mask...]
	ρranM  = m3disM[:d][mask...]
	vranM  = m3disM[:uz][mask...] ./1e5

	hpM = scaleheight(mean(pgranM), mean(ρranM), loggM)
	ctM = convective_turnover(hpM, mean(vranM))

	@info "Pressure scale height [km]: $(hpM)"
	@info "Convective turnover time [s]: $(ctM)"
	@info "Convective turnover time [min]: $(ctM ./60)"
	@info "________________________________________________________"
	@info "Convective turnover time [code units]: $(ctM ./100)"

	hpM, ctM
end

# ╔═╡ f6736f9c-f76e-40ae-bb85-63c1aef6605c
granularstatistics(modelM)

# ╔═╡ 9cbc4ead-5f1d-41b9-a7ac-dea24f4b91cb


# ╔═╡ f33f10e7-506d-4d55-a7c2-1b573ee6e13e
md"We can also use the optical surface for this statistics. 
For this we add the pressure to the model and compute the scale height and turn over time at the optical surface using the pre-defined functions."

# ╔═╡ c11b167f-ba1f-4ae6-8e38-e3c011f6f036
function optical_turnover(model; logg)
	m3disM = pick_snapshot(out_folder[model], -1) |> first
	pgM = exp.(lookup(eos[model], :lnPg, log.(m3disM[:d]), log.(m3disM[:ee])))
	
	add!(m3disM, Pg=pgM)

	hpM = MUST.pressurescaleheight(m3disM, τ_ross=5.0, logspace=true, logg=logg)
	ctM = MUST.convectiveturnovertime(m3disM, τ_ross=5.0, logspace=true, logg=logg)

	@info "Pressure scale height [km]: $(hpM/100000.0)"
	@info "Convective turnover time [s]: $(ctM)"
	@info "Convective turnover time [min]: $(ctM ./60)"

	hpM, ctM
end

# ╔═╡ 3f194d77-08cf-4d51-81ee-f811ed587967
optical_turnover(models["best"], logg=4.44)

# ╔═╡ 95849fd6-e8bf-48e4-a0d9-2dd86cee58bb
optical_turnover(models["t55g45m00"], logg=4.5)

# ╔═╡ a15e914f-13f6-456f-a74f-b9edd63a7ab9
optical_turnover(models["t60g45m00"], logg=4.5)

# ╔═╡ 059978b8-a8b0-4f1c-a681-ea41d865a488
optical_turnover(models["t65g45m00"], logg=4.5)

# ╔═╡ f3a0a1b6-af99-46ec-8355-ab9fe7e5881d


# ╔═╡ 2587d8e9-7612-42f6-8a0c-01c72b1d32d3
md"### (N) Multiple vertical slices"

# ╔═╡ 5003a0ce-7aa0-4a5e-9765-6129d7660e1b
granularstatistics(models["t55g45m00"])

# ╔═╡ 6ce71e68-3e46-4ee7-9777-7a04342b5e2e
granularstatistics(models["t60g45m00"])

# ╔═╡ 09c3002d-1839-4adf-9d47-1811f1c81bcd
granularstatistics(models["t65g45m00"])

# ╔═╡ e8da9c4e-3474-41ea-a441-a22b089d20d6
rellim(arr) = begin
	(minimum(arr) + 0.005*(maximum(arr)-minimum(arr)),
	 maximum(arr) - 0.005*(maximum(arr)-minimum(arr)))
end

# ╔═╡ 0f3f07fa-660e-4117-a46e-792150a707b9
function vertical_velocity_slice_with_contours(
	mN, y0i, label, marcsT, marcsZ, levels; 
	pickevery=1, quantity=:T, scale=identity,
	cbar_label="")

	ipN = MUST.linear_interpolation(marcsT, marcsZ)
	
	plt.close()
	fN, axN = plt.subplots(1, 1, figsize=(10, 6))		
	
	m3disNi = pick_snapshot(out_folder[mN], -1) |> first
	#m3disNi = gresample(
	#	m3disNi, 
	#	nz=ceil(Int, size(m3disNi, 1)+1), 
	#	nx=ceil(Int, size(m3disNi, 1)+1), 
		#ny=ceil(Int, size(m3disNi, 1) /2)
	#)
	
	y_y0N = MUST.axis(m3disNi, :y)[y0i] ./1e8
	uxN = m3disNi[:ux][1:pickevery:end, y0i, 1:pickevery:end] ./1e5
	uzN = m3disNi[:uz][1:pickevery:end, y0i, 1:pickevery:end] ./1e5

	xxN = m3disNi.x[1:pickevery:end, y0i, 1:pickevery:end] ./1e8
	zzN = m3disNi.z[1:pickevery:end, y0i, 1:pickevery:end] ./1e8

	TN = scale.(m3disNi[quantity][1:pickevery:end, y0i, 1:pickevery:end])
	TN .= MUST.pyconvert(Array, (ndimage.gaussian_filter(TN, sigma=2.8)))

	csN = axN.contourf(
		xxN, zzN, TN, 
		cmap="coolwarm", levels=levels, alpha=0.99, #linewidths=3,
	)
	axN.contour(
		xxN, zzN, TN, colors="k", #cmap="terrain",
		levels=levels, alpha=0.3, linewidths=1.5, linestyles="solid",
	)
	#plt.clabel(csN, inline=1)

	#axN.set_aspect("equal")
	
	normN = matplotlib.colors.Normalize(
		vmin=csN.cvalues.min(), vmax=csN.cvalues.max()
	)
	smN = plt.cm.ScalarMappable(norm=normN, cmap=csN.cmap)
	fN.colorbar(smN, ax=axN, label=cbar_label).ax.invert_yaxis()
	#fraction=visual.cbar_fraction, pad=visual.cbar_pad)

	
	#axN.quiver(
	#	xxN, zzN, uxN, uzN, color="k", scale=220, zorder=100, headwidth=3.2
	#)

	
	xxN = np.linspace(first(xxN), last(xxN), size(xxN, 1))
	zzN = np.linspace(first(zzN), last(zzN), size(zzN, 2))	

	lw = sqrt.(uxN.^2 .+ uzN .^2)
	a = 0.05
	b = 4.5
	min_val = minimum(lw)
    max_val = maximum(lw)
    lw = a .+ (lw .- min_val) .* ((b - a) / (max_val - min_val))

	color_binary = similar(uzN)
	a = -1.0
	b = 1.0
	min_val = minimum(uzN)
    max_val = maximum(uzN)
    color_binary = a .+ (uzN .- min_val) .* ((b - a) / (max_val - min_val))

	
	uxNs = MUST.pyconvert(Array, (ndimage.gaussian_filter(uxN, sigma=2.8)))
	uzNs = MUST.pyconvert(Array, (ndimage.gaussian_filter(uzN, sigma=2.8)))
	
	axN.streamplot(
		xxN, 
		zzN,
		MUST.Py(uxNs).to_numpy().T, 
		MUST.Py(uzNs).to_numpy().T, 
		color="k",#MUST.Py(TN).to_numpy().T, #MUST.Py(color_binary).to_numpy().T,
		zorder=100,
		cmap="terrain",
		density=0.5, 
		arrowstyle="fancy",
		#integration_direction="forward",
		broken_streamlines=false,
		linewidth=MUST.Py(lw).to_numpy().T, 
		arrowsize=1.3
	)
	
	axN.axhline(0.0, color="k", ls="-", lw=2)

	axN.set_xlim(rellim(xxN)...)
	axN.set_ylim(rellim(zzN)...)


	# add a 1D model on top
	#if mN == models["best"]
		to_the_rightN = MUST.pyconvert(Any, axN.get_xlim()[1] - axN.get_xlim()[0])
		offsetN = MUST.pyconvert(Any, axN.get_xlim()[0])
		byN = -0.13
		y_top = MUST.pyconvert(Any, axN.get_ylim()[1])
	
		levelsN = MUST.pyconvert(Vector, csN.levels)
		cmapN = csN.get_cmap()
		colorsN = MUST.pyconvert(
			Vector, cmapN(range(0, 1, length=length(levelsN)) |> collect)
		)
	
		xmarcsN = [l for l in levelsN if marcsT[1] < l < marcsT[end]]
		colorsN = [
			c for (j, c) in enumerate(colorsN) 
				if marcsT[1] < levelsN[j] < marcsT[end]
		]
		ymarcsN = ipN.(xmarcsN)
	
		for point in eachindex(xmarcsN)
			if ymarcsN[point]/1e8 > y_top
				continue
			end
			#@show point xmarcsN[point] ymarcsN[point]/1e8
			axN.text(
				offsetN + (byN - 0.07)*to_the_rightN, 
				ymarcsN[point]/1e8, 
				"$(xmarcsN[point]) ",
				ha="right", va="center", fontsize="small",
				color="k" #, zorder=110+point
			)
			axN.text(
				offsetN + byN*to_the_rightN, 
				ymarcsN[point]/1e8, 
				L"-",
				ha="right", va="center", fontsize=25,
				color=colorsN[point]#, zorder=110+point
			)
		end

		axN.text(
				offsetN + byN*to_the_rightN, 
				min(maximum(marcsZ)/1e8, y_top), 
				"MARCS ",
				color="k",
				ha="right", va="bottom"
			)
	
		lineN = matplotlib.lines.Line2D(
			[offsetN + byN*to_the_rightN, offsetN + byN*to_the_rightN],
			[minimum(marcsZ)/1e8, min(maximum(marcsZ)/1e8, y_top)], 
			color="k", lw=2,
			zorder=150, clip_on=false
		)
		axN.add_line(lineN)
	#end
	

	#=τcN = granularstatistics(
		modelsN[i], 
		x_limit=x_limitN[i], y_limit=y_y0N.+y_limitN[i], z_limit=z_limitN[i]
	)

	labi = labelsN[i]*L"\rm \tau_{c}="*"$(τcN./60) min"=#

	#=axN.text(
		0.97, 0.05, label, 
		ha="right", va="bottom", 
		transform=axN.transAxes,
		color="white", backgroundcolor="k", zorder=150
	)=#
	axN.set_title(
		label
	)

	axN.set_ylabel("\n\n"*L"\rm Z\ [Mm]")
	axN.set_ylabel("\n\n"*L"\rm Z\ [Mm]")
	
	axN.set_xlabel(L"\rm X\ [Mm]")
	axN.set_xlabel(L"\rm X\ [Mm]")
	
	quantity = quantity==:τ_ross ? "tross" : quantity
	fN.savefig(
		"vertical_slice_velocity-$(labels[mN])_$(quantity).png", dpi=600#, bbox_inches="tight"
	)
	fN.savefig(
		"vertical_slice_velocity-$(labels[mN])_$(quantity).pdf"#, bbox_inches="tight"
	)

	fN, axN
end

# ╔═╡ 22bc6765-400a-470b-9f9d-ff1365d97e33
begin
	modelsN = [
		models["best"],
		models["t55g45m00"], 
		#models["t60g45m00"], 
		models["t65g45m00"],
		models["t45g40m00"]
	]
	
	y0N = [
		200,
		20,
		#140,
		210,
		200
	]

	x_limitN = []
	y_limitN = []
	z_limitN = []

	labelsN = [
		L"\rm T_{eff}\ 5777\ K, log(g)\ 4.44\ dex",
		L"\rm T_{eff}\ 5500\ K, log(g)\ 4.5\ dex",
		#L"\rm 6000\ K, 4.5\ dex", 
		L"\rm T_{eff}\ 6500\ K, log(g)\ 4.5\ dex",
		L"\rm T_{eff}\ 4500\ K, log(g)\ 4.0\ dex"
	]

	marcs_TN = [
		marcs_models[5].structure["T"],
		marcs_models[1].structure["T"],
		marcs_models[3].structure["T"],
		marcs_models[4].structure["T"]
	]

	marcs_TauN = [
		marcs_models[5].structure["lgTauR"],
		marcs_models[1].structure["lgTauR"],
		marcs_models[3].structure["lgTauR"],
		marcs_models[4].structure["lgTauR"]
	]

	marcs_zN = [
		-marcs_models[5].structure["Depth"],
		-marcs_models[1].structure["Depth"],
		-marcs_models[3].structure["Depth"],
		-marcs_models[4].structure["Depth"]
	]


	levelsGlobN = [[
		3000, 3500, 4000, 4500, 5000, 6000, 8500, 11500, 12500, 13500, 14500, 15500, 16500, 17500, 18500
	] for _ in eachindex(labelsN)]

	#=levelsGlobN[3] = [
		5000, 5500, 6000, 8500, 11500, 12500, 13500, 14500, 15500, 16500, 17500
	]
	levelsGlobN[4] = [
		4000, 4500, 5000, 6000, 8500, 11500, 12500, 13500, 14000
	]=#
	#=levelsGlobN = [
		exp.(range(log(4000), log(16000), length=20)) |> collect
		for _ in eachindex(modelsN)
	]=#
	
	pickeveryN = 1
	
	for (i, mN) in enumerate(modelsN)
		vertical_velocity_slice_with_contours(
			mN, y0N[i], labelsN[i], marcs_TN[i], marcs_zN[i], levelsGlobN[i],
			pickevery=pickeveryN, cbar_label=L"\rm T\ [K]"
		)
	end
	
	gcf()
end

# ╔═╡ 4874c7a5-6dd4-457a-876d-431e31b5875d
begin
	levelsGlobTauN = [
		[-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7] for _ in eachindex(labelsN)
	]

	levelsGlobTauN[4] = [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
		
	for (i, mN) in enumerate(modelsN)
		vertical_velocity_slice_with_contours(
			mN, y0N[i], labelsN[i], marcs_TauN[i], marcs_zN[i], levelsGlobTauN[i],
			pickevery=pickeveryN, quantity=:τ_ross, scale=log10, 
			cbar_label=L"\rm \log \tau_{ross}"
		)
	end
	
	gcf()
end

# ╔═╡ 64a75b01-b653-45d5-9c87-0f61391ca7a6
md"### (O) Actual 3D profiles compared to MARCS"

# ╔═╡ 8db25fe1-0e7d-4004-98f1-8cae9c057ac9
begin
	modelsO = [
		models["t55g45m00"], 
		models["best"], 
		models["t65g45m00"],
		models["t45g40m00"]
	]

	marcs_TO = [
		marcs_models[1].structure["T"],
		marcs_models[5].structure["T"],
		marcs_models[3].structure["T"],
		marcs_models[4].structure["T"]
	]

	marcs_trO = [
		marcs_models[1].structure["lgTauR"],
		marcs_models[5].structure["lgTauR"],
		marcs_models[3].structure["lgTauR"],
		marcs_models[4].structure["lgTauR"]
	]
	
	lwO = [1.5, 1.5, 1.5, 1.5] .*2
	
	labelsO = [
		L"\rm T_{eff}\ 5500\ K, log(g)\ 4.5\ dex", 
		L"\rm T_{eff}\ 5777\ K, log(g)\ 4.44\ dex", 
		L"\rm T_{eff}\ 6500\ K, log(g)\ 4.5\ dex",
		L"\rm T_{eff}\ 4500\ K, log(g)\ 4.0\ dex"
	]

	colorsO = ["k", "k", "k", "k"]
end

# ╔═╡ ff449d80-1610-45ef-9b0b-907d6de848f8
begin
	for (i, mO) in enumerate(modelsO)
		plt.close()
		fO, axO = plt.subplots(1, 1, figsize=(5, 6))
		#visual.basic_plot!(axO)
		
		m3disOi, m3disOi_τ = pick_snapshot(out_folder[mO], -1) 

		x = reshape(log10.(m3disOi[:τ_ross]), :) #[1:100:end]
		y = reshape(m3disOi[:T], :) #[1:100:end]
		axO.hist2d(
			x, y, bins=500, cmin=1, cmap="Blues",
			norm=matplotlib.colors.LogNorm()
		)

		axO.plot(
			profile(
				mean,
				m3disOi_τ,
				:log10τ_ross,
				:T
			)...,
			lw=lwO[i],
			label="M3DIS (D)",
			color=colorsO[i]
		)

		axO.plot(
			marcs_trO[i], 
			marcs_TO[i],
			ls=":", 
			color="k", 
			label=L"\rm MARCS",
			lw=2.5
		)
		
		if mO == models["best"]
			#=axO.plot(
				profile(
					mean,
					stagger_τ,
					:log10τ_ross,
					:T
				)...,
				lw=1.5,
				label=L"\rm Stagger",
				color="k"
			)
			axO.plot(
				co5boldmean_τ[:, 1], co5boldmean_τ[:, 2],
				lw=2,
				color="k",
				label=L"\rm CO5BOLD",
				ls="--"
			)=#
		end


		#=axO.text(
			0.97, 0.05, labelsO[i], 
			ha="right", va="bottom", 
			transform=axO.transAxes,
			color="white", backgroundcolor="k", zorder=150
		)=#
		axO.set_title(
			labelsO[i]
		)

		
		axO.set_xlabel(
			L"\rm \log \tau_{ross}"
		)
		axO.set_ylabel(
			L"\rm T \ [K]"
		)

		axO.set_xlim(-4.1, 2.1)

		xmi, xma = MUST.pyconvert.(Any, axO.get_xlim())
		ydata = y[xmi .< x .< xma]
		axO.set_ylim(minimum(ydata), maximum(ydata))
		axO.legend(
			framealpha=0, ncol=1, labelspacing=0.01, handlelength=hl, loc="upper left"
		)	

		fO.savefig("full3D_T-$(labels[mO]).png", dpi=600, bbox_inches="tight")
		fO.savefig("full3D_T-$(labels[mO]).pdf", bbox_inches="tight")
	end

	gcf()
end

# ╔═╡ Cell order:
# ╟─485a693d-4872-47d4-975d-e91a7bbc0be7
# ╠═7be90c6e-260a-11ee-12d4-93fd873d1015
# ╟─a61a48b3-cfe9-4d0b-a0f7-7008add4e1d9
# ╟─feefb5e4-c4aa-428f-87ea-82b32072fb26
# ╟─b20211ed-1d69-49f0-9af3-002affd58ff2
# ╠═c4846fea-8b7a-4338-bfa9-df9b9c0aff6f
# ╟─5e79333c-94c1-47cd-86d5-a6ee22b4a62e
# ╟─55bbe1f3-3ced-438b-9b57-777eb2fbc41d
# ╟─0c1ca595-6c09-48f8-ac15-9a499b9072c4
# ╟─22bacb57-9f7f-446b-93c6-4cf0fca5f731
# ╠═0b02d772-a51b-4aa5-9d73-23c85bac5a6a
# ╟─911bed42-e814-4a42-9be4-148334315fe2
# ╠═eb51ab2b-d6e5-45a8-9131-5383545d02f0
# ╟─1e9b3054-c949-4a94-9db5-cba24602a470
# ╟─152c58ee-6568-4fd3-ae3d-2f8e397bbd04
# ╟─ec5e9d64-e461-44ec-9753-c306cc8d29dd
# ╟─33690f11-20d9-48ee-82b9-bc0df0f9e560
# ╟─8e521645-e67e-463d-ad8b-9a9a8dd9056e
# ╟─b8bd4aa7-8e9e-40de-8188-0ddb8b233b5e
# ╟─f6edb03c-5177-4777-b3f2-fc69513a7249
# ╟─f0327365-b8ff-44ac-9aa1-d42ecd5a0a3f
# ╟─b67ca362-66e0-4e6b-ac1e-094515498379
# ╠═9a504c51-4e79-4ef4-8a4f-570f549422c0
# ╠═0a2286e2-29e4-454b-b3a7-6a33a6828760
# ╟─2e44bee0-da4c-47cf-8c7e-e8c46e696bd2
# ╟─9465b4d5-08e2-453a-b758-6d4d6744c20b
# ╠═d87a8635-7527-4584-81f6-eef0da1101d6
# ╟─a22bef7c-6409-4eb5-a17e-66ea1a42ab73
# ╟─0016626d-bb90-47eb-83f8-ec1f9d905a03
# ╠═dfb2ddf1-10ec-426d-8907-eec2b694fe5a
# ╟─4eecac61-40f9-45d4-b789-d017aef6f149
# ╟─2270321b-6581-4d45-a00f-0e728022bc56
# ╠═5b42ec94-a9fc-4ada-aa7e-5fb6bac656b1
# ╟─0f310354-ca88-468f-a760-8194322855f3
# ╟─c9a52157-81c4-4b8c-8dcd-b728a0890b00
# ╠═515e5910-4bf1-4769-84b9-cff2a08751a0
# ╠═fa715a70-036e-420a-a0fa-06cefc3ef3c0
# ╟─87cc7333-4559-467a-be0f-0f2ec1dc0611
# ╟─3babe5e7-82ee-4d5b-a871-060b6050a53f
# ╠═c3b8db0b-c9ee-4f4d-a1a4-75e3551d7675
# ╠═6c3227a7-993b-45bc-809e-be6a0907f384
# ╠═239ce3e3-5522-4d03-96d7-607c69418781
# ╠═f78e1d7f-6756-47d9-bb6c-5b6c623dc899
# ╠═8b046999-9214-468b-8b6f-90df998ecebd
# ╠═f6c1c537-9a8d-4746-abcd-17fd0a712353
# ╠═cff24989-efe0-4c50-8163-7260869895d2
# ╠═4fb7c54f-7f7a-48b7-a6a6-f80cb9c31360
# ╠═1005536b-f9c3-4526-be24-17ad9b063524
# ╟─68e9c5c9-34df-41f6-a4e6-7f5deb935d59
# ╟─dd5a0ce9-2daa-4d99-9dbb-443c113e3af5
# ╟─089e4d19-22ad-4da4-8886-6980d70b5f31
# ╠═d9272f17-8f18-457e-80d3-aa7dd61831a9
# ╟─418c34ff-1c08-4ae6-82bf-7ac2fced7244
# ╟─fd9d2f80-3b17-4c3a-b8e8-37de365c0b2a
# ╟─23728e99-134a-445f-802b-f01467a03e10
# ╠═114947c5-43c4-4d82-9031-8f27903d0415
# ╟─26ad98c3-8c6a-42f9-944b-035691fa9ad0
# ╟─c5e1dc32-de92-4deb-92a9-06818d898307
# ╟─574fe0c0-fdcc-493f-8aad-bbfdc5b55c6a
# ╟─d1163018-162d-4489-aa43-c423bcfd64f7
# ╟─01f1b55d-ba4b-4d74-8eba-12ede107cfbf
# ╟─3fb509e5-c611-45e2-b475-f43bd28ac53a
# ╠═84fc6b08-0825-41ab-b214-5fef7e2a2dc0
# ╠═35d2d2c5-6b70-4bb9-acdc-216557859282
# ╟─ae765c38-1c13-4914-85a1-41c41bff4f3f
# ╠═0ef33323-1d6a-484c-9292-3485b4871b38
# ╟─381be17e-b257-4e51-aa79-941db153f398
# ╠═e1517fd6-523f-4083-bb74-ebf666813002
# ╟─e73c5a0e-ae3b-4bbd-b098-d8e6ff215c2e
# ╟─f6f45939-55f8-42c3-b0d2-79c1cefb645c
# ╟─d77216a7-038d-491f-b7c2-ba74d627e2a2
# ╟─518657e8-1865-4ecf-bb57-73cb502b5b19
# ╟─65a010db-d74f-4326-b909-8586c4d033a7
# ╠═13b57468-6d36-4413-95c6-0743e2633924
# ╟─58bb3285-544f-4e8a-a27e-9b7e90804fe8
# ╟─dd667e0e-e554-4dfe-9493-aa38a505210f
# ╟─6de3b5df-9766-41ca-8d51-f18023a419a3
# ╟─8743067b-00ff-43c6-a6fd-800063b1d9de
# ╟─6b63bde1-9a9a-4178-ac06-630ae811964b
# ╟─792cd263-3322-4a1f-9c1c-babcdcd88fe1
# ╟─e02162d4-6345-43e2-bb3f-2a453cc1b9eb
# ╟─b9890085-fcca-472d-b192-da7abf349b45
# ╠═b14da0d1-63da-44aa-b1fd-4340ec64528e
# ╟─516fe726-9d67-44dd-bbdb-f0ec17841185
# ╟─dc93e1ef-c947-465c-ab37-b1fd251b695d
# ╟─64213279-f47a-4b46-b1ac-730945d5b8f1
# ╠═da45e33e-9013-443e-be60-e0b53d23d305
# ╠═d59e82db-89f3-46fd-8af4-6bf96fab14cf
# ╠═c9004d63-d6e5-415a-87ef-3db717a3f788
# ╠═b214c18c-58ee-474f-809d-d07582ab5647
# ╠═4935697c-974b-4de8-8182-b5c33f3628f4
# ╟─acf85ab7-3d8b-4e83-8a73-1ed00598882f
# ╟─4ae6b55e-952a-4b70-937f-c81a2f790e83
# ╠═8324d6c9-f46a-4724-be76-5eb92f753e21
# ╟─faccedff-2a3a-40b1-ae88-c42a69112d16
# ╟─1634883c-2a93-4b31-bc3a-662a894733c4
# ╟─f409f3e8-ef97-4fb9-a8e6-7a8f1e2b2d22
# ╟─759b0406-1a86-42d9-b8d7-8784a1574c02
# ╟─82b1bd73-d81a-404c-aefd-643b7008d2b7
# ╟─306dd6bc-c27d-4af4-9188-6e96981541fb
# ╟─521f69d2-c668-4e95-b422-bc95a858286c
# ╟─bec4d310-281f-4a74-bd61-cebbfd1872f5
# ╟─3c3d3319-1e82-4027-8ac4-a460acb15518
# ╟─04543c5d-d519-4ed2-9d98-050fbd1361a1
# ╟─6dc706bc-0fbd-4277-948f-df7b4306bf4d
# ╟─6a2d0807-a437-40db-9ed3-934b4a1018b9
# ╟─848900a2-3a56-48e2-b771-45ec8dc3b392
# ╟─ba8f433c-400a-4829-b534-1b88e35808f9
# ╟─1977a610-41f0-4288-9c5e-58a140acdab3
# ╟─1f458a9f-a0dd-459b-8a60-aded8d8e4617
# ╟─d602d302-5964-4525-b62e-d7eb7aec5182
# ╠═78ce2d34-9898-4140-a7fd-8fe2e807d972
# ╠═94c55529-523a-4a44-9ed8-1caaf34e4da5
# ╟─3466e7cb-17ab-4a1a-ac4f-7b3f1264584d
# ╟─cb14997c-ed1f-4be0-80ce-bba1d94765d4
# ╠═f6736f9c-f76e-40ae-bb85-63c1aef6605c
# ╟─9cbc4ead-5f1d-41b9-a7ac-dea24f4b91cb
# ╟─f33f10e7-506d-4d55-a7c2-1b573ee6e13e
# ╠═c11b167f-ba1f-4ae6-8e38-e3c011f6f036
# ╠═3f194d77-08cf-4d51-81ee-f811ed587967
# ╠═95849fd6-e8bf-48e4-a0d9-2dd86cee58bb
# ╠═a15e914f-13f6-456f-a74f-b9edd63a7ab9
# ╠═059978b8-a8b0-4f1c-a681-ea41d865a488
# ╟─f3a0a1b6-af99-46ec-8355-ab9fe7e5881d
# ╟─2587d8e9-7612-42f6-8a0c-01c72b1d32d3
# ╠═5003a0ce-7aa0-4a5e-9765-6129d7660e1b
# ╠═6ce71e68-3e46-4ee7-9777-7a04342b5e2e
# ╠═09c3002d-1839-4adf-9d47-1811f1c81bcd
# ╠═e8da9c4e-3474-41ea-a441-a22b089d20d6
# ╟─0f3f07fa-660e-4117-a46e-792150a707b9
# ╟─22bc6765-400a-470b-9f9d-ff1365d97e33
# ╟─4874c7a5-6dd4-457a-876d-431e31b5875d
# ╟─64a75b01-b653-45d5-9c87-0f61391ca7a6
# ╠═8db25fe1-0e7d-4004-98f1-8cae9c057ac9
# ╟─ff449d80-1610-45ef-9b0b-907d6de848f8
