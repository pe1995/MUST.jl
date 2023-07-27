### A Pluto.jl notebook ###
# v0.19.26

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

	plt = matplotlib.pyplot
end;

# ╔═╡ 485a693d-4872-47d4-975d-e91a7bbc0be7
md"# Paper I: Validation & The Sun (v2)"

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
	PythonPlot.matplotlib.rcParams["font.size"] = 12
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

# ╔═╡ 9465b4d5-08e2-453a-b758-6d4d6744c20b
md"### MARCS"

# ╔═╡ d87a8635-7527-4584-81f6-eef0da1101d6
marcs_model = readdlm(
	"/u/peitner/DISPATCH/examples/from_stagger/marcs_sun.txt", 
	skipstart=1
)

# ╔═╡ 3babe5e7-82ee-4d5b-a871-060b6050a53f
md"## Dispatch models"

# ╔═╡ c3b8db0b-c9ee-4f4d-a1a4-75e3551d7675
names = [
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1"
]

# ╔═╡ 6c3227a7-993b-45bc-809e-be6a0907f384
out_folder = [
	"models/sun_magg_150x300",
	"models/sun_magg_90x90",
	"models/sun_magg_90x180"
]

# ╔═╡ 239ce3e3-5522-4d03-96d7-607c69418781
in_folder = [
	MUST.@in_dispatch "input_data/grd/$(name)" for name in names
]

# ╔═╡ f78e1d7f-6756-47d9-bb6c-5b6c623dc899
eos_folder = [
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3"),
	MUST.@in_dispatch("input_data/binned/DIS_MARCS_E_v1.6.3")
]

# ╔═╡ f6c1c537-9a8d-4746-abcd-17fd0a712353
colors = [
	"tomato",
	"cyan",
	"magenta"
]

# ╔═╡ cff24989-efe0-4c50-8163-7260869895d2
ls = [
	"-",
	"-",
	"-"
]

# ╔═╡ 4fb7c54f-7f7a-48b7-a6a6-f80cb9c31360
labels = [
	"M3DIS",
	"M3DIS - lres",
	"M3DIS - ires"
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
	"best" => 1,
	"lres" => 2,
	"ires" => 3
)

# ╔═╡ 418c34ff-1c08-4ae6-82bf-7ac2fced7244
begin
	eos = [reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5")) 
				for i in eachindex(eos_folder)]
	opa = [reload(SqOpacity, joinpath(eos_folder[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder)]
end

# ╔═╡ 23728e99-134a-445f-802b-f01467a03e10
for i in eachindex(eos)
	@info "Opacity table size: $(size(opa[i].κ))"
end

# ╔═╡ 114947c5-43c4-4d82-9031-8f27903d0415
resolution(snap) = @sprintf "%.1f km" first(diff(MUST.axis(snap, :z) ./1e5))

# ╔═╡ c5e1dc32-de92-4deb-92a9-06818d898307
for (i,snap) in enumerate(snapshots)
	@info "Resolution $(out_folder[i]): $(resolution(snap))"
end

# ╔═╡ 381be17e-b257-4e51-aa79-941db153f398
md"## Figures
### (A) Internal Comparison"

# ╔═╡ e1517fd6-523f-4083-bb74-ebf666813002
begin
	modelA = models["best"]
	plt.close()
	
	fA, axA = plt.subplots(2, 1, sharex=true)
	plt.subplots_adjust(wspace=0, hspace=0)
	visual.basic_plot!.(axA)

	
	# logτ vs. T
	axA[0].plot(
		MUST.profile(
			mean, 
			pick_snapshot(out_folder[modelA], -1) |> last, :log10τ_ross, :T
		)...,
		color="k",
		label=labels[modelA]
	)

	# logτ vs. logρ
	axA[1].plot(
		MUST.profile(
			mean,
			pick_snapshot(out_folder[modelA], -1) |> last, :log10τ_ross, :log10d
		)...,
		color="k"
	)

	axA[0].legend(framealpha=0, fontsize="large")
	
	axA[0].set_ylabel(L"\rm T\ [K]", fontsize="large")
	axA[1].set_ylabel(L"\rm \log \rho\ [g \times cm^{-3}]", fontsize="large")
	axA[1].set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")

	gcf()
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
		common_tauB, ip(xA, yA) .- ip(xB, yB)
	end
	absolute_difference(f, snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(f, snapA, x, y)
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, ip(xA, yA) .- ip(xB, yB)
	end
	absolute_difference(f, xA, yA, snapB, x, y) = begin
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, ip(xA, yA) .- ip(xB, yB)
	end
	absolute_difference(f, xA, yA, xB, yB, x, y) = begin
		common_tauB, ip(xA, yA) .- ip(xB, yB)
	end
end

# ╔═╡ 6de3b5df-9766-41ca-8d51-f18023a419a3
begin
	relative_difference(snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(mean, snapA, x, y)
		xB, yB = MUST.profile(mean, snapB, x, y)
		common_tauB, (ip(xA, yA) .- ip(xB, yB)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, snapA, snapB, x, y) = begin
		xA, yA = MUST.profile(f, snapA, x, y)
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, (ip(xA, yA) .- ip(xB, yB)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, xA, yA, snapB, x, y) = begin
		xB, yB = MUST.profile(f, snapB, x, y)
		common_tauB, (ip(xA, yA) .- ip(xB, yB)) ./ ip(xA, yA) * 100.0
	end
	relative_difference(f, xA, yA, xB, yB, x, y) = begin
		common_tauB, (ip(xA, yA) .- ip(xB, yB)) ./ ip(xA, yA) * 100.0
	end
end

# ╔═╡ 8743067b-00ff-43c6-a6fd-800063b1d9de
begin
	modelB = models["best"]
	
	# Figure
	plt.close()
	
	fB, axB = plt.subplots(3, 2, sharex=true, figsize=(12, 10))
	plt.subplots_adjust(wspace=0.25, hspace=0.0)
	for j in 0:1
		for i in 0:2
			visual.basic_plot!(axB[i, j])
		end
	end
	
	

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

	

	lwD = 3
	# Absolute comparison
	begin
		## logτ vs. T
		axB[0, 0].plot(
			MUST.profile(mean, stagger_τ, :log10τ_ross, :T)...,
			color="k",
			label="Stagger",
			lw=lwD
		)
		#=axB[0, 0].plot(
			MUST.profile(mean, muram_τ, :log10τ_ross, :T)...,
			color="k",
			label="MURaM",
			ls="--",
			lw=lwD
		)=#
		axB[0, 0].plot(
			marcs_model[:, 2], marcs_model[:, 5],
			color="k",
			ls=":",
			label="MARCS",
			lw=lwD
		)
		axB[0, 0].plot(
			m3disB_x, m3disB_y,
			color="cornflowerblue",
			label=labels[modelB],
			lw=lwD
		)
		
		## logτ vs. logρ
		axB[1, 0].plot(
			MUST.profile(mean, stagger_τ, :log10τ_ross, :log10d)...,
			color="k",
			label="Stagger",
			lw=lwD
		)
		#=axB[1, 0].plot(
			MUST.profile(mean, muram_τ, :log10τ_ross, :log10d)...,
			color="k",
			label="MURaM",
			ls="--",
			lw=lwD
		)=#
		axB[1, 0].plot(
			marcs_model[:, 2], log10.(marcs_model[:, 11]),
			color="k",
			ls=":",
			label="MARCS",
			lw=lwD
		)
		axB[1, 0].plot(
			MUST.profile(mean, m3disB, :log10τ_ross, :log10d)...,
			color="cornflowerblue",
			lw=lwD
		)
		
		## logτ vs. rms velocity
		axB[2, 0].plot(
			MUST.profile(rms5, stagger_τ, :log10τ_ross, :uz)...,
			color="k",
			label="Stagger",
			lw=lwD
		)
		#=axB[2, 0].plot(
			MUST.profile(rms5, muram_τ, :log10τ_ross, :uz)...,
			color="k",
			label="MURaM",
			ls="--",
			lw=lwD
		)=#
		axB[2, 0].plot(
			MUST.profile(rms5, m3disB, :log10τ_ross, :uz)...,
			color="cornflowerblue",
			lw=lwD
		)

		
		axB[0, 0].set_ylim(3800, 11000)
		axB[1, 0].set_ylim(-9, -5.85)
		axB[2, 0].set_ylim(0.5, 3.75)
		
		
		axB[0, 0].legend(framealpha=0, fontsize="large")
		
		axB[0, 0].set_ylabel(L"\rm T\ [K]", fontsize="large")
		axB[1, 0].set_ylabel(L"\rm \log \rho\ [g \times cm^{-3}]", fontsize="large")
		axB[2, 0].set_ylabel(
			L"\rm U_{z}\ [km\ \times s^{-1}]", 
			fontsize="large"
		)
	end


	# relative comparison
	begin
		## logτ vs. T
		axB[0, 1].axhline(0.0, color="0.5", alpha=0.2)
		axB[0, 1].plot(
				relative_difference(
				mean,
				m3disB_x, 
				m3disB_y,
				stagger_τ,
				:log10τ_ross, :T
			)...,
			color="k",
			label="M3DIS - Stagger",
			lw=lwD
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
		axB[0, 1].plot(
			common_tauB, 
			(ip(m3disB_x, m3disB_y) .- 
				ip(marcs_model[:, 2], marcs_model[:, 5])) ./ 
				ip(m3disB_x, m3disB_y) *100.0,
			color="k",
			label="M3DIS - MARCS",
			ls=":",
			lw=lwD
		)
		
		## logτ vs. logρ
		axB[1, 1].axhline(0.0, color="0.5", alpha=0.2)
		axB[1, 1].plot(
			absolute_difference(
				mean,
				m3disB_xr, m3disB_yr,
				stagger_τ,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label="M3DIS - Stagger",
			lw=lwD
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
		axB[1, 1].plot(
			common_tauB, 
			ip(m3disB_xr, m3disB_yr) .- 
				ip(marcs_model[:, 2], log10.(marcs_model[:, 11])),
			color="k",
			label="M3DIS - MARCS",
			ls=":",
			lw=lwD
		)

		## logτ vs. rms
		axB[2, 1].axhline(0.0, color="0.5", alpha=0.2)
		axB[2, 1].plot(
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
	

		axB[0, 1].set_xlim(-4, 2)
		axB[0, 1].set_ylim(-6.25, 6.25)
		axB[1, 1].set_ylim(-0.17, 0.17)
		axB[2, 1].set_ylim(-0.47, 0.47)
	
		
		axB[0, 1].legend(framealpha=0, fontsize="large")
		axB[0, 1].set_ylabel(L"\rm \Delta\ T\ /\ T\ [\%]", fontsize="large")
		axB[1, 1].set_ylabel(
			L"\rm \Delta\ \log \rho\ [g \times cm^{-3}]", fontsize="large"
		)
		axB[2, 1].set_ylabel(
			L"\rm \Delta\ U_{z}\ [km\ \times s^{-1}]", 
			fontsize="large"
		)
	end

	axB[2, 0].set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")
	axB[2, 1].set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")
		

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
	
	lresC = pick_snapshot(out_folder[lresmodelC], :recent) |> last
	iresC = pick_snapshot(out_folder[iresmodelC], :recent) |> last
	bestC = pick_snapshot(out_folder[bestmodelC], :recent) |> last
	
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
end

# ╔═╡ e02162d4-6345-43e2-bb3f-2a453cc1b9eb
begin
	# Figure
	plt.close()
	
	fC, axC = plt.subplots(3, 2, sharex=true, figsize=(12, 10))
	plt.subplots_adjust(wspace=0.25, hspace=0.0)
	for j in 0:1
		for i in 0:2
			visual.basic_plot!(axC[i, j])
		end
	end
	

	lwC = 2.5
	# Absolute comparison
	begin
		## logτ vs. T
		axC[0, 0].plot(
			bestC_x, bestC_y,
			color="cornflowerblue",
			label=@sprintf("%i km", 2.3/300 *1000),
			lw=lwC*2
		)
		axC[0, 0].plot(
			iresC_x, iresC_y,
			color="k",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[0, 0].plot(
			lresC_x, lresC_y,
			color="k",
			label=@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)
		
		
		## logτ vs. logρ
		axC[1, 0].plot(
			bestC_xr, bestC_yr,
			color="cornflowerblue",
			label=@sprintf("%i km", 2.3/300 *1000),
			lw=lwC
		)
		axC[1, 0].plot(
			iresC_xr, iresC_yr,
			color="k",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[1, 0].plot(
			lresC_xr, lresC_yr,
			color="k",
			label=@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)

		
		## logτ vs. uz
		axC[2, 0].plot(
			bestC_xu, bestC_yu,
			color="cornflowerblue",
			label=@sprintf("%i km", 2.3/300 *1000),
			lw=lwC
		)
		axC[2, 0].plot(
			iresC_xu, iresC_yu,
			color="k",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[2, 0].plot(
			lresC_xu, lresC_yu,
			color="k",
			label=@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)

		
		axC[0, 0].set_ylim(3800, 11000)
		axC[1, 0].set_ylim(-9, -5.85)
		axC[2, 0].set_ylim(0.5, 3.75)
		
		axC[0, 0].legend(framealpha=0, fontsize="large")
		
		axC[0, 0].set_ylabel(L"\rm T\ [K]", fontsize="large")
		axC[1, 0].set_ylabel(L"\rm \log \rho\ [g \times cm^{-3}]", fontsize="large")
		axC[2, 0].set_ylabel(
			L"\rm U_{z}\ [km\ \times s^{-1}]", 
			fontsize="large"
		)
	end


	# relative comparison
	begin
		## logτ vs. T
		axC[0, 1].axhline(0.0, color="0.5", alpha=0.2)
			axC[0, 1].plot(
			relative_difference(
				mean,
				bestC_x, bestC_y,
				iresC_x, iresC_y,
				:log10τ_ross, :T
			)...,
			color="k",
			label=@sprintf("best - %i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[0, 1].plot(
			relative_difference(
				mean,
				bestC_x, bestC_y,
				lresC_x, lresC_y,
				:log10τ_ross, :T
			)...,
			color="k",
			label=@sprintf("best - %i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)
		
		## logτ vs. logρ
		axC[1, 1].axhline(0.0, color="0.5", alpha=0.2)
		axC[1, 1].plot(
			absolute_difference(
				mean,
				bestC_xr, bestC_yr,
				iresC_xr, iresC_yr,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label=@sprintf("%i km", 2.3/180 *1000),
			lw=lwC,
			ls=":"
		)
		axC[1, 1].plot(
			absolute_difference(
				mean,
				bestC_xr, bestC_yr,
				lresC_xr, lresC_yr,
				:log10τ_ross, :log10d
			)...,
			color="k",
			label=@sprintf("%i km", 2.3/90 *1000),
			lw=lwC,
			ls="--"
		)

		## logτ vs. rms
		axC[2, 1].axhline(0.0, color="0.5", alpha=0.2)
		axC[2, 1].plot(
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
		axC[2, 1].axhline(0.0, color="0.5", alpha=0.2)
		axC[2, 1].plot(
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
			

		axC[0, 1].set_xlim(-4, 2)
		axC[0, 1].set_ylim(-6.25, 6.25)
		axC[1, 1].set_ylim(-0.17, 0.17)
		axC[2, 1].set_ylim(-0.47, 0.47)
	
		
		axC[0, 1].legend(framealpha=0, fontsize="large")
		axC[0, 1].set_ylabel(L"\rm \Delta\ T\ / \ T\ [\%]", fontsize="large")
		axC[1, 1].set_ylabel(
			L"\rm \Delta\ \log \rho\ [g \times cm^{-3}]", fontsize="large"
		)
		axC[2, 1].set_ylabel(
			L"\rm \Delta\ U_{z}\ [km\ \times s^{-1}]", 
			fontsize="large"
		)
	end

	axC[2, 0].set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")
	axC[2, 1].set_xlabel(L"\rm \log \tau_{ross}", fontsize="large")
		

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
	fD, axD = plt.subplots(1, 1, figsize=(6,7))
	visual.basic_plot!(axD)
	
	im = axD.scatter(
		log10.(λ), -log10.(fopaD.κ_ross), 
		c=binassignment,
		s=1,
		cmap=cD,
		norm=normD,
	)
	cbarD = fD.colorbar(
		im, 
		ax=axD, 
		ticks=range(1, length(opaD.λ), step=1) |> collect
	)


	axD.set_ylabel(
		L"\rm - \log \tau_{ross}\ (\log \tau_{\lambda}=1)",
		fontsize="large"
	)
	axD.set_xlabel(
		L"\rm \log \lambda\ [\AA]",
		fontsize="large"
	)
	cbarD.set_label(
		L"\rm opacity\ bins", 
		fontsize="large"
	)
	
	fD.savefig("formation_opacities.pdf", bbox_inches="tight")
	fD.savefig("formation_opacities.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ 516fe726-9d67-44dd-bbdb-f0ec17841185
md"### (E) Effective Temperature"

# ╔═╡ dc93e1ef-c947-465c-ab37-b1fd251b695d
md"### (F) Radiative Heating (Binning)"

# ╔═╡ 64213279-f47a-4b46-b1ac-730945d5b8f1
md"### (G) 2D surface"

# ╔═╡ 9c01846f-d4cc-43a8-9c0f-93f7e08bcefb
uz_surface_optical(snap) = begin
	isurf = MUST.closest(log10.(MUST.axis(snap, :τ_ross, 3)), 0)
	snap[:uz][:, :, isurf] #./1e5
end

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
	
	fG, axG = plt.subplots(1, 1, figsize=(6, 7))

	# limits for color bar
	uzG  = uz_surface_optical(m3disG) ./1e5
	vmin = minimum([minimum(u) for u in uzG])
	vmax = maximum([maximum(u) for u in uzG])

	# plot the surface in uz
	imG = axG.imshow(
		uzG,
		origin="lower",
		vmin=vmin, vmax=vmax,
		cmap="hot",
		extent=extent(m3disG)
	)

	# setup colorbar
	cG = fG.colorbar(imG, ax=axG, fraction=0.046, pad=0.04)
	cG.set_label(L"\rm U_z\ [km \times s^{-1}]", fontsize="large")

	# Other labels
	axG.set_xlabel("x [Mm]", fontsize="large")
	axG.set_ylabel("y [Mm]", fontsize="large")

	#fG.savefig("vertical_velocity_surface.pdf", bbox_inches="tight")
	fG.savefig("vertical_velocity_surface.png", bbox_inches="tight", dpi=600)
	
	gcf()
end

# ╔═╡ faccedff-2a3a-40b1-ae88-c42a69112d16
md"### (H) 3D cube with velocities"

# ╔═╡ 1634883c-2a93-4b31-bc3a-662a894733c4
begin
	plt.close()

	modelH = models["best"]
	m3disH = pick_snapshot(out_folder[modelH], -1) |> first
	
	fH, axH = visual.cube_with_velocities(m3disH, vmax_3d=17000, cmap="hot")


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
	snap_m3disI_τ = MUST.height_scale(snap_m3disI, :τ_ross)
	

	# We look at the surface and average profile of this snapshot for different 
	# resamplings
	samplings = reshape(collect(Iterators.product((5, 10, 20, 80), (299,))), :)

	fI, axI = plt.subplots(4, 2, figsize=(14,21))
	for j in 0:1
		for i in 0:3
			visual.basic_plot!(axI[i, j])
		end
	end
	plt.subplots_adjust(wspace=0.2, hspace=0.1)

	vminI = 5500
	vmaxI = 6900
	i0I_full = argmin(abs.(log.(MUST.axis(m3disI_τ, :τ_ross, 3))))
	TI_full = m3disI[:T][:, :, i0I_full]
	
	for (i, resolution) in enumerate(samplings)
		snap = gresample(
			m3disI, 
			nx=resolution[1], 
			ny=resolution[1],
			nz=resolution[2]
		)

		snap_r = MUST.height_scale(snap, :τ_ross)
		
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

		xI_i , yI_i = profile(MUST.mean, snap_r, :log10τ_ross, :T)		
		axI[i-1, 1].plot(
			xI_i, yI_i,
			color="k", 
			marker="o", 
			markersize=5,
			ls="",
			label="$(resolution[1])x$(resolution[1])x$(resolution[2])"
		)
		
		xI_i , yI_i = profile(MUST.mean, snap_m3disI_τ, :log10τ_ross, :T)		
		axI[i-1, 1].plot(
			xI_i, yI_i,
			color="tomato", 
			marker="",
			ls="-",
			lw=3,
			label="original"
		)

		axI[i-1, 1].set_ylim(4000, 10000)
		axI[i-1, 1].set_xlim(-4.0, 1.75)
		axI[i-1, 1].legend(framealpha=0, fontsize="x-large")
		axI[i-1, 1].set_ylabel(L"\rm T\ [K]", fontsize="x-large")
	end

	for (i, resolution) in enumerate(samplings)
		axI[i-1, 0].set_ylabel(L"\rm Y\ [Mm]", fontsize="x-large")	
	end
	axI[3, 0].set_xlabel(L"\rm X\ [Mm]", fontsize="x-large")
	axI[3, 1].set_xlabel(L"\rm \log \tau_{ross}", fontsize="x-large")	
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─485a693d-4872-47d4-975d-e91a7bbc0be7
# ╠═7be90c6e-260a-11ee-12d4-93fd873d1015
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
# ╟─f6edb03c-5177-4777-b3f2-fc69513a7249
# ╟─f0327365-b8ff-44ac-9aa1-d42ecd5a0a3f
# ╟─b67ca362-66e0-4e6b-ac1e-094515498379
# ╠═9a504c51-4e79-4ef4-8a4f-570f549422c0
# ╠═0a2286e2-29e4-454b-b3a7-6a33a6828760
# ╟─9465b4d5-08e2-453a-b758-6d4d6744c20b
# ╟─d87a8635-7527-4584-81f6-eef0da1101d6
# ╟─3babe5e7-82ee-4d5b-a871-060b6050a53f
# ╠═c3b8db0b-c9ee-4f4d-a1a4-75e3551d7675
# ╠═6c3227a7-993b-45bc-809e-be6a0907f384
# ╠═239ce3e3-5522-4d03-96d7-607c69418781
# ╠═f78e1d7f-6756-47d9-bb6c-5b6c623dc899
# ╠═f6c1c537-9a8d-4746-abcd-17fd0a712353
# ╠═cff24989-efe0-4c50-8163-7260869895d2
# ╠═4fb7c54f-7f7a-48b7-a6a6-f80cb9c31360
# ╟─68e9c5c9-34df-41f6-a4e6-7f5deb935d59
# ╟─dd5a0ce9-2daa-4d99-9dbb-443c113e3af5
# ╟─089e4d19-22ad-4da4-8886-6980d70b5f31
# ╠═d9272f17-8f18-457e-80d3-aa7dd61831a9
# ╟─418c34ff-1c08-4ae6-82bf-7ac2fced7244
# ╟─23728e99-134a-445f-802b-f01467a03e10
# ╟─114947c5-43c4-4d82-9031-8f27903d0415
# ╟─c5e1dc32-de92-4deb-92a9-06818d898307
# ╟─381be17e-b257-4e51-aa79-941db153f398
# ╟─e1517fd6-523f-4083-bb74-ebf666813002
# ╟─58bb3285-544f-4e8a-a27e-9b7e90804fe8
# ╟─dd667e0e-e554-4dfe-9493-aa38a505210f
# ╟─6de3b5df-9766-41ca-8d51-f18023a419a3
# ╟─8743067b-00ff-43c6-a6fd-800063b1d9de
# ╟─6b63bde1-9a9a-4178-ac06-630ae811964b
# ╟─792cd263-3322-4a1f-9c1c-babcdcd88fe1
# ╟─e02162d4-6345-43e2-bb3f-2a453cc1b9eb
# ╟─b9890085-fcca-472d-b192-da7abf349b45
# ╟─b14da0d1-63da-44aa-b1fd-4340ec64528e
# ╟─516fe726-9d67-44dd-bbdb-f0ec17841185
# ╟─dc93e1ef-c947-465c-ab37-b1fd251b695d
# ╟─64213279-f47a-4b46-b1ac-730945d5b8f1
# ╟─9c01846f-d4cc-43a8-9c0f-93f7e08bcefb
# ╟─acf85ab7-3d8b-4e83-8a73-1ed00598882f
# ╟─4ae6b55e-952a-4b70-937f-c81a2f790e83
# ╟─faccedff-2a3a-40b1-ae88-c42a69112d16
# ╟─1634883c-2a93-4b31-bc3a-662a894733c4
# ╟─f409f3e8-ef97-4fb9-a8e6-7a8f1e2b2d22
# ╟─759b0406-1a86-42d9-b8d7-8784a1574c02
