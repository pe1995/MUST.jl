### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 32b14b00-8ae8-11ef-2aa2-937a4b8d3a0c
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("..")
	using TSO
	using MUST
	using PythonPlot
	using LaTeXStrings 
	using PlutoUI
	using KernelDensity
end

# ╔═╡ bbb5fb60-6669-440c-84d6-5e0edbd90e3d
begin
	@import_dispatch "../../../../dispatch2"
	@import_m3dis "../../../../Multi3D"
end

# ╔═╡ 024bd32f-a8c3-43eb-bba3-05d590f9632f
auxfuncs = MUST.pyimport("m3dis.auxfuncs")

# ╔═╡ 9967494d-13ad-4997-95c9-341ca3402ded
begin
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Hoppe2024.mplstyle"))
	matplotlib.rcParams.update(Dict(
		"font.size"=> 14,
		"axes.linewidth"=> 2,
		"xtick.major.width"=> 2,
		"xtick.minor.width"=> 1.8,
		"ytick.major.width"=> 2,
		"ytick.minor.width"=> 1.8,
		"xtick.major.size"=> 6.5,
		"ytick.major.size"=> 6.5,
		"xtick.minor.size"=> 4,
		"ytick.minor.size"=> 4,
		
	))
	TableOfContents()
end

# ╔═╡ bd2ecf16-79e4-48e9-947a-0a323a1d3dd7
md"# Models"

# ╔═╡ 8b84192a-765a-43a2-a7af-006e74be01fb
md"To cover the relevant parameterspace in the Yoon-Beers diagram, we pick specific models that represent Group II and III CEMP-no stars. Those models will be inspected with respect to their atmospheric structure and the direct influence on the visible spectrum."

# ╔═╡ 1ae7a204-5a76-4705-94f2-ab248c3ecc20


# ╔═╡ af7646aa-a4f8-4eb2-8121-1e99559a16f0
md"## 3D Models"

# ╔═╡ b5b1d8c0-b67a-45d2-b21e-e0e04555a9fd
dispatchPath(runname, dir="CEMP_models2") = MUST.@in_dispatch(joinpath(dir, runname))

# ╔═╡ c7a6739a-8ae7-45c0-ab4f-e800a9ab1d2e
runs3D = Dict(
	"t50.00g25.00m-5.000 CEMP"    => dispatchPath("CEMP_t50.00g25.00m-5.000_v1.0"),
	"t57.50g45.00m-5.000 CEMP"    => dispatchPath("CEMP_t57.50g45.00m-5.000_v1.0"),
	"t52.50g30.00m-5.000 CEMP"    => dispatchPath("CEMP_t52.50g30.00m-5.000_v1.0"),
	"t50.00g25.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-5.000_v1.0"),
	"t57.50g45.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-5.000_v1.0"),
	"t52.50g30.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-5.000_v1.0"),
	"t50.00g25.00m-4.000 CEMP"    => dispatchPath("CEMP_t50.00g25.00m-4.000_v1.0"),
	"t57.50g45.00m-4.000 CEMP"    => dispatchPath("CEMP_t57.50g45.00m-4.000_v1.0"),
	"t52.50g30.00m-4.000 CEMP"    => dispatchPath("CEMP_t52.50g30.00m-4.000_v1.0"),
	"t50.00g25.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-4.000_v1.0"),
	"t57.50g45.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-4.000_v1.0"),
	"t52.50g30.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-4.000_v1.0"),
	"t50.00g25.00m-3.000 CEMP"    => dispatchPath("CEMP_t50.00g25.00m-3.000_v1.0"),
	"t57.50g45.00m-3.000 CEMP"    => dispatchPath("CEMP_t57.50g45.00m-3.000_v1.0"),
	"t52.50g30.00m-3.000 CEMP"    => dispatchPath("CEMP_t52.50g30.00m-3.000_v1.0"),
	"t50.00g25.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-3.000_v1.0"),
	"t57.50g45.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-3.000_v1.0"),
	"t52.50g30.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-3.000_v1.0"),
	"t57.50g45.00m-2.000 CEMP"    => dispatchPath("CEMP_t57.50g45.00m-2.000_v1.0"),
	"t52.50g30.00m-2.000 CEMP"    => dispatchPath("CEMP_t52.50g30.00m-2.000_v1.0"),
	"t57.50g45.00m-2.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-2.000_v1.0"),
	"t52.50g30.00m-2.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-2.000_v1.0"),
)

# ╔═╡ a343cf41-e425-4a59-81bd-87f8765251c0
snapshots3D = Dict(
	k=>-1 for k in keys(runs3D)
)

# ╔═╡ cd8b84b0-6891-4d3e-ba04-22fd4890e09d
begin
	models3D = Dict()
	for k in keys(runs3D)
		models3D[k] = pick_snapshot(runs3D[k], snapshots3D[k])
	end
end

# ╔═╡ d2884a6e-b1e5-4234-ae78-05f682ac6908


# ╔═╡ 863499a2-1e33-4ba4-8dad-935a9a38a4ec
md"## 1D Models"

# ╔═╡ 9d8761fd-4880-4ff3-855f-f3713a55a6e5
runs1D = Dict(
	"MARCS t50.00g25.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-5.000_v1.0/s5000_g+2.5_m1.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t57.50g45.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-5.000_v1.0/p5750_g+4.5_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t52.50g30.00m-5.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-5.000_v1.0/p5250_g+3.0_m0.0_t01_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t50.00g25.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-4.000_v1.0/s5000_g+2.5_m1.0_t05_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t57.50g45.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-4.000_v1.0/p5750_g+4.5_m0.0_t01_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t52.50g30.00m-4.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-4.000_v1.0/p5250_g+3.0_m0.0_t01_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t50.00g25.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t50.00g25.00m-3.000_v1.0/s5000_g+2.5_m1.0_t02_st_z-3.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t57.50g45.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-3.000_v1.0/p5750_g+4.5_m0.0_t01_st_z-3.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t52.50g30.00m-3.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-3.000_v1.0/p5250_g+3.0_m0.0_t01_st_z-3.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t57.50g45.00m-2.000 no-CEMP" => dispatchPath("ScaledSolar_t57.50g45.00m-2.000_v1.0/p5750_g+4.5_m0.0_t01_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
	"MARCS t52.50g30.00m-2.000 no-CEMP" => dispatchPath("ScaledSolar_t52.50g30.00m-2.000_v1.0/p5250_g+3.0_m0.0_t01_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
)

# ╔═╡ b26dbd2d-85ce-41c2-a32a-c9a488f9469f
find_eos(dir) = begin
	all_dir = MUST.glob("*/", dir)
	all_dir[findfirst([isfile(joinpath(d, "eos.hdf5")) for d in all_dir])]
end

# ╔═╡ 5ac9c28b-be99-4c27-9bdd-3aff553f9bc3
begin
	read_simple(name, dir) = begin
		eos_path = find_eos(dir)
		eos = reload(SqEoS, joinpath(eos_path, "eos_T.hdf5"))
		
		# read as actual MARCS model
		marcsB = MUST.marcsBox(joinpath(dir, name))
		
		m = TSO.flip!(Model1D(
			z=marcsB.z[1, 1, :], 
			lnρ=log.(marcsB[:d][1, 1, :]),
			lnT=log.(marcsB[:T][1, 1, :])
		), depth=true)
		
		m_tau = @optical m eos
	end
	read_simple(key) = read_simple(basename(runs1D[key]), dirname(runs1D[key]))
end

# ╔═╡ e4eef27a-65c4-4ff5-a3e7-b0119a7e33c5
models1D = Dict(k=>read_simple(k) for (k,v) in runs1D)

# ╔═╡ 67c8479d-7c41-4baa-b2ec-13da36a48900
models1D["MARCS t57.50g45.00m-5.000 no-CEMP"]

# ╔═╡ 4bc7b154-2fdd-4057-bf40-e44c7851bc7d


# ╔═╡ 0a0779a4-7d87-40d2-9dfc-454b366c2dff
md"# Spectra"

# ╔═╡ 5beedaf9-b12a-46db-a576-7913eb5f2caf
begin
	spectraPath(rundir, snapshot, extension, dimension="m3dis_$(snapshot)") = begin
		if snapshot < 0
			allsnaps = MUST.glob("spectra_m3dis_*_lam_4297-4303$(extension)", rundir)
			@show rundir snapshot extension
			allsnaps[end + snapshot + 1]
		else
			joinpath(
				rundir, "spectra_$(dimension)_lam_4297-4303$(extension)"
			)
		end
	end
	spectraPath(k, extension) = spectraPath(runs3D[k], snapshots3D[k], extension)
	spectraPathSnap(k, snapshot, extension) = spectraPath(runs3D[k], snapshot, extension)
	spectraPath1D(k, extension) = spectraPath(dirname(runs1D[k]), 0, extension, basename(runs1D[k]))
end

# ╔═╡ a0cdfa21-bc3f-4462-aa4d-a84cc795bece
md"## 3D spectra"

# ╔═╡ 6613157c-d199-4825-b779-a477ac4f3f30
spectrapath3D = Dict(
	"t50.00g25.00m-5.000 CEMP - C3.0"    => spectraPath("t50.00g25.00m-5.000 CEMP","C_3.0"),
	"t50.00g25.00m-5.000 no-CEMP - C3.0" => spectraPath("t50.00g25.00m-5.000 no-CEMP","C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0"    => spectraPath("t57.50g45.00m-5.000 CEMP","C_3.0"),
	"t57.50g45.00m-5.000 no-CEMP - C3.0" => spectraPath("t57.50g45.00m-5.000 no-CEMP","C_3.0"),
	"t52.50g30.00m-5.000 CEMP - C3.0"    => spectraPath("t52.50g30.00m-5.000 CEMP","C_3.0"),
	"t52.50g30.00m-5.000 no-CEMP - C3.0" => spectraPath("t52.50g30.00m-5.000 no-CEMP","C_3.0"),
	# --> 3D effect fitting
	"t57.50g45.00m-5.000 CEMP - C3.2" => spectraPath("t57.50g45.00m-5.000 CEMP","C_3.2"), 	
	"t57.50g45.00m-5.000 CEMP - C2.8" => spectraPath("t57.50g45.00m-5.000 CEMP","C_2.8"),
	"t57.50g45.00m-5.000 CEMP - C2.6" => spectraPath("t57.50g45.00m-5.000 CEMP","C_2.6"),
	"t57.50g45.00m-5.000 CEMP - C2.4" => spectraPath("t57.50g45.00m-5.000 CEMP","C_2.4"),
	"t57.50g45.00m-5.000 CEMP - C2.2" => spectraPath("t57.50g45.00m-5.000 CEMP","C_2.2"),
	"t57.50g45.00m-5.000 CEMP - C2.0" => spectraPath("t57.50g45.00m-5.000 CEMP","C_2.0"),
	"t52.50g30.00m-5.000 CEMP - C3.4" => spectraPath("t52.50g30.00m-5.000 CEMP","C_3.4"),
	"t52.50g30.00m-5.000 CEMP - C3.2" => spectraPath("t52.50g30.00m-5.000 CEMP","C_3.2"),
	"t52.50g30.00m-5.000 CEMP - C2.8" => spectraPath("t52.50g30.00m-5.000 CEMP","C_2.8"),
	"t52.50g30.00m-5.000 CEMP - C2.6" => spectraPath("t52.50g30.00m-5.000 CEMP","C_2.6"),
	"t52.50g30.00m-5.000 CEMP - C2.4" => spectraPath("t52.50g30.00m-5.000 CEMP","C_2.4"),
	"t52.50g30.00m-5.000 CEMP - C2.2" => spectraPath("t52.50g30.00m-5.000 CEMP","C_2.2"),
	"t52.50g30.00m-5.000 CEMP - C2.0" => spectraPath("t52.50g30.00m-5.000 CEMP","C_2.0"),
	"t50.00g25.00m-5.000 CEMP - C3.4" => spectraPath("t50.00g25.00m-5.000 CEMP","C_3.4"),
	"t50.00g25.00m-5.000 CEMP - C3.2" => spectraPath("t50.00g25.00m-5.000 CEMP","C_3.2"),
	"t50.00g25.00m-5.000 CEMP - C2.8" => spectraPath("t50.00g25.00m-5.000 CEMP","C_2.8"),
	"t50.00g25.00m-5.000 CEMP - C2.6" => spectraPath("t50.00g25.00m-5.000 CEMP","C_2.6"),
	"t50.00g25.00m-5.000 CEMP - C2.4" => spectraPath("t50.00g25.00m-5.000 CEMP","C_2.4"),
	"t50.00g25.00m-5.000 CEMP - C2.2" => spectraPath("t50.00g25.00m-5.000 CEMP","C_2.2"),
	"t50.00g25.00m-5.000 CEMP - C2.0" => spectraPath("t50.00g25.00m-5.000 CEMP","C_2.0"),
	
	# FeH = -4
	"t50.00g25.00m-4.000 CEMP - C2.0"    => spectraPath("t50.00g25.00m-4.000 CEMP","C_2.0"),
	"t50.00g25.00m-4.000 no-CEMP - C2.0" => spectraPath("t50.00g25.00m-4.000 no-CEMP","C_2.0"),
	"t57.50g45.00m-4.000 CEMP - C2.0"    => spectraPath("t57.50g45.00m-4.000 CEMP","C_2.0"),
	"t57.50g45.00m-4.000 no-CEMP - C2.0" => spectraPath("t57.50g45.00m-4.000 no-CEMP","C_2.0"),
	"t52.50g30.00m-4.000 CEMP - C2.0"    => spectraPath("t52.50g30.00m-4.000 CEMP","C_2.0"),
	"t52.50g30.00m-4.000 no-CEMP - C2.0" => spectraPath("t52.50g30.00m-4.000 no-CEMP","C_2.0"),
	# --> 3D effect fitting
	"t57.50g45.00m-4.000 CEMP - C2.2" => spectraPath("t57.50g45.00m-4.000 CEMP","C_2.2"),
	"t57.50g45.00m-4.000 CEMP - C1.8" => spectraPath("t57.50g45.00m-4.000 CEMP","C_1.8"),
	"t57.50g45.00m-4.000 CEMP - C1.6" => spectraPath("t57.50g45.00m-4.000 CEMP","C_1.6"),
	"t57.50g45.00m-4.000 CEMP - C1.4" => spectraPath("t57.50g45.00m-4.000 CEMP","C_1.4"),
	"t57.50g45.00m-4.000 CEMP - C1.2" => spectraPath("t57.50g45.00m-4.000 CEMP","C_1.2"),
	"t57.50g45.00m-4.000 CEMP - C1.0" => spectraPath("t57.50g45.00m-4.000 CEMP","C_1.0"),
	"t52.50g30.00m-4.000 CEMP - C2.4" => spectraPath("t52.50g30.00m-4.000 CEMP","C_2.4"),
	"t52.50g30.00m-4.000 CEMP - C2.2" => spectraPath("t52.50g30.00m-4.000 CEMP","C_2.2"),
	"t52.50g30.00m-4.000 CEMP - C1.8" => spectraPath("t52.50g30.00m-4.000 CEMP","C_1.8"),
	"t52.50g30.00m-4.000 CEMP - C1.6" => spectraPath("t52.50g30.00m-4.000 CEMP","C_1.6"),
	"t52.50g30.00m-4.000 CEMP - C1.4" => spectraPath("t52.50g30.00m-4.000 CEMP","C_1.4"),
	"t52.50g30.00m-4.000 CEMP - C1.2" => spectraPath("t52.50g30.00m-4.000 CEMP","C_1.2"),
	#"t52.50g30.00m-4.000 CEMP - C1.0" => 
	#spectraPath("t52.50g30.00m-4.000 CEMP","C_1.0"),
	"t50.00g25.00m-4.000 CEMP - C2.4" => spectraPath("t50.00g25.00m-4.000 CEMP","C_2.4"),
	"t50.00g25.00m-4.000 CEMP - C2.2" => spectraPath("t50.00g25.00m-4.000 CEMP","C_2.2"),
	"t50.00g25.00m-4.000 CEMP - C1.8" => spectraPath("t50.00g25.00m-4.000 CEMP","C_1.8"),
	"t50.00g25.00m-4.000 CEMP - C1.6" => spectraPath("t50.00g25.00m-4.000 CEMP","C_1.6"),
	"t50.00g25.00m-4.000 CEMP - C1.4" => spectraPath("t50.00g25.00m-4.000 CEMP","C_1.4"),
	"t50.00g25.00m-4.000 CEMP - C1.2" => spectraPath("t50.00g25.00m-4.000 CEMP","C_1.2"),
	

	# FeH = -3
	"t50.00g25.00m-3.000 CEMP - C1.0"    => spectraPath("t50.00g25.00m-3.000 CEMP","C_1.0"),
	"t50.00g25.00m-3.000 no-CEMP - C1.0" => spectraPath("t50.00g25.00m-3.000 no-CEMP","C_1.0"),
	"t57.50g45.00m-3.000 CEMP - C1.0"    => spectraPath("t57.50g45.00m-3.000 CEMP","C_1.0"),
	"t57.50g45.00m-3.000 no-CEMP - C1.0" => spectraPath("t57.50g45.00m-3.000 no-CEMP","C_1.0"),
	"t52.50g30.00m-3.000 CEMP - C1.0"    => spectraPath("t52.50g30.00m-3.000 CEMP","C_1.0"),
	"t52.50g30.00m-3.000 no-CEMP - C1.0" => spectraPath("t52.50g30.00m-3.000 no-CEMP","C_1.0"),
	# --> 3D effect fitting
	"t57.50g45.00m-3.000 CEMP - C1.2" => spectraPath("t57.50g45.00m-3.000 CEMP","C_1.2"),
	"t57.50g45.00m-3.000 CEMP - C0.8" => spectraPath("t57.50g45.00m-3.000 CEMP","C_0.8"),
	"t57.50g45.00m-3.000 CEMP - C0.6" => spectraPath("t57.50g45.00m-3.000 CEMP","C_0.6"),
	"t57.50g45.00m-3.000 CEMP - C0.4" => spectraPath("t57.50g45.00m-3.000 CEMP","C_0.4"),
	"t57.50g45.00m-3.000 CEMP - C0.2" => spectraPath("t57.50g45.00m-3.000 CEMP","C_0.2"),
	"t57.50g45.00m-3.000 CEMP - C0.0" => spectraPath("t57.50g45.00m-3.000 CEMP","C_0.0"),
	"t52.50g30.00m-3.000 CEMP - C1.2" => spectraPath("t52.50g30.00m-3.000 CEMP","C_1.2"),
	"t52.50g30.00m-3.000 CEMP - C0.8" => spectraPath("t52.50g30.00m-3.000 CEMP","C_0.8"),
	"t52.50g30.00m-3.000 CEMP - C0.6" => spectraPath("t52.50g30.00m-3.000 CEMP","C_0.6"),
	"t52.50g30.00m-3.000 CEMP - C0.4" => spectraPath("t52.50g30.00m-3.000 CEMP","C_0.4"),
	"t52.50g30.00m-3.000 CEMP - C0.2" => spectraPath("t52.50g30.00m-3.000 CEMP","C_0.2"),
	"t52.50g30.00m-3.000 CEMP - C0.0" => spectraPath("t52.50g30.00m-3.000 CEMP","C_0.0"),
	"t50.00g25.00m-3.000 CEMP - C1.2" => spectraPath("t50.00g25.00m-3.000 CEMP","C_1.2"),
	"t50.00g25.00m-3.000 CEMP - C0.8" => spectraPath("t50.00g25.00m-3.000 CEMP","C_0.8"),
	"t50.00g25.00m-3.000 CEMP - C0.6" => spectraPath("t50.00g25.00m-3.000 CEMP","C_0.6"),
	"t50.00g25.00m-3.000 CEMP - C0.4" => spectraPath("t50.00g25.00m-3.000 CEMP","C_0.4"),
	"t50.00g25.00m-3.000 CEMP - C0.2" => spectraPath("t50.00g25.00m-3.000 CEMP","C_0.2"),
	"t50.00g25.00m-3.000 CEMP - C0.0" => spectraPath("t50.00g25.00m-3.000 CEMP","C_0.0"),

	
	# FeH = -2
	"t57.50g45.00m-2.000 no-CEMP - C1.0" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_1.0"),
	"t52.50g30.00m-2.000 no-CEMP - C1.0" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_1.0"),
	# --> 3D effect fitting
	"t57.50g45.00m-2.000 no-CEMP - C1.2" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_1.2"),
	"t57.50g45.00m-2.000 no-CEMP - C0.8" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_0.8"),
	"t57.50g45.00m-2.000 no-CEMP - C0.6" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_0.6"),
	"t57.50g45.00m-2.000 no-CEMP - C0.4" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_0.4"),
	"t57.50g45.00m-2.000 no-CEMP - C0.2" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_0.2"),
	"t57.50g45.00m-2.000 no-CEMP - C0.0" => spectraPath("t57.50g45.00m-2.000 no-CEMP","C_0.0"),
	"t52.50g30.00m-2.000 no-CEMP - C1.2" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_1.2"),
	"t52.50g30.00m-2.000 no-CEMP - C0.8" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_0.8"),
	"t52.50g30.00m-2.000 no-CEMP - C0.6" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_0.6"),
	"t52.50g30.00m-2.000 no-CEMP - C0.4" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_0.4"),
	"t52.50g30.00m-2.000 no-CEMP - C0.2" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_0.2"),
	"t52.50g30.00m-2.000 no-CEMP - C0.0" => spectraPath("t52.50g30.00m-2.000 no-CEMP","C_0.0"),


	# Contribution functions
	"t50.00g25.00m-5.000 CEMP - contr" => spectraPath("t50.00g25.00m-5.000 CEMP" ,"contr"),
	"t57.50g45.00m-5.000 CEMP - contr" => spectraPath("t57.50g45.00m-5.000 CEMP" ,"contr"),
	"t52.50g30.00m-5.000 CEMP - contr" => spectraPath("t52.50g30.00m-5.000 CEMP" ,"contr"),
	"t57.50g45.00m-4.000 CEMP - contr" => spectraPath("t57.50g45.00m-4.000 CEMP" ,"contr"),
	"t52.50g30.00m-4.000 CEMP - contr" => spectraPath("t52.50g30.00m-4.000 CEMP" ,"contr"),
	"t57.50g45.00m-3.000 CEMP - contr" => spectraPath("t57.50g45.00m-3.000 CEMP" ,"contr"),
	"t52.50g30.00m-3.000 CEMP - contr" => spectraPath("t52.50g30.00m-3.000 CEMP" ,"contr"),
)

# ╔═╡ e67fed9e-3a82-45ed-9865-1c5862ebba2c
spectra3D = Dict(
	k=>M3DISRun(v) for (k, v) in spectrapath3D
)

# ╔═╡ f3493218-6013-4567-b1c4-442d0a48c9d0


# ╔═╡ 6ba576df-1a11-4be4-811b-212ed1f23347
md"## 1D spectra"

# ╔═╡ 63ed5788-011d-4a31-ac82-efcd135f609c
spectrapath1D = Dict(
	#=
	# Feh = -5
	"MARCS t57.50g45.00m-5.000 no-CEMP - C3.0" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_3.0",),
	"MARCS t57.50g45.00m-5.000 no-CEMP - C0.0" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_0.0"),

	# feh = -4
	"MARCS t57.50g45.00m-4.000 no-CEMP - C0.75" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_0.75",),
	"MARCS t57.50g45.00m-4.000 no-CEMP - C0.0" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_0.0"),

	# contribution functions
	"MARCS t57.50g45.00m-5.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_3.0"),
	"MARCS t57.50g45.00m-4.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_0.75")
	=#

	# After fixing molecular opacities
	# Feh = -5
	"MARCS t50.00g25.00m-5.000 no-CEMP - C3.0" => spectraPath1D("MARCS t50.00g25.00m-5.000 no-CEMP","C_3.0"),
	"MARCS t57.50g45.00m-5.000 no-CEMP - C3.0" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_3.0"),
	"MARCS t52.50g30.00m-5.000 no-CEMP - C3.0" => spectraPath1D("MARCS t52.50g30.00m-5.000 no-CEMP","C_3.0"),

	# Feh = -4
	"MARCS t50.00g25.00m-4.000 no-CEMP - C2.0" => spectraPath1D("MARCS t50.00g25.00m-4.000 no-CEMP","C_2.0"),
	"MARCS t57.50g45.00m-4.000 no-CEMP - C2.0" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_2.0"),
	"MARCS t52.50g30.00m-4.000 no-CEMP - C2.0" => spectraPath1D("MARCS t52.50g30.00m-4.000 no-CEMP","C_2.0"),

	# Feh = -3
	"MARCS t50.00g25.00m-3.000 no-CEMP - C1.0" => spectraPath1D("MARCS t50.00g25.00m-3.000 no-CEMP","C_1.0"),
	"MARCS t57.50g45.00m-3.000 no-CEMP - C1.0" => spectraPath1D("MARCS t57.50g45.00m-3.000 no-CEMP","C_1.0"),
	"MARCS t52.50g30.00m-3.000 no-CEMP - C1.0" => spectraPath1D("MARCS t52.50g30.00m-3.000 no-CEMP","C_1.0"),

	# Feh = -2
	"MARCS t57.50g45.00m-2.000 no-CEMP - C1.0" => spectraPath1D("MARCS t57.50g45.00m-2.000 no-CEMP","C_1.0"),
	"MARCS t52.50g30.00m-2.000 no-CEMP - C1.0" => spectraPath1D("MARCS t52.50g30.00m-2.000 no-CEMP","C_1.0"),

	# contribution functions
	"MARCS t50.00g25.00m-5.000 no-CEMP - contr" => spectraPath1D("MARCS t50.00g25.00m-5.000 no-CEMP","contr"),
	"MARCS t57.50g45.00m-5.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","contr"),
	"MARCS t52.50g30.00m-5.000 no-CEMP - contr" => spectraPath1D("MARCS t52.50g30.00m-5.000 no-CEMP","contr"),
	"MARCS t57.50g45.00m-4.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","contr"),
	"MARCS t52.50g30.00m-4.000 no-CEMP - contr" => spectraPath1D("MARCS t52.50g30.00m-4.000 no-CEMP","contr"),
	"MARCS t57.50g45.00m-3.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-3.000 no-CEMP","contr"),
	"MARCS t52.50g30.00m-3.000 no-CEMP - contr" => spectraPath1D("MARCS t52.50g30.00m-3.000 no-CEMP","contr"),
)

# ╔═╡ 4bf4013c-ee0d-42fd-9dd2-4f26f090c15b
begin
	spectra1D = Dict()
	for (k, v) in spectrapath1D 
		!isdir(v) && @show(v)
		spectra1D[k] = M3DISRun(v)
	end	
end

# ╔═╡ 4ca89537-5cc7-42e6-b131-362ad7ceb0b0


# ╔═╡ 0dee6625-70b1-489d-be3c-436062ca992b
md"# Design of Figures"

# ╔═╡ b231118c-bfed-4c9f-9522-cf7de097a12b
md"It is useful to assign design elements to a common interface."

# ╔═╡ 5d5df14a-3b77-4781-9a53-eb6031b71c95
struct PlotDesign
	key
	label
	ls
	lw
	color
end

# ╔═╡ 4cf9ee24-a75b-45c2-81b6-7012879a2217


# ╔═╡ c1c36104-a84d-4894-aa8c-33d3fdafc25b
md"# Structure comparison"

# ╔═╡ 8d7b17da-ead8-445c-9fb0-4a6d98ca5dda
md"## Carbon effect"

# ╔═╡ 4b636278-0e13-44d0-80c0-3c4fd9e1f0b6
design_p1 = [
	PlotDesign("t57.50g45.00m-5.000 CEMP", "carbon enhanced", "-", 2.5, "tomato"),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP", "not carbon enhanced", "--", 2, "k")
]

# ╔═╡ d786a454-b820-4c2f-b968-f7dec0604598
name1_p1 = "CEMP_vs_noCEMP_T_z_t57.50g45.00m-5.000.pdf"

# ╔═╡ 52ffd5ea-6215-4644-8e27-5777e594567a
name1_p1_1 = "CEMP_vs_noCEMP_T_tau_t57.50g45.00m-5.000.pdf"

# ╔═╡ cbb98c8a-8ec3-46fa-a38f-4afc3d292f41
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>first, :z, :T)
		x ./= 1e5
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-20, 280)
	ax.set_ylim(3000, 6000)
	ax.set_yscale("log")
	ax.set_xlabel("geometrical height [km]")
	ax.set_ylabel("temperature [K]")

	ax.legend()

	f.savefig(name1_p1)

	f
end

# ╔═╡ 2aadb41f-a64e-4769-a202-aa9d2d700104
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>last, :log10τ_ross, :T)
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(3200, 7900)
	#ax.set_yscale("log")
	ax.set_title(L"\rm T_{eff}=5\,750\ K,\ log(g)=4.5,\ [Fe/H] = -5,\ [C/Fe] = 3")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{ross}]")
	ax.set_ylabel(L"\rm temperature\ [K]")

	ax.legend()

	f.savefig(name1_p1_1)

	f
end

# ╔═╡ 56d11775-333b-44ac-ba5f-39e1b372564d
name2_p1 = "CEMP_vs_noCEMP_d_z_t57.50g45.00m-5.000.pdf"

# ╔═╡ 58ff1560-c4ad-4ac5-bf29-33977b211b11
name2_p1_1 = "CEMP_vs_noCEMP_d_tau_t57.50g45.00m-5.000.pdf"

# ╔═╡ c421ceca-ee21-455c-bfec-371c6c260ea3
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>first, :z, :log10d)
		x ./= 1e5
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-50, 230)
	ax.set_ylim(-7, -6.1)
	ax.set_xlabel("geometrical height [km]")
	ax.set_ylabel("log density [g cm-3]")

	ax.legend()

	f.savefig(name2_p1)

	f
end

# ╔═╡ 41671152-f0cb-4d1f-8aec-740a486c69b4
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>last, :log10τ_ross, :log10d)
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(-7, -5.9)
	#ax.set_yscale("log")
	ax.set_title(L"\rm T_{eff}=5\,750\ K,\ log(g)=4.5,\ [Fe/H] = -5,\ [C/Fe] = 3")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{ross}]")
	ax.set_ylabel(L"\rm density\ [g\ cm^{-3}]")

	ax.legend()

	f.savefig(name2_p1_1)

	f
end

# ╔═╡ c437b8d0-3da5-4d74-b796-0accda21a1f6


# ╔═╡ 13bdade4-3eb8-41cb-aa64-463124ccacf1
design_p1B = [
	PlotDesign("t52.50g30.00m-5.000 CEMP", "carbon enhanced", "-", 2.5, "tomato"),
	PlotDesign("t52.50g30.00m-5.000 no-CEMP", "not carbon enhanced", "--", 2, "k")
]

# ╔═╡ f3398799-5ad0-4487-9ae5-cc9d7b679f3f
name1_p1B = "CEMP_vs_noCEMP_T_z_t52.50g30.00m-5.000.pdf"

# ╔═╡ 01720034-cc2c-4bc5-9d99-f2073b5d0bcc
name1_p1B_1 = "CEMP_vs_noCEMP_T_tau_t52.50g30.00m-5.000.pdf"

# ╔═╡ 455ca2ac-6fe2-457a-941d-4f2b8c2f58b3
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1B)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>first, :z, :T)
		x ./= 1e5
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-1000, 9000)
	ax.set_ylim(3000, 7200)
	ax.set_yscale("log")
	ax.set_xlabel("geometrical height [km]")
	ax.set_ylabel("temperature [K]")

	ax.legend()

	f.savefig(name1_p1B)

	f
end

# ╔═╡ 124e57c4-82e6-4abf-bb48-7f882da71c5b
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1B)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>last, :log10τ_ross, :T)
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(3200, 7200)
	#ax.set_yscale("log")
	ax.set_title(L"\rm T_{eff}=5\,250\ K,\ log(g)=3.0,\ [Fe/H] = -5,\ [C/Fe] = 3")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{ross}]")
	ax.set_ylabel(L"\rm temperature\ [K]")

	ax.legend()

	f.savefig(name1_p1B_1)

	f
end

# ╔═╡ 4d2a898f-8bc4-4c88-b35d-cb7542d55932
name2_p1B = "CEMP_vs_noCEMP_d_z_t52.50g30.00m-5.000.pdf"

# ╔═╡ ae9d6151-45cd-4301-bf3f-df319baa0870
name2_p1B_1 = "CEMP_vs_noCEMP_d_tau_t52.50g30.00m-5.000.pdf"

# ╔═╡ f8c93a98-ad4f-4d9c-b0cd-c66d432753ed
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1B)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>first, :z, :log10d)
		x ./= 1e5
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-1000, 10000)
	ax.set_ylim(-8, -6.8)
	ax.set_xlabel("geometrical height [km]")
	ax.set_ylabel("log density [g cm-3]")

	ax.legend()

	f.savefig(name2_p1B)

	f
end

# ╔═╡ 0d0234bc-c6b5-4bcf-8409-01cf02c8c843
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1B)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>last, :log10τ_ross, :log10d)
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-4.2, 1)
	ax.set_ylim(-7.78, -6.65)
	#ax.set_yscale("log")
	ax.set_title(L"\rm T_{eff}=5\,250\ K,\ log(g)=3.0,\ [Fe/H] = -5,\ [C/Fe] = 3")
	ax.set_xlabel(L"\rm optical\ depth\ [\tau_{ross}]")
	ax.set_ylabel(L"\rm density\ [g\ cm^{-3}]")

	ax.legend()

	f.savefig(name2_p1B_1)

	f
end

# ╔═╡ 773a7896-a492-46a3-9813-fa739647142d


# ╔═╡ 244f1b18-80fe-419a-8fd0-892aad803ca5
md"## 3D metallicity effect"

# ╔═╡ accdaadc-8912-48d2-8a0c-3c19e5b4ba19
design_p2 = [
	PlotDesign("t57.50g45.00m-5.000 no-CEMP", "DISPATCH <3D>", "-", 2.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP", "MARCS", "--", 2., "k"),
]

# ╔═╡ ffb5a590-77e5-41c5-8964-eafe1153a5dc
name1_p2 = "vs_1D_t_tau_t57.50g45.00m-5.000.pdf"

# ╔═╡ 778a7f04-0841-4801-8aa8-b5686e1e498c
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	for (i, p) in enumerate(design_p2)
		k = p.key
		if haskey(models3D, k)
			b = models3D[k][2]
			#=x, y = log10.(b[:τ_ross]), b[:T]
			
			H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=50)

			xe = MUST.pyconvert(Array, xedges)
			ye = MUST.pyconvert(Array, yedges)
			H = MUST.pyconvert(Array, H.T)
			
			ax.imshow(
				H, interpolation="bicubic", origin="lower", aspect="auto",
				extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
			)=#

			x, y = profile(MUST.mean, b, :log10τ_ross, :T)
			ax.plot(x, y, color=p.color, label=p.label, ls=p.ls, lw=p.lw)
		elseif haskey(models1D, k)
			b = models1D[k]
			x = log10.(b.τ)
			y = exp.(b.lnT)
			ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
		end
	end

	ax.set_xlim(-5.5, 4)
	ax.set_ylim(3000, 12000)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{ross}"*"]")
	ax.set_ylabel("temperature [K]")

	ax.legend(loc="upper left")
	f.savefig(name1_p2)

	f
end

# ╔═╡ 983f1c2d-1532-491d-be04-e52e64c42775


# ╔═╡ 5a5704c2-2479-4bbd-87c5-cc015d546612
design_p2B = [
	PlotDesign("t52.50g30.00m-4.000 no-CEMP", "DISPATCH <3D>", "-", 2.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-4.000 no-CEMP", "MARCS", "--", 2., "k"),
]

# ╔═╡ 785a473c-1a46-4ef3-8c2f-e7170452064b
name1_p2B = "vs_1D_t_tau_t52.50g30.00m-4.000.pdf"

# ╔═╡ 51f58ec5-f1bc-4b66-8ef1-c65ad08f15c9
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	for (i, p) in enumerate(design_p2B)
		k = p.key
		if haskey(models3D, k)
			b = models3D[k][2]
			#=x, y = log10.(b[:τ_ross]), b[:T]
			
			H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=50)

			xe = MUST.pyconvert(Array, xedges)
			ye = MUST.pyconvert(Array, yedges)
			H = MUST.pyconvert(Array, H.T)
			
			ax.imshow(
				H, interpolation="bicubic", origin="lower", aspect="auto",
				extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
			)=#

			x, y = profile(MUST.mean, b, :log10τ_ross, :T)
			ax.plot(x, y, color=p.color, label=p.label, ls=p.ls, lw=p.lw)
		elseif haskey(models1D, k)
			b = models1D[k]
			x = log10.(b.τ)
			y = exp.(b.lnT)
			ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
		end
	end

	ax.set_xlim(-5.0, 4)
	ax.set_ylim(2800, 12000)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{ross}"*"]")
	ax.set_ylabel("temperature [K]")

	ax.legend(loc="upper left")
	f.savefig(name1_p2B)

	f
end

# ╔═╡ f61fb319-518f-4b56-9e19-a172e499307e


# ╔═╡ 0487f738-1a0e-49f9-96c3-357634226467
design_p2C = [
	PlotDesign("t52.50g30.00m-3.000 no-CEMP", "DISPATCH <3D>", "-", 2.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-3.000 no-CEMP", "MARCS", "--", 2., "k"),
]

# ╔═╡ cb4ca3c9-6c8a-4876-b927-a1419e403f43
name1_p2C = "vs_1D_t_tau_t52.50g30.00m-3.000.pdf"

# ╔═╡ 56e5777e-e80d-458d-b6a0-cde6a8f811ea
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	for (i, p) in enumerate(design_p2C)
		k = p.key
		if haskey(models3D, k)
			b = models3D[k][2]
			#=x, y = log10.(b[:τ_ross]), b[:T]
			
			H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=50)

			xe = MUST.pyconvert(Array, xedges)
			ye = MUST.pyconvert(Array, yedges)
			H = MUST.pyconvert(Array, H.T)
			
			ax.imshow(
				H, interpolation="bicubic", origin="lower", aspect="auto",
				extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
			)=#

			x, y = profile(MUST.mean, b, :log10τ_ross, :T)
			ax.plot(x, y, color=p.color, label=p.label, ls=p.ls, lw=p.lw)
		elseif haskey(models1D, k)
			b = models1D[k]
			x = log10.(b.τ)
			y = exp.(b.lnT)
			ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
		end
	end

	ax.set_xlim(-5.0, 4)
	ax.set_ylim(2800, 12000)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{ross}"*"]")
	ax.set_ylabel("temperature [K]")

	ax.legend(loc="upper left")
	f.savefig(name1_p2C)

	f
end

# ╔═╡ 89d5c1d8-733b-4f11-b306-a7410ae48e75


# ╔═╡ 85d716f2-542a-4793-8fa7-a53383fb37ec
design_p2D = [
	PlotDesign("t52.50g30.00m-2.000 no-CEMP", "DISPATCH <3D>", "-", 2.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-2.000 no-CEMP", "MARCS", "--", 2., "k"),
]

# ╔═╡ 8992c977-8832-4102-89b7-430b4abd7e46
name1_p2D = "vs_1D_t_tau_t52.50g30.00m-2.000.pdf"

# ╔═╡ 43c0fc35-8336-43ba-a2c1-c4ae2ce33d43
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	for (i, p) in enumerate(design_p2D)
		k = p.key
		if haskey(models3D, k)
			b = models3D[k][2]
			#=x, y = log10.(b[:τ_ross]), b[:T]
			
			H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=50)

			xe = MUST.pyconvert(Array, xedges)
			ye = MUST.pyconvert(Array, yedges)
			H = MUST.pyconvert(Array, H.T)
			
			ax.imshow(
				H, interpolation="bicubic", origin="lower", aspect="auto",
				extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
			)=#

			x, y = profile(MUST.mean, b, :log10τ_ross, :T)
			ax.plot(x, y, color=p.color, label=p.label, ls=p.ls, lw=p.lw)
		elseif haskey(models1D, k)
			b = models1D[k]
			x = log10.(b.τ)
			y = exp.(b.lnT)
			ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
		end
	end

	ax.set_xlim(-5.0, 4)
	ax.set_ylim(2800, 12000)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{ross}"*"]")
	ax.set_ylabel("temperature [K]")

	ax.legend(loc="upper left")
	f.savefig(name1_p2D)

	f
end

# ╔═╡ b0d6cf77-7197-4c1a-a59a-679dca8ddf7e


# ╔═╡ d40dd522-128d-4407-a53e-62977f43b27d
md"## Optical surface slices"

# ╔═╡ 0b060c44-9da6-46d5-9c8d-4a6d0285c2a2
md"### [Fe/H] = -5"

# ╔═╡ 07c86b08-0181-491e-952b-2e468ea8896c
design_p2E = PlotDesign("t50.00g25.00m-5.000 CEMP", L"\rm T_{eff}=5\,000\ K, log(g)=2.5,\ [Fe/H]=-5", "-", 2.5, "tomato")

# ╔═╡ 87581ad4-4fe8-4e46-9529-9255448f22fd
name1_p2E = "optical_surface_t50.00g25.00m-5.000.pdf" 

# ╔═╡ 1ef702b2-1a94-4794-acd7-4b5516595203
let
	b = models3D[design_p2E.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2E.label)

	f.savefig(name1_p2E)
	f
end

# ╔═╡ 26e9ea76-9c43-471e-a195-73150208ae19
design_p2F = PlotDesign("t52.50g30.00m-5.000 CEMP", L"\rm T_{eff}=5\,250\ K, log(g)=3.0,\ [Fe/H]=-5 ", "-", 2.5, "tomato")

# ╔═╡ f8461a99-8d27-4293-bdc4-9afe19728a4e
name1_p2F = "optical_surface_t52.50g30.00m-5.000.pdf"

# ╔═╡ 86e48b71-a6a9-4c3a-9227-7f013324c813
let
	b = models3D[design_p2F.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2F.label)

	f.savefig(name1_p2F)
	
	f
end

# ╔═╡ 228cda78-3c92-41fa-bd8d-74b8714ce094
design_p2G = PlotDesign("t57.50g45.00m-5.000 CEMP", L"\rm T_{eff}=5\,750\ K, log(g)=4.5,\ [Fe/H]=-5 ", "-", 2.5, "tomato")

# ╔═╡ 2158c40c-a464-45a3-a0c9-2fc69fe635da
name1_p2G = "optical_surface_t57.50g45.00m-5.000.pdf"

# ╔═╡ 326ee706-2b03-4fd5-bc16-e396e2ec4a1d
let
	b = models3D[design_p2G.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2G.label)

	f.savefig(name1_p2G)
	
	f
end

# ╔═╡ 96e30c0f-9493-42b0-9c38-888b653723d1


# ╔═╡ 59566cb5-45c3-41a6-977c-b6af33e15c1f
md"### [Fe/H] = -3"

# ╔═╡ 30fa8414-9468-49ac-9c9e-8e61823a0263
design_p2E_2 = PlotDesign("t50.00g25.00m-3.000 CEMP", L"\rm T_{eff}=5\,000\ K, log(g)=2.5,\ [Fe/H]=-3", "-", 2.5, "tomato")

# ╔═╡ d95e2320-31bb-4e50-8292-d7c713d58a87
name1_p2E_2 = "optical_surface_t50.00g25.00m-3.000.pdf"

# ╔═╡ 1ed777ae-6b23-4bdf-8a6b-5336e88fddde
let
	b = models3D[design_p2E_2.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2E_2.label)

	f.savefig(name1_p2E_2)
	f
end

# ╔═╡ 5c7630d5-e305-467c-87d2-1caebd62c1f2
design_p2F_2 = PlotDesign("t52.50g30.00m-3.000 CEMP", L"\rm T_{eff}=5\,250\ K, log(g)=3.0,\ [Fe/H]=-3 ", "-", 2.5, "tomato")

# ╔═╡ 234cfaf6-099c-4fc1-a776-c44de65620c0
name1_p2F_2 = "optical_surface_t52.50g30.00m-3.000.pdf"

# ╔═╡ bfdcd5c1-621f-49b5-b597-fef0d87d30d9
let
	b = models3D[design_p2F_2.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2F_2.label)

	f.savefig(name1_p2F_2)
	
	f
end

# ╔═╡ c2c504e6-8f4d-4220-9a2a-099c17a06565
design_p2G_2 = PlotDesign("t57.50g45.00m-3.000 CEMP", L"\rm T_{eff}=5\,750\ K, log(g)=4.5,\ [Fe/H]=-3 ", "-", 2.5, "tomato")

# ╔═╡ 0265aaca-4836-4082-9d21-c11a5f0fde86
name1_p2G_2 = "optical_surface_t57.50g45.00m-3.000.pdf"

# ╔═╡ 11a5eed7-b234-4d42-81b8-b415c8679ff5
let
	b = models3D[design_p2G_2.key][2]
	Tplane = MUST.interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
	
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	extent = [minimum(b.x), maximum(b.x), minimum(b.y), maximum(b.y)] ./1e8
	im = ax.imshow(
		Tplane, 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
	cbar.set_label(L"\rm temperature\ [K]")

	ax.set_xlabel(L"\rm x\ [Mm]")
	ax.set_ylabel(L"\rm y\ [Mm]")
	ax.set_title(design_p2G_2.label)

	f.savefig(name1_p2G_2)
	
	f
end

# ╔═╡ 83729ef7-9c05-4c0b-8bcf-0ecfccc48d66


# ╔═╡ 0d4a0a2f-758f-44a0-9692-23c082befac8
md"# Spectra"

# ╔═╡ 1f6dfdaf-b82e-46ee-989f-dbeb685d72ff
md"## Carbon effect"

# ╔═╡ 391e98d1-344e-4ec3-837b-58327c02c74b
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 10d7d97f-5c54-4066-9618-c78983c6c18c
design_p3 = [
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "carbon enhanced", "-", 1.5, "tomato"),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "not carbon enhanced", "--", 1.5, "k")
]

# ╔═╡ b97665fe-1536-4047-b3e9-d6176d47f3c1
other_spectra_p3 = Dict(
	#="t57.50g45.00m-5.000 CEMP - C3.0"=>[
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 375", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 374", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 373", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 372", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 369", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 368", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 367", "", "-", 1.5, "0.5"), 
	],=#
)

# ╔═╡ a04336bf-9cdf-401e-b237-a41f44987320
name1_p3 = "CEMP_vs_noCEMP_spectrum_t57.50g45.00m-5.000.pdf"

# ╔═╡ 1deeb995-4fda-4071-8d1f-fe9a9bc8a09c
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3
		k = p.key

		if k in keys(other_spectra_p3)
			for p2 in other_spectra_p3[k]
				k2 = p2.key
				s = spectra3D[k2]
				λ, F = flux(s, norm=true)
				ax.plot(λ, F, color=p2.color, ls=p2.ls, lw=p2.lw, label=p2.label)
			end
		end
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	#ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	ax.legend(loc="lower center", ncol=2)

	ax.set_title(L"\rm [Fe/H]=-5,\ [C/Fe]=-3")

	f.savefig(name1_p3)
	
	f
end

# ╔═╡ 379e05ed-20db-499a-a0bd-b93f0e3e3642


# ╔═╡ 0b18ead3-5444-4737-99a8-cd86404343d3
design_p3_2 = [
	PlotDesign("t52.50g30.00m-5.000 CEMP - C3.0", "carbon enhanced", "-", 1.5, "tomato"),
	PlotDesign("t52.50g30.00m-5.000 no-CEMP - C3.0", "not carbon enhanced", "--", 1.5, "k")
]

# ╔═╡ 9e9ce7c5-2860-4cb0-9ddb-52323a4beb1d
other_spectra_p3_2 = Dict(
	#="t57.50g45.00m-5.000 CEMP - C3.0"=>[
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 375", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 374", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 373", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 372", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 369", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 368", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 367", "", "-", 1.5, "0.5"), 
	],=#
)

# ╔═╡ b28b1962-34c6-49f8-973a-3f63eed6c4d0
name1_p3_2 = "CEMP_vs_noCEMP_spectrum_t52.50g30.00m-5.000.pdf"

# ╔═╡ c0855af1-d407-4ec2-83fc-4397330ac35c
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3_2
		k = p.key

		if k in keys(other_spectra_p3)
			for p2 in other_spectra_p3[k]
				k2 = p2.key
				s = spectra3D[k2]
				λ, F = flux(s, norm=true)
				ax.plot(λ, F, color=p2.color, ls=p2.ls, lw=p2.lw, label=p2.label)
			end
		end
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p3_2)
	
	f
end

# ╔═╡ 41367bd1-d61c-4b45-8781-6b5b71c7a511


# ╔═╡ 1943029c-2f0c-4590-9619-89870040932c
design_p3_3 = [
	PlotDesign("t50.00g25.00m-5.000 CEMP - C3.0", "carbon enhanced", "-", 1.5, "tomato"),
	PlotDesign("t50.00g25.00m-5.000 no-CEMP - C3.0", "not carbon enhanced", "--", 1.5, "k")
]

# ╔═╡ 5f7c6598-73d9-4e78-90ff-a53fa8ccf6c7
other_spectra_p3_3 = Dict(
	#="t57.50g45.00m-5.000 CEMP - C3.0"=>[
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 375", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 374", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 373", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 372", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 369", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 368", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 367", "", "-", 1.5, "0.5"), 
	],=#
)

# ╔═╡ ff48212b-40ae-4bcf-8576-034cd46f4738
name1_p3_3 = "CEMP_vs_noCEMP_spectrum_t50.00g25.00m-5.000.pdf"

# ╔═╡ 7df2c30a-445e-498e-9ebd-1a206ece821e
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3_3
		k = p.key

		if k in keys(other_spectra_p3_3)
			for p2 in other_spectra_p3_3[k]
				k2 = p2.key
				s = spectra3D[k2]
				λ, F = flux(s, norm=true)
				ax.plot(λ, F, color=p2.color, ls=p2.ls, lw=p2.lw, label=p2.label)
			end
		end
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	#ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	ax.legend(loc="lower center", ncol=2)

	ax.set_title(L"\rm [Fe/H]=-5,\ [C/Fe]=-3")

	f.savefig(name1_p3_3)
	
	f
end

# ╔═╡ 2aa33bb2-3cd3-4104-a6a2-0f1b12fe2960


# ╔═╡ c905fd9a-e2f6-485c-9096-c6a588b31b4d
md"### [Fe/H]=-4, [C/Fe]=2.0"

# ╔═╡ 91699cc6-cb81-4382-bc5a-01931914343b
design_p3B = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - C2.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t57.50g45.00m-4.000 no-CEMP - C2.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ aeccdb4b-63f9-47aa-830e-21aee1093dc4
name2_p3 = "CEMP_vs_noCEMP_spectrum_t57.50g45.00m-4.000.pdf"

# ╔═╡ a55d5b04-a41b-47b6-81e8-39073bf538e2
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3B
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p3)
	
	f
end

# ╔═╡ a056ba9a-f392-4945-a2c3-a64fc1849831


# ╔═╡ 77bc95ff-f456-46da-a81d-5412e6be032e
design_p3B_2 = [
	PlotDesign("t52.50g30.00m-4.000 CEMP - C2.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t52.50g30.00m-4.000 no-CEMP - C2.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ baede960-e93b-490a-a027-1caf661428c5
name2_p3_2 = "CEMP_vs_noCEMP_spectrum_t52.50g30.00m-4.000.pdf"

# ╔═╡ 3537629b-0766-43e3-9506-30cfe0cc4a30
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3B_2
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p3_2)
	
	f
end

# ╔═╡ 136c2f2a-8902-4e05-920d-1d28c918a3f1


# ╔═╡ dddd50fb-8044-4755-8180-647311268a91
design_p3B_3 = [
	PlotDesign("t50.00g25.00m-4.000 CEMP - C2.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t50.00g25.00m-4.000 no-CEMP - C2.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ e551c818-ed62-4a2e-b2ac-74bfe6fb9b5f
name2_p3_3 = "CEMP_vs_noCEMP_spectrum_t50.00g25.00m-4.000.pdf"

# ╔═╡ 809c6511-ff14-47e4-80d7-4c984b3e8d3f
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3B_3
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p3_3)
	
	f
end

# ╔═╡ f0ff834f-203c-460c-af2e-ee1fed8be671


# ╔═╡ 64073a6c-f9b1-45ad-8494-b948f6ec08d2
md"### [Fe/H]=-3, [C/Fe]=1.0"

# ╔═╡ c7631a99-5344-4a84-b1e1-5d9d96ddfeeb
design_p3C = [
	PlotDesign("t57.50g45.00m-3.000 CEMP - C1.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t57.50g45.00m-3.000 no-CEMP - C1.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ fe6d9b1c-a5c3-4b03-b66a-980191b21ba8
name3_p3 = "CEMP_vs_noCEMP_spectrum_t57.50g45.00m-3.000.pdf"

# ╔═╡ 00327ee0-ef59-43aa-bfe9-768f0f154f2c
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3C
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p3)
	
	f
end

# ╔═╡ b0d32b4a-a107-47be-b12c-b4ff6ce46b53


# ╔═╡ 825031cd-5e40-49b6-936d-4a263a0c0378
design_p3C_2 = [
	PlotDesign("t52.50g30.00m-3.000 CEMP - C1.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t52.50g30.00m-3.000 no-CEMP - C1.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ 948ee564-7a5f-4581-9ca5-b0787ac41823
name3_p3_2 = "CEMP_vs_noCEMP_spectrum_t52.50g30.00m-3.000.pdf"

# ╔═╡ 18d9db1c-67b4-40f4-8548-057c7039020f
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3C_2
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p3_2)
	
	f
end

# ╔═╡ f42278a0-ed19-47a8-bd2e-d78d1ef83d3a


# ╔═╡ d12753f2-9c57-4af0-9623-194955c17cb7
design_p3C_3 = [
	PlotDesign("t50.00g25.00m-3.000 CEMP - C1.0", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t50.00g25.00m-3.000 no-CEMP - C1.0", "scaled-solar model", "--", 1.5, "k")
]

# ╔═╡ 77ac5fd5-959a-452c-af88-73225291a890
name3_p3_3 = "CEMP_vs_noCEMP_spectrum_t50.00g25.00m-3.000.pdf"

# ╔═╡ aa555374-910d-4202-9862-6742cb49f0af
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p3C_3
		k = p.key
		
		s = spectra3D[k]
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p3_3)
	
	f
end

# ╔═╡ 1397cb72-51db-45bf-8a67-cf1352f7f990


# ╔═╡ bd3ee310-3f68-4e08-8ba8-d881aeecb76a
md"## 3D metallicity effect"

# ╔═╡ f064d1a4-be5b-4cfc-9bf9-67002a333d2d
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 9987b27b-2437-45a3-9e83-116a01f13d1d
design_p4 = [
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "CEMP", "-", 1.7, "k"),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "non-CEMP", "--", 1.7, "steelblue"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", "MARCS", "-", 1.7, "tomato")
]

# ╔═╡ 17ced5c3-080f-4d65-aa1c-b24a2950e783
name1_p4 = "vs_1D_spectrum_t57.50g45.00m-5.000.pdf"

# ╔═╡ 75600b27-d640-438b-bd00-d35e18ff46f6
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(12, 4))
	
	for p in design_p4
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.11, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	#ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05))

	ax.set_title(L"\rm T_{eff}=5\,750\ K,\ log(g)=4.5,\ [Fe/H] = -5")
	
	f.savefig(name1_p4)
	
	f
end

# ╔═╡ 14ee9206-f322-4e06-a65e-271c19d4c792


# ╔═╡ 6c428722-39e9-4f1c-851a-ee90cb617377
design_p4_2 = [
	PlotDesign("t52.50g30.00m-5.000 CEMP - C3.0", "CEMP", "-", 1.7, "k"),
	PlotDesign("t52.50g30.00m-5.000 no-CEMP - C3.0", "non-CEMP", "--", 1.5, "steelblue"),
	PlotDesign("MARCS t52.50g30.00m-5.000 no-CEMP - C3.0", "MARCS", "-", 1.7, "tomato")
]

# ╔═╡ 059c7c73-acaa-468b-b7a4-e63a4748d358
name1_p4_2 = "vs_1D_spectrum_t52.50g30.00m-5.000.pdf"

# ╔═╡ df0a0a2d-52ea-4852-b450-6fc1baf5838d
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(12, 4))
	
	for p in design_p4_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.11, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	#ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05))
	ax.set_title(L"\rm T_{eff}=5\,250\ K,\ log(g)=3.0,\ [Fe/H] = -5")

	f.savefig(name1_p4_2)
	
	f
end

# ╔═╡ 9669b2e7-396d-4e09-827b-3138c3d8f70a


# ╔═╡ 20e610c2-9054-4dc6-a551-6f63909715ee
design_p4_3 = [
	PlotDesign("t50.00g25.00m-5.000 CEMP - C3.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t50.00g25.00m-5.000 no-CEMP - C3.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 5d9fc015-4d8b-4811-8f48-c4197eb95709
name1_p4_3 = "vs_1D_spectrum_t50.00g25.00m-5.000.pdf"

# ╔═╡ e09757e6-29c2-414b-9ac8-2a5fe237e0e6
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	#ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	ax.legend(loc="lower center", ncol=2)
	ax.set_title(L"\rm T_{eff}=5\,000\ K,\ log(g)=2.5,\ [Fe/H] = -5")

	f.savefig(name1_p4_3)
	
	f
end

# ╔═╡ 59ab2b58-2e38-4cce-9ba4-bc6899616a7d


# ╔═╡ a61f7ec3-88c9-42f7-bf5c-acc792a05fb5
md"### [Fe/H]=-4, [C/Fe]=2.0"

# ╔═╡ 2aa8d303-dae3-4a23-b253-c6091e2a638d
design_p4B = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - C2.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C2.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ cec8b45b-e4b6-4ad0-a83d-e9afe62b486a
name2_p4 = "vs_1D_spectrum_t57.50g45.00m-4.000.pdf" 

# ╔═╡ ff0787a1-b1cb-4eb0-85d7-b33bd39b2533
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4B
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p4)
	f
end

# ╔═╡ f8bf1f6c-5522-463f-9ff0-cb59c1a950a2


# ╔═╡ fda2a05e-48b0-471a-a5a5-4ada3d0d5d78
design_p4B_2 = [
	PlotDesign("t52.50g30.00m-4.000 CEMP - C2.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-4.000 no-CEMP - C2.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 9d7f44dc-a938-4f7e-8355-b8dc31e707a8
name2_p4_2 = "vs_1D_spectrum_t52.50g30.00m-4.000.pdf" 

# ╔═╡ 334059d2-ea76-4c0d-84d9-7b1b9ef1ea87
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4B_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p4_2)
	f
end

# ╔═╡ d3530325-728c-4392-b25d-18f4e5651719


# ╔═╡ 745b79ab-48f7-424b-ba20-9d7761c1cd98
md"### [Fe/H]=-3, [C/Fe]=1.0"

# ╔═╡ 916bfaf4-10da-4206-ad4e-a77e8274d0c5
design_p4C = [
	PlotDesign("t57.50g45.00m-3.000 CEMP - C1.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-3.000 no-CEMP - C1.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ f9a45d78-cba7-461b-9237-eeac330cf6d1
name3_p4 = "vs_1D_spectrum_t57.50g45.00m-3.000.pdf"

# ╔═╡ 5af5f222-d401-4dfb-aeea-0531d9a160ff
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4C
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p4)
	f
end

# ╔═╡ fd2cb9f8-59da-467b-9ac0-81002549ce51


# ╔═╡ f5e18ac4-8058-4b25-b1f6-86ad30f68dbe
design_p4C_2 = [
	PlotDesign("t52.50g30.00m-3.000 CEMP - C1.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-3.000 no-CEMP - C1.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 16a3524c-1367-4f62-a657-b0c6ed85159d
name3_p4_2 = "vs_1D_spectrum_t52.50g30.00m-3.000.pdf"

# ╔═╡ 681fc74f-cecd-4771-a967-f28864520323
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4C_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p4_2)
	f
end

# ╔═╡ 5924f436-2453-4b92-8b83-8e8033b86c7f


# ╔═╡ a8ea82ac-7516-4223-9fa9-8ec93b63416d
design_p4C_3 = [
	PlotDesign("t50.00g25.00m-3.000 CEMP - C1.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t50.00g25.00m-3.000 no-CEMP - C1.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 2a23cc14-6593-4f54-a42b-65ded01323af
name3_p4_3 = "vs_1D_spectrum_t50.00g25.00m-3.000.pdf"

# ╔═╡ 28e040c9-a4fb-4e58-9e00-2427c69ba268
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4C_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p4_3)
	f
end

# ╔═╡ c09bfd6d-4801-4d8d-8755-147af007c93b


# ╔═╡ f928a5ce-3022-4dc8-a70c-75971734398c
md"### [Fe/H]=-2, [C/Fe]=1.0"

# ╔═╡ 0d68e17d-da25-4cb4-a89b-7465644e0f20
design_p4D = [
	PlotDesign("t57.50g45.00m-2.000 no-CEMP - C1.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-2.000 no-CEMP - C1.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 721f3a04-aa21-4ff5-b566-4968d486438a
name4_p4 = "vs_1D_spectrum_t57.50g45.00m-2.000.pdf"

# ╔═╡ 335f1b86-a585-4ffa-a1d8-6790739f1bbb
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4D
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name4_p4)
	f
end

# ╔═╡ 6ec5f942-94d4-4dfa-9eee-3baa7636476d
design_p4D_2 = [
	PlotDesign("t52.50g30.00m-2.000 no-CEMP - C1.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t52.50g30.00m-2.000 no-CEMP - C1.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 20e07853-44ed-47d9-8937-b923bc570df9
name4_p4_2 = "vs_1D_spectrum_t52.50g30.00m-4.000.pdf"

# ╔═╡ 45ef7372-5e7c-4f1a-a837-584c5dc16b6e
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
	for p in design_p4D_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW
	
		#name = replace(k, ' '=>'_')
		#name = name * "_4297_4303.txt"
		#MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end
	
	
	#ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	#ax.set_ylim(0.94, 1.0)
	
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name4_p4_2)
	f
end

# ╔═╡ da1e169c-1c79-4d84-ac61-e1fab544f696


# ╔═╡ 76464024-a571-4c72-bb59-f881b595c124
md"# CH Abundance"

# ╔═╡ e0e4254a-e9f3-491c-9cac-6ed8a349e1ae
"""
	plot_abund_3D(s, name; ax)

Plot abundances as a function of optical depth (500).
"""
plot_abund_3D(s, name; ax, cmap="Greys") = begin
	totn = MUST.pyconvert(
		Array, 
		s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0]
	)
	toth = MUST.pyconvert(
		Array, s.run.get_toth()
	)
	totn[totn .<= 0.0] .= minimum(totn[totn .> 0.0])
	C1 = log10.(totn ./ toth) .+ 12

	#= Plot =#
	x = MUST.pyconvert(Array, s.run.ltau)
	y = C1
	
	H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=300)

	xe = MUST.pyconvert(Array, xedges)
	ye = MUST.pyconvert(Array, yedges)
	H = MUST.pyconvert(Array, H.T)

	xx, yy = meshgrid(xe, ye)
	
	#=ax.imshow(
		H, interpolation="bicubic", origin="lower", aspect="auto",
		extent=[xe[1], xe[end], ye[1], ye[end]], cmap=cmap, norm=matplotlib.colors.LogNorm()
	)=#
	ax.pcolormesh(xx, yy, H', cmap=cmap, norm=matplotlib.colors.LogNorm(), rasterized=true)
end

# ╔═╡ 36c4a557-54d8-4787-a0c9-4942b9dff4ef
"""
	plot_abund_av_3D(s, name; ax, kwargs...)

Plot abundances (averages) as a function of optical depth (500).
"""
plot_abund_av_3D(s, name; ax, kwargs...) = begin
	try
		totn = MUST.pyconvert(
			Array, 
			s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0]
		)
		toth = MUST.pyconvert(Array, s.run.get_toth())
		C1 = log10.(totn ./ toth) .+ 12

		#= Plot =#
		x = MUST.pyconvert(Array, s.run.ltau)
		y = C1

		x = [MUST.mean(x[:, :, z]) for z in axes(x, 3)]
		y = [MUST.mean(y[:, :, z]) for z in axes(y, 3)]
		
		ax.plot(x, y; kwargs...)

		x, y
	catch
		totn = MUST.pyconvert(
			Array, 
			s.run.read_patch_save(
				name, lazy=false, concat=true, fdim=0, zdim=2
			)[0][0][0]
		)
		toth = MUST.pyconvert(Array, s.run.get_toth())
		C1 = log10.(totn ./ toth) .+ 12

		#= Plot =#
		x = MUST.pyconvert(Array, s.run.ltau)
		y = C1

		ax.plot(x, y; kwargs...)

		x, y
	end
end

# ╔═╡ 1ddb0459-d832-49a4-9bf2-af1627edd2e5


# ╔═╡ 6d6b3b86-0bec-427a-b795-3c31c96f9968
md"## Carbon effect"

# ╔═╡ ce9723af-dbc8-440c-aeaf-c88dd902d57b
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ bf3b3496-3d18-40af-a89c-dc3dbc4d54c2
design_p5A = [
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "CEMP", "-", 2, ""),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "no CEMP", "--", 2, "")
]

# ╔═╡ b2f9818f-f434-4a2f-9db3-67ea2f229105
name1_p5 = "CEMP_vs_noCEMP_molec_abund_t57.50g45.00m-5.000.pdf"

# ╔═╡ bb3097f0-048c-477a-ad54-890c705437ff
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(design_p5A)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_abund_3D(
				s, "CH", 
				ax=ax
			)
			
			x, y = plot_abund_av_3D(
				s, "C_I", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ I", 
				color="tomato", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "C_II", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ II", 
				color="steelblue", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "CH", 
				ax=ax, 
				label=p.label*L"\rm\ -\ CH", 
				color="k", ls=p.ls, lw=p.lw*1.5
			)
		else
			s = spectra1D[k]
			x, y = plot_abund_av_3D(
				s, "C_I", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ I", 
				color="tomato", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "C_II", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ II", 
				color="steelblue", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "CH", 
				ax=ax, 
				label=p.label*L"\rm\ -\ CH", 
				color="k", ls=p.ls, lw=p.lw*1.5
			)
		end
		
	end

	ax.set_xlim(-5, 4)
	ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("[X/H]")

	ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p5)
	f
end

# ╔═╡ 59457f84-e696-4cb1-914f-65305102a4fb


# ╔═╡ 804fbe02-c7d1-4bf9-91da-9cd3ef441e46
md"## 3D metallicty effect"

# ╔═╡ 452e4faa-57be-4d74-a7c5-36673c21aa9b
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ fdcf50ee-af07-466f-9382-2a3af5de5a25
design_p5 = [
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "3D", "-", 2, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", "1D", "--", 2, "steelblue")
]

# ╔═╡ 3a5bad17-f447-43bb-9d9a-1679b4acf9e6
name2_p5 = "vs_1D_molec_abund_t57.50g45.00m-5.000.pdf"

# ╔═╡ 764090a5-8132-4310-9932-fe80e25dbce5
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(design_p5)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_abund_3D(
				s, "CH", 
				ax=ax
			)
			
			x, y = plot_abund_av_3D(
				s, "C_I", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ I", 
				color="tomato", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "C_II", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ II", 
				color="steelblue", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "CH", 
				ax=ax, 
				label=p.label*L"\rm\ -\ CH", 
				color="k", ls=p.ls, lw=p.lw*1.5
			)
		else
			s = spectra1D[k]
			x, y = plot_abund_av_3D(
				s, "C_I", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ I", 
				color="tomato", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "C_II", 
				ax=ax, 
				label=p.label*L"\rm\ -\ C\ II", 
				color="steelblue", ls=p.ls, lw=p.lw
			)
			x, y = plot_abund_av_3D(
				s, "CH", 
				ax=ax, 
				label=p.label*L"\rm\ -\ CH", 
				color="k", ls=p.ls, lw=p.lw*1.5
			)
		end
		
	end

	ax.set_xlim(-5, 4)
	ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("[X/H]")

	ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p5)

	f
end

# ╔═╡ 92442a47-82f9-465a-9dc0-89ca9f1dddea


# ╔═╡ f7d263c5-ddd6-4d94-945d-303a8c81b8c6


# ╔═╡ ec5b18a4-70f9-4078-8d64-7d3f52dbb7b8
md"# Curve of growth"

# ╔═╡ 6fd82808-3881-4cdf-83d1-f2095cb553ac
abundance_corrections = Dict()

# ╔═╡ ec97726a-ad08-454a-874c-e5ddf54a3e85


# ╔═╡ 483897d7-fbf3-4a3d-a32b-2992ce5fccd2
md"## [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 7021b190-4a8b-4ae4-bd71-06b3840f20c6
abundances_p6 = Dict(
	2.0 => PlotDesign("t57.50g45.00m-5.000 CEMP - C2.0", "", "-", 1.5, nothing),
	2.2 => PlotDesign("t57.50g45.00m-5.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.4 => PlotDesign("t57.50g45.00m-5.000 CEMP - C2.4", "", "-", 1.5, nothing),
	2.6 => PlotDesign("t57.50g45.00m-5.000 CEMP - C2.6", "", "-", 1.5, nothing),
	2.8 => PlotDesign("t57.50g45.00m-5.000 CEMP - C2.8", "", "-", 1.5, nothing),
	#3.0 => PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "", "-", 1.5, nothing),
	3.2 => PlotDesign("t57.50g45.00m-5.000 CEMP - C3.2", "", "-", 1.5, nothing)
)

# ╔═╡ b7a0779f-bd55-4ec3-afc8-ea8caf9db6fc
reference_p6 = [
	(3.0, PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(3.0, PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ d9818823-8afb-4ad2-876d-24ff473a93e7
name1_p6 = "cog_t57.50g45.00m-5.000.pdf"

# ╔═╡ f3f50db8-66b7-4d48-b067-bea39c9c934c
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	abund = abund[sortmask]
	ew = ew[sortmask]
	f_ab = MUST.linear_interpolation(ew, abund, extrapolation_bc=MUST.Line())
	f_ew = MUST.linear_interpolation(abund, ew, extrapolation_bc=MUST.Line())

	x = range(1.95, 3.35, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.07*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color, 
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5750.0, 4.5, -5.0)] = corrs
	
	#ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(1.96, 3.25)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p6)
	
	f
end

# ╔═╡ e4d71c15-0b34-4e99-a738-d41892f93e00


# ╔═╡ 6deaf7e4-7d5a-41d2-9f8d-e024d8150506
abundances_p6_2 = Dict(
	2.0 => PlotDesign("t52.50g30.00m-5.000 CEMP - C2.0", "", "-", 1.5, nothing),
	2.2 => PlotDesign("t52.50g30.00m-5.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.4 => PlotDesign("t52.50g30.00m-5.000 CEMP - C2.4", "", "-", 1.5, nothing),
	2.6 => PlotDesign("t52.50g30.00m-5.000 CEMP - C2.6", "", "-", 1.5, nothing),
	2.8 => PlotDesign("t52.50g30.00m-5.000 CEMP - C2.8", "", "-", 1.5, nothing),
	3.0 => PlotDesign("t52.50g30.00m-5.000 CEMP - C3.0", "", "-", 1.5, nothing),
	3.2 => PlotDesign("t52.50g30.00m-5.000 CEMP - C3.2", "", "-", 1.5, nothing),
	3.4 => PlotDesign("t52.50g30.00m-5.000 CEMP - C3.4", "", "-", 1.5, nothing)
)

# ╔═╡ 82013231-d110-4791-acc4-46bf165bec09
reference_p6_2 = [
	(3.0, PlotDesign("t52.50g30.00m-5.000 no-CEMP - C3.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(3.0, PlotDesign("MARCS t52.50g30.00m-5.000 no-CEMP - C3.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ f2d61607-6728-4520-a3ef-e1b4b29b1fbb
name1_p6_2 = "cog_t52.50g30.00m-5.000.pdf"

# ╔═╡ a8abb2a7-687c-45a4-9c47-2273a181b6fb
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	abund = abund[sortmask]
	ew = ew[sortmask]
	f_ab = MUST.linear_interpolation(ew, abund, extrapolation_bc=MUST.Line())
	f_ew = MUST.linear_interpolation(abund, ew, extrapolation_bc=MUST.Line())

	x = range(1.95, 3.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.07*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color, 
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5250.0, 3.0, -5.0)] = corrs
	
	
	#ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(1.96, 3.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p6_2)
	
	f
end

# ╔═╡ fc260ef0-2b8b-415f-9e16-de551b281d43


# ╔═╡ eb6f1896-6945-459f-859f-2594c76de6d8
abundances_p6_3 = Dict(
	2.0 => PlotDesign("t50.00g25.00m-5.000 CEMP - C2.0", "", "-", 1.5, nothing),
	2.2 => PlotDesign("t50.00g25.00m-5.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.4 => PlotDesign("t50.00g25.00m-5.000 CEMP - C2.4", "", "-", 1.5, nothing),
	2.6 => PlotDesign("t50.00g25.00m-5.000 CEMP - C2.6", "", "-", 1.5, nothing),
	2.8 => PlotDesign("t50.00g25.00m-5.000 CEMP - C2.8", "", "-", 1.5, nothing),
	3.0 => PlotDesign("t50.00g25.00m-5.000 CEMP - C3.0", "", "-", 1.5, nothing),
	3.2 => PlotDesign("t50.00g25.00m-5.000 CEMP - C3.2", "", "-", 1.5, nothing),
	3.4 => PlotDesign("t50.00g25.00m-5.000 CEMP - C3.4", "", "-", 1.5, nothing)
)

# ╔═╡ b742f64a-b6fe-4c6e-9368-9c4e60fc702c
reference_p6_3 = [
	(3.0, PlotDesign("t50.00g25.00m-5.000 no-CEMP - C3.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(3.0, PlotDesign("MARCS t50.00g25.00m-5.000 no-CEMP - C3.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ 4ddd121e-8b2f-4da4-9474-5cc9894fc43c
name1_p6_3 = "cog_t50.00g25.00m-5.000.pdf"

# ╔═╡ 2dae804b-ef8e-4f39-b6f6-011bb92f0384
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	abund = abund[sortmask]
	ew = ew[sortmask]
	f_ab = MUST.linear_interpolation(ew, abund, extrapolation_bc=MUST.Line())
	f_ew = MUST.linear_interpolation(abund, ew, extrapolation_bc=MUST.Line())

	x = range(1.95, 3.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.07*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color, 
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5000.0, 2.5, -5.0)] = corrs
	
	
	#ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(1.96, 3.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p6_3)
	
	f
end

# ╔═╡ 38ce6bb7-f6f2-44d3-8685-8049c42ffeee


# ╔═╡ 2cd6317a-97b4-45bb-8e87-b09b40a470f7
md"## [Fe/H]=-4, [C/Fe]=2.0"

# ╔═╡ de419d2f-eb6d-47cd-a69c-99d5cb3161a3
abundances_p6B = Dict(
	2.2 => PlotDesign("t57.50g45.00m-4.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.0 => PlotDesign("t57.50g45.00m-4.000 CEMP - C2.0", "", "-", 1.5, nothing),
	1.8 => PlotDesign("t57.50g45.00m-4.000 CEMP - C1.8", "", "-", 1.5, nothing),
	1.6 => PlotDesign("t57.50g45.00m-4.000 CEMP - C1.6", "", "-", 1.5, nothing),
	1.4 => PlotDesign("t57.50g45.00m-4.000 CEMP - C1.4", "", "-", 1.5, nothing),
	1.2 => PlotDesign("t57.50g45.00m-4.000 CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t57.50g45.00m-4.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ 135f931b-387f-420d-8a0e-d7f391c1569f
reference_p6B = [
	(2.0, PlotDesign("t57.50g45.00m-4.000 no-CEMP - C2.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(2.0, PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C2.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ 5f9b9dc4-f43b-4853-85d4-ae51247c13d0
name2_p6 = "cog_t57.50g45.00m-4.000.pdf"

# ╔═╡ b0dfbe89-c719-4622-8f7b-55e5a09590df
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6B
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.7, 2.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6B
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5750.0, 4.5, -4.0)] = corrs
	
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.7, 2.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p6)
	
	f
end

# ╔═╡ d165ca55-95ab-4702-b194-83d123f9d333


# ╔═╡ b15fc564-a58d-49f6-9cf6-6587aeb61e1c
abundances_p6B_2 = Dict(
	2.4 => PlotDesign("t52.50g30.00m-4.000 CEMP - C2.4", "", "-", 1.5, nothing),
	2.2 => PlotDesign("t52.50g30.00m-4.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.0 => PlotDesign("t52.50g30.00m-4.000 CEMP - C2.0", "", "-", 1.5, nothing),
	1.8 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.8", "", "-", 1.5, nothing),
	1.6 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.6", "", "-", 1.5, nothing),
	1.4 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.4", "", "-", 1.5, nothing),
	1.2 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.2", "", "-", 1.5, nothing),
	#1.0 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ 99cc258f-19df-4ba8-8947-60a730db5f39
reference_p6B_2 = [
	(2.0, PlotDesign("t52.50g30.00m-4.000 no-CEMP - C2.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(2.0, PlotDesign("MARCS t52.50g30.00m-4.000 no-CEMP - C2.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ c1eb200f-dd1b-40b2-8f1a-ea53908778db
name2_p6_2 = "cog_t52.50g30.00m-4.000.pdf"

# ╔═╡ bf484955-9cdb-4e12-99a0-bb744e977e9b
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6B_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(1.1, 2.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6B_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5250.0, 3.0, -4.0)] = corrs
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(1.1, 2.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p6_2)
	
	f
end

# ╔═╡ ee74de1f-4f7f-4c35-afb8-1186813b8a7b


# ╔═╡ 5ce8c3a7-d12e-4a3b-a3ed-81dbcd3a6118
abundances_p6B_3 = Dict(
	2.4 => PlotDesign("t50.00g25.00m-4.000 CEMP - C2.4", "", "-", 1.5, nothing),
	2.2 => PlotDesign("t50.00g25.00m-4.000 CEMP - C2.2", "", "-", 1.5, nothing),
	2.0 => PlotDesign("t50.00g25.00m-4.000 CEMP - C2.0", "", "-", 1.5, nothing),
	1.8 => PlotDesign("t50.00g25.00m-4.000 CEMP - C1.8", "", "-", 1.5, nothing),
	1.6 => PlotDesign("t50.00g25.00m-4.000 CEMP - C1.6", "", "-", 1.5, nothing),
	1.4 => PlotDesign("t50.00g25.00m-4.000 CEMP - C1.4", "", "-", 1.5, nothing),
	1.2 => PlotDesign("t50.00g25.00m-4.000 CEMP - C1.2", "", "-", 1.5, nothing),
	#1.0 => PlotDesign("t50.00g25.00m-4.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ 4df7f92c-db30-4261-b757-86f9af54c6b3
reference_p6B_3 = [
	(2.0, PlotDesign("t50.00g25.00m-4.000 no-CEMP - C2.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(2.0, PlotDesign("MARCS t50.00g25.00m-4.000 no-CEMP - C2.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ ad4dcf0b-89a9-44f9-9de9-185fce79a510
name2_p6_3 = "cog_t50.00g25.00m-4.000.pdf"

# ╔═╡ f7f9e7fd-dafc-40f3-8ce6-4af55d2a89ac
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6B_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(1.1, 2.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6B_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5000.0, 2.5, -4.0)] = [corrs[1], NaN]
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(1.1, 2.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p6_3)
	
	f
end

# ╔═╡ e7f0e3ea-3204-4c2f-a464-e270baaf8db1


# ╔═╡ 1ead0a28-b327-476e-8be8-0e84ac546d2f
md"## [Fe/H]=-3, [C/Fe]=1.0"

# ╔═╡ 44717756-9575-47b7-a3e5-59c271f79384
abundances_p6C = Dict(
	1.2 => PlotDesign("t57.50g45.00m-3.000 CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t57.50g45.00m-3.000 CEMP - C1.0", "", "-", 1.5, nothing),
	0.8 => PlotDesign("t57.50g45.00m-3.000 CEMP - C0.8", "", "-", 1.5, nothing),
	0.6 => PlotDesign("t57.50g45.00m-3.000 CEMP - C0.6", "", "-", 1.5, nothing),
	0.4 => PlotDesign("t57.50g45.00m-3.000 CEMP - C0.4", "", "-", 1.5, nothing),
	0.2 => PlotDesign("t57.50g45.00m-3.000 CEMP - C0.2", "", "-", 1.5, nothing),
	0.0 => PlotDesign("t57.50g45.00m-3.000 CEMP - C0.0", "", "-", 1.5, nothing),
)

# ╔═╡ bc158515-fd34-412f-b335-3934ebea6d02
reference_p6C = [
	(1.0, PlotDesign("t57.50g45.00m-3.000 no-CEMP - C1.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(1.0, PlotDesign("MARCS t57.50g45.00m-3.000 no-CEMP - C1.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ 21f23414-0891-4aee-ab2b-c2f3c156dd2f
name3_p6 = "cog_t57.50g45.00m-3.000.pdf"

# ╔═╡ 2bdfa4d9-3c0e-40c2-950e-61b7860d529d
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6C
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.0, 1.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6C
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5750.0, 4.5, -3.0)] = corrs
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.0, 1.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p6)
	
	f
end

# ╔═╡ 060d0d6b-f52e-4a5a-acec-858f0ea1384d
abundances_p6C_2 = Dict(
	1.2 => PlotDesign("t52.50g30.00m-3.000 CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t52.50g30.00m-3.000 CEMP - C1.0", "", "-", 1.5, nothing),
	#0.8 => PlotDesign("t52.50g30.00m-3.000 CEMP - C0.8", "", "-", 1.5, nothing),
	0.6 => PlotDesign("t52.50g30.00m-3.000 CEMP - C0.6", "", "-", 1.5, nothing),
	0.4 => PlotDesign("t52.50g30.00m-3.000 CEMP - C0.4", "", "-", 1.5, nothing),
	0.2 => PlotDesign("t52.50g30.00m-3.000 CEMP - C0.2", "", "-", 1.5, nothing),
	#1.0 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ 859af4d4-1e06-4567-ae2d-2d4083da2891
reference_p6C_2 = [
	(1.0, PlotDesign("t52.50g30.00m-3.000 no-CEMP - C1.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(1.0, PlotDesign("MARCS t52.50g30.00m-3.000 no-CEMP - C1.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ de377555-e268-4185-82d8-588b0ee36b38
name3_p6_2 = "cog_t52.50g30.00m-3.000.pdf"

# ╔═╡ 0140b523-2d02-4522-aab4-201c0014ae05
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6C_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.0, 1.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6C_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5250.0, 3.0, -3.0)] = corrs
	
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.0, 1.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p6_2)
	
	f
end

# ╔═╡ 659cb51a-06ea-432a-91e2-9bd3f5006eb4


# ╔═╡ c2e9d134-6bab-454c-b472-9b33965e031b
abundances_p6C_3 = Dict(
	1.2 => PlotDesign("t50.00g25.00m-3.000 CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t50.00g25.00m-3.000 CEMP - C1.0", "", "-", 1.5, nothing),
	#0.8 => PlotDesign("t50.00g25.00m-3.000 CEMP - C0.8", "", "-", 1.5, nothing),
	0.6 => PlotDesign("t50.00g25.00m-3.000 CEMP - C0.6", "", "-", 1.5, nothing),
	0.4 => PlotDesign("t50.00g25.00m-3.000 CEMP - C0.4", "", "-", 1.5, nothing),
	0.2 => PlotDesign("t50.00g25.00m-3.000 CEMP - C0.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t50.00g25.00m-3.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ b0ae8644-36f4-4f09-8602-9a36ab5ea08f
reference_p6C_3 = [
	(1.0, PlotDesign("t50.00g25.00m-3.000 no-CEMP - C1.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(1.0, PlotDesign("MARCS t50.00g25.00m-3.000 no-CEMP - C1.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ d779ed13-aa80-4844-adcb-4662f59ad0c6
name3_p6_3 = "cog_t50.00g25.00m-3.000.pdf"

# ╔═╡ 298ab49e-f1a3-431f-a9c1-6c5bcce6771b
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6C_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.0, 1.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6C_3
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5000.0, 2.5, -3.0)] = corrs
	
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.0, 1.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name3_p6_3)
	
	f
end

# ╔═╡ 2c44e5cc-43e9-4181-ab36-d629f4c61080


# ╔═╡ 098be297-ba34-4648-a08b-935db2bf5b7f
md"## [Fe/H]=-2, [C/Fe]=1.0"

# ╔═╡ 63a51e50-4cb5-4df2-92ab-c24eeb23e482
abundances_p6D = Dict(
	1.2 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C1.0", "", "-", 1.5, nothing),
	0.8 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C0.8", "", "-", 1.5, nothing),
	0.6 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C0.6", "", "-", 1.5, nothing),
	0.4 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C0.4", "", "-", 1.5, nothing),
	0.2 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C0.2", "", "-", 1.5, nothing),
	0.0 => PlotDesign("t57.50g45.00m-2.000 no-CEMP - C0.0", "", "-", 1.5, nothing),
)

# ╔═╡ 96980f1d-4e4c-48c9-8fe6-080cb9eb3ef8
reference_p6D = [
	(1.0, PlotDesign("t57.50g45.00m-2.000 no-CEMP - C1.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(1.0, PlotDesign("MARCS t57.50g45.00m-2.000 no-CEMP - C1.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ 5c6e6453-5911-45ac-9fc0-1d69d92516fd
name4_p6 = "cog_t57.50g45.00m-2.000.pdf"

# ╔═╡ 693a142f-f26f-4543-8751-54dd8425a6b9
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6D
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.0, 1.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6D
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5750.0, 4.5, -2.0)] = corrs
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.0, 1.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name4_p6)
	
	f
end

# ╔═╡ ca64a419-dc88-45e0-904e-535cfcd94745
abundances_p6D_2 = Dict(
	1.2 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C1.2", "", "-", 1.5, nothing),
	1.0 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C1.0", "", "-", 1.5, nothing),
	0.8 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C0.8", "", "-", 1.5, nothing),
	0.6 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C0.6", "", "-", 1.5, nothing),
	0.4 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C0.4", "", "-", 1.5, nothing),
	0.2 => PlotDesign("t52.50g30.00m-2.000 no-CEMP - C0.2", "", "-", 1.5, nothing),
	#1.0 => PlotDesign("t52.50g30.00m-4.000 CEMP - C1.0", "", "-", 1.5, nothing),
)

# ╔═╡ b050a849-317f-4c2e-bc05-f79ff513dea4
reference_p6D_2 = [
	(1.0, PlotDesign("t52.50g30.00m-2.000 no-CEMP - C1.0", "non-CEMP (3D)", "s", 8.5, "tomato")),
	(1.0, PlotDesign("MARCS t52.50g30.00m-2.000 no-CEMP - C1.0", "MARCS", "s", 8.5, "steelblue"))
]

# ╔═╡ 0d028a71-18f0-4c6a-bd1b-75c230a1cd08
name4_p6_2 = "cog_t52.50g30.00m-2.000.pdf"

# ╔═╡ 02f3d4b0-8435-4bab-9218-0ba803df0faf
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	ew = []
	abund = []
	for (ab, p) in abundances_p6D_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		append!(ew, [sum(w .* (1 .- F))])
		append!(abund, [ab])
	end

	sortmask = sortperm(abund)
	sortmask2 = sortperm(ew)
	f_ab = MUST.linear_interpolation(
		ew[sortmask2], abund[sortmask2], extrapolation_bc=MUST.Line()
	)
	f_ew = MUST.linear_interpolation(
		abund[sortmask], ew[sortmask], extrapolation_bc=MUST.Line()
	)

	abund = abund[sortmask]
	ew = ew[sortmask]

	x = range(0.0, 1.4, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")

	corrs = []
	for (aref, p) in reference_p6D_2
		k = p.key

		s = if haskey(spectra3D, k)
			spectra3D[k]
		else
			spectra1D[k]
		end
		λ, F = flux(s, norm=true)
	
		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))

		ΔC = round(f_ab(EW) - aref, sigdigits=2)
		append!(corrs, [ΔC])
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="right", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end

	abundance_corrections[(5250.0, 3.0, -2.0)] = corrs
	
	
	
	#ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.0, 1.4)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name4_p6_2)
	
	f
end

# ╔═╡ 26108373-70d1-42ea-afea-dbd1586c0798


# ╔═╡ 27a3ef9a-d39a-4df5-8644-794a9196b1be
md"## Summary"

# ╔═╡ f8ca6c4c-0de3-4ae8-b860-9ef17e04f2ad
let
	@info "Abundance corrections"
	ks = keys(abundance_corrections) |> collect
	sortmask = sortperm(ks, by=x->(x[3], x[2], x[1]))
	for k in ks[sortmask]
		Δ_CEMP_nonCEMP = first(abundance_corrections[k])
		Δ_3D_1D = last(abundance_corrections[k])
		@info "Teff: $(k[1]), log(g): $(k[2]), [Fe/H]:$(k[3])" Δ_CEMP_nonCEMP Δ_3D_1D
	end
end

# ╔═╡ d8117e26-d1e4-48ed-a543-dbabc7b7cfc8


# ╔═╡ 579ec7d6-9535-4b71-9c9f-5c19d21d3da9


# ╔═╡ 517d0552-7f2a-4f70-8e1f-fae0cd8cd05e
md"## Correction function
We define a correction function based on Teff and log(g) of the models in the grid. This function simply takes whatever correctios are available and picks the closest one in metallicity. In Teff/logg it interpolates, or picks the closes whatever is more convenient."

# ╔═╡ 7714d6d3-fd29-47ef-a52f-621ec02800ff
"""
	Finds the closest point in a set of points to a given test point.

Args:
	points: An array of points, where each point is a tuple (teff, logg).
    test_point: The test point, as a tuple (teff, logg).

Returns:
    The closest point to the test point.
"""
function closest_point(points, test_point)
  # Normalize the points to account for different magnitudes
  teff_max = maximum(p[1] for p in points)
  logg_max = maximum(p[2] for p in points)
  normalized_points = [(p[1]/teff_max, p[2]/logg_max) for p in points]
  normalized_test_point = (test_point[1]/teff_max, test_point[2]/logg_max)

  # Calculate the distances between the normalized points and the normalized test point
  distances = [sqrt((p[1] - normalized_test_point[1])^2 + (p[2] - normalized_test_point[2])^2) for p in normalized_points]

  # Find the index of the closest point
  closest_index = argmin(distances)

  # Return the original (unnormalized) closest point
  return points[closest_index]
end

# ╔═╡ 641f1109-17a5-4401-bc27-70fc4db7583e
"""
	metallicity_interpolator(corrections=abundance_corrections)

Compute the metallicity interpolation for every point in the grid of models.
"""
function metallicity_interpolator(corrections=abundance_corrections; which=2, skip=nothing)
	ks = keys(corrections) |> collect
	
	# for all the points in the grid, compute the metallicity interpolator
	interpolators = Dict()

	for key in ks
		knew = (key[1], key[2])
		if !isnothing(skip)
			(knew in skip) && continue
		end
		if !(knew in keys(interpolators))
			metalicities = sort([m[3] for m in ks if (m[1]==key[1])&(m[2]==key[2])])
			corrs = [corrections[(knew...,m)][which] for m in metalicities]
			nanfilter = .!isnan.(corrs)
			@info "Interpolate $(knew):" metalicities corrs
			interpolators[knew] = MUST.linear_interpolation(
				metalicities[nanfilter],
				corrs[nanfilter],
				extrapolation_bc=MUST.Flat()
			)
		end
	end

	interpolators
end

# ╔═╡ eb865dfb-b206-4daa-a1e7-d8c706ba2485
interpolators = metallicity_interpolator(skip=[(5000.0, 2.5)])

# ╔═╡ a206fa91-09cd-4faa-b8ae-ccda70d36c16
interpolators_carbon = metallicity_interpolator(which=1, skip=[(5000.0, 2.5)])

# ╔═╡ 4b8de50b-b68a-458b-8408-1ded9b10976c


# ╔═╡ 2d3423ee-7a2f-4cc9-8800-267caf418034
"""
	interpolate_corrections(teff, logg, feh; feh_interpolators=interpolators)

Interpolate the corrections in metallicity based on the closest point in the Teff, logg plane.
"""
function interpolate_corrections(teff, logg, feh; feh_interpolators=interpolators)
	ks = keys(feh_interpolators) |> collect
	logg_gr = getindex.(ks, 2)
	teff_gr = getindex.(ks, 1)

	# closest point teff-logg plane
	tefflogg = closest_point(zip(teff_gr,logg_gr)|>collect, (teff, logg))

	# interpolate the corrections in metallicity
	feh_interpolators[tefflogg](feh)
end

# ╔═╡ ef746dad-2b45-42c8-beb4-02ca0f4d983d


# ╔═╡ ae657bc5-ad83-4675-bdda-f578d7054af7
design_pCorr = [
	PlotDesign((5750.0, 4.5), L"\rm T_{eff} = 5750K"*"\n"*L"\rm log(g)=4.5", "-", 1.5, "k"),
	PlotDesign((5250.0, 3.0), L"\rm T_{eff} = 5250K"*"\n"*L"\rm log(g)=3.0", "-", 1.5, "tomato"),
	#PlotDesign((5000.0, 2.5), L"\rm T_{eff} = 5000K,\ log(g)=2.5", "-", 1.5, "steelblue"),
]

# ╔═╡ 24c0eb4f-5dcf-48bf-b305-2e7737c1828b
design_pCorr_c = [
	PlotDesign((5750.0, 4.5),  L"\rm T_{eff} = 5750K"*"\n"*L"\rm log(g)=4.5", "-", 1.5, "k"),
	PlotDesign((5250.0, 3.0), L"\rm T_{eff} = 5250K"*"\n"*L"\rm log(g)=3.0", "-", 1.5, "tomato"),
	#PlotDesign((5000.0, 2.5), L"\rm T_{eff} = 5000K,\ log(g)=2.5", "--", 1.5, "steelblue"),
]

# ╔═╡ 6a495bfa-fc58-4cbd-89a7-1a8e198a2e91
name1_p62 = "abundance_correction.pdf"

# ╔═╡ b768e044-d7fc-49a9-bc5e-1859c359f0d1
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(7, 5), sharex=true, sharey=true)
	plt.subplots_adjust(wspace=0)
	
	ax[0].axhline(0.0, color="0.5", lw=1, ls=":")
	ax[1].axhline(0.0, color="0.5", lw=1, ls=":")

	for p in design_pCorr
		key = p.key
		metalicities = sort([
			m[3] for m in keys(abundance_corrections) if (m[1]==key[1])&(m[2]==key[2])
		])
		corrs = [abundance_corrections[(key...,m)][2] for m in metalicities]
		nanfilter = .!isnan.(corrs)
		#=ax[0].plot(
			metalicities[nanfilter], corrs[nanfilter], marker="s", ls=p.ls, lw=p.lw, markersize=9, color=p.color, label=p.label, markerfacecolor="w", markeredgecolor=p.color, markeredgewidth=2.5
		)=#
		ax[0].plot(
			metalicities[nanfilter], corrs[nanfilter], marker="s", ls=p.ls, lw=p.lw, markersize=9, color=p.color, label=p.label
		)
	end

	for p in design_pCorr_c
		key = p.key
		metalicities = sort([
			m[3] for m in keys(abundance_corrections) if (m[1]==key[1])&(m[2]==key[2])
		])
		corrs = [abundance_corrections[(key...,m)][1] for m in metalicities]

		ax[1].plot(
			metalicities, corrs, marker="s", ls=p.ls, lw=p.lw, markersize=9, color=p.color, label=p.label
		)
	end

	#ax[0].legend(ncol=1, loc="upper center")
	ax[1].legend(ncol=1, loc="lower center")
	ax[0].set_xlabel(L"\rm [Fe/H]\ [dex]")
	ax[1].set_xlabel(L"\rm [Fe/H]\ [dex]")
	ax[0].set_ylabel(L"\rm \Delta\ A(C)\ [dex]")
	#ax[1].set_ylabel(L"\rm \Delta\ A(C)\ [dex]")
	ax[0].set_ylim(-0.77, 0.47)

	ax[0].set_title(L"\rm 3D\ CEMP - 1D")
	ax[1].set_title(L"\rm 3D\ CEMP - 3D\ non\ CEMP")

	f.savefig(name1_p62)

	f
end

# ╔═╡ 05c3da3f-6654-4ff4-953e-9f6d6c315f16


# ╔═╡ 349d1a5b-aed0-4315-b7b4-cc766e2d0bc3
md"# Yoon-Beers 2019"

# ╔═╡ 8cc66430-8635-4729-b595-a7df6a28cba5
yb19 = MUST.readdlm("../yoon_beers19.txt", ';')

# ╔═╡ 3bea6988-24ec-46bf-b20d-1f9817e7cc0b
feh = yb19[4:end, 6]

# ╔═╡ ac78dfbc-1423-4741-b0b4-dcd8c809787f
c = yb19[4:end, 8]

# ╔═╡ 90f2ca39-e67f-477b-ad1b-ff48c578174f


# ╔═╡ 009e81c3-7b8c-433a-a931-982f9a406d94
md"Use the corrections computed above to correct abundances."

# ╔═╡ b5896895-1bd7-40de-b86b-b326d502e587
begin
	models_t57g45 = sort([k for k in keys(abundance_corrections) if (k[1]==5750.0) && (k[2]==4.5)], by=x->(x[3], x[2], x[1]))

	models_t52g30 = sort([k for k in keys(abundance_corrections) if (k[1]==5250.0) && (k[2]==3.0)], by=x->(x[3], x[2], x[1]))
	
	ip_t57g45 = MUST.linear_interpolation(
		last.(models_t57g45),
		[abundance_corrections[m][2] for m in models_t57g45],
		extrapolation_bc=MUST.Flat()
	)

	ip_t52g30 = MUST.linear_interpolation(
		last.(models_t52g30),
		[abundance_corrections[m][2] for m in models_t52g30],
		extrapolation_bc=MUST.Flat()
	)

	corrections = ip_t57g45.(feh)
end

# ╔═╡ 5c5bf31a-b1c4-478c-b2ae-6f42957d0c19
name1_p7 = "yoon_beers_2019_corrected.pdf"

# ╔═╡ 103de55b-aea6-44eb-b3c1-e733a8c4c6e3
let
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	
	ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.8), -6, 0.8, color="orange", alpha=0.3
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.7, 6.0), -4, 0.8, color="green", alpha=0.3, angle=45
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.6, 7.9), -3, 2, color="blue", alpha=0.3, angle=10
	)
	ax.add_patch(ellipse1)
	ax.add_patch(ellipse2)
	ax.add_patch(ellipse3)

	ax.scatter(feh, c, s=10, marker="o", c="k", label="Yoon, Beers at al. 2019")
	#ax.scatter(feh, c+corrections, s=10, marker="o", c="tomato", label="corrected")

	cemp_line(x) = 0.75 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="k", alpha=0.5, ls="--")

	ax.set_xlim(-8.5, -1)
	ax.set_ylim(4, 10)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(loc="upper left")
	
	f.savefig(name1_p7)

	f
end

# ╔═╡ d0c31252-9b07-47a8-a22c-25dcfca07d2a


# ╔═╡ 5f76e7c0-7366-4d1b-b88e-e9e3639d44de
md"# Contribution functions"

# ╔═╡ dc0509a9-a1bb-438f-8250-e52a56ffd1e4
begin
	"""
		plot_contr_av(s; ax, kwargs...)
	
	Plot contribution function (averages) as a function of optical depth (500).
	"""
	plot_contr_av(s; ax, kwargs...) = begin
		krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
		line = s.run.atom.trans[krs[1]]
		x, y = line.get_cntrbf(gmean=true, df=0, norm=true)
	
		x = MUST.pyconvert(Array, x)
		y = MUST.pyconvert(Array, y)
	
		ax.plot(x, y; kwargs...)
	end
	
	"""
		plot_contr_1D(s; ax, kwargs...)
	
	Plot contribution function (averages) as a function of optical depth (500).
	"""
	plot_contr_1D(s; ax, kwargs...) = begin
		krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
		line = s.run.atom.trans[krs[1]]
		x, y = line.get_cntrbf(df=0, norm=true)

		x = MUST.pyconvert(Array, x)
		y = MUST.pyconvert(Array, y)
	
		ax.plot(x, y; kwargs...)
	end
end

# ╔═╡ 1773e7cb-bdff-4a97-be64-e0aa8284225a
"""
	plot_abund_av_3D(s, name; ax, kwargs...)

Plot contribution function at the CH surface.
"""
plot_contr_3D(s; ax, kwargs...) = begin
	krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
	line = s.run.atom.trans[krs[1]]
	xtau, ycontr = line.get_cntrbf(df=0)

	x = MUST.pyconvert(Array, s.run.xx)
	y = MUST.pyconvert(Array, s.run.yy)
	z = MUST.pyconvert(Array, s.run.zz)
	data = Dict(
		:τ500 => exp10.(MUST.pyconvert(Array, s.run.ltau)),
		:cont => MUST.pyconvert(Array, ycontr),
		:T => MUST.pyconvert(Array, s.run.temp),
		:d => MUST.pyconvert(Array, s.run.rho),
	)

	data[:cont] .= data[:cont] ./ maximum(data[:cont])

	xx, yy, zz = MUST.meshgrid(x, y, z)
	b = MUST.Box(xx, yy, zz, data, MUST.AtmosphericParameters())
	bsurf = MUST.interpolate_to(b, :cont, τ500=-3.5, logspace=true)

	cplot = bsurf[:cont][:,:,1]

	im = ax.pcolormesh(xx[:,:,1], yy[:,:,1], cplot, rasterized=true, cmap="gist_heat")
end

# ╔═╡ 95ad457e-9d59-45a4-b983-5e2e2b0f6c76
"""
	plot_contr_yslice(s, yslice; ax, kwargs...)

Plot contribution function in yslices.
"""
plot_contr_yslice(s, yslice; ax, kwargs...) = begin
	krs = MUST.pyconvert(Array, s.run.nml["line_mask"]["cntrbf_lines"]) .-1
	line = s.run.atom.trans[krs[1]]
	xtau, ycontr = line.get_cntrbf(df=0)

	x = MUST.pyconvert(Array, s.run.xx)
	y = MUST.pyconvert(Array, s.run.yy)
	z = MUST.pyconvert(Array, s.run.zz)
	data = Dict(
		:τ500 => exp10.(MUST.pyconvert(Array, s.run.ltau)),
		:cont => MUST.pyconvert(Array, ycontr),
		:T => MUST.pyconvert(Array, s.run.temp),
		:d => MUST.pyconvert(Array, s.run.rho),
	)

	data[:cont] .= data[:cont] ./ maximum(data[:cont])
	xx, yy, zz = MUST.meshgrid(x, y, z)
	b = MUST.Box(xx, yy, zz, data, MUST.AtmosphericParameters())
	bt = MUST.height_scale_fast(b, :τ500, logspace=true, dont_log=[:T, :cont])

	
	cplot = bt[:cont][:,yslice,:] 
	im = ax.pcolormesh(xx[:,yslice,:], log10.(bt[:τ500][:,yslice,:]), cplot, rasterized=true, cmap="gist_heat")
end

# ╔═╡ fe710591-6d00-4dad-9173-d8cf17ffefc1
"""
	plot_profile(s, q; ax, kwargs...)

Plot profile function of optical depth cube (500).
"""
plot_profile(s, q; ax=nothing, kwargs...) = begin
	x = MUST.pyconvert(Array, s.run.xx)
	y = MUST.pyconvert(Array, s.run.yy)
	z = MUST.pyconvert(Array, s.run.zz)
	data = Dict(
		:τ500 => exp10.(MUST.pyconvert(Array, s.run.ltau)),
		:T => MUST.pyconvert(Array, s.run.temp),
		:d => MUST.pyconvert(Array, s.run.rho),
	)

	xx, yy, zz = MUST.meshgrid(x, y, z)
	b = MUST.Box(xx, yy, zz, data, MUST.AtmosphericParameters())
	bt = MUST.height_scale_fast(b, :τ500, logspace=true)

	qs, logq = MUST.is_log(q)
	
	# Background distribution
	h, x, y, p = plt.hist2d(
		reshape(log10.(b[:τ500]), :), reshape(logq.(b[qs]), :), bins=200
	)
	if !isnothing(ax)
		ax.imshow(
			h.T, 
			origin="lower",
			interpolation = "bicubic", 
			extent=[minimum(x), maximum(x), minimum(y), maximum(y)],
			cmap="GnBu", norm=matplotlib.colors.LogNorm(vmin=1), aspect="auto"
		)
	end
	
	#ax.hist2d(reshape(log10.(b[:τ500]), :), reshape(b[:T], :), rasterized=true, cmap="YlGnBu", bins=300, norm=matplotlib.colors.LogNorm())

	x, y = profile(MUST.mean, bt, :log10τ500, q)
end

# ╔═╡ 294fa362-da67-4a22-8a68-8619e185df6c


# ╔═╡ 6addc6c2-ad44-4310-a39a-292acf6b0dac
md"## [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 9962676e-40cb-4fe8-ab50-8abb5a854612
spectra_p8 = [
	PlotDesign("t57.50g45.00m-5.000 CEMP - contr", "3D", "-", 2, "k"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - contr", "1D", "-", 5, "tomato"),
]

# ╔═╡ d5af99b0-8e95-4037-af01-1437f9aa0e04
differences_p8 = [
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ 19c42123-ae32-468b-b6a1-4fd7457aa526
inset_p8 = PlotDesign("t57.50g45.00m-5.000 CEMP", "", "-", 2, "k")

# ╔═╡ 43b4cf33-7efd-48b8-af90-66555c827b52
name1_p8 = "contribution_tau_t57.50g45.00m-5.000.pdf"

# ╔═╡ 71545110-8314-45c7-9e83-e15035c3f60d
name2_p8 = "contribution_surface_t57.50g45.00m-5.000.pdf"

# ╔═╡ 266c3319-78c1-4a2b-b06c-386f54907ecc
name3_p8 = "contribution_height_t57.50g45.00m-5.000.pdf"

# ╔═╡ d1c307dd-ab6f-44bc-a0d5-a879960231a5
name4_p8 = "profile_tau500_T_t57.50g45.00m-5.000.pdf"

# ╔═╡ 57c74ea4-a4ff-4018-9507-2753fc2a440c
name5_p8 = "profile_tau500_rho_t57.50g45.00m-5.000.pdf"

# ╔═╡ 2765daec-fb70-4b55-889f-0ceb0808e84e
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p8)

	f
end

# ╔═╡ 515b8fb8-2568-4200-b10f-6ef8e2318615
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p8)

	f
end

# ╔═╡ d08caaf5-e763-498d-a594-8915fc8eaaf9
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -5.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p8)

	f
end

# ╔═╡ 1dcf4d03-7a20-4ff7-b891-958003c23d61
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)
	ax2.tick_params(labelsize=12)

	k = differences_p8[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p8)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend(fontsize=11, loc="upper center")

	#ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	#ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.hlines(3500, -3.5, 0.0, ls="-", color="k", alpha=0.2, lw=4)
	ax.hlines(4800, -1.2, 0.0, ls="-", color="tomato", alpha=0.2, lw=4)

	#ax.axvline(log10(2/3))
	
	ax.set_xlim(-4.5, 2)
	ax2.set_xlim(-4.5, 2)
	ax.set_ylim(2800, 9900)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm T_{eff}=5\,750\ K,\ log(g)=4.5,\ [Fe/H] = -5,\ [C/Fe] = 3")
	

	f.savefig(name4_p8)

	f
end

# ╔═╡ c8a6af9e-625a-4864-9ca2-1efad14fccd4
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
			
			#ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, ls=p.ls)
		end
	end

	
	#=k = differences_p8[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p8)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end=#

	left, bottom, width, height = [0.11, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	b_in = models3D[inset_p8.key][2]
	Tplane = MUST.interpolate_to(b_in, :T; logspace=true, τ_ross=0.0)[:T][:,:,1]
	extent = [minimum(b_in.x), maximum(b_in.x), minimum(b_in.y), maximum(b_in.y)] ./1e8
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")

	ax2.set_xlabel(L"\rm x\ [Mm]")
	#ax2.set_ylabel(L"\rm y\ [Mm]")
	#ax2.set_title(L"\rm \Delta T\ [K]")
	#ax2.legend(fontsize=11, loc="upper center")

	#ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	#ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	#ax.hlines(3500, -3.5, 0.0, ls="-", color="k", alpha=0.2, lw=4)
	#ax.hlines(4800, -1.2, 0.0, ls="-", color="tomato", alpha=0.2, lw=4)
	ax.hlines(3100, -3.5, 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
	ax.hlines(3450, -1.2, 0.0, ls="-", color="tomato", alpha=0.8, lw=6)

	ax.text(0.1, 3100, "3D", ha="left", va="center", color="steelblue")
	ax.text(0.1, 3450, "1D", ha="left", va="center", color="tomato")

	
	ax.set_xlim(-4.5, 1.7)
	#ax2.set_xlim(-4.5, 2)
	ax.set_ylim(2700, 9900)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm T_{eff}=5\,750\ K,\ log(g)=4.5,\ [Fe/H] = -5,\ [C/Fe] = 3")
	

	f.savefig(name4_p8)

	f
end

# ╔═╡ a21d6f31-f766-432f-856c-7726e989d42f
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p8)

	f
end

# ╔═╡ 67aab606-2b0e-44f8-90ab-99d63830742d


# ╔═╡ 29ff7a48-9110-4da3-b024-2a9c7af59477
spectra_p8B = [
	PlotDesign("t52.50g30.00m-5.000 CEMP - contr", "3D", "-", 2, "k"),
	PlotDesign("MARCS t52.50g30.00m-5.000 no-CEMP - contr", "1D", "-", 5, "tomato"),
]

# ╔═╡ a007184b-f10c-4592-ad0c-aca561a48a45
differences_p8B = [
	PlotDesign("t52.50g30.00m-5.000 CEMP - C3.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t52.50g30.00m-5.000 no-CEMP - C3.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t52.50g30.00m-5.000 no-CEMP - C3.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ 932b8ec8-efcf-4851-8b3c-46fc18d0a674
inset_p8B = PlotDesign("t52.50g30.00m-5.000 CEMP", "", "-", 2, "k")

# ╔═╡ 9ad8b431-b875-4ebd-844a-e70177aa6fde
name1_p8B = "contribution_tau_t52.50g30.00m-5.000.pdf"

# ╔═╡ 4b07347d-8e66-4087-84b7-5efb4076139f
name2_p8B = "contribution_surface_t52.50g30.00m-5.000.pdf"

# ╔═╡ b6ba9f13-2f1d-407a-a496-661f89a86205
name3_p8B = "contribution_height_t52.50g30.00m-5.000.pdf"

# ╔═╡ 4afb0820-64bb-4ef3-aa4f-820621da1452
name4_p8B = "profile_tau500_T_t52.50g30.00m-5.000.pdf"

# ╔═╡ b14074bc-e4ef-4669-8c42-5f6194306029
name5_p8B = "profile_tau500_rho_t52.50g30.00m-5.000.pdf"

# ╔═╡ 5b154765-c22f-45e0-9ebc-b7dfa384ebf7
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			@show k
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p8B)

	f
end

# ╔═╡ ec4b9098-48f3-4967-84fe-fd2d7cc4020e
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p8B)

	f
end

# ╔═╡ 31c72e56-4f33-4328-be62-817a827d2161
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -4.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p8B)

	f
end

# ╔═╡ d8c47e25-67c5-414f-8416-8d26030c9cde
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)
	ax2.tick_params(labelsize=12)
	

	k = differences_p8B[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p8B)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend(fontsize=10, loc="upper center")

	#ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	#ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.hlines(3350, -3.5, 0.0, ls="-", color="k", alpha=0.2, lw=4)
	ax.hlines(4400, -1.2, 0.0, ls="-", color="tomato", alpha=0.2, lw=4)
	
	ax.set_xlim(-4.5, 2)
	ax2.set_xlim(-4.5, 2)
	ax2.set_ylim(-1200, 1500)
	ax.set_ylim(2300, 9900)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm T_{eff}=5\,250\ K,\ log(g)=3.0,\ [Fe/H] = -5,\ [C/Fe] = 3")
	

	f.savefig(name4_p8B)

	f
end

# ╔═╡ bf028bf1-7480-4020-8760-6622d0057188
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
			
			#ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.12, 0.51, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height])
	ax2.tick_params(labelsize=12)
	b_in = models3D[inset_p8B.key][2]
	Tplane = MUST.interpolate_to(b_in, :T; logspace=true, τ_ross=0.0)[:T][:,:,1]
	extent = [minimum(b_in.x), maximum(b_in.x), minimum(b_in.y), maximum(b_in.y)] ./1e8
	im = ax2.imshow(
		Tplane', 
		rasterized=true, 
		origin="lower",
		extent=extent,
		cmap="gist_heat",
		aspect="equal"
	)
	cbar = f.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
	ax2.set_title(L"\rm temperature\ [K]")

	ax2.set_xlabel(L"\rm x\ [Mm]")
	#ax2.set_ylabel(L"\rm y\ [Mm]")
	#ax2.set_title(L"\rm \Delta T\ [K]")
	#ax2.legend(fontsize=11, loc="upper center")

	#ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	#ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	#ax.hlines(3500, -3.5, 0.0, ls="-", color="k", alpha=0.2, lw=4)
	#ax.hlines(4800, -1.2, 0.0, ls="-", color="tomato", alpha=0.2, lw=4)
	ax.hlines(2550, -3.0, 0.0, ls="-", color="steelblue", alpha=0.8, lw=6)
	ax.hlines(2900, -1.3, 0.0, ls="-", color="tomato", alpha=0.8, lw=6)

	ax.text(0.1, 2550, "3D", ha="left", va="center", color="steelblue")
	ax.text(0.1, 2900, "1D", ha="left", va="center", color="tomato")

	
	ax.set_xlim(-4.2, 1.7)
	#ax2.set_xlim(-4.5, 2)
	ax.set_ylim(2100, 9900)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm T_{eff}=5\,250\ K,\ log(g)=3.0,\ [Fe/H] = -5,\ [C/Fe] = 3")	

	f.savefig(name4_p8B)

	f
end

# ╔═╡ 896797bd-8fb0-4127-a2d6-57438ee03fae
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p8B)

	f
end

# ╔═╡ 75480722-5bf5-453d-9f13-690cda6cdc3f


# ╔═╡ 009eabf7-0385-4ed7-b648-1e0b8c4cfc47
spectra_p8C = [
	PlotDesign("t50.00g25.00m-5.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t50.00g25.00m-5.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ aee45896-0f70-4bc8-8c8a-e685d92cb60a
differences_p8C = [
	PlotDesign("t50.00g25.00m-5.000 CEMP - C3.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t50.00g25.00m-5.000 no-CEMP - C3.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t50.00g25.00m-5.000 no-CEMP - C3.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ 9f452aaa-e239-4201-ab9f-21136f8768be
name1_p8C = "contribution_tau_t50.00g25.00m-5.000.pdf"

# ╔═╡ e1d754a5-cb5e-4fa1-9df2-fb6c806b17ea
name2_p8C = "contribution_surface_t50.00g25.00m-5.000.pdf"

# ╔═╡ 1882ff98-7c04-4909-8010-3ca0ed5f6e51
name3_p8C = "contribution_height_t50.00g25.00m-5.000.pdf"

# ╔═╡ 8e757a6a-cfc2-458f-bc10-f31e6606e3ee
name4_p8C = "profile_tau500_T_t50.00g25.00m-5.000.pdf"

# ╔═╡ b1e90b50-628c-49ef-85ea-5691cebd2243
name5_p8C = "profile_tau500_rho_t50.00g25.00m-5.000.pdf"

# ╔═╡ 6d528aab-ecc9-4907-966b-e76250f08a2d
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8C)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			@show k
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p8C)

	f
end

# ╔═╡ de50c004-56b1-4719-9a93-6543d0f75f9b
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8C)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p8C)

	f
end

# ╔═╡ 02e275ca-adf2-4497-a269-c5279bcbc5e9
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8C)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -4.1)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p8C)

	f
end

# ╔═╡ c90d990c-b17a-46d7-a145-00b5eedc0682
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p8C)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

	k = differences_p8C[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p8C)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend()

	ax.axvline(-2.5, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.5, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.1, 3)
	ax.set_ylim(2800, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm [Fe/H] = -5,\ [C/Fe] = 3")
	ax.set_title(L"\rm T_{eff}=5\,000\ K,\ log(g)=2.5,\ [Fe/H] = -5,\ [C/Fe] = 3")
	

	f.savefig(name4_p8C)

	f
end

# ╔═╡ 31d31e74-c494-4ea6-a28c-d2bb7afd57a4
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p8C)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p8C)

	f
end

# ╔═╡ b5cda9d6-480b-42a1-8091-7ceab841c3df


# ╔═╡ 6225b614-dc48-4158-8abe-55a1c07b651e
md"## [Fe/H]=-4, [C/Fe]=2.0"

# ╔═╡ b04303d1-bd59-4104-bc7c-8fdef9477fee
spectra_p9 = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ 901d430e-254b-4774-b7cc-4dcbd8ea1772
differences_p9 = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - C2.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t57.50g45.00m-4.000 no-CEMP - C2.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C2.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ 479b4bca-166c-43c2-8199-04218a24fb12
name1_p9 = "contribution_tau_t57.50g45.00m-4.000.pdf"

# ╔═╡ 7a734c2e-175e-4300-8152-718669ce8137
name2_p9 = "contribution_surface_t57.50g45.00m-4.000.pdf"

# ╔═╡ 77bb0169-c777-425a-a8b5-041cb01279f1
name3_p9 = "contribution_height_t57.50g45.00m-4.000.pdf"

# ╔═╡ 865f7e7a-91a9-4d16-b12d-43a84479e096
name4_p9 = "profile_tau500_T_t57.50g45.00m-4.000.pdf"

# ╔═╡ 1fd35559-c9d6-4be2-bee0-6d43417e4617
name5_p9 = "profile_tau500_rho_t57.50g45.00m-4.000.pdf"

# ╔═╡ 5e7ef044-c782-4ee5-b408-c394118a8f72
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p9)

	f
end

# ╔═╡ 5e23e8f0-cc1a-4912-b841-c9d9f01e9681
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p9)

	f
end

# ╔═╡ 3d1174e7-01ba-4e1d-bbce-cf8d28642ca2
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -5.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p9)

	f
end

# ╔═╡ 576ddad7-10a6-4549-80cc-18d0c6793024
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p9)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

	k = differences_p9[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p9)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend()

	ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(3000, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.set_title(L"\rm [Fe/H] = -4,\ [C/Fe] = 2")
	

	ax.legend(loc="lower right")

	f.savefig(name4_p9)

	f
end

# ╔═╡ fad7f786-7f00-4c99-a208-dcafcb50ef56
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p9)

	f
end

# ╔═╡ 21304a7d-22a0-4ef4-aa59-7d66930b850f


# ╔═╡ d5f1af06-31b6-4ec9-a788-df58f2ee892f
spectra_p9B = [
	PlotDesign("t52.50g30.00m-4.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t52.50g30.00m-4.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ 60772fa6-fdc9-4a02-bf58-1eb1fbf5d156
differences_p9B = [
	PlotDesign("t52.50g30.00m-4.000 CEMP - C2.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t52.50g30.00m-4.000 no-CEMP - C2.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t52.50g30.00m-4.000 no-CEMP - C2.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ bc997051-e470-4ce0-8ff9-dd09ca5d229d
name1_p9B = "contribution_tau_t52.50g30.00m-4.000.pdf"

# ╔═╡ e339d118-74bb-4681-bb12-b94dc6a502a9
name2_p9B = "contribution_surface_t52.50g30.00m-4.000.pdf"

# ╔═╡ bf0a5a0c-32e1-47e6-bc00-66cdafd2dfde
name3_p9B = "contribution_height_t52.50g30.00m-4.000.pdf"

# ╔═╡ b3e9b8d9-91b1-47bf-9fde-b480561948e4
name4_p9B = "profile_tau500_T_t52.50g30.00m-4.000.pdf"

# ╔═╡ fd3849bd-80b1-43ea-9319-e36341f2f8c7
name5_p9B = "profile_tau500_rho_t52.50g30.00m-4.000.pdf"

# ╔═╡ 13ef8a7e-7749-4102-8f0b-2ce0265d0906
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p9B)

	f
end

# ╔═╡ 3fcc4340-02d5-4988-9d5c-1e1e232b5dfa
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p9B)

	f
end

# ╔═╡ ecd51ba0-5871-46a3-9b24-2c608b27af54
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -4.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p9B)

	f
end

# ╔═╡ d29863ce-5612-424c-a563-b8dcc94cf303
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p9B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

	k = differences_p9B[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p9B)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend()

	ax.axvline(-3.0, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(2700, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm [Fe/H] = -4,\ [C/Fe] = 2")
	
	

	f.savefig(name4_p9B)

	f
end

# ╔═╡ ce5b1ffa-b789-445f-8792-ae3c8d1aa60d
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p9B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p9B)

	f
end

# ╔═╡ efedcfaa-d90a-49f1-ae13-cb6c93bfcd7b


# ╔═╡ 8b3b6194-3448-4605-a366-fc51aa2dc453
md"## [Fe/H]=-3, [C/Fe]=1.0"

# ╔═╡ 7e467805-7bcf-4221-be05-27a1e4d89249
spectra_p92 = [
	PlotDesign("t57.50g45.00m-3.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t57.50g45.00m-3.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ ad048639-aebb-4d1c-a09f-3e9a3c9decb5
differences_p92 = [
	PlotDesign("t57.50g45.00m-3.000 CEMP - C1.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t57.50g45.00m-3.000 no-CEMP - C1.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t57.50g45.00m-3.000 no-CEMP - C1.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ eb2f9c36-0920-40ae-bc92-45333ab5385e
name1_p92 = "contribution_tau_t57.50g45.00m-3.000.pdf"

# ╔═╡ 4d3eeb0f-b981-4979-a220-620cd5d7bef1
name2_p92 = "contribution_surface_t57.50g45.00m-3.000.pdf"

# ╔═╡ df6f5bdd-b442-4501-8ab7-05a520dd016f
name3_p92 = "contribution_height_t57.50g45.00m-3.000.pdf"

# ╔═╡ 4c3fdce3-be2b-4196-907d-e7600c74303e
name4_p92 = "profile_tau500_T_t57.50g45.00m-3.000.pdf"

# ╔═╡ b24ea01e-77c7-4d14-a712-ec5e3f37a0a4
name5_p92 = "profile_tau500_rho_t57.50g45.00m-3.000.pdf"

# ╔═╡ 3ec7ed5c-4369-44d4-9c0c-7527415b2de1
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p92)

	f
end

# ╔═╡ 07562017-cbdd-40b2-bb1e-592e8abb8feb
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p92)

	f
end

# ╔═╡ 16d4de6b-4dc9-44c4-9a6b-bfd1252df846
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -5.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p92)

	f
end

# ╔═╡ 6af041e4-7476-44a6-af15-87bf4a1decbc
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p92)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

	k = differences_p92[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p92)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend()

	ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(2700, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm [Fe/H] = -3,\ [C/Fe] = 1")
	
	

	f.savefig(name4_p92)

	f
end

# ╔═╡ b94c008d-4a41-468d-8e65-a42e586f83ce
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p92)

	f
end

# ╔═╡ 9195232e-556c-47d5-820b-c492916eed54
spectra_p92B = [
	PlotDesign("t52.50g30.00m-3.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t52.50g30.00m-3.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ 77514283-cf30-4362-be95-d3f8e8018833
differences_p92B = [
	PlotDesign("t52.50g30.00m-3.000 CEMP - C1.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t52.50g30.00m-3.000 no-CEMP - C1.0", L"\rm CEMP - non\ CEMP", "-", 2., "k"),
	PlotDesign("MARCS t52.50g30.00m-3.000 no-CEMP - C1.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

# ╔═╡ 9d771e98-4e0b-418d-acc5-ca9e0dafa06a
name1_p92B = "contribution_tau_t52.50g30.00m-3.000.pdf"

# ╔═╡ 2cd63880-9111-4ce2-97cf-4e615990b618
name2_p92B = "contribution_surface_t52.50g30.00m-3.000.pdf"

# ╔═╡ 5c376243-3c38-449d-8ae1-d9ecdcc44aaf
name3_p92B = "contribution_height_t52.50g30.00m-3.000.pdf"

# ╔═╡ 5e932f76-ca1a-4305-9a1c-a8538fdef710
name4_p92B = "profile_tau500_T_t52.50g30.00m-3.000.pdf"

# ╔═╡ c60880f2-3073-4a7e-8c2c-0979a3e254d3
name5_p92B = "profile_tau500_rho_t52.50g30.00m-3.000.pdf"

# ╔═╡ b8a2c430-0c1a-42c0-9a38-5b1af3f12ba6
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_contr_av(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			plot_contr_1D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		end
		
	end

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{500}"*"]")
	ax.set_ylabel("contribution function [normalized]")

	ax.legend(loc="upper right")

	f.savefig(name1_p92B)

	f
end

# ╔═╡ 513b431d-5246-4c48-adce-afcd34b957a4
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_3D(s, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation surface", fontsize="x-large")

	#ax.set_xlim(-5, 4)
	#ax.set_ylim(0.1, 7)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel("Y [Mm]")

	f.savefig(name2_p92B)

	f
end

# ╔═╡ a2f277f0-a017-45bd-b573-88cb06bd43a8
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_contr_yslice(s, 10, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			cbar = f.colorbar(im, ax=ax)
			cbar.set_label("contribution function [normalized]")
		else
			s = spectra1D[k]
			
		end
		
	end

	ax.set_title("CH formation Height", fontsize="x-large")

	ax.set_ylim(4, -4.5)
	ax.set_xlabel("X [Mm]")
	ax.set_ylabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	f.savefig(name3_p92B)

	f
end

# ╔═╡ 2dc6c4db-0cb4-4f0c-b700-b990d48704b1
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()

	for (i, p) in enumerate(spectra_p92B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			x, y = plot_profile(s, :T, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)

			ax.plot(x, y; ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = s.run.ltau
			y = s.run.temp
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
	end

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

	k = differences_p92B[1].key
	xref, yref = if haskey(spectra3D, k)
		s = spectra3D[k]
		plot_profile(s, :T)
	else
		s = spectra1D[k]
		s.run.ltau, s.run.temp
	end
	
	for (i, p) in enumerate(differences_p92B)
		(i == 1) && continue 
		k = p.key
		x, y = if haskey(spectra3D, k)
			s = spectra3D[k]
			plot_profile(s, :T)
		else
			s = spectra1D[k]
			MUST.pyconvert(Array, s.run.ltau), MUST.pyconvert(Array, s.run.temp)
		end
		mask = sortperm(x)
		fref = MUST.linear_interpolation(x[mask], y[mask], extrapolation_bc=MUST.Line())

		ax2.plot(
			xref, yref .- fref.(xref), 
			color=p.color, lw=p.lw, label=p.label, ls=p.ls
		)
	end

	ax2.set_title(L"\rm \Delta T\ [K]")
	ax2.legend()

	ax.axvline(-2.7, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.3, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(2700, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

	ax.set_title(L"\rm [Fe/H] = -3,\ [C/Fe] = 1")
	
	

	f.savefig(name4_p92B)

	f
end

# ╔═╡ 3e437ffb-0c97-473a-bff6-929f8c720abd
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))
	plt.close()
	
	for (i, p) in enumerate(spectra_p92B)
		k = p.key
		if haskey(spectra3D, k)
			s = spectra3D[k]
			im = plot_profile(s, :log10d, ax=ax, ls=p.ls, lw=p.lw, label=p.label, color=p.color)
		else
			s = spectra1D[k]
			x = MUST.pyconvert(Array, s.run.ltau)
			y = log10.(MUST.pyconvert(Array, s.run.rho))
			ax.plot(x, y, color=p.color, lw=p.lw, label=p.label, ls=p.ls)
		end
		
	end

	ax.set_xlim(-4.5, 3)
	ax.set_ylim(-8.5, -5.7)
	ax.set_ylabel(L"\rm log\ density\ [g\ cm^{-3}]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="upper left")

	f.savefig(name5_p92B)

	f
end

# ╔═╡ 3fdeaaaf-2f00-4ef8-b6bd-cc1d29dbcb74


# ╔═╡ 04ea0005-24c6-4c39-816e-94a142fe65b8
md"# SAGA Data"

# ╔═╡ a215c938-11e3-495d-92cb-dc0559ed86f1
saga_paths = Dict(
	"CEMP" => "../cgisess_486fcb28413b27ba518e2836eea91aa5_CEMP_MS.dat",
	"CEMP-no" => "../cgisess_486fcb28413b27ba518e2836eea91aa5_CEMP-no_MS.dat",
	"CEMP-s" => "../cgisess_486fcb28413b27ba518e2836eea91aa5_CEMP-s_MS.dat",
	"EMP" => "../cgisess_486fcb28413b27ba518e2836eea91aa5_EMP_MS.dat",
	"MP" => "../cgisess_486fcb28413b27ba518e2836eea91aa5_MP.dat"
)

# ╔═╡ e45f8e09-b588-4e84-98c0-602f96f2c649
saga_all_path = "../cgisess_32aefcdc405d0a0c9fb7936f6f1324aa_data.tsv"

# ╔═╡ 284b014d-1cfe-45c6-9c58-2a7144eee417


# ╔═╡ b669a28e-642f-44ad-94bf-8e61af14d2be


# ╔═╡ fa1bdf35-c482-4869-b195-10b6e4861a2f
md"We can also download the un-split dataset"

# ╔═╡ fd0357ee-086b-4394-add3-fdff0de4cd39
saga_data = Dict(k=>MUST.readdlm(p, skipstart=1) for (k, p) in saga_paths)

# ╔═╡ 2c6fbe92-5a73-4360-8ed2-f51f59dc5e3f
begin
	saga_all_data = MUST.readdlm(saga_all_path, '\t', String, skipstart=1)

	teff_all, logg_all, z_all = saga_all_data[:, 5], saga_all_data[:, 6], saga_all_data[:, 7]

	cfe_all, bafe_all, eufe_all, lafe_all, feh_all = saga_all_data[:, 9], saga_all_data[:, 10], saga_all_data[:, 11], saga_all_data[:, 12], saga_all_data[:, 13]

	binarity_all = saga_all_data[:, 12]

	parameters_all = Dict{Any, Any}(
		"teff" => teff_all, 
		"logg" => logg_all, 
		"z"    => z_all, 
		"cfe"  => cfe_all, 
		"eufe" => eufe_all, 
		"bafe" => bafe_all, 
		"feh"  => feh_all,
		"lafe" => lafe_all,
		"b"    => binarity_all
	)

	for (i, p) in parameters_all
		d = zeros(length(p))

		for (j, pj) in enumerate(p)
			d[j] = if pj == ""
				NaN
			elseif (pj[1] == '<') |  (pj[1] == '>')
				parse(Float64, pj[2:end])
			else
				try
					parse(Float64, pj)
				catch
					NaN
				end
			end
		end
		
		parameters_all[i] = d
	end

	mask = (parameters_all["feh"].>0.0) 
	parameters_all["feh"][mask] .= NaN

	parameters_all
end

# ╔═╡ c1d54755-4b3f-4ad8-aac5-af7b365500c5


# ╔═╡ bc8ada92-8975-403e-a94e-ddc2d9ead536
md"To isolate the giant stars, one has to define a limiting log(g)."

# ╔═╡ 9e88185c-9c1f-4485-b465-30b3140b0e58
logg_limit = 2.7

# ╔═╡ e63719e6-d832-490f-ad1f-e903ec979c18


# ╔═╡ 49da5db2-1ff5-41c2-b0b8-15b2444e0522
md"We can furthermore define multiple sub-limits for different CEMP groups"

# ╔═╡ 7fa5a22a-49be-455d-b37d-f2683b29d5ed
# For CEMP in general
cfe_limit = 0.7

# ╔═╡ 9453fe67-731e-4832-918a-34a638d24042
# For CEMP-s and CEMP-r/s
ba_limit = 1.0

# ╔═╡ 1404ac62-c99c-4ee4-99b8-2e9b1e64608e
# For CEMP-r
eu_limit = 1.0

# ╔═╡ cc251a0f-c9bb-46f8-9640-d8670bff4526


# ╔═╡ 0575e66a-2008-42a4-9e05-3d8e3b2d6561
md"Apply r/s cut from Yoon et al. (2016) at A(C)=7.1"

# ╔═╡ bd5a7a2c-75cd-4a67-88a2-e09cf152f6a9
c_limit = 7.1

# ╔═╡ 25822bc1-6e1b-41e8-928f-4d088f0eede5
parameters_all["r/s"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit

# ╔═╡ 0aaeeb0d-2740-4817-83a4-fd56d9f2c9ce
parameters_all["no"] = (parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit

# ╔═╡ 92293f19-535a-4b35-8d58-bf99cddabcd3
selection_mask = parameters_all["logg"] .>= logg_limit;

# ╔═╡ 60b1f1ac-4fa5-48f6-b1fa-1c50a00447b3
let
	f, ax = plt.subplots(figsize=(6, 5))

	x = parameters_all["teff"]
	y = parameters_all["logg"]

	nanmask = .!isnan.(x) .& .!isnan.(y) 
	
	ax.scatter(
		parameters_all["teff"][.!selection_mask], 
		parameters_all["logg"][.!selection_mask], 
		color="0.5", alpha=0.1,
		s=10, rasterized=true
	)
	
	im = ax.scatter(
		parameters_all["teff"][selection_mask], 
		parameters_all["logg"][selection_mask], 
		c=parameters_all["feh"][selection_mask], 
		s=5, cmap="rainbow",
		rasterized=true
	)
	cbar = f.colorbar(im, ax=ax)

	xlim = MUST.pyconvert(Array, ax.get_xlim())
	ylim = MUST.pyconvert(Array, ax.get_ylim())
	ax.set_xlim(xlim[end:-1:1])
	ax.set_ylim(ylim[end:-1:1])

	ax.axhline(logg_limit, color="0.5", ls="--", alpha=0.5, lw=2)

	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	ax.set_ylabel(L"\rm log(g)\ [dex]")
	cbar.set_label(L"\rm [Fe/H]")
	
	f
end

# ╔═╡ 65d7b527-baa6-4af2-972f-1960c0072d10
begin
	@info "Number of stars in the Main Sequence: $(count(selection_mask))"
	@info "Number of stars outside the Main Sequence: $(count(.!selection_mask))"
end

# ╔═╡ 29d58861-8237-4cc9-a05a-c37b8228efaf


# ╔═╡ 6f3f7e16-5611-43a6-a5ca-5d25bd3f29f5
md"## Correction"

# ╔═╡ 6b2af110-702a-4572-9ead-24d2b90db0f6
begin
	p1_saga = (-4.0, 0.75-4+8.560, -0.55)
	p2_saga = (-5.0, 3.0-5+8.560, -0.79+0.15)
	
	scatter_int_saga(x, y) = MUST.pyconvert(typeof(x),
		MUST.scipy_interpolate.griddata(
		([p1_saga[1], p2_saga[1]], 
		[p1_saga[2], p2_saga[2]]), 
		[p1_saga[3], p2_saga[3]], 
		(x, y), method="nearest")
	)
	corrections_saga = Dict(
		k=>scatter_int_saga(saga_data[k][:, 1], saga_data[k][:, 2])
		for k in keys(saga_data)
	)
end

# ╔═╡ b120b287-060f-417f-902b-3ef4bc26148b


# ╔═╡ 80b491c3-7c54-4a8c-8460-c76a82ee7b86
design_p10 = [
	#PlotDesign("MP", "MP", "s", 20, "nothing"),
	#PlotDesign("EMP", "EMP", "s", 20, "k"),
	PlotDesign("CEMP", "CEMP", "s", 20, "tomato"),
	PlotDesign("CEMP-no", "CEMP-no", "s", 20, "k"),
	PlotDesign("CEMP-s", "CEMP-s", "s", 20, "lime"),
]

# ╔═╡ bd6060f6-0e09-4ccd-a9ed-12d7b63b54f2
name1_p10 = "saga_corrected.pdf"

# ╔═╡ dfe42b5b-9c74-426e-bfdc-b21748179f59
let
	f, ax = plt.subplots(1, 2, figsize=(8, 5), sharey=true, sharex=true)
	plt.subplots_adjust(wspace=0)

	for p in design_p10
		d = saga_data[p.key]
		x = d[:, 1]
		y = d[:, 2] 
		ax[0].scatter(x, y, label=p.label, s=p.lw, marker=p.ls, color=p.color)
		y = d[:, 2] .+ corrections_saga[p.key]
		ax[1].scatter(x, y, label=p.label, s=p.lw, marker=p.ls, color=p.color)
	end

	ax[0].axhline(0.7, ls=":", color="0.2", lw=2, alpha=0.7)
	#ax[0].text(-5.5, 0.76, "CEMP", color="0.2", alpha=0.7, ha="left", va="bottom")
	#ax[0].text(-5.5, 0.72, "not-CEMP", color="0.2", alpha=0.7, ha="left", va="top")

	ax[1].axhline(0.7, ls=":", color="0.2", lw=2, alpha=0.7)
	ax[1].text(-5.7, 0.71, "CEMP", color="0.2", alpha=0.7, ha="left", va="bottom")
	ax[1].text(-5.7, 0.67, "not-CEMP", color="0.2", alpha=0.7, ha="left", va="top")

	ax[0].set_ylabel("[C/Fe]")
	ax[0].set_xlabel("[Fe/H]")
	ax[1].set_xlabel("[Fe/H]")
	
	ax[0].legend(loc="upper right")

	f.savefig(name1_p10)

	f
end

# ╔═╡ 8b807632-3dd5-4b03-a539-1f348af4852c


# ╔═╡ 740ebde2-5ae2-476d-a7c8-b6541be938f5
md"Correction of the big dataset"

# ╔═╡ 99598049-9494-4a72-afd4-13dd12a09ff4
begin
	corrections_saga_all = fill!(similar(parameters_all["cfe"]), 0.0)

	interpolators_logg = [3.0, 4.5]
	interpolators_f = [ip_t52g30, ip_t57g45]
	
	for i in eachindex(parameters_all["cfe"])
		logg = parameters_all["logg"][i]
		teff = parameters_all["teff"][i]
		feh = parameters_all["feh"][i]
		#iint = argmin(abs.(logg .- interpolators_logg))
		#iint = logg > 3.4 ? 2 : 1
		#corrections_saga_all[i] = interpolators_f[iint](parameters_all["feh"][i])
		corrections_saga_all[i] = interpolate_corrections(
			teff, logg, feh, feh_interpolators=interpolators
		)
	end
end

# ╔═╡ d9b1abef-4cae-4171-904f-6e4a046a9bf0


# ╔═╡ 210a4e3c-6998-400f-9e0a-3da3e0efe557
md"## Galactic CEMP Distribution"

# ╔═╡ 0e5c7037-aa27-44f2-b563-94aca1da290b
#metallicity_bin_edges = [-6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]
metallicity_bin_edges = [-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.0]

# ╔═╡ f3442668-a393-4136-8d43-d597af23279b
metallicity_bin_centers = (metallicity_bin_edges[2:end] .+ metallicity_bin_edges[1:end-1]) ./ 2

# ╔═╡ 0f92fa2a-7f2f-4a04-bcec-5eaf307ed3f3
function bin_parameter(bins, general_mask=trues(size(saga_all_data, 1)); name="feh")
	bin_centers = (bins[2:end] .+ bins[1:end-1]) ./ 2
	bin_masks = falses(size(saga_all_data, 1), length(bin_centers))
	for i in eachindex(bin_centers)
		bin_masks[:, i] .= (bins[i] .<= parameters_all[name]) .&
		(parameters_all[name] .< bins[i+1]) .& (.!isnan.(parameters_all[name]))

		bin_masks[:, i] .= bin_masks[:, i] .& general_mask .& selection_mask
	end

	bin_masks
end

# ╔═╡ ab9353b8-4b9d-4e93-b169-aafb0a44247b
(253-189)/253

# ╔═╡ 0b489280-582d-47af-aed0-ac54a79a060c
begin
	bins_general = bin_parameter(
		metallicity_bin_edges, 
		.!isnan.(parameters_all["cfe"])
	)
	count_bins_general = reshape(count(bins_general, dims=1), :)

	# binarity
	bins_binary = bin_parameter(
		metallicity_bin_edges, 
		abs.(parameters_all["b"]) .>= 0.5
	)
	count_bins_binary = reshape(count(bins_binary, dims=1), :)

	bins_nonbinary = bin_parameter(
		metallicity_bin_edges, 
		abs.(parameters_all["b"]) .< 0.5
	)
	count_bins_nonbinary = reshape(count(bins_nonbinary, dims=1), :)
	
	# CEMP
	bins_CEMP = bin_parameter(
		metallicity_bin_edges, 
		parameters_all["cfe"] .>= cfe_limit
	)
	count_bins_CEMP = reshape(count(bins_CEMP, dims=1), :)

	# CEMP (corrected)
	bins_CEMP_corr = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit
	)
	count_bins_CEMP_corr = reshape(count(bins_CEMP_corr, dims=1), :)
	

	# CEMP-r/s
	bins_CEMP_rs = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		#(parameters_all["bafe"] .>= ba_limit) .& 
		#(parameters_all["eufe"] .>= eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .> c_limit)
	)
	count_bins_CEMP_rs = reshape(count(bins_CEMP_rs, dims=1), :)

	# CEMP-r/s (corrected)
	bins_CEMP_rs_corr = bin_parameter(
		metallicity_bin_edges, 
		((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& 
		#(parameters_all["bafe"] .>= ba_limit) .& 
		#(parameters_all["eufe"] .>= eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .> c_limit)
	)
	count_bins_CEMP_rs_corr = reshape(count(bins_CEMP_rs_corr, dims=1), :)
	

	# CEMP-s
	bins_CEMP_s = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		(parameters_all["bafe"] .>= ba_limit) .& 
		(parameters_all["eufe"] .< eu_limit)
	)
	count_bins_CEMP_s = reshape(count(bins_CEMP_s, dims=1), :)
	

	# CEMP-no
	bins_CEMP_no = bin_parameter(
		metallicity_bin_edges, 
		(parameters_all["cfe"] .>= cfe_limit) .& 
		#(parameters_all["bafe"] .< ba_limit) .& 
		#(parameters_all["eufe"] .< eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560) .<= c_limit)
	)
	count_bins_CEMP_no = reshape(count(bins_CEMP_no, dims=1), :)

	# CEMP-no (corrected)
	bins_CEMP_no_corr = bin_parameter(
		metallicity_bin_edges, 
		((parameters_all["cfe"] .+ corrections_saga_all) .>= cfe_limit) .& 
		#(parameters_all["bafe"] .< ba_limit) .& 
		#(parameters_all["eufe"] .< eu_limit)
		((parameters_all["cfe"] .+ parameters_all["feh"] .+ 8.560 .+ corrections_saga_all) .<= c_limit)
	)
	count_bins_CEMP_no_corr = reshape(count(bins_CEMP_no_corr, dims=1), :)
	

	# CEMP-unclassified
	bins_CEMP_nothing = bins_CEMP .& (.!bins_CEMP_rs) .& (.!bins_CEMP_s) .& (.!bins_CEMP_no) 
	count_bins_CEMP_nothing = reshape(count(bins_CEMP_nothing, dims=1), :)
	

	@info "Bin counts (all): $(count_bins_general)"
	@info "Bin counts (CEMP): $(count_bins_CEMP)"
	@info "Bin counts (CEMP corrected): $(count_bins_CEMP_corr)"
	@info "Bin counts (CEMP-r/s): $(count_bins_CEMP_rs)"
	@info "Bin counts (CEMP-r/s corrected): $(count_bins_CEMP_rs_corr)"
	@info "Bin counts (CEMP-s): $(count_bins_CEMP_s)"
	@info "Bin counts (CEMP-no): $(count_bins_CEMP_no)"
	@info "Bin counts (CEMP-no corrected): $(count_bins_CEMP_no_corr)"
	@info "Bin counts (CEMP unclassified): $(count_bins_CEMP_nothing)"
	@info "Bin counts (binarity): $(count_bins_binary)"
	@info "Bin counts (non-binarity): $(count_bins_binary)"
end

# ╔═╡ 6f03a4bb-63e2-46fc-a39a-680c67dc82b1
begin
	@info "cumsum of those bins:"
	@info cumsum(count_bins_general)
	@info cumsum(count_bins_CEMP)
	@info cumsum(count_bins_CEMP_corr)
end

# ╔═╡ c58fef48-f544-4377-97f4-7daf9ba218ac
name1_p11 = "saga_binned_cumsum.pdf"

# ╔═╡ aa05228a-9c86-4276-bdbf-c47479ccc7dc
name1_2_p11 = "saga_binned_cumsum_rsno.pdf"

# ╔═╡ 921e3b89-4baf-4984-85b1-087af8e8fc2c
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ax.plot(metallicity_bin_centers, count_bins_CEMP ./ count_bins_general, zorder=10, color="k", lw=2.0, label=L"\rm SAGA\ (MS)", marker="s")
	ax.plot(metallicity_bin_centers, count_bins_CEMP_corr ./ count_bins_general, zorder=10, color="tomato", lw=2.0, label=L"\rm SAGA\ (MS)\ -\ corrected", marker="s")

	ax2 = ax.twinx()
	ax.set_zorder(10)
	ax.patch.set_visible(false)
	
	ax2.scatter(parameters_all["feh"][selection_mask], parameters_all["cfe"][selection_mask], s=10, rasterized=true, alpha=0.2, marker="o", color="0.5", zorder=0)
	
	ax2.axhline(0.7, ls="--", color="k", alpha=0.3, lw=1.5)
	ax2.text(-0.2, 0.71, "CEMP", color="k", alpha=1, ha="right", va="bottom", fontsize=11)
	ax2.text(-0.2, 0.67, "not-CEMP", color="k", alpha=1, ha="right", va="top", fontsize=11)

	ax.set_xlabel(L"\rm [Fe/H]")
	ax2.set_ylabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm N_{\Delta[Fe/H], CEMP}\ /\ N_{\Delta[Fe/H]}")
	ax.legend(loc="lower left")
	#ax.set_yscale("log")

	ax.set_xlim(-6.8, 0.1)
	ax2.set_ylim(-0.75, 5.)
	f
end

# ╔═╡ 646546e9-8e80-428d-af73-02e3b033e092
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=true, sharey=true)

	plt.subplots_adjust(wspace=0)

	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, marker="", alpha=0.3, ls=":")
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, marker="", alpha=0.3, ls=":")
	
	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm CEMP-no,\ 1D", marker="s", markersize=7)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_no_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm CEMP-no,\ 3D", marker="s", markersize=7)

	ax[0].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm CEMP-r/s,\ 1D", marker="o", ls="--", markeredgecolor="k", markerfacecolor="w", markersize=8)
	ax[1].plot(metallicity_bin_centers, cumsum(count_bins_CEMP_rs_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm CEMP-r/s,\ 3D", marker="o", ls="--", markeredgecolor="tomato", markerfacecolor="w", markersize=8)

	
	
	#ax2.scatter(parameters_all["feh"][selection_mask], parameters_all["cfe"][selection_mask], s=10, rasterized=true, alpha=0.2, marker="o", color="0.5", zorder=0)
	
	#ax2.axhline(0.7, ls="--", color="k", alpha=0.3, lw=1.5)
	#ax2.text(-0.2, 0.71, "CEMP", color="k", alpha=1, ha="right", va="bottom", fontsize=11)
	#ax2.text(-0.2, 0.67, "not-CEMP", color="k", alpha=1, ha="right", va="top", fontsize=11)

	ax[0].set_xlabel(L"\rm [Fe/H]")
	ax[1].set_xlabel(L"\rm [Fe/H]")
	ax[0].set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax[0].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.04))
	ax[1].legend(loc="upper right", ncol=1, bbox_to_anchor=(1.05, 1.04))
	#ax.set_yscale("log")

	ax[0].set_xlim(-6.2, -1.5)
	#ax.set_ylim(0.24*100, 1.04*100)
	#ax2.set_ylim(-0.75, 5.)

	f.savefig(name1_2_p11)
	
	f
end

# ╔═╡ 823716ae-f1fd-4805-be6e-a9bdb3827940
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP) ./ cumsum(count_bins_general)*100, zorder=10, color="k", lw=2.0, label=L"\rm 1D", marker="s")
	ax.plot(metallicity_bin_centers, cumsum(count_bins_CEMP_corr) ./ cumsum(count_bins_general)*100, zorder=10, color="tomato", lw=2.0, label=L"\rm 3D", marker="s")

	#ax2 = ax.twinx()
	ax.set_zorder(10)
	ax.patch.set_visible(false)
	
	#ax2.scatter(parameters_all["feh"][selection_mask], parameters_all["cfe"][selection_mask], s=10, rasterized=true, alpha=0.2, marker="o", color="0.5", zorder=0)
	
	#ax2.axhline(0.7, ls="--", color="k", alpha=0.3, lw=1.5)
	#ax2.text(-0.2, 0.71, "CEMP", color="k", alpha=1, ha="right", va="bottom", fontsize=11)
	#ax2.text(-0.2, 0.67, "not-CEMP", color="k", alpha=1, ha="right", va="top", fontsize=11)

	ax.set_xlabel(L"\rm [Fe/H]")
	#ax2.set_ylabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm N_{\leq [Fe/H], CEMP}\ /\ N_{\leq [Fe/H]}\ [\%]")
	ax.legend(loc="lower left")
	#ax.set_yscale("log")

	ax.set_xlim(-6.2, -1.8)
	ax.set_ylim(0.24*100, 1.04*100)
	#ax2.set_ylim(-0.75, 5.)

	f.savefig(name1_p11)
	
	f
end

# ╔═╡ e6f3d47b-fb81-4100-8343-a8c822c67137


# ╔═╡ 5cc34af8-b49e-401c-8320-3a18bd8bcbd3
feh_p11 = -1.5

# ╔═╡ d17f4ee1-b3ae-49e6-80ee-26a322b2d661
relative_p11 = false

# ╔═╡ c9e781bb-1f8d-4a67-9db1-e142aab69e26
name2_p11 = "saga_feh4_hist_$(feh_p11).pdf"

# ╔═╡ 67ed94c0-acf1-4dfa-b0ce-165d0b8e9a9d
let
	gaussian(μ, σ) = x -> exp(-(x-μ)^2/(2*σ^2))
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ibin = argmin(abs.(metallicity_bin_centers .- (feh_p11)))
	
	binmask_all = bins_general[:, ibin] .& (.!isnan.(parameters_all["cfe"]))
	binmask = bins_CEMP[:, ibin] .& (.!isnan.(parameters_all["cfe"]))
	x = parameters_all["cfe"][binmask_all]
	xm = mean(x)
	xs = MUST.std(x)
	xr = range(minimum(x)-2.5*xs, maximum(x)+2.5*xs, length=100) |> collect
	g = gaussian(xm, xs)
	k = kde(x)
	c, b = MUST.numpy.histogram(x, bins=10)
	norm = 1 / maximum(k.density) 
	#ax.plot(xr, g.(xr) .*norm , label="1D", color="k")
	ax.plot(k.x, k.density .*norm , label="1D", color="k")
	#ax.hist(x)
	
	binmask = bins_CEMP_corr[:, ibin] .& (.!isnan.(parameters_all["cfe"]))
	x = parameters_all["cfe"][binmask_all] .+ corrections_saga_all[binmask_all]
	xm = mean(x)
	xs = MUST.std(x)
	xr = range(minimum(x)-2.5*xs, maximum(x)+2.5*xs, length=100) |> collect
	g = gaussian(xm, xs)
	k = kde(x)
	c2, b2 = MUST.numpy.histogram(x, bins=10)
	norm = maximum(c2) / maximum(k.density) / maximum(c)
	#ax.plot(xr, g.(xr) .*norm , label="3D", color="tomato")
	ax.plot(k.x, k.density .*norm , label="3D", color="tomato")
	#ax.hist(x)

	ax.axvline(0.7, color="0.5", alpha=0.5, ls="--", lw=2)
	
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm N_{\Delta [Fe/H]}\ [normalized]")
	ax.set_title("$(metallicity_bin_edges[ibin])"*L"\rm \leq [Fe/H] \leq"*"$(metallicity_bin_edges[ibin+1])")
	ax.set_xlim(-0.2, 4)
	ax.legend()

	f.savefig(name2_p11)

	f
end

# ╔═╡ a32770f6-8fc2-482d-9326-8476465b7ce3


# ╔═╡ 8d599135-649c-4418-8937-511dc02bb97a
md"## SAGA in Yoon-Beers diagram"

# ╔═╡ 1f039bf1-85ae-489c-b614-a2b4e24a73a0
name1_p12 = "saga_yoon_beers.pdf" 

# ╔═╡ 67219339-4d5d-4bb7-a76e-0465e8aebe62
name2_p12 = "saga_yoon_beers_reclassify.pdf"

# ╔═╡ 37bd6ab5-45d9-464a-874d-2ff12d869285
name3_p12 = "saga_cfe_feh_binned.pdf"

# ╔═╡ ac58c4eb-b800-48cc-8a95-96f7b7917f85
name4_p12 = "saga_cfe_feh.pdf"

# ╔═╡ 2ec865f9-0cc1-4147-827a-f4e4c38c44e9
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.8), -6, 0.8, color="orange", alpha=0.1
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.7, 6.0), -4, 0.8, color="green", alpha=0.1, angle=45
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.6, 7.9), -3, 2, color="blue", alpha=0.1, angle=10
	)
	#ax.add_patch(ellipse1)
	#ax.add_patch(ellipse2)
	#ax.add_patch(ellipse3)

	CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe .+ feh .+ 8.560
	ax.scatter(feh, c, s=20, marker="s", c="k", alpha=0.05)
	ax.scatter(feh, c+corrections_saga_all[CEMP_mask], s=20, marker="o", c="tomato", alpha=0.4)

	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="0.5", alpha=0.15, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="k", lw=3.5, label=L"\rm 1D")

	
	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="tomato", alpha=0.15, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="tomato", label=L"\rm 3D", lw=3.5)
	
	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=1.0, ls="--", lw=2.5)

	ax.set_xlim(-6.1, -1.1)
	ax.set_ylim(4, 9.2)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(loc="upper left")
	
	f.savefig(name1_p12)

	f
end

# ╔═╡ 254cb2e0-4580-4858-a1fa-098038930b37
let
	plt.close()
	f = plt.figure(figsize=(8, 8))

	gs = matplotlib.gridspec.GridSpec(
		3, 3, 
		width_ratios=[3, 1, 0.2], 
		height_ratios=[1, 3, 0.2],
        wspace=0.00, hspace=0.00
	)
	
	# Create main scatter plot
	ax = plt.subplot(gs[1, 0])
	
	# Create top histogram
	ax_xhist = plt.subplot(gs[0, 0], sharex=ax)
	
	# Create right histogram
	ax_yhist = plt.subplot(gs[1, 1], sharey=ax)


	
	ellipse3 = matplotlib.patches.Ellipse(
		(-5.5, 6.8), -4.5, 0.8, color="orange", alpha=0.25, ec="0.5"
	)
	ellipse2 = matplotlib.patches.Ellipse(
		(-3.7, 5.9), -3, 0.8, color="green", alpha=0.25, angle=45, ec="0.5"
	)
	ellipse1 = matplotlib.patches.Ellipse(
		(-2.6, 7.9), -2, 3, color="blue", alpha=0.25, angle=-30, ec="0.5"
	)
	ax.add_patch(ellipse1)
	ax.add_patch(ellipse2)
	ax.add_patch(ellipse3)

	ax.text(-6.7, 7.4, L"\rm \mathbf{Group\ III}", color="0.3", fontsize="small")
	ax.text(-3.2, 5.3, L"\rm \mathbf{Group\ II}", color="0.3", fontsize="small")
	ax.text(-1.4, 9.4, L"\rm \mathbf{Group\ I}", color="0.3", fontsize="small", ha="right")
	
	ax.text(-5.5, 3.8, L"carbon"*"\n"*L"enhanced", color="0.5", fontsize="x-small", ha="right", va="bottom")
	ax.text(-4.6, 3.8, L"carbon"*"\n"*L"normal", color="0.5", fontsize="x-small", ha="left", va="bottom")

	text_group3 = """
	$(L"\rm \mathbf{CEMP-no}")
	Faint SNe
	Pop III
	"""
	ax.text(-6.7, 7.5, text_group3, color="0.5", fontsize="x-small", ha="left", va="bottom")

	
	text_group2 = """
	$(L"\rm \mathbf{CEMP-no}")
	Standard CCSNe
	Pop II
	"""
	ax.text(-3.2, 5.15, text_group2, color="0.5", fontsize="x-small", ha="left", va="top")


	text_group1 = """
	$(L"\rm \mathbf{CEMP-s\ +\ CEMP-r/s}")
	MRSNe, NS-NS or NS-BH merger
	AGB contribution
	"""
	ax.text(-1.4, 9.5, text_group1, color="0.5", fontsize="x-small", ha="right", va="bottom")
	

	
	CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe .+ feh .+ 8.560
	#ax.scatter(feh, c, s=20, marker="s", c="None", alpha=0.2, linewidth=1, edgecolor="k")
	ax.scatter(feh, c+corrections_saga_all[CEMP_mask], s=15, marker="o", c="None", alpha=0.6, edgecolor="tomato", rasterized=true)

	
	# add histograms
	ax_xhist.hist(feh, bins=30, color="0.5")
	ax_yhist.hist(c, orientation="horizontal", bins=30, label="1D", color="0.5")
	ax_yhist.hist(c+corrections_saga_all[CEMP_mask], orientation="horizontal", bins=30, label="3D", histtype="step", lw=3, color="tomato")
	ax_yhist.legend()
	

	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP[:, i]] .+ 
		parameters_all["feh"][bins_CEMP[:, i]] .+ 
		8.560
	)  for i in eachindex(metallicity_bin_centers)]
	#ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="0.5", alpha=0.3, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="k", lw=3, label=L"\rm 1D")

	
	c_mean =  [mean(
		parameters_all["cfe"][bins_CEMP_corr[:, i]] .+ 
		parameters_all["feh"][bins_CEMP_corr[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP_corr[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][bins_CEMP_corr[:, i]] .+ 
		parameters_all["feh"][bins_CEMP_corr[:, i]] .+ 
		8.560 .+corrections_saga_all[bins_CEMP_corr[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	#ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="cyan", alpha=0.3, ls="")
	ax.plot(metallicity_bin_centers, c_mean, color="tomato", label=L"\rm 3D", lw=3)
	
	
	
	cemp_line(x) = 0.7 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=0.9, ls="--", lw=3)

	

	ax.set_xlim(-6.9, -1.1)
	ax.set_ylim(3.6, 10.9)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(loc="upper left")

	# Hide labels on shared axes
	plt.setp(ax_xhist.get_xticklabels(), visible=false)
	plt.setp(ax_yhist.get_yticklabels(), visible=false)
	
	# Adjust spines
	#=for ax in [ax_xhist, ax_yhist]
	    ax.spines["top"].set_visible(false)
	    ax.spines["right"].set_visible(false)
	end=#

	f.savefig(name2_p12)

	f
end

# ╔═╡ 4ee57a5a-be38-4a56-a686-95f9582c0910


# ╔═╡ 35e906c9-d6a2-4424-baed-186dd4f85f8d
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))


	#CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	CEMP_mask = selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe 

	group = bins_CEMP

	c_mean =  [mean(
		parameters_all["cfe"][group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="0.5", alpha=0.1, ls="")
	ax.errorbar(
		metallicity_bin_centers, c_mean, yerr=c_sigma,
		color="k", lw=2.2, label=L"\rm 1D", marker="s", capsize=5,
		markersize=8, ls="--"
	)
	@show collect(zip(metallicity_bin_centers, c_sigma))

	
	c_mean =  [mean(
		parameters_all["cfe"][group[:, i]]
		.+ corrections_saga_all[group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	c_sigma =  [MUST.std(
		parameters_all["cfe"][group[:, i]]
		.+corrections_saga_all[group[:, i]]
	)  for i in eachindex(metallicity_bin_centers)]
	ax.fill_between(metallicity_bin_centers, c_mean .- c_sigma, c_mean .+ c_sigma, color="tomato", alpha=0.1, ls="")
	ax.errorbar(
		metallicity_bin_centers, c_mean, yerr=c_sigma,
		color="tomato", label=L"\rm 3D", lw=2.2, marker="s", capsize=5,
		markersize=8
	)
	@show collect(zip(metallicity_bin_centers, c_sigma))
	
	
	cemp_line(x) = 0.7 
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="0.5", alpha=1.0, ls=":", lw=2.5)

	ax.set_xlim(-6.1, -1.5)
	#ax.set_ylim(4, 9.2)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("[C/Fe]")
	ax.legend(loc="upper right")
	
	f.savefig(name3_p12)

	f
end

# ╔═╡ ca1c8ba9-ef43-4ca1-8ea4-f2e08415ed75
let
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=true, sharex=true)
	plt.subplots_adjust(wspace=0)


	CEMP_mask = (parameters_all["cfe"] .>= cfe_limit) .& selection_mask
	#CEMP_mask = selection_mask
	feh = parameters_all["feh"][CEMP_mask]
	cfe = parameters_all["cfe"][CEMP_mask]
	c = cfe 


	ax[0].scatter(feh, cfe, label="1D", color="k", s=15, rasterized=true)
	ax[1].scatter(feh, cfe.+ corrections_saga_all[CEMP_mask], label=L"\rm 3D", color="tomato", s=15, rasterized=true)
	
	
	cemp_line(x) = 0.7 
	x = range(-9, 0, length=100)|>collect
	ax[0].plot(x, cemp_line.(x), color="0.5", alpha=1.0, ls="--", lw=2.5)
	ax[1].plot(x, cemp_line.(x), color="0.5", alpha=1.0, ls="--", lw=2.5)

	ax[0].set_xlim(-6.1, -1.5)
	#ax.set_ylim(4, 9.2)
	ax[0].set_xlabel("[Fe/H]")
	ax[1].set_xlabel("[Fe/H]")
	ax[0].set_ylabel("[C/Fe]")
	ax[0].legend(loc="upper right")
	ax[1].legend(loc="upper right")
	
	f.savefig(name4_p12)

	f
end

# ╔═╡ c615896c-cdb9-4cfa-af81-6ee3ad6b3e02


# ╔═╡ Cell order:
# ╠═32b14b00-8ae8-11ef-2aa2-937a4b8d3a0c
# ╠═bbb5fb60-6669-440c-84d6-5e0edbd90e3d
# ╠═024bd32f-a8c3-43eb-bba3-05d590f9632f
# ╠═9967494d-13ad-4997-95c9-341ca3402ded
# ╟─bd2ecf16-79e4-48e9-947a-0a323a1d3dd7
# ╟─8b84192a-765a-43a2-a7af-006e74be01fb
# ╟─1ae7a204-5a76-4705-94f2-ab248c3ecc20
# ╟─af7646aa-a4f8-4eb2-8121-1e99559a16f0
# ╠═b5b1d8c0-b67a-45d2-b21e-e0e04555a9fd
# ╠═c7a6739a-8ae7-45c0-ab4f-e800a9ab1d2e
# ╠═a343cf41-e425-4a59-81bd-87f8765251c0
# ╠═cd8b84b0-6891-4d3e-ba04-22fd4890e09d
# ╟─d2884a6e-b1e5-4234-ae78-05f682ac6908
# ╟─863499a2-1e33-4ba4-8dad-935a9a38a4ec
# ╠═9d8761fd-4880-4ff3-855f-f3713a55a6e5
# ╟─b26dbd2d-85ce-41c2-a32a-c9a488f9469f
# ╠═5ac9c28b-be99-4c27-9bdd-3aff553f9bc3
# ╠═e4eef27a-65c4-4ff5-a3e7-b0119a7e33c5
# ╠═67c8479d-7c41-4baa-b2ec-13da36a48900
# ╟─4bc7b154-2fdd-4057-bf40-e44c7851bc7d
# ╟─0a0779a4-7d87-40d2-9dfc-454b366c2dff
# ╟─5beedaf9-b12a-46db-a576-7913eb5f2caf
# ╟─a0cdfa21-bc3f-4462-aa4d-a84cc795bece
# ╠═6613157c-d199-4825-b779-a477ac4f3f30
# ╠═e67fed9e-3a82-45ed-9865-1c5862ebba2c
# ╟─f3493218-6013-4567-b1c4-442d0a48c9d0
# ╟─6ba576df-1a11-4be4-811b-212ed1f23347
# ╠═63ed5788-011d-4a31-ac82-efcd135f609c
# ╠═4bf4013c-ee0d-42fd-9dd2-4f26f090c15b
# ╟─4ca89537-5cc7-42e6-b131-362ad7ceb0b0
# ╟─0dee6625-70b1-489d-be3c-436062ca992b
# ╟─b231118c-bfed-4c9f-9522-cf7de097a12b
# ╠═5d5df14a-3b77-4781-9a53-eb6031b71c95
# ╟─4cf9ee24-a75b-45c2-81b6-7012879a2217
# ╟─c1c36104-a84d-4894-aa8c-33d3fdafc25b
# ╟─8d7b17da-ead8-445c-9fb0-4a6d98ca5dda
# ╠═4b636278-0e13-44d0-80c0-3c4fd9e1f0b6
# ╠═d786a454-b820-4c2f-b968-f7dec0604598
# ╠═52ffd5ea-6215-4644-8e27-5777e594567a
# ╟─cbb98c8a-8ec3-46fa-a38f-4afc3d292f41
# ╟─2aadb41f-a64e-4769-a202-aa9d2d700104
# ╠═56d11775-333b-44ac-ba5f-39e1b372564d
# ╠═58ff1560-c4ad-4ac5-bf29-33977b211b11
# ╟─c421ceca-ee21-455c-bfec-371c6c260ea3
# ╟─41671152-f0cb-4d1f-8aec-740a486c69b4
# ╟─c437b8d0-3da5-4d74-b796-0accda21a1f6
# ╠═13bdade4-3eb8-41cb-aa64-463124ccacf1
# ╠═f3398799-5ad0-4487-9ae5-cc9d7b679f3f
# ╠═01720034-cc2c-4bc5-9d99-f2073b5d0bcc
# ╟─455ca2ac-6fe2-457a-941d-4f2b8c2f58b3
# ╟─124e57c4-82e6-4abf-bb48-7f882da71c5b
# ╠═4d2a898f-8bc4-4c88-b35d-cb7542d55932
# ╠═ae9d6151-45cd-4301-bf3f-df319baa0870
# ╟─f8c93a98-ad4f-4d9c-b0cd-c66d432753ed
# ╟─0d0234bc-c6b5-4bcf-8409-01cf02c8c843
# ╟─773a7896-a492-46a3-9813-fa739647142d
# ╟─244f1b18-80fe-419a-8fd0-892aad803ca5
# ╠═accdaadc-8912-48d2-8a0c-3c19e5b4ba19
# ╠═ffb5a590-77e5-41c5-8964-eafe1153a5dc
# ╟─778a7f04-0841-4801-8aa8-b5686e1e498c
# ╟─983f1c2d-1532-491d-be04-e52e64c42775
# ╠═5a5704c2-2479-4bbd-87c5-cc015d546612
# ╠═785a473c-1a46-4ef3-8c2f-e7170452064b
# ╟─51f58ec5-f1bc-4b66-8ef1-c65ad08f15c9
# ╟─f61fb319-518f-4b56-9e19-a172e499307e
# ╠═0487f738-1a0e-49f9-96c3-357634226467
# ╠═cb4ca3c9-6c8a-4876-b927-a1419e403f43
# ╟─56e5777e-e80d-458d-b6a0-cde6a8f811ea
# ╟─89d5c1d8-733b-4f11-b306-a7410ae48e75
# ╠═85d716f2-542a-4793-8fa7-a53383fb37ec
# ╠═8992c977-8832-4102-89b7-430b4abd7e46
# ╟─43c0fc35-8336-43ba-a2c1-c4ae2ce33d43
# ╟─b0d6cf77-7197-4c1a-a59a-679dca8ddf7e
# ╟─d40dd522-128d-4407-a53e-62977f43b27d
# ╟─0b060c44-9da6-46d5-9c8d-4a6d0285c2a2
# ╠═07c86b08-0181-491e-952b-2e468ea8896c
# ╠═87581ad4-4fe8-4e46-9529-9255448f22fd
# ╟─1ef702b2-1a94-4794-acd7-4b5516595203
# ╠═26e9ea76-9c43-471e-a195-73150208ae19
# ╠═f8461a99-8d27-4293-bdc4-9afe19728a4e
# ╟─86e48b71-a6a9-4c3a-9227-7f013324c813
# ╠═228cda78-3c92-41fa-bd8d-74b8714ce094
# ╠═2158c40c-a464-45a3-a0c9-2fc69fe635da
# ╟─326ee706-2b03-4fd5-bc16-e396e2ec4a1d
# ╟─96e30c0f-9493-42b0-9c38-888b653723d1
# ╟─59566cb5-45c3-41a6-977c-b6af33e15c1f
# ╠═30fa8414-9468-49ac-9c9e-8e61823a0263
# ╠═d95e2320-31bb-4e50-8292-d7c713d58a87
# ╟─1ed777ae-6b23-4bdf-8a6b-5336e88fddde
# ╠═5c7630d5-e305-467c-87d2-1caebd62c1f2
# ╠═234cfaf6-099c-4fc1-a776-c44de65620c0
# ╟─bfdcd5c1-621f-49b5-b597-fef0d87d30d9
# ╠═c2c504e6-8f4d-4220-9a2a-099c17a06565
# ╠═0265aaca-4836-4082-9d21-c11a5f0fde86
# ╟─11a5eed7-b234-4d42-81b8-b415c8679ff5
# ╟─83729ef7-9c05-4c0b-8bcf-0ecfccc48d66
# ╟─0d4a0a2f-758f-44a0-9692-23c082befac8
# ╟─1f6dfdaf-b82e-46ee-989f-dbeb685d72ff
# ╟─391e98d1-344e-4ec3-837b-58327c02c74b
# ╟─10d7d97f-5c54-4066-9618-c78983c6c18c
# ╟─b97665fe-1536-4047-b3e9-d6176d47f3c1
# ╟─a04336bf-9cdf-401e-b237-a41f44987320
# ╟─1deeb995-4fda-4071-8d1f-fe9a9bc8a09c
# ╟─379e05ed-20db-499a-a0bd-b93f0e3e3642
# ╠═0b18ead3-5444-4737-99a8-cd86404343d3
# ╠═9e9ce7c5-2860-4cb0-9ddb-52323a4beb1d
# ╠═b28b1962-34c6-49f8-973a-3f63eed6c4d0
# ╟─c0855af1-d407-4ec2-83fc-4397330ac35c
# ╟─41367bd1-d61c-4b45-8781-6b5b71c7a511
# ╠═1943029c-2f0c-4590-9619-89870040932c
# ╠═5f7c6598-73d9-4e78-90ff-a53fa8ccf6c7
# ╠═ff48212b-40ae-4bcf-8576-034cd46f4738
# ╟─7df2c30a-445e-498e-9ebd-1a206ece821e
# ╟─2aa33bb2-3cd3-4104-a6a2-0f1b12fe2960
# ╟─c905fd9a-e2f6-485c-9096-c6a588b31b4d
# ╠═91699cc6-cb81-4382-bc5a-01931914343b
# ╠═aeccdb4b-63f9-47aa-830e-21aee1093dc4
# ╟─a55d5b04-a41b-47b6-81e8-39073bf538e2
# ╟─a056ba9a-f392-4945-a2c3-a64fc1849831
# ╠═77bc95ff-f456-46da-a81d-5412e6be032e
# ╠═baede960-e93b-490a-a027-1caf661428c5
# ╟─3537629b-0766-43e3-9506-30cfe0cc4a30
# ╟─136c2f2a-8902-4e05-920d-1d28c918a3f1
# ╠═dddd50fb-8044-4755-8180-647311268a91
# ╠═e551c818-ed62-4a2e-b2ac-74bfe6fb9b5f
# ╟─809c6511-ff14-47e4-80d7-4c984b3e8d3f
# ╟─f0ff834f-203c-460c-af2e-ee1fed8be671
# ╟─64073a6c-f9b1-45ad-8494-b948f6ec08d2
# ╠═c7631a99-5344-4a84-b1e1-5d9d96ddfeeb
# ╠═fe6d9b1c-a5c3-4b03-b66a-980191b21ba8
# ╟─00327ee0-ef59-43aa-bfe9-768f0f154f2c
# ╟─b0d32b4a-a107-47be-b12c-b4ff6ce46b53
# ╠═825031cd-5e40-49b6-936d-4a263a0c0378
# ╠═948ee564-7a5f-4581-9ca5-b0787ac41823
# ╟─18d9db1c-67b4-40f4-8548-057c7039020f
# ╟─f42278a0-ed19-47a8-bd2e-d78d1ef83d3a
# ╠═d12753f2-9c57-4af0-9623-194955c17cb7
# ╠═77ac5fd5-959a-452c-af88-73225291a890
# ╟─aa555374-910d-4202-9862-6742cb49f0af
# ╟─1397cb72-51db-45bf-8a67-cf1352f7f990
# ╟─bd3ee310-3f68-4e08-8ba8-d881aeecb76a
# ╟─f064d1a4-be5b-4cfc-9bf9-67002a333d2d
# ╠═9987b27b-2437-45a3-9e83-116a01f13d1d
# ╠═17ced5c3-080f-4d65-aa1c-b24a2950e783
# ╠═75600b27-d640-438b-bd00-d35e18ff46f6
# ╟─14ee9206-f322-4e06-a65e-271c19d4c792
# ╠═6c428722-39e9-4f1c-851a-ee90cb617377
# ╠═059c7c73-acaa-468b-b7a4-e63a4748d358
# ╟─df0a0a2d-52ea-4852-b450-6fc1baf5838d
# ╟─9669b2e7-396d-4e09-827b-3138c3d8f70a
# ╠═20e610c2-9054-4dc6-a551-6f63909715ee
# ╠═5d9fc015-4d8b-4811-8f48-c4197eb95709
# ╟─e09757e6-29c2-414b-9ac8-2a5fe237e0e6
# ╟─59ab2b58-2e38-4cce-9ba4-bc6899616a7d
# ╟─a61f7ec3-88c9-42f7-bf5c-acc792a05fb5
# ╠═2aa8d303-dae3-4a23-b253-c6091e2a638d
# ╠═cec8b45b-e4b6-4ad0-a83d-e9afe62b486a
# ╟─ff0787a1-b1cb-4eb0-85d7-b33bd39b2533
# ╟─f8bf1f6c-5522-463f-9ff0-cb59c1a950a2
# ╠═fda2a05e-48b0-471a-a5a5-4ada3d0d5d78
# ╠═9d7f44dc-a938-4f7e-8355-b8dc31e707a8
# ╟─334059d2-ea76-4c0d-84d9-7b1b9ef1ea87
# ╟─d3530325-728c-4392-b25d-18f4e5651719
# ╟─745b79ab-48f7-424b-ba20-9d7761c1cd98
# ╠═916bfaf4-10da-4206-ad4e-a77e8274d0c5
# ╠═f9a45d78-cba7-461b-9237-eeac330cf6d1
# ╟─5af5f222-d401-4dfb-aeea-0531d9a160ff
# ╟─fd2cb9f8-59da-467b-9ac0-81002549ce51
# ╠═f5e18ac4-8058-4b25-b1f6-86ad30f68dbe
# ╠═16a3524c-1367-4f62-a657-b0c6ed85159d
# ╟─681fc74f-cecd-4771-a967-f28864520323
# ╟─5924f436-2453-4b92-8b83-8e8033b86c7f
# ╠═a8ea82ac-7516-4223-9fa9-8ec93b63416d
# ╠═2a23cc14-6593-4f54-a42b-65ded01323af
# ╟─28e040c9-a4fb-4e58-9e00-2427c69ba268
# ╟─c09bfd6d-4801-4d8d-8755-147af007c93b
# ╟─f928a5ce-3022-4dc8-a70c-75971734398c
# ╠═0d68e17d-da25-4cb4-a89b-7465644e0f20
# ╠═721f3a04-aa21-4ff5-b566-4968d486438a
# ╟─335f1b86-a585-4ffa-a1d8-6790739f1bbb
# ╠═6ec5f942-94d4-4dfa-9eee-3baa7636476d
# ╠═20e07853-44ed-47d9-8937-b923bc570df9
# ╟─45ef7372-5e7c-4f1a-a837-584c5dc16b6e
# ╟─da1e169c-1c79-4d84-ac61-e1fab544f696
# ╟─76464024-a571-4c72-bb59-f881b595c124
# ╟─e0e4254a-e9f3-491c-9cac-6ed8a349e1ae
# ╟─36c4a557-54d8-4787-a0c9-4942b9dff4ef
# ╟─1ddb0459-d832-49a4-9bf2-af1627edd2e5
# ╟─6d6b3b86-0bec-427a-b795-3c31c96f9968
# ╟─ce9723af-dbc8-440c-aeaf-c88dd902d57b
# ╠═bf3b3496-3d18-40af-a89c-dc3dbc4d54c2
# ╠═b2f9818f-f434-4a2f-9db3-67ea2f229105
# ╟─bb3097f0-048c-477a-ad54-890c705437ff
# ╟─59457f84-e696-4cb1-914f-65305102a4fb
# ╟─804fbe02-c7d1-4bf9-91da-9cd3ef441e46
# ╟─452e4faa-57be-4d74-a7c5-36673c21aa9b
# ╠═fdcf50ee-af07-466f-9382-2a3af5de5a25
# ╠═3a5bad17-f447-43bb-9d9a-1679b4acf9e6
# ╟─764090a5-8132-4310-9932-fe80e25dbce5
# ╟─92442a47-82f9-465a-9dc0-89ca9f1dddea
# ╟─f7d263c5-ddd6-4d94-945d-303a8c81b8c6
# ╟─ec5b18a4-70f9-4078-8d64-7d3f52dbb7b8
# ╠═6fd82808-3881-4cdf-83d1-f2095cb553ac
# ╟─ec97726a-ad08-454a-874c-e5ddf54a3e85
# ╟─483897d7-fbf3-4a3d-a32b-2992ce5fccd2
# ╠═7021b190-4a8b-4ae4-bd71-06b3840f20c6
# ╠═b7a0779f-bd55-4ec3-afc8-ea8caf9db6fc
# ╠═d9818823-8afb-4ad2-876d-24ff473a93e7
# ╟─f3f50db8-66b7-4d48-b067-bea39c9c934c
# ╟─e4d71c15-0b34-4e99-a738-d41892f93e00
# ╠═6deaf7e4-7d5a-41d2-9f8d-e024d8150506
# ╠═82013231-d110-4791-acc4-46bf165bec09
# ╠═f2d61607-6728-4520-a3ef-e1b4b29b1fbb
# ╟─a8abb2a7-687c-45a4-9c47-2273a181b6fb
# ╟─fc260ef0-2b8b-415f-9e16-de551b281d43
# ╠═eb6f1896-6945-459f-859f-2594c76de6d8
# ╠═b742f64a-b6fe-4c6e-9368-9c4e60fc702c
# ╠═4ddd121e-8b2f-4da4-9474-5cc9894fc43c
# ╟─2dae804b-ef8e-4f39-b6f6-011bb92f0384
# ╟─38ce6bb7-f6f2-44d3-8685-8049c42ffeee
# ╟─2cd6317a-97b4-45bb-8e87-b09b40a470f7
# ╠═de419d2f-eb6d-47cd-a69c-99d5cb3161a3
# ╠═135f931b-387f-420d-8a0e-d7f391c1569f
# ╠═5f9b9dc4-f43b-4853-85d4-ae51247c13d0
# ╟─b0dfbe89-c719-4622-8f7b-55e5a09590df
# ╟─d165ca55-95ab-4702-b194-83d123f9d333
# ╠═b15fc564-a58d-49f6-9cf6-6587aeb61e1c
# ╠═99cc258f-19df-4ba8-8947-60a730db5f39
# ╠═c1eb200f-dd1b-40b2-8f1a-ea53908778db
# ╟─bf484955-9cdb-4e12-99a0-bb744e977e9b
# ╟─ee74de1f-4f7f-4c35-afb8-1186813b8a7b
# ╠═5ce8c3a7-d12e-4a3b-a3ed-81dbcd3a6118
# ╠═4df7f92c-db30-4261-b757-86f9af54c6b3
# ╠═ad4dcf0b-89a9-44f9-9de9-185fce79a510
# ╟─f7f9e7fd-dafc-40f3-8ce6-4af55d2a89ac
# ╟─e7f0e3ea-3204-4c2f-a464-e270baaf8db1
# ╟─1ead0a28-b327-476e-8be8-0e84ac546d2f
# ╠═44717756-9575-47b7-a3e5-59c271f79384
# ╠═bc158515-fd34-412f-b335-3934ebea6d02
# ╠═21f23414-0891-4aee-ab2b-c2f3c156dd2f
# ╟─2bdfa4d9-3c0e-40c2-950e-61b7860d529d
# ╠═060d0d6b-f52e-4a5a-acec-858f0ea1384d
# ╠═859af4d4-1e06-4567-ae2d-2d4083da2891
# ╠═de377555-e268-4185-82d8-588b0ee36b38
# ╟─0140b523-2d02-4522-aab4-201c0014ae05
# ╟─659cb51a-06ea-432a-91e2-9bd3f5006eb4
# ╠═c2e9d134-6bab-454c-b472-9b33965e031b
# ╠═b0ae8644-36f4-4f09-8602-9a36ab5ea08f
# ╠═d779ed13-aa80-4844-adcb-4662f59ad0c6
# ╟─298ab49e-f1a3-431f-a9c1-6c5bcce6771b
# ╟─2c44e5cc-43e9-4181-ab36-d629f4c61080
# ╟─098be297-ba34-4648-a08b-935db2bf5b7f
# ╠═63a51e50-4cb5-4df2-92ab-c24eeb23e482
# ╠═96980f1d-4e4c-48c9-8fe6-080cb9eb3ef8
# ╠═5c6e6453-5911-45ac-9fc0-1d69d92516fd
# ╟─693a142f-f26f-4543-8751-54dd8425a6b9
# ╠═ca64a419-dc88-45e0-904e-535cfcd94745
# ╠═b050a849-317f-4c2e-bc05-f79ff513dea4
# ╠═0d028a71-18f0-4c6a-bd1b-75c230a1cd08
# ╟─02f3d4b0-8435-4bab-9218-0ba803df0faf
# ╟─26108373-70d1-42ea-afea-dbd1586c0798
# ╟─27a3ef9a-d39a-4df5-8644-794a9196b1be
# ╟─f8ca6c4c-0de3-4ae8-b860-9ef17e04f2ad
# ╟─d8117e26-d1e4-48ed-a543-dbabc7b7cfc8
# ╟─579ec7d6-9535-4b71-9c9f-5c19d21d3da9
# ╟─517d0552-7f2a-4f70-8e1f-fae0cd8cd05e
# ╟─7714d6d3-fd29-47ef-a52f-621ec02800ff
# ╟─641f1109-17a5-4401-bc27-70fc4db7583e
# ╠═eb865dfb-b206-4daa-a1e7-d8c706ba2485
# ╠═a206fa91-09cd-4faa-b8ae-ccda70d36c16
# ╟─4b8de50b-b68a-458b-8408-1ded9b10976c
# ╟─2d3423ee-7a2f-4cc9-8800-267caf418034
# ╟─ef746dad-2b45-42c8-beb4-02ca0f4d983d
# ╠═ae657bc5-ad83-4675-bdda-f578d7054af7
# ╠═24c0eb4f-5dcf-48bf-b305-2e7737c1828b
# ╠═6a495bfa-fc58-4cbd-89a7-1a8e198a2e91
# ╟─b768e044-d7fc-49a9-bc5e-1859c359f0d1
# ╟─05c3da3f-6654-4ff4-953e-9f6d6c315f16
# ╟─349d1a5b-aed0-4315-b7b4-cc766e2d0bc3
# ╠═8cc66430-8635-4729-b595-a7df6a28cba5
# ╠═3bea6988-24ec-46bf-b20d-1f9817e7cc0b
# ╠═ac78dfbc-1423-4741-b0b4-dcd8c809787f
# ╟─90f2ca39-e67f-477b-ad1b-ff48c578174f
# ╟─009e81c3-7b8c-433a-a931-982f9a406d94
# ╠═b5896895-1bd7-40de-b86b-b326d502e587
# ╠═5c5bf31a-b1c4-478c-b2ae-6f42957d0c19
# ╟─103de55b-aea6-44eb-b3c1-e733a8c4c6e3
# ╟─d0c31252-9b07-47a8-a22c-25dcfca07d2a
# ╟─5f76e7c0-7366-4d1b-b88e-e9e3639d44de
# ╟─dc0509a9-a1bb-438f-8250-e52a56ffd1e4
# ╟─1773e7cb-bdff-4a97-be64-e0aa8284225a
# ╟─95ad457e-9d59-45a4-b983-5e2e2b0f6c76
# ╟─fe710591-6d00-4dad-9173-d8cf17ffefc1
# ╟─294fa362-da67-4a22-8a68-8619e185df6c
# ╟─6addc6c2-ad44-4310-a39a-292acf6b0dac
# ╠═9962676e-40cb-4fe8-ab50-8abb5a854612
# ╠═d5af99b0-8e95-4037-af01-1437f9aa0e04
# ╠═19c42123-ae32-468b-b6a1-4fd7457aa526
# ╠═43b4cf33-7efd-48b8-af90-66555c827b52
# ╠═71545110-8314-45c7-9e83-e15035c3f60d
# ╠═266c3319-78c1-4a2b-b06c-386f54907ecc
# ╠═d1c307dd-ab6f-44bc-a0d5-a879960231a5
# ╠═57c74ea4-a4ff-4018-9507-2753fc2a440c
# ╟─2765daec-fb70-4b55-889f-0ceb0808e84e
# ╟─515b8fb8-2568-4200-b10f-6ef8e2318615
# ╟─d08caaf5-e763-498d-a594-8915fc8eaaf9
# ╟─1dcf4d03-7a20-4ff7-b891-958003c23d61
# ╠═c8a6af9e-625a-4864-9ca2-1efad14fccd4
# ╟─a21d6f31-f766-432f-856c-7726e989d42f
# ╟─67aab606-2b0e-44f8-90ab-99d63830742d
# ╠═29ff7a48-9110-4da3-b024-2a9c7af59477
# ╠═a007184b-f10c-4592-ad0c-aca561a48a45
# ╠═932b8ec8-efcf-4851-8b3c-46fc18d0a674
# ╠═9ad8b431-b875-4ebd-844a-e70177aa6fde
# ╠═4b07347d-8e66-4087-84b7-5efb4076139f
# ╠═b6ba9f13-2f1d-407a-a496-661f89a86205
# ╠═4afb0820-64bb-4ef3-aa4f-820621da1452
# ╠═b14074bc-e4ef-4669-8c42-5f6194306029
# ╟─5b154765-c22f-45e0-9ebc-b7dfa384ebf7
# ╟─ec4b9098-48f3-4967-84fe-fd2d7cc4020e
# ╟─31c72e56-4f33-4328-be62-817a827d2161
# ╟─d8c47e25-67c5-414f-8416-8d26030c9cde
# ╟─bf028bf1-7480-4020-8760-6622d0057188
# ╟─896797bd-8fb0-4127-a2d6-57438ee03fae
# ╟─75480722-5bf5-453d-9f13-690cda6cdc3f
# ╠═009eabf7-0385-4ed7-b648-1e0b8c4cfc47
# ╠═aee45896-0f70-4bc8-8c8a-e685d92cb60a
# ╠═9f452aaa-e239-4201-ab9f-21136f8768be
# ╠═e1d754a5-cb5e-4fa1-9df2-fb6c806b17ea
# ╠═1882ff98-7c04-4909-8010-3ca0ed5f6e51
# ╠═8e757a6a-cfc2-458f-bc10-f31e6606e3ee
# ╠═b1e90b50-628c-49ef-85ea-5691cebd2243
# ╟─6d528aab-ecc9-4907-966b-e76250f08a2d
# ╟─de50c004-56b1-4719-9a93-6543d0f75f9b
# ╟─02e275ca-adf2-4497-a269-c5279bcbc5e9
# ╟─c90d990c-b17a-46d7-a145-00b5eedc0682
# ╟─31d31e74-c494-4ea6-a28c-d2bb7afd57a4
# ╟─b5cda9d6-480b-42a1-8091-7ceab841c3df
# ╟─6225b614-dc48-4158-8abe-55a1c07b651e
# ╠═b04303d1-bd59-4104-bc7c-8fdef9477fee
# ╠═901d430e-254b-4774-b7cc-4dcbd8ea1772
# ╠═479b4bca-166c-43c2-8199-04218a24fb12
# ╠═7a734c2e-175e-4300-8152-718669ce8137
# ╠═77bb0169-c777-425a-a8b5-041cb01279f1
# ╠═865f7e7a-91a9-4d16-b12d-43a84479e096
# ╠═1fd35559-c9d6-4be2-bee0-6d43417e4617
# ╟─5e7ef044-c782-4ee5-b408-c394118a8f72
# ╟─5e23e8f0-cc1a-4912-b841-c9d9f01e9681
# ╟─3d1174e7-01ba-4e1d-bbce-cf8d28642ca2
# ╟─576ddad7-10a6-4549-80cc-18d0c6793024
# ╟─fad7f786-7f00-4c99-a208-dcafcb50ef56
# ╟─21304a7d-22a0-4ef4-aa59-7d66930b850f
# ╠═d5f1af06-31b6-4ec9-a788-df58f2ee892f
# ╠═60772fa6-fdc9-4a02-bf58-1eb1fbf5d156
# ╠═bc997051-e470-4ce0-8ff9-dd09ca5d229d
# ╠═e339d118-74bb-4681-bb12-b94dc6a502a9
# ╠═bf0a5a0c-32e1-47e6-bc00-66cdafd2dfde
# ╠═b3e9b8d9-91b1-47bf-9fde-b480561948e4
# ╠═fd3849bd-80b1-43ea-9319-e36341f2f8c7
# ╟─13ef8a7e-7749-4102-8f0b-2ce0265d0906
# ╟─3fcc4340-02d5-4988-9d5c-1e1e232b5dfa
# ╟─ecd51ba0-5871-46a3-9b24-2c608b27af54
# ╟─d29863ce-5612-424c-a563-b8dcc94cf303
# ╟─ce5b1ffa-b789-445f-8792-ae3c8d1aa60d
# ╟─efedcfaa-d90a-49f1-ae13-cb6c93bfcd7b
# ╟─8b3b6194-3448-4605-a366-fc51aa2dc453
# ╠═7e467805-7bcf-4221-be05-27a1e4d89249
# ╠═ad048639-aebb-4d1c-a09f-3e9a3c9decb5
# ╠═eb2f9c36-0920-40ae-bc92-45333ab5385e
# ╠═4d3eeb0f-b981-4979-a220-620cd5d7bef1
# ╠═df6f5bdd-b442-4501-8ab7-05a520dd016f
# ╠═4c3fdce3-be2b-4196-907d-e7600c74303e
# ╠═b24ea01e-77c7-4d14-a712-ec5e3f37a0a4
# ╟─3ec7ed5c-4369-44d4-9c0c-7527415b2de1
# ╟─07562017-cbdd-40b2-bb1e-592e8abb8feb
# ╟─16d4de6b-4dc9-44c4-9a6b-bfd1252df846
# ╟─6af041e4-7476-44a6-af15-87bf4a1decbc
# ╟─b94c008d-4a41-468d-8e65-a42e586f83ce
# ╠═9195232e-556c-47d5-820b-c492916eed54
# ╠═77514283-cf30-4362-be95-d3f8e8018833
# ╠═9d771e98-4e0b-418d-acc5-ca9e0dafa06a
# ╠═2cd63880-9111-4ce2-97cf-4e615990b618
# ╠═5c376243-3c38-449d-8ae1-d9ecdcc44aaf
# ╠═5e932f76-ca1a-4305-9a1c-a8538fdef710
# ╠═c60880f2-3073-4a7e-8c2c-0979a3e254d3
# ╟─b8a2c430-0c1a-42c0-9a38-5b1af3f12ba6
# ╟─513b431d-5246-4c48-adce-afcd34b957a4
# ╟─a2f277f0-a017-45bd-b573-88cb06bd43a8
# ╟─2dc6c4db-0cb4-4f0c-b700-b990d48704b1
# ╟─3e437ffb-0c97-473a-bff6-929f8c720abd
# ╟─3fdeaaaf-2f00-4ef8-b6bd-cc1d29dbcb74
# ╟─04ea0005-24c6-4c39-816e-94a142fe65b8
# ╠═a215c938-11e3-495d-92cb-dc0559ed86f1
# ╠═e45f8e09-b588-4e84-98c0-602f96f2c649
# ╟─284b014d-1cfe-45c6-9c58-2a7144eee417
# ╟─b669a28e-642f-44ad-94bf-8e61af14d2be
# ╟─fa1bdf35-c482-4869-b195-10b6e4861a2f
# ╠═fd0357ee-086b-4394-add3-fdff0de4cd39
# ╟─2c6fbe92-5a73-4360-8ed2-f51f59dc5e3f
# ╟─c1d54755-4b3f-4ad8-aac5-af7b365500c5
# ╟─bc8ada92-8975-403e-a94e-ddc2d9ead536
# ╠═9e88185c-9c1f-4485-b465-30b3140b0e58
# ╟─e63719e6-d832-490f-ad1f-e903ec979c18
# ╟─49da5db2-1ff5-41c2-b0b8-15b2444e0522
# ╠═7fa5a22a-49be-455d-b37d-f2683b29d5ed
# ╠═9453fe67-731e-4832-918a-34a638d24042
# ╠═1404ac62-c99c-4ee4-99b8-2e9b1e64608e
# ╟─cc251a0f-c9bb-46f8-9640-d8670bff4526
# ╟─0575e66a-2008-42a4-9e05-3d8e3b2d6561
# ╠═bd5a7a2c-75cd-4a67-88a2-e09cf152f6a9
# ╠═25822bc1-6e1b-41e8-928f-4d088f0eede5
# ╠═0aaeeb0d-2740-4817-83a4-fd56d9f2c9ce
# ╠═92293f19-535a-4b35-8d58-bf99cddabcd3
# ╟─60b1f1ac-4fa5-48f6-b1fa-1c50a00447b3
# ╟─65d7b527-baa6-4af2-972f-1960c0072d10
# ╟─29d58861-8237-4cc9-a05a-c37b8228efaf
# ╟─6f3f7e16-5611-43a6-a5ca-5d25bd3f29f5
# ╠═6b2af110-702a-4572-9ead-24d2b90db0f6
# ╟─b120b287-060f-417f-902b-3ef4bc26148b
# ╠═80b491c3-7c54-4a8c-8460-c76a82ee7b86
# ╠═bd6060f6-0e09-4ccd-a9ed-12d7b63b54f2
# ╟─dfe42b5b-9c74-426e-bfdc-b21748179f59
# ╟─8b807632-3dd5-4b03-a539-1f348af4852c
# ╟─740ebde2-5ae2-476d-a7c8-b6541be938f5
# ╠═99598049-9494-4a72-afd4-13dd12a09ff4
# ╟─d9b1abef-4cae-4171-904f-6e4a046a9bf0
# ╟─210a4e3c-6998-400f-9e0a-3da3e0efe557
# ╠═0e5c7037-aa27-44f2-b563-94aca1da290b
# ╠═f3442668-a393-4136-8d43-d597af23279b
# ╠═0f92fa2a-7f2f-4a04-bcec-5eaf307ed3f3
# ╠═ab9353b8-4b9d-4e93-b169-aafb0a44247b
# ╟─0b489280-582d-47af-aed0-ac54a79a060c
# ╟─6f03a4bb-63e2-46fc-a39a-680c67dc82b1
# ╠═c58fef48-f544-4377-97f4-7daf9ba218ac
# ╠═aa05228a-9c86-4276-bdbf-c47479ccc7dc
# ╟─921e3b89-4baf-4984-85b1-087af8e8fc2c
# ╟─646546e9-8e80-428d-af73-02e3b033e092
# ╟─823716ae-f1fd-4805-be6e-a9bdb3827940
# ╟─e6f3d47b-fb81-4100-8343-a8c822c67137
# ╠═5cc34af8-b49e-401c-8320-3a18bd8bcbd3
# ╠═d17f4ee1-b3ae-49e6-80ee-26a322b2d661
# ╠═c9e781bb-1f8d-4a67-9db1-e142aab69e26
# ╟─67ed94c0-acf1-4dfa-b0ce-165d0b8e9a9d
# ╟─a32770f6-8fc2-482d-9326-8476465b7ce3
# ╟─8d599135-649c-4418-8937-511dc02bb97a
# ╠═1f039bf1-85ae-489c-b614-a2b4e24a73a0
# ╠═67219339-4d5d-4bb7-a76e-0465e8aebe62
# ╠═37bd6ab5-45d9-464a-874d-2ff12d869285
# ╠═ac58c4eb-b800-48cc-8a95-96f7b7917f85
# ╟─2ec865f9-0cc1-4147-827a-f4e4c38c44e9
# ╟─254cb2e0-4580-4858-a1fa-098038930b37
# ╟─4ee57a5a-be38-4a56-a686-95f9582c0910
# ╟─35e906c9-d6a2-4424-baed-186dd4f85f8d
# ╟─ca1c8ba9-ef43-4ca1-8ea4-f2e08415ed75
# ╟─c615896c-cdb9-4cfa-af81-6ee3ad6b3e02
