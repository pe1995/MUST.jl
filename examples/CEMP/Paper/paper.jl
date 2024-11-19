### A Pluto.jl notebook ###
# v0.20.0

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
dispatchPath(runname) = MUST.@in_dispatch(joinpath("shipped_models", runname))

# ╔═╡ c7a6739a-8ae7-45c0-ab4f-e800a9ab1d2e
runs3D = Dict(
	"t57.50g45.00m-5.000 CEMP" => dispatchPath("shipped_P3G_E_t57.50g45.00m-5.000_v1.0"),
	"t57.50g45.00m-5.000 no-CEMP" => dispatchPath("shipped_P4G_E_t57.50g45.00m-5.000_v1.0"),
	"t57.50g45.00m-4.000 CEMP" => dispatchPath("shipped_P5G_E_t57.50g45.00m-4.000_v1.0"),
	"t57.50g45.00m-4.000 no-CEMP" => dispatchPath("shipped_P6G_E_t57.50g45.00m-4.000_v1.0"),
)

# ╔═╡ a343cf41-e425-4a59-81bd-87f8765251c0
snapshots3D = Dict(
	"t57.50g45.00m-5.000 CEMP" => 371,
	"t57.50g45.00m-5.000 no-CEMP" => 371,
	"t57.50g45.00m-4.000 CEMP" => 371,
	"t57.50g45.00m-4.000 no-CEMP" => 371
)

# ╔═╡ cd8b84b0-6891-4d3e-ba04-22fd4890e09d
models3D = Dict(k=>pick_snapshot(runs3D[k], snapshots3D[k]) for k in keys(runs3D))

# ╔═╡ d2884a6e-b1e5-4234-ae78-05f682ac6908


# ╔═╡ 863499a2-1e33-4ba4-8dad-935a9a38a4ec
md"## 1D Models"

# ╔═╡ 9d8761fd-4880-4ff3-855f-f3713a55a6e5
runs1D = Dict(
	# feh = -5
	"MARCS t57.50g45.00m-5.000 no-CEMP" => dispatchPath("shipped_P4G_E_t57.50g45.00m-5.000_v1.0/p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),

	# feh = -4
	"MARCS t57.50g45.00m-4.000 no-CEMP" => dispatchPath("shipped_P6G_E_t57.50g45.00m-4.000_v1.0/p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"),
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
		model = MUST.readdlm(joinpath(dir, name), skipstart=2)

		m = TSO.flip!(Model1D(
			model[:, 1], 
			log.(model[:, 4]), 
			log.(model[:, 2]),
			fill!(similar(model[:, 1]), 0.0),
			-99.
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
	spectraPath(rundir, snapshot, extension, dimension="m3dis_$(snapshot)") = joinpath(
		rundir, "spectra_$(dimension)_lam_4297-4303$(extension)"
	)
	spectraPath(k, extension) = spectraPath(runs3D[k], snapshots3D[k], extension)
	spectraPathSnap(k, snapshot, extension) = spectraPath(runs3D[k], snapshot, extension)
	spectraPath1D(k, extension) = spectraPath(dirname(runs1D[k]), "", extension, basename(runs1D[k]))
end

# ╔═╡ 6613157c-d199-4825-b779-a477ac4f3f30
spectrapath3D = Dict(
	# feh=-5 models
	"t57.50g45.00m-5.000 CEMP - C3.0" => 
	spectraPath("t57.50g45.00m-5.000 CEMP","C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C0.0" => 
	spectraPath("t57.50g45.00m-5.000 CEMP","C_0.0"),
	"t57.50g45.00m-5.000 no-CEMP - C3.0" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_3.0"),
	"t57.50g45.00m-5.000 no-CEMP - C0.0" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_0.0"),

	# feh=-5 models - 3D effect fitting
	"t57.50g45.00m-5.000 no-CEMP - C2.8" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_2.8"),
	"t57.50g45.00m-5.000 no-CEMP - C2.6" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_2.6"),
	"t57.50g45.00m-5.000 no-CEMP - C2.4" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_2.4"),
	"t57.50g45.00m-5.000 no-CEMP - C2.2" => 
	spectraPath("t57.50g45.00m-5.000 no-CEMP" ,"C_2.2"),
	
	# feh=-4 models
	"t57.50g45.00m-4.000 CEMP - C0.75" => 
	spectraPath("t57.50g45.00m-4.000 CEMP","C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.0" => 
	spectraPath("t57.50g45.00m-4.000 CEMP","C_0.0"),
	"t57.50g45.00m-4.000 no-CEMP - C0.75" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.75"),
	"t57.50g45.00m-4.000 no-CEMP - C0.0" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.0"),

	# feh=-4 models - 3D effect fitting
	"t57.50g45.00m-4.000 no-CEMP - C0.40" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.40"),
	"t57.50g45.00m-4.000 no-CEMP - C0.35" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.35"),
	"t57.50g45.00m-4.000 no-CEMP - C0.30" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.30"),
	"t57.50g45.00m-4.000 no-CEMP - C0.25" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.25"),
	"t57.50g45.00m-4.000 no-CEMP - C0.20" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.20"),
	"t57.50g45.00m-4.000 no-CEMP - C0.15" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.15"),

	"t57.50g45.00m-4.000 no-CEMP - C0.60" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.60"),
	"t57.50g45.00m-4.000 no-CEMP - C0.68" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.68"),
	"t57.50g45.00m-4.000 no-CEMP - C0.70" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.70"),
	"t57.50g45.00m-4.000 no-CEMP - C0.72" => 
	spectraPath("t57.50g45.00m-4.000 no-CEMP" ,"C_0.72"),

	# Contribution functions
	"t57.50g45.00m-5.000 CEMP - contr" => 
	spectraPath("t57.50g45.00m-5.000 CEMP" ,"contr"),
	"t57.50g45.00m-4.000 CEMP - contr" => 
	spectraPath("t57.50g45.00m-4.000 CEMP" ,"contr"),

	# Spectra of other snapshots
	"t57.50g45.00m-5.000 CEMP - C3.0 - 375" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",375,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 374" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",374,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 373" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",373,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 372" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",372,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 369" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",369,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 368" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",368,"C_3.0"),
	"t57.50g45.00m-5.000 CEMP - C3.0 - 367" =>spectraPathSnap("t57.50g45.00m-5.000 CEMP",367,"C_3.0"),

	"t57.50g45.00m-4.000 CEMP - C0.75 - 373" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",373,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 372" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",372,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 370" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",370,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 369" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",369,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 368" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",368,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 367" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",367,"C_0.75"),
	"t57.50g45.00m-4.000 CEMP - C0.75 - 366" =>spectraPathSnap("t57.50g45.00m-4.000 CEMP",366,"C_0.75"),
	
)

# ╔═╡ e67fed9e-3a82-45ed-9865-1c5862ebba2c
spectra3D = Dict(
	k=>M3DISRun(v) for (k, v) in spectrapath3D
)

# ╔═╡ f3493218-6013-4567-b1c4-442d0a48c9d0


# ╔═╡ 63ed5788-011d-4a31-ac82-efcd135f609c
spectrapath1D = Dict(
	# Feh = -5
	"MARCS t57.50g45.00m-5.000 no-CEMP - C3.0" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_3.0",),
	"MARCS t57.50g45.00m-5.000 no-CEMP - C0.0" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","C_0.0"),

	# feh = -4
	"MARCS t57.50g45.00m-4.000 no-CEMP - C0.75" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_0.75",),
	"MARCS t57.50g45.00m-4.000 no-CEMP - C0.0" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","C_0.0"),

	# contribution functions
	"MARCS t57.50g45.00m-5.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-5.000 no-CEMP","contr"),
	"MARCS t57.50g45.00m-4.000 no-CEMP - contr" => spectraPath1D("MARCS t57.50g45.00m-4.000 no-CEMP","contr")
)

# ╔═╡ 4bf4013c-ee0d-42fd-9dd2-4f26f090c15b
spectra1D = Dict(
	k=>M3DISRun(v) for (k, v) in spectrapath1D
)

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

# ╔═╡ cbb98c8a-8ec3-46fa-a38f-4afc3d292f41
let
	f, ax = plt.subplots(1, 1, figsize=(6,5))

	for (i, p) in enumerate(design_p1)
		k = p.key
		x, y = profile(MUST.mean, models3D[k]|>first, :z, :T)
		x ./= 1e5
		ax.plot(x, y, color=p.color, ls=p.ls, lw=p.lw, label=p.label)
	end

	ax.set_xlim(-50, 230)
	ax.set_ylim(3000, 7000)
	ax.set_xlabel("geometrical height [km]")
	ax.set_ylabel("temperature [K]")

	ax.legend()

	f.savefig(name1_p1)

	f
end

# ╔═╡ 56d11775-333b-44ac-ba5f-39e1b372564d
name2_p1 = "CEMP_vs_noCEMP_d_z_t57.50g45.00m-5.000.pdf"

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

# ╔═╡ 773a7896-a492-46a3-9813-fa739647142d


# ╔═╡ 244f1b18-80fe-419a-8fd0-892aad803ca5
md"## 3D metallicity effect"

# ╔═╡ accdaadc-8912-48d2-8a0c-3c19e5b4ba19
design_p2 = [
	PlotDesign("t57.50g45.00m-4.000 no-CEMP", "DISPATCH <3D>", "-", 2.5, "tomato"),
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

	ax.set_xlim(-4, 4)
	ax.set_ylim(3200, 12000)
	ax.set_xlabel("optical depth ["*L"\rm\tau_{ross}"*"]")
	ax.set_ylabel("temperature [K]")

	ax.legend(loc="upper left")
	f.savefig(name1_p2)

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
	"t57.50g45.00m-5.000 CEMP - C3.0"=>[
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 375", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 374", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 373", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 372", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 369", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 368", "", "-", 1.5, "0.5"), 
		PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0 - 367", "", "-", 1.5, "0.5"), 
	],
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
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p3)
	
	f
end

# ╔═╡ 2aa33bb2-3cd3-4104-a6a2-0f1b12fe2960


# ╔═╡ c905fd9a-e2f6-485c-9096-c6a588b31b4d
md"### [Fe/H]=-4, [C/Fe]=0.75"

# ╔═╡ 91699cc6-cb81-4382-bc5a-01931914343b
design_p3B = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - C0.75", "CEMP model", "-", 1.5, "tomato"),
	PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.75", "scaled-solar model", "--", 1.5, "k")
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
	
	
	ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p3)
	
	f
end

# ╔═╡ 1397cb72-51db-45bf-8a67-cf1352f7f990


# ╔═╡ bd3ee310-3f68-4e08-8ba8-d881aeecb76a
md"## 3D metallicity effect"

# ╔═╡ f064d1a4-be5b-4cfc-9bf9-67002a333d2d
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 9987b27b-2437-45a3-9e83-116a01f13d1d
design_p4 = [
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", "MARCS", "--", 1.5, "k")
]

# ╔═╡ 17ced5c3-080f-4d65-aa1c-b24a2950e783
name1_p4 = "vs_1D_spectrum_t57.50g45.00m-5.000.pdf"

# ╔═╡ 75600b27-d640-438b-bd00-d35e18ff46f6
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(10, 4))
	
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
	
	
	ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)
	f.savefig(name1_p4)
	
	f
end

# ╔═╡ 59ab2b58-2e38-4cce-9ba4-bc6899616a7d


# ╔═╡ a61f7ec3-88c9-42f7-bf5c-acc792a05fb5
md"### [Fe/H]=-4, [C/Fe]=0.75"

# ╔═╡ 2aa8d303-dae3-4a23-b253-c6091e2a638d
design_p4B = [
	PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.75", "DISPATCH 3D", "-", 1.5, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C0.75", "MARCS", "--", 1.5, "k")
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
	
	
	ax.set_ylim(0.68, 1.01)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel("flux [normalized]")
	ax.set_xlabel("wavelength [Å]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p4)
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
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "CEMP model", "-", 2, ""),
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "scaled-solar model", "--", 2, "")
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
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "DISPATCH 3D", "-", 2, "tomato"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", "MARCS", "--", 2, "steelblue")
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

# ╔═╡ 483897d7-fbf3-4a3d-a32b-2992ce5fccd2
md"### [Fe/H]=-5, [C/Fe]=3.0"

# ╔═╡ 7021b190-4a8b-4ae4-bd71-06b3840f20c6
abundances_p6 = Dict(
	2.2 => PlotDesign("t57.50g45.00m-5.000 no-CEMP - C2.2", "", "-", 1.5, nothing),
	2.4 => PlotDesign("t57.50g45.00m-5.000 no-CEMP - C2.4", "", "-", 1.5, nothing),
	2.6 => PlotDesign("t57.50g45.00m-5.000 no-CEMP - C2.6", "", "-", 1.5, nothing),
	2.8 => PlotDesign("t57.50g45.00m-5.000 no-CEMP - C2.8", "", "-", 1.5, nothing),
	3.0 => PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "", "-", 1.5, nothing)
)

# ╔═╡ b7a0779f-bd55-4ec3-afc8-ea8caf9db6fc
reference_p6 = [
	(3.0, PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", "DISPATCH 3D - CEMP", "s", 8.5, "tomato")),
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

	x = range(2.05, 3.15, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")
	
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
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.07*(maximum(ew)-minimum(ew)), clab, ha="center", va="bottom", color=p.color, 
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end
	
	
	#ax.set_ylim(-0.05, 1.01)
	ax.set_xlim(2.05, 3.15)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name1_p6)
	
	f
end

# ╔═╡ 38ce6bb7-f6f2-44d3-8685-8049c42ffeee


# ╔═╡ 2cd6317a-97b4-45bb-8e87-b09b40a470f7
md"### [Fe/H]=-4, [C/Fe]=0.75"

# ╔═╡ de419d2f-eb6d-47cd-a69c-99d5cb3161a3
abundances_p6B = Dict(
	0.75 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.75", "", "-", 1.5, nothing),
	0.60 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.60", "", "-", 1.5, nothing),
	0.68 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.68", "", "-", 1.5, nothing),
	0.70 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.70", "", "-", 1.5, nothing),
	0.72 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.72", "", "-", 1.5, nothing),
	0.40 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.40", "", "-", 1.5, nothing),
	0.35 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.35", "", "-", 1.5, nothing),
	0.30 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.30", "", "-", 1.5, nothing),
	0.25 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.25", "", "-", 1.5, nothing),
	0.20 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.20", "", "-", 1.5, nothing),
	0.15 => PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.15", "", "-", 1.5, nothing)
)

# ╔═╡ 135f931b-387f-420d-8a0e-d7f391c1569f
reference_p6B = [
	(0.75, PlotDesign("t57.50g45.00m-4.000 CEMP - C0.75", "DISPATCH 3D - CEMP", "s", 8.5, "tomato")),
	(0.75, PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C0.75", "MARCS", "s", 8.5, "steelblue"))
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

	x = range(0.07, 0.78, length=100) |> collect
	y = f_ew.(x)
	ax.plot(x, y, lw=2, color="k", ls="-")
	ax.plot(abund, ew, lw=2, marker="s", color="k", markersize=10, ls="", markerfacecolor="white")
	
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
		clab = L"\rm \Delta_{[C/Fe]}\ =\ "*"$(ΔC)" *L"\rm\ dex"
	
		ax.plot([aref], [EW], color=p.color, marker=p.ls, markersize=p.lw, label=p.label)
		ax.axhline(EW, color=p.color, lw=1.5, ls="--")
		ax.axvline(aref + ΔC, color=p.color, lw=1.5, ls="--")
		ax.text(
			aref + ΔC, EW + 0.08*(maximum(ew)-minimum(ew)), clab, ha="center", va="bottom", color=p.color,
			bbox=Dict("facecolor"=>"white", "edgecolor"=>p.color, "alpha"=>0.7)
		)
	end
	
	
	ax.set_ylim(0.08, 0.55)
	ax.set_xlim(0.07, 0.9)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm [C/Fe]")
	ax.set_ylabel(L"\rm equivalent\ width\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.95))
	#ax.legend(loc="lower center", ncol=2)

	f.savefig(name2_p6)
	
	f
end

# ╔═╡ 2c44e5cc-43e9-4181-ab36-d629f4c61080


# ╔═╡ 349d1a5b-aed0-4315-b7b4-cc766e2d0bc3
md"# Yoon-Beers 2019"

# ╔═╡ 8cc66430-8635-4729-b595-a7df6a28cba5
yb19 = MUST.readdlm("../yoon_beers19.txt", ';')

# ╔═╡ 3bea6988-24ec-46bf-b20d-1f9817e7cc0b
feh = yb19[4:end, 6]

# ╔═╡ ac78dfbc-1423-4741-b0b4-dcd8c809787f
c = yb19[4:end, 8]

# ╔═╡ b5896895-1bd7-40de-b86b-b326d502e587
begin
	p1 = (-4.0, 0.75-4+8.560, -0.55)
	p2 = (-5.0, 3.0-5+8.560, -0.79+0.15)
	
	#=scatter_int(x, y) = MUST.pyconvert(typeof(x),
		MUST.scipy_interpolate.griddata(
		([p1[1], p2[1]], 
		[p1[2], p2[2]]), 
		[p1[3], p2[3]], 
		(x, y), method="linear")
	)
	corrections = scatter_int(feh, c)
	mask = isnan.(corrections)=#
	scatter_int2(x, y) = MUST.pyconvert(typeof(x),
		MUST.scipy_interpolate.griddata(
		([p1[1], p2[1]], 
		[p1[2], p2[2]]), 
		[p1[3], p2[3]], 
		(x, y), method="nearest")
	)
	corrections = scatter_int2(feh, c)
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
	ax.scatter(feh, c+corrections, s=10, marker="o", c="tomato", label="corrected")

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
	PlotDesign("t57.50g45.00m-5.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ d5af99b0-8e95-4037-af01-1437f9aa0e04
differences_p8 = [
	PlotDesign("t57.50g45.00m-5.000 no-CEMP - C3.0", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t57.50g45.00m-5.000 CEMP - C3.0", L"\rm not\ CEMP - CEMP", "-", 2., "k"),
	PlotDesign("MARCS t57.50g45.00m-5.000 no-CEMP - C3.0", L"\rm 3D - 1D", "--", 2., "tomato"),
]

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

	left, bottom, width, height = [0.22, 0.52, 0.35, 0.3]
	ax2 = f.add_axes([left, bottom, width, height], sharex=ax)

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
	ax2.legend()

	ax.axvline(-3.5, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(3000, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

	ax.legend(loc="lower right")

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


# ╔═╡ 6225b614-dc48-4158-8abe-55a1c07b651e
md"## [Fe/H]=-4, [C/Fe]=0.75"

# ╔═╡ b04303d1-bd59-4104-bc7c-8fdef9477fee
spectra_p9 = [
	PlotDesign("t57.50g45.00m-4.000 CEMP - contr", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - contr", "MARCS", "--", 2.5, "tomato"),
]

# ╔═╡ 901d430e-254b-4774-b7cc-4dcbd8ea1772
differences_p9 = [
	PlotDesign("t57.50g45.00m-4.000 no-CEMP - C0.75", "DISPATCH <3D>", "-", 2, "k"),
	PlotDesign("t57.50g45.00m-4.000 CEMP - C0.75", L"\rm not\ CEMP - CEMP", "-", 2., "k"),
	PlotDesign("MARCS t57.50g45.00m-4.000 no-CEMP - C0.75", L"\rm 3D - 1D", "--", 2., "tomato"),
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

	ax.axvline(-3.0, ls="-", color="k", alpha=0.1, lw=4)
	ax.axvline(-1.2, ls="-", color="tomato", alpha=0.1, lw=4)
	ax.set_xlim(-4.5, 3)
	ax2.set_xlim(-4.5, 3)
	ax.set_ylim(3000, 11000)
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.set_xlabel(L"\rm optical\ depth\ [log\ \tau_{500}]")

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

	cfe_all, eufe_all, bafe_all, feh_all = saga_all_data[:, 9], saga_all_data[:, 10], saga_all_data[:, 11], saga_all_data[:, 12]

	parameters_all = Dict{Any, Any}(
		"teff" => teff_all, 
		"logg" => logg_all, 
		"z"    => z_all, 
		"cfe"  => cfe_all, 
		"eufe" => eufe_all, 
		"bafe" => bafe_all, 
		"feh"  => feh_all
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
logg_limit = 3.5

# ╔═╡ 92293f19-535a-4b35-8d58-bf99cddabcd3
selection_mask = parameters_all["logg"] .>= logg_limit;

# ╔═╡ 60b1f1ac-4fa5-48f6-b1fa-1c50a00447b3
let
	f, ax = plt.subplots(figsize=(6, 5))

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
		s=5, cmap="rainbow_r",
		rasterized=true
	)
	cbar = f.colorbar(im, ax=ax)

	xlim = MUST.pyconvert(Array, ax.get_xlim())
	ylim = MUST.pyconvert(Array, ax.get_ylim())
	ax.set_xlim(xlim[end:-1:1])
	ax.set_ylim(ylim[end:-1:1])

	ax.axhline(logg_limit, color="0.5", ls="--", alpha=0.5, lw=2)

	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	ax.set_ylabel(L"\rm log(g)\ [K]")
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
	p1_saga_all = (-4.0, -0.55)
	p2_saga_all = (-5.0, -0.79+0.15)
	
	ip_saga_all = MUST.linear_interpolation(
		[p2_saga_all[1],p1_saga_all[1]],
		[p2_saga_all[2],p1_saga_all[2]],
		extrapolation_bc=MUST.Flat()
	)
	
	corrections_saga_all = ip_saga_all.(parameters_all["cfe"])
end

# ╔═╡ d9b1abef-4cae-4171-904f-6e4a046a9bf0


# ╔═╡ 210a4e3c-6998-400f-9e0a-3da3e0efe557
md"## Galactic CEMP Distribution"

# ╔═╡ 0e5c7037-aa27-44f2-b563-94aca1da290b
metallicity_bin_edges = [-7, -6, -5, -4, -3.5, -3, -2.5, -2., -1.5, -1., -0.5, 0.0]

# ╔═╡ 775ac5e9-6933-4d6d-a1a2-e40ec840cae4
begin
	metallicity_bin_centers = (metallicity_bin_edges[2:end] .+ metallicity_bin_edges[1:end-1]) ./ 2
	bin_masks = falses(size(saga_all_data, 1), length(metallicity_bin_centers))
	for i in eachindex(metallicity_bin_centers)
		bin_masks[:, i] .= (metallicity_bin_edges[i] .< parameters_all["feh"]) .&
		(parameters_all["feh"] .< metallicity_bin_edges[i+1])

		bin_masks[:, i] .= bin_masks[:, i] .& selection_mask
	end
	count_bins = count.([bin_masks[:, i] for i in axes(bin_masks, 2)])
	@info "Bin counts (all): $(count_bins)"

	# now the same for the CEMP stars
	bin_masks_CEMP = falses(size(saga_all_data, 1), length(metallicity_bin_centers))
	for i in eachindex(metallicity_bin_centers)
		bin_masks_CEMP[:, i] .= bin_masks[:, i] .& (parameters_all["cfe"] .> 0.7)
	end

	count_bins_CEMP = count.([bin_masks_CEMP[:, i] for i in axes(bin_masks_CEMP, 2)])
	@info "Bin counts (CEMP): $(count_bins_CEMP)"


	# now the same for the CEMP stars (after correction)
	bin_masks_CEMP2 = falses(size(saga_all_data, 1), length(metallicity_bin_centers))
	for i in eachindex(metallicity_bin_centers)
		bin_masks_CEMP2[:, i] .= bin_masks[:, i] .& (parameters_all["cfe"] .+ corrections_saga_all .> 0.7)
	end

	count_bins_CEMP2 = count.([bin_masks_CEMP2[:, i] for i in axes(bin_masks_CEMP2, 2)])
	@info "Bin counts (CEMP corrected): $(count_bins_CEMP2)"
end

# ╔═╡ 921e3b89-4baf-4984-85b1-087af8e8fc2c
let
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	ax.plot(metallicity_bin_centers, count_bins_CEMP ./ count_bins, zorder=10, color="k", lw=2.0, label=L"\rm SAGA\ (MS)", marker="s")
	ax.plot(metallicity_bin_centers, count_bins_CEMP2 ./ count_bins, zorder=10, color="tomato", lw=2.0, label=L"\rm SAGA\ (MS)\ -\ corrected", marker="s")

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

	ax.set_xlim(-6.8, 0.1)
	ax2.set_ylim(-0.75, 5.)
	f
end

# ╔═╡ Cell order:
# ╠═32b14b00-8ae8-11ef-2aa2-937a4b8d3a0c
# ╠═bbb5fb60-6669-440c-84d6-5e0edbd90e3d
# ╠═024bd32f-a8c3-43eb-bba3-05d590f9632f
# ╟─9967494d-13ad-4997-95c9-341ca3402ded
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
# ╟─5ac9c28b-be99-4c27-9bdd-3aff553f9bc3
# ╠═e4eef27a-65c4-4ff5-a3e7-b0119a7e33c5
# ╠═67c8479d-7c41-4baa-b2ec-13da36a48900
# ╟─4bc7b154-2fdd-4057-bf40-e44c7851bc7d
# ╟─0a0779a4-7d87-40d2-9dfc-454b366c2dff
# ╠═5beedaf9-b12a-46db-a576-7913eb5f2caf
# ╠═6613157c-d199-4825-b779-a477ac4f3f30
# ╠═e67fed9e-3a82-45ed-9865-1c5862ebba2c
# ╟─f3493218-6013-4567-b1c4-442d0a48c9d0
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
# ╟─cbb98c8a-8ec3-46fa-a38f-4afc3d292f41
# ╠═56d11775-333b-44ac-ba5f-39e1b372564d
# ╟─c421ceca-ee21-455c-bfec-371c6c260ea3
# ╟─773a7896-a492-46a3-9813-fa739647142d
# ╟─244f1b18-80fe-419a-8fd0-892aad803ca5
# ╠═accdaadc-8912-48d2-8a0c-3c19e5b4ba19
# ╠═ffb5a590-77e5-41c5-8964-eafe1153a5dc
# ╟─778a7f04-0841-4801-8aa8-b5686e1e498c
# ╟─83729ef7-9c05-4c0b-8bcf-0ecfccc48d66
# ╟─0d4a0a2f-758f-44a0-9692-23c082befac8
# ╟─1f6dfdaf-b82e-46ee-989f-dbeb685d72ff
# ╟─391e98d1-344e-4ec3-837b-58327c02c74b
# ╠═10d7d97f-5c54-4066-9618-c78983c6c18c
# ╠═b97665fe-1536-4047-b3e9-d6176d47f3c1
# ╠═a04336bf-9cdf-401e-b237-a41f44987320
# ╟─1deeb995-4fda-4071-8d1f-fe9a9bc8a09c
# ╟─2aa33bb2-3cd3-4104-a6a2-0f1b12fe2960
# ╟─c905fd9a-e2f6-485c-9096-c6a588b31b4d
# ╠═91699cc6-cb81-4382-bc5a-01931914343b
# ╠═aeccdb4b-63f9-47aa-830e-21aee1093dc4
# ╟─a55d5b04-a41b-47b6-81e8-39073bf538e2
# ╟─1397cb72-51db-45bf-8a67-cf1352f7f990
# ╟─bd3ee310-3f68-4e08-8ba8-d881aeecb76a
# ╟─f064d1a4-be5b-4cfc-9bf9-67002a333d2d
# ╠═9987b27b-2437-45a3-9e83-116a01f13d1d
# ╠═17ced5c3-080f-4d65-aa1c-b24a2950e783
# ╟─75600b27-d640-438b-bd00-d35e18ff46f6
# ╟─59ab2b58-2e38-4cce-9ba4-bc6899616a7d
# ╟─a61f7ec3-88c9-42f7-bf5c-acc792a05fb5
# ╠═2aa8d303-dae3-4a23-b253-c6091e2a638d
# ╠═cec8b45b-e4b6-4ad0-a83d-e9afe62b486a
# ╟─ff0787a1-b1cb-4eb0-85d7-b33bd39b2533
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
# ╟─483897d7-fbf3-4a3d-a32b-2992ce5fccd2
# ╠═7021b190-4a8b-4ae4-bd71-06b3840f20c6
# ╠═b7a0779f-bd55-4ec3-afc8-ea8caf9db6fc
# ╠═d9818823-8afb-4ad2-876d-24ff473a93e7
# ╟─f3f50db8-66b7-4d48-b067-bea39c9c934c
# ╟─38ce6bb7-f6f2-44d3-8685-8049c42ffeee
# ╟─2cd6317a-97b4-45bb-8e87-b09b40a470f7
# ╠═de419d2f-eb6d-47cd-a69c-99d5cb3161a3
# ╠═135f931b-387f-420d-8a0e-d7f391c1569f
# ╠═5f9b9dc4-f43b-4853-85d4-ae51247c13d0
# ╟─b0dfbe89-c719-4622-8f7b-55e5a09590df
# ╟─2c44e5cc-43e9-4181-ab36-d629f4c61080
# ╟─349d1a5b-aed0-4315-b7b4-cc766e2d0bc3
# ╠═8cc66430-8635-4729-b595-a7df6a28cba5
# ╠═3bea6988-24ec-46bf-b20d-1f9817e7cc0b
# ╠═ac78dfbc-1423-4741-b0b4-dcd8c809787f
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
# ╠═43b4cf33-7efd-48b8-af90-66555c827b52
# ╠═71545110-8314-45c7-9e83-e15035c3f60d
# ╠═266c3319-78c1-4a2b-b06c-386f54907ecc
# ╠═d1c307dd-ab6f-44bc-a0d5-a879960231a5
# ╠═57c74ea4-a4ff-4018-9507-2753fc2a440c
# ╟─2765daec-fb70-4b55-889f-0ceb0808e84e
# ╟─515b8fb8-2568-4200-b10f-6ef8e2318615
# ╟─d08caaf5-e763-498d-a594-8915fc8eaaf9
# ╟─1dcf4d03-7a20-4ff7-b891-958003c23d61
# ╟─a21d6f31-f766-432f-856c-7726e989d42f
# ╟─67aab606-2b0e-44f8-90ab-99d63830742d
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
# ╟─92293f19-535a-4b35-8d58-bf99cddabcd3
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
# ╟─775ac5e9-6933-4d6d-a1a2-e40ec840cae4
# ╟─921e3b89-4baf-4984-85b1-087af8e8fc2c
