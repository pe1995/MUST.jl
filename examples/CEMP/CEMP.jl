### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ cc6eeb7e-5974-11ef-0436-fd682756392a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using PythonPlot
	using LaTeXStrings
end

# ╔═╡ a73ba74d-30ca-4716-9a0c-455aa24cf975
begin
	@import_dispatch "../../../dispatch2"
	@import_m3dis "../../../Multi3D"
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
end;

# ╔═╡ 39703f8e-40b2-476a-8ccb-cf5d208866fa
md"# DISPATCH Models
Load a couple of Dispatch models that can be analysed."

# ╔═╡ 5e8a2581-4cb0-4337-a913-4c4df3121455
models = Dict(
	"metal-poor sun, CEMP" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP C0.0" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP C2.0" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP C2.5" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP C2.25" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, av CEMP C3.0" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, av CEMP C3.1" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, av CEMP C2.9" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun C3.1" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun C2.9" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun C2.8" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun C2.7" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c0.75" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c0.85" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c0.95" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c1.05" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c1.25" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c1.50" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c1.75" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-4 c2.75" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, CEMP m-4 c0.75" => "M8_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, CEMP" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, MARCS m-7 c5.0" => "M11.2_E_t57.50g45.00m-7.000_v1.0", 
	"T5750 g45, CEMP m-7 c5.0" => "M11.2_E_t57.50g45.00m-7.000_v1.0_hres", 
	"T5750 g45, CEMP m-7 c4.25" => "M11.2_E_t57.50g45.00m-7.000_v1.0",
	"T5750 g45, CEMP m-7 c4.5" => "M11.2_E_t57.50g45.00m-7.000_v1.0",
	"T5750 g45, CEMP m-7 c4.75" => "M11.2_E_t57.50g45.00m-7.000_v1.0",
	"T5750 g45, CEMP C2.1" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, CEMP C2.5" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, CEMP C3.2" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, CEMP C3.1" => "M3_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45, c0.4" => "M4_E_t57.50g45.00m-4.000_v1.0",
	"T5750 g45" => "M4_E_t57.50g45.00m-4.000_v1.0",
	"metal-poor sun MARCS 5500" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun MARCS 5750" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun MARCS 6000" => "M2_E_t57.77g44.40m-4.000_v1.0",
	"metal-poor sun, CEMP MARCS 5600" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5700" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5800" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772 C3.5" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772 C4.1" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772 C4.3" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772 C4.5" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
	"metal-poor sun, CEMP MARCS 5772 C4.7" => "M1_E_t57.77g44.40m-4.000_v1.0_long",
)

# ╔═╡ d0815aec-c693-4b7f-ab3d-2ba6b8f7176c
spectra_4297_4303_name(isnap, ext) = "spectra_sn$(isnap)_lam_4297-4303$(ext)"

# ╔═╡ cea7a484-48f5-4c38-9966-00e0f4409b50
spectra_4297_4303_name_new(isnap, ext) = "spectra_m3dis_$(isnap)_lam_4297-4303$(ext)"

# ╔═╡ f17e97d8-3239-46d5-9883-8124cb0170a3
snapshots = Dict(
	"metal-poor sun, CEMP" => 399,
	"metal-poor sun, CEMP C0.0" => 399,
	"metal-poor sun, CEMP C2.0" => 399,
	"metal-poor sun, CEMP C2.5" => 399,
	"metal-poor sun, CEMP C2.25" => 399,
	"metal-poor sun" => 399,
	"metal-poor sun C3.1" => 399,
	"metal-poor sun C2.9" => 399,
	"metal-poor sun C2.8" => 399,
	"metal-poor sun C2.7" => 399,
	"T5750 g45, CEMP m-4 c0.75" => 378,
	"T5750 g45, CEMP" => 398,
	"T5750 g45, CEMP C2.1" => 398,
	"T5750 g45, CEMP C2.5" => 398,
	"T5750 g45, CEMP C3.2" => 398,
	"T5750 g45, CEMP C3.1" => 398,
	"T5750 g45, c0.4" => 398,
	"T5750 g45" =>398,
	"T5750 g45, CEMP m-7 c5.0" => 597,
	"T5750 g45, CEMP m-7 c4.25" => 567,
	"T5750 g45, CEMP m-7 c4.5" => 567,
	"T5750 g45, CEMP m-7 c4.75" => 567,
)

# ╔═╡ 2d84b97b-b835-421b-8b8a-2800da47e130
extension = Dict(
	"metal-poor sun, CEMP" => "CEMP-CEMP",
	"metal-poor sun, CEMP C0.0" => "_C0.0",
	"metal-poor sun, CEMP C2.0" => "_C2.0",
	"metal-poor sun, CEMP C2.5" => "_C2.5",
	"metal-poor sun, CEMP C2.25" => "_C2.25",
	"metal-poor sun" => "MP-CEMP",
	"metal-poor sun C3.1" => "_C3.1",
	"metal-poor sun C2.9" => "_C2.9",
	"metal-poor sun C2.8" => "_C2.8",
	"metal-poor sun C2.7" => "_C2.7",
	"T5750 g45, CEMP m-4 c0.75"  => "",
	"T5750 g45, MARCS m-4 c0.85" => "_C_0.85",
	"T5750 g45, MARCS m-4 c0.95" => "_C_0.95",
	"T5750 g45, MARCS m-4 c1.05" => "_C_1.05",
	"T5750 g45, MARCS m-4 c1.25" => "_C_1.25",
	"T5750 g45, MARCS m-4 c1.50" => "_C_1.50",
	"T5750 g45, MARCS m-4 c1.75" => "_C_1.75",
	"T5750 g45, MARCS m-4 c2.75" => "_C_2.75",
	"T5750 g45, c0.4" => "_C_0.4",
	"T5750 g45, CEMP"  => "",
	"T5750 g45, CEMP C2.1"  => "_C2.1",
	"T5750 g45, CEMP C2.5"  => "_C2.5",
	"T5750 g45, CEMP C3.2"  => "_C3.2",
	"T5750 g45, CEMP C3.1"  => "_C3.1",
	"T5750 g45" => "",
	"T5750 g45, CEMP m-7 c5.0" => "",
	"T5750 g45, CEMP m-7 c4.25" => "_C_4.25",
	"T5750 g45, CEMP m-7 c4.5" => "_C_4.5",
	"T5750 g45, CEMP m-7 c4.75" => "_C_4.75",
)

# ╔═╡ b6f2f6e8-67e6-4bc8-936d-6d0f4a6a5111
models1Dnames = Dict(
	"T5750 g45, MARCS m-4 c0.75" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS" => "p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c0.85" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c0.95" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c1.05" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c1.25" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c1.50" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c1.75" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-4 c2.75" => "p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"T5750 g45, MARCS m-7 c5.0" => "p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"metal-poor sun MARCS 5500" => "p5500_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"metal-poor sun MARCS 5750" => "p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"metal-poor sun MARCS 6000" => "p6000_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
	"metal-poor sun, CEMP MARCS 5600" => "5600g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5700" => "5700g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5800" => "5800g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772 C4.1" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772 C4.3" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772 C4.5" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772 C4.7" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod",
	"metal-poor sun, CEMP MARCS 5772 C3.5" => "5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod"
)

# ╔═╡ a7ee3a6e-5855-4dc1-a683-7f32f1bccb18
begin
	spectranames = Dict(
		(k=>spectra_4297_4303_name_new(v, extension[k]) for (k, v) in snapshots)...,
		"metal-poor sun, av CEMP C3.0" => "spectra_av_sn399CEMP-CEMP_lam_4297-4303",
		"metal-poor sun, av CEMP C3.1" => "spectra_av_sn399CEMP-CEMP_lam_4297-4303C3.1",
		"metal-poor sun, av CEMP C2.9" => "spectra_av_sn399CEMP-CEMP_lam_4297-4303C2.9",
		"metal-poor sun MARCS 5500" => "spectra_p5500_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"metal-poor sun MARCS 5750" => "spectra_p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"metal-poor sun MARCS 6000" => "spectra_p6000_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"metal-poor sun, CEMP MARCS 5600" => "spectra_5600g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303",
		"metal-poor sun, CEMP MARCS 5700" => "spectra_5700g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303",
		"metal-poor sun, CEMP MARCS 5772" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303",
		"metal-poor sun, CEMP MARCS 5800" => "spectra_5800g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303",
		"metal-poor sun, CEMP MARCS 5772 C3.5" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303_C3.5",
		"metal-poor sun, CEMP MARCS 5772 C4.1" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303_C4.1",
		"metal-poor sun, CEMP MARCS 5772 C4.3" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303_C4.3",
		"metal-poor sun, CEMP MARCS 5772 C4.5" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303_C4.5",
		"metal-poor sun, CEMP MARCS 5772 C4.7" => "spectra_5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod_lam_4297-4303_C4.7",
		"T5750 g45, MARCS m-7 c5.0" => "spectra_p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"T5750 g45, MARCS" => "spectra_p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"T5750 g45, MARCS m-4 c0.75" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303",
		"T5750 g45, MARCS m-4 c0.85" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_0.85",
		"T5750 g45, MARCS m-4 c0.95" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_0.95",
		"T5750 g45, MARCS m-4 c1.05" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_1.05",
		"T5750 g45, MARCS m-4 c1.75" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_1.75",
		"T5750 g45, MARCS m-4 c2.75" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_2.75",
		"T5750 g45, MARCS m-4 c1.25" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_1.25",
		"T5750 g45, MARCS m-4 c1.50" => "spectra_p5750_g+4.5_m0.0_t02_st_z-4.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod_lam_4297-4303_C_1.50",
	)
end

# ╔═╡ 307e6d9a-d7ac-435c-b4eb-c1c9e7825110


# ╔═╡ 05895238-a947-429d-9ecb-751ebe199ca9


# ╔═╡ b737400b-0e76-40ae-98d4-01c00e8f2dd3
spectra_4297_4303 = Dict(
	k=>M3DISRun(joinpath(@in_dispatch(joinpath("data", models[k], v)))) for (k, v) in spectranames
)

# ╔═╡ 4c85536c-87d8-4d0b-87d1-b2fc32ca5118
cubes = Dict(
	k=>pick_snapshot(@in_dispatch(joinpath("data", models[k])), snapshots[k]) 
	for (k, m) in snapshots
)

# ╔═╡ 4ad62be3-8d9f-446e-bf61-7d77d0475d25


# ╔═╡ f00cc719-1090-46fc-ba84-31bb1351c726
md"Average the cubes for <3D> test:"

# ╔═╡ 44429920-277b-4a09-a8a4-a6e143056628
begin
	#=for (m, boxes) in cubes
		model_name = models[m]
		isnap = snapshots[m]
		newname = @in_dispatch(
			joinpath("data",model_name,"av_sn$(isnap)$(extension[m])")
		)
		short_name = "av_sn$(isnap)$(extension[m])"

		z, T = profile(MUST.mean, boxes[1], :z, :T)
		z, logr = profile(MUST.mean, boxes[1], :z, :log10d)

		Z, T, logr = -reverse(z), reverse(T), reverse(logr)
		open(newname, "w") do f
			write(f, short_name*"\n")
			write(f, "$(length(z))\n")
			for i in eachindex(z)
				write(f, "$(z[i]) $(T[i]) 1.0 $(exp10(logr[i])) 2.0 \n")
			end
		end
	end=#
end

# ╔═╡ 37bcf775-4ccc-44a7-af80-b10ffa8a26c7


# ╔═╡ 1f6071d5-ab91-4eef-9dd2-dea9237f1ec8
md"## Style"

# ╔═╡ 21e58231-22b5-48e6-a8aa-c4bbe30a8c4c
colors = Dict(
	"metal-poor sun, CEMP" => "k",
	"metal-poor sun, CEMP C0.0" => "k",
	"metal-poor sun, CEMP C2.0" => "k",
	"metal-poor sun, CEMP C2.5" => "k",
	"metal-poor sun, CEMP C2.25" => "k",
	"metal-poor sun, av CEMP C3.0" => "0.5",
	"metal-poor sun, av CEMP C3.1" => "r",
	"metal-poor sun, av CEMP C2.9" => "lime",
	"metal-poor sun" => "k",
	"metal-poor sun C3.1" => "lime",
	"metal-poor sun C2.9" => "magenta",
	"metal-poor sun C2.8" => "cyan",
	"metal-poor sun C2.7" => "lime",
	"T5750 g45, MARCS m-4 c0.75" => "r",
	"T5750 g45, MARCS m-4 c0.85" => "r",
	"T5750 g45, MARCS m-4 c0.95" => "r",
	"T5750 g45, MARCS m-4 c1.05" => "r",
	"T5750 g45, MARCS m-4 c1.25" => "tomato",
	"T5750 g45, MARCS m-4 c1.50" => "royalblue",
	"T5750 g45, MARCS m-4 c1.75" => "r",
	"T5750 g45, MARCS m-4 c2.75" => "r",
	"T5750 g45, CEMP m-4 c0.75"  => "k",
	"T5750 g45, CEMP"  => "k",
	"T5750 g45, CEMP C2.1"  => "tomato",
	"T5750 g45, CEMP C2.5"  => "royalblue",
	"T5750 g45, CEMP C3.2"  => "cyan",
	"T5750 g45, CEMP C3.1"  => "tomato",
	"T5750 g45, c0.4"  => "magenta",
	"T5750 g45"  => "royalblue",
	"T5750 g45, MARCS"  => "red",
	"T5750 g45, MARCS m-7 c5.0" => "r",
	"T5750 g45, CEMP m-7 c5.0" => "k",
	"T5750 g45, CEMP m-7 c4.25" => "magenta",
	"T5750 g45, CEMP m-7 c4.5" => "cyan",
	"T5750 g45, CEMP m-7 c4.75" => "lime",
	"metal-poor sun MARCS 5500" => "r",
	"metal-poor sun MARCS 5750" => "r",
	"metal-poor sun MARCS 6000" => "r",
	"metal-poor sun, CEMP MARCS 5600" => "r",
	"metal-poor sun, CEMP MARCS 5700" => "r",
	"metal-poor sun, CEMP MARCS 5772" => "r",
	"metal-poor sun, CEMP MARCS 5800" => "r",
	"metal-poor sun, CEMP MARCS 5772 C3.5" => "0.1",
	"metal-poor sun, CEMP MARCS 5772 C4.1" => "0.2",
	"metal-poor sun, CEMP MARCS 5772 C4.3" => "0.3",
	"metal-poor sun, CEMP MARCS 5772 C4.5" => "0.4",
	"metal-poor sun, CEMP MARCS 5772 C4.7" => "0.5",
)

# ╔═╡ 91bb18f2-3f87-4c10-a722-02232bc71ea5
ls = Dict(
	"metal-poor sun, CEMP" => "-",
	"metal-poor sun, CEMP C0.0" => "--",
	"metal-poor sun, CEMP C2.0" => "--",
	"metal-poor sun, CEMP C2.5" => "--",
	"metal-poor sun, CEMP C2.25" => "--",
	"metal-poor sun, av CEMP C3.0" => "--",
	"metal-poor sun, av CEMP C3.1" => "-",
	"metal-poor sun, av CEMP C2.9" => "-",
	"metal-poor sun" => "--",
	"metal-poor sun C3.1" => "-",
	"metal-poor sun C2.9" => "-",
	"metal-poor sun C2.8" => "-",
	"metal-poor sun C2.7" => "-",
	"T5750 g45, MARCS m-4 c0.75" => "-",
	"T5750 g45, MARCS m-4 c0.85" => "--",
	"T5750 g45, MARCS m-4 c0.95" => "--",
	"T5750 g45, MARCS m-4 c1.05" => "--",
	"T5750 g45, MARCS m-4 c1.25" => "--",
	"T5750 g45, MARCS m-4 c1.50" => "--",
	"T5750 g45, MARCS m-4 c1.75" => "--",
	"T5750 g45, MARCS m-4 c2.75" => "--",
	"T5750 g45, CEMP m-4 c0.75"  => "-",
	"T5750 g45, CEMP"  => "-",
	"T5750 g45, CEMP C2.1"  => "--",
	"T5750 g45, CEMP C2.5"  => "--",
	"T5750 g45, CEMP C3.2"  => "-",
	"T5750 g45, CEMP C3.1"  => "-",
	"T5750 g45, c0.4"  => "-",
	"T5750 g45, MARCS"  => "-",
	"T5750 g45"  => "--",
	"T5750 g45, MARCS m-7 c5.0" => "-",
	"T5750 g45, CEMP m-7 c5.0" => "-",
	"T5750 g45, CEMP m-7 c4.25" => "--",
	"T5750 g45, CEMP m-7 c4.5" => "--",
	"T5750 g45, CEMP m-7 c4.75" => "--",
	"metal-poor sun MARCS 5500" => "-",
	"metal-poor sun MARCS 5750" => "-",
	"metal-poor sun MARCS 6000" => "-",
	"metal-poor sun, CEMP MARCS 5600" => "-",
	"metal-poor sun, CEMP MARCS 5700" => "-",
	"metal-poor sun, CEMP MARCS 5772" => "-",
	"metal-poor sun, CEMP MARCS 5800" => "-",
	"metal-poor sun, CEMP MARCS 5772 C4.1" => "--",
	"metal-poor sun, CEMP MARCS 5772 C4.3" => "--",
	"metal-poor sun, CEMP MARCS 5772 C4.5" => "--",
	"metal-poor sun, CEMP MARCS 5772 C4.7" => "--",
	"metal-poor sun, CEMP MARCS 5772 C3.5" => "--"
)

# ╔═╡ 7f14c7c4-b39e-4fef-8e30-0381ac27cb8c
lw = Dict(
	"metal-poor sun, CEMP" => 2.5,
	"metal-poor sun, CEMP C0.0" => 1.5,
	"metal-poor sun, CEMP C2.0" => 1.5,
	"metal-poor sun, CEMP C2.5" => 1.0,
	"metal-poor sun, CEMP C2.25" => 1.0,
	"metal-poor sun, av CEMP C3.0" => 1.5,
	"metal-poor sun, av CEMP C3.1" => 1.5,
	"metal-poor sun, av CEMP C2.9" => 1.5,
	"metal-poor sun" => 1.5,
	"metal-poor sun C3.1" => 1.5,
	"metal-poor sun C2.9" => 1.5,
	"metal-poor sun C2.8" => 1.5,
	"metal-poor sun C2.7" => 1.5,
	"T5750 g45, MARCS m-4 c0.75" => 1.5,
	"T5750 g45, MARCS m-4 c0.85" => 1.5,
	"T5750 g45, MARCS m-4 c0.95" => 1.5,
	"T5750 g45, MARCS m-4 c1.05" => 1.5,
	"T5750 g45, MARCS m-4 c1.25" => 1.5,
	"T5750 g45, MARCS m-4 c1.50" => 1.5,
	"T5750 g45, MARCS m-4 c1.75" => 1.5,
	"T5750 g45, MARCS m-4 c2.75" => 1.5,
	"T5750 g45, CEMP m-4 c0.75"  => 1.5,
	"T5750 g45, CEMP"  => 1.5,
	"T5750 g45, CEMP C2.1"  => 1.5,
	"T5750 g45, CEMP C2.5"  => 1.5,
	"T5750 g45, CEMP C3.2"  => 1.5,
	"T5750 g45, CEMP C3.1"  => 1.5,
	"T5750 g45, c0.4"  => 1.5,
	"T5750 g45, MARCS"  => 1.5,
	"T5750 g45"  => 1.5,
	"T5750 g45, MARCS m-7 c5.0" => 1.5,
	"T5750 g45, CEMP m-7 c5.0" => 1.5,
	"T5750 g45, CEMP m-7 c4.25" => 1.5,
	"T5750 g45, CEMP m-7 c4.5" => 1.5,
	"T5750 g45, CEMP m-7 c4.75" => 1.5,
	"metal-poor sun MARCS 5500" => 1.0,
	"metal-poor sun MARCS 5750" => 1.5,
	"metal-poor sun MARCS 6000" => 1.0,
	"metal-poor sun, CEMP MARCS 5600" => 1.5,
	"metal-poor sun, CEMP MARCS 5700" => 2.5,
	"metal-poor sun, CEMP MARCS 5772" => 1.5,
	"metal-poor sun, CEMP MARCS 5800" => 1.5,
	"metal-poor sun, CEMP MARCS 5772 C4.1" => 1.5,
	"metal-poor sun, CEMP MARCS 5772 C4.3" => 1.5,
	"metal-poor sun, CEMP MARCS 5772 C4.5" => 1.5,
	"metal-poor sun, CEMP MARCS 5772 C4.7" => 1.5,
	"metal-poor sun, CEMP MARCS 5772 C3.5" => 1.5
)

# ╔═╡ 6ca44b98-6cbd-44ab-9ca6-4356591b4fac


# ╔═╡ 2e6e1727-29c6-4608-852b-527a76b4cec6
md"Plotting the following data (with label, comment out to remove from plot)"

# ╔═╡ 448e6a35-04da-4d51-bbb7-a60e5c3c00f1
labels = Dict(
	#"metal-poor sun, CEMP" => L"\rm [Fe/H]=-5.0,\ [C/Fe] = +3.0",
	#"metal-poor sun" => L"\rm [Fe/H]=-5.0,\ [C/Fe] = +0.4",
	"T5750 g45, CEMP m-4 c0.75"  => L"\rm [Fe/H]=-4.0,\ [C/Fe] = +0.75",
	"T5750 g45, CEMP"  => L"\rm [Fe/H]=-5.0,\ [C/Fe] = +3.0",
	"T5750 g45"  => L"\rm [Fe/H]=-5.0,\ [C/Fe] = +0.4",
	"T5750 g45, CEMP m-7 c5.0" =>  L"\rm [Fe/H]=-7.0,\ [C/Fe] = +5.0"
);

# ╔═╡ 437ab9e3-1eea-4a2a-95db-0ad9317e654f
labels_4297_4303 = Dict(
	"metal-poor sun, CEMP" => L"\rm CEMP\ DISPATCH\ (3D)",
	"metal-poor sun, CEMP C0.0" =>  L"\rm \Delta_{A(C)} = -3.0",
	"metal-poor sun, CEMP C2.0" =>  L"\rm \Delta_{A(C)} = -1.0",
	"metal-poor sun, CEMP C2.5" =>  L"\rm \Delta_{A(C)} = -0.5",
	"metal-poor sun, CEMP C2.25" =>  L"\rm \Delta_{A(C)} = -0.75",
	"metal-poor sun, av CEMP C3.0" => L"\rm <3D>",
	"metal-poor sun, av CEMP C3.1" => L"\rm <3D>, \Delta_{[C/Fe]} = +0.1",
	"metal-poor sun, av CEMP C2.9" => L"\rm <3D>, \Delta_{[C/Fe]} = -0.1",
	"metal-poor sun" => L"\rm scaled-solar\ atmosphere", 
	"metal-poor sun C3.1" => L"\rm [C/Fe]_{spectra} = 3.1",
	"metal-poor sun C2.9" => L"\rm \Delta_{A(C)} = +0.1",
	"metal-poor sun C2.8" => L"\rm \Delta_{A(C)} = +0.2",
	"metal-poor sun C2.7" => L"\rm \Delta_{A(C)} = +0.3",
	"T5750 g45, MARCS m-4 c0.75" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=0.75",
	"T5750 g45, MARCS m-4 c0.85" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=0.85",
	"T5750 g45, MARCS m-4 c0.95" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=0.95",
	"T5750 g45, MARCS m-4 c1.05" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=1.05",
	"T5750 g45, MARCS m-4 c1.25" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=1.25",
	"T5750 g45, MARCS m-4 c1.50" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=1.5",
	"T5750 g45, MARCS m-4 c1.75" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=1.75",
	"T5750 g45, MARCS m-4 c2.75" => L"\rm MARCS,\ [Fe/H]=-4,\ [C/Fe]=2.75",
	"T5750 g45, CEMP m-4 c0.75"  => L"\rm DISPATCH\ (3D),\ [Fe/H]=-4,\ [C/Fe]=0.75",
	"T5750 g45, CEMP"  => L"\rm CEMP\ DISPATCH\ (3D),\ [Fe/H]=-5,\ [C/Fe]=3.0",
	"T5750 g45, CEMP C2.1"  => L"\rm CEMP\ DISPATCH\ (3D),\ [Fe/H]=-5,\ [C/Fe]=2.1",
	"T5750 g45, CEMP C2.5"  => L"\rm CEMP\ DISPATCH\ (3D),\ [Fe/H]=-5,\ [C/Fe]=2.5",
	"T5750 g45, CEMP C3.2"  => L"\rm CEMP\ DISPATCH\ (3D),\ [Fe/H]=-5,\ [C/Fe]=3.2",
	"T5750 g45, CEMP C3.1"  => L"\rm CEMP\ DISPATCH\ (3D),\ [Fe/H]=-5,\ [C/Fe]=3.1",
	"T5750 g45, c0.4"  => L"\rm scaled-solar\ DISPATCH\ (3D), [Fe/H]=-5,\ [C/Fe]=3.05",
	"T5750 g45"  => L"\rm scaled-solar\ DISPATCH\ (3D)",
	"T5750 g45, MARCS"  => L"\rm MARCS,\ [Fe/H]=-5,\ [C/Fe]=3.0",
	"T5750 g45, MARCS m-7 c5.0" =>  L"\rm MARCS,\ [Fe/H]=-7,\ [C/Fe]=5",
	"T5750 g45, CEMP m-7 c5.0" =>  L"\rm DISPATCH\ (3D),\ [Fe/H]=-7,\ [C/Fe]=5",
	"T5750 g45, CEMP m-7 c4.25" => L"\rm DISPATCH\ (3D),\ [Fe/H]=-7,\ [C/Fe]=4.25",
	"T5750 g45, CEMP m-7 c4.5" => L"\rm DISPATCH\ (3D),\ [Fe/H]=-7,\ [C/Fe]=4.5",
	"T5750 g45, CEMP m-7 c4.75" => L"\rm DISPATCH\ (3D),\ [Fe/H]=-7,\ [C/Fe]=4.75",
	"metal-poor sun MARCS 5500" => L"\rm scaled-solar\ MARCS",
	"metal-poor sun MARCS 5750" => L"\rm MARCS,\ [Fe/H]=-5,\ [C/Fe]=3.0",
	"metal-poor sun MARCS 6000" => L"\rm scaled-solar\ MARCS",
	"metal-poor sun, CEMP MARCS 5600" => "CEMP MARCS",
	"metal-poor sun, CEMP MARCS 5700" => "CEMP MARCS",
	"metal-poor sun, CEMP MARCS 5772" => "CEMP MARCS",
	"metal-poor sun, CEMP MARCS 5800" => "CEMP MARCS",
	"metal-poor sun, CEMP MARCS 5772 C4.1" => L"\Delta_{A(C)} = -1.1",
	"metal-poor sun, CEMP MARCS 5772 C4.3" => L"\Delta_{A(C)} = -1.2",
	"metal-poor sun, CEMP MARCS 5772 C4.5" => L"\Delta_{A(C)} = -1.3",
	"metal-poor sun, CEMP MARCS 5772 C4.7" => L"\Delta_{A(C)} = -1.4",
	"metal-poor sun, CEMP MARCS 5772 C3.5" => L"\Delta_{A(C)} = -0.5"
);

# ╔═╡ 59bab38c-5a3f-46a1-b006-2c3a466cf10f
# solar like
#=plot_order = [
	"metal-poor sun, CEMP", "metal-poor sun"
]=#

# 3D effect
plot_order = [
	 "T5750 g45, CEMP", "T5750 g45"
]

# CEMP group
#plot_order = [
#	 "T5750 g45, CEMP", "T5750 g45, CEMP m-4 c0.75", "T5750 g45, CEMP m-7 c5.0"
#]

# ╔═╡ 24f1cca3-59d3-4d01-823b-75c63dea045b
# comparison between CEMP and non-CEMP models
#=plot_order_4297_4303 = [
	"T5750 g45, CEMP", "T5750 g45", "T5750 g45, CEMP C3.2", "T5750 g45, CEMP C3.1"
]=#

# comparuson 1D MP MARCS - 3D CEMP
plot_order_4297_4303 = [
	#"T5750 g45, MARCS", 
	#"T5750 g45",
	#"T5750 g45, CEMP", 
	#"T5750 g45, CEMP C2.1", 
	#"T5750 g45, CEMP C2.5", 
	#
	#"T5750 g45, MARCS m-4 c0.75",
	#"T5750 g45, MARCS m-4 c0.85",
	#"T5750 g45, MARCS m-4 c0.95",
	#"T5750 g45, MARCS m-4 c1.05",
	#"T5750 g45, MARCS m-4 c1.25",
	#"T5750 g45, MARCS m-4 c1.50",
	#"T5750 g45, MARCS m-4 c1.75",
	#"T5750 g45, MARCS m-4 c2.75",
	#
	#"T5750 g45, CEMP m-4 c0.75", 
	#
	"T5750 g45, MARCS m-7 c5.0", 
	"T5750 g45, CEMP m-7 c5.0", 
	#"T5750 g45, CEMP m-7 c4.25", 
	#"T5750 g45, CEMP m-7 c4.5", 
	#"T5750 g45, CEMP m-7 c4.75", 
	#
	#"metal-poor sun MARCS 5500",
	#"metal-poor sun MARCS 6000",
	#"metal-poor sun, CEMP MARCS 5700",
	#"metal-poor sun, CEMP MARCS 5772 C4.1",
	#"metal-poor sun, CEMP MARCS 5772 C4.3",
	#"metal-poor sun, CEMP MARCS 5772 C4.5",
	#"metal-poor sun, CEMP MARCS 5772 C4.7",
	#"metal-poor sun, CEMP MARCS 5772 C3.5",	
	#"metal-poor sun, CEMP MARCS 5772",
	#"metal-poor sun, CEMP C2.5",
	#"metal-poor sun, CEMP C2.25",
	#"metal-poor sun, CEMP C2.0",
	#"metal-poor sun, av CEMP C3.0",
	#"metal-poor sun", 
	#"metal-poor sun, CEMP"
]

# comparison between CEMP 3D and <3D>
#=plot_order_4297_4303 = [
	"metal-poor sun, av CEMP C3.0",
	"metal-poor sun, av CEMP C3.1",
	"metal-poor sun, av CEMP C2.9",
	"metal-poor sun, CEMP"
]=#

# ╔═╡ fce9aeb5-8c56-4774-b686-1d93b5857276
plot_order_molecules = [
	"T5750 g45",
	"T5750 g45, CEMP",
	#"T5750 g45, MARCS"
]

# ╔═╡ 5bca5b7c-5261-458d-92de-958f907d9531


# ╔═╡ 218bdee1-929f-413b-bdab-bcc0886c4208


# ╔═╡ 150d3049-ab7f-4d59-8e3f-3931fc0bff90
md"# Metal fraction"

# ╔═╡ ee9da77b-cfdb-45c5-86de-8475b2c26bb4
abund = MUST.readdlm(@in_m3dis("input_multi3d/abund_magg_a0.4_c3.0"))

# ╔═╡ 5291439a-7d01-4361-912d-1fe9f48eccfd
md"The abundances are scaled to the metallicity."

# ╔═╡ e058beb7-be6e-4ce8-afef-93e2af62d55a
begin
	FeH = -5.0
	abundances = Dict(
		e => e in ["H", "He", "Li"] ? a : a + FeH 
		for (e, a) in zip(abund[:,1], abund[:,2])
	)
end;

# ╔═╡ c4982179-9b76-424e-8428-b0fbbdcf9196


# ╔═╡ f75b8038-5a97-4173-b11f-4d0a339f9630
md"""
Abundance ratios are in units of H, such that

$\rm A(X) = log(N_X/N_H) + 12 \rightarrow N_X/N_H = 10^{A(X) - 12}$

That means that the metal fraction Z, with Hydorgen fraction X and Helium fraction Y such that 

$\rm X+Y+Z\equiv1$

can be conputed with

$\rm \frac{N_{He}}{N_H} + \sum_i \frac{Z_i}{N_H} = \frac{N_{He} + \sum_i Z_i}{N_H} = \frac{Y+Z}{X}$

This means that 

$\rm Y+Z = X\cdot \left(\frac{N_{He}}{N_H} + \sum_i \frac{Z_i}{N_H} \right)$

and hence

$\rm X+Y+Z = X\cdot \left(1 + \frac{N_{He}}{N_H} + \sum_i \frac{Z_i}{N_H} \right) \equiv 1$

which means that the Hydrogen fraction can be computes as

$\rm X = \frac{1}{1 + \frac{N_{He}}{N_H} + \sum_i \frac{Z_i}{N_H}}$
"""

# ╔═╡ c4763ca6-ee61-4cc0-a201-d6efb73a9c15


# ╔═╡ b6aa69ca-36ad-4ff3-93fa-bce19c6d0dda
abundance(a) = exp10(a - 12)

# ╔═╡ cd57629b-9337-49b3-ba3c-74b2eccca8be
metalfractions(a) = begin
	YX = abundance(a["He"])
	ZX = sum(abundance(a[e]) for e in keys(a) if !(e in ["H", "He"]))
	X = 1.0 / (1.0 + YX + ZX)
	Y = YX * X
	Z = ZX * X

	X, Y, Z
end

# ╔═╡ 2a885218-86b8-4286-9346-1ba4eefc5471
begin
	X, Y, Z = metalfractions(abundances)

	@info "This results in a composition of:" X Y Z
end

# ╔═╡ d3059747-d6ab-4690-9751-0e50bc6869ce


# ╔═╡ e80482ee-abb4-404d-95e8-2a9ff7d00c76


# ╔═╡ 8ee15ed7-4582-4967-ba7f-4f67f121a606
md"# Average profiles"

# ╔═╡ cf8bfd7b-c571-444f-b57a-8bc35a549a6e
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	for name in plot_order
		boxes = cubes[name]
		!(name in keys(labels)) && continue
		
		n = replace(name, ' '=>'_')
		n = n * "_T"
		x, y = profile(mean, boxes[2], :log10τ_ross, :T)
		MUST.writedlm(n, [x y])
		
		ax.plot(x, y, color=colors[name], ls=ls[name], lw=lw[name], label=labels[name])
	end

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	#ax.set_xlabel(L"\rm z\ [km]")
	ax.set_ylabel(L"\rm\ \log\ temperature\ [K]")

	ax.set_xlim(-4, 4)
	ax.set_ylim(3500, 12000)
	ax.legend()
	f.savefig("temperature_3Deffect_z.pdf")

	
	f
end

# ╔═╡ 1bb13217-808d-43f8-9705-0792d968f86d
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 5))

	for name in plot_order
		boxes=cubes[name]
		!(name in keys(labels)) && continue
		
		n = replace(name, ' '=>'_')
		n = n * "_rho"
		x, y = profile(mean, boxes[2], :log10τ_ross, :log10d)
		MUST.writedlm(n, [x y])
		
		ax.plot(x, y, color=colors[name], ls=ls[name], lw=lw[name], label=labels[name])
	end

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	ax.set_ylabel(L"\rm\ \log\ density\ [g\ cm^{-3}]")

	ax.set_xlim(-4, 4)
	ax.set_ylim(-6.9, -5.6)
	ax.legend()

	#f.savefig("density_3Deffect.pdf")

	f
end

# ╔═╡ da40ed04-2c17-4057-9803-3a04175a80ec


# ╔═╡ b99fc273-2918-4f4b-bacd-3b897291638e
md"# Spectra"

# ╔═╡ 202a1cb6-b7b3-446c-8c44-da438b645180
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(10, 5))
	
	for k in plot_order_4297_4303
		!haskey(labels_4297_4303, k) && continue
		s = spectra_4297_4303[k]
		λ, F = flux(s, norm=true)

		w = TSO.ω_midpoint(λ) 
		EW = sum(w .* (1 .- F))
		@info k*" [Å]" EW

		name = replace(k, ' '=>'_')
		name = name * "_4297_4303.txt"
		MUST.writedlm(name, [λ F])

		ax.plot(λ, F, color=colors[k], ls=ls[k], lw=lw[k], label=labels_4297_4303[k])
	end

	
	#ax.set_ylim(-0.05, 1.1)
	ax.set_xlim(4297, 4303)
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_ylabel(L"\rm flux\ [normalized]")
	ax.set_xlabel(L"\rm wavelength\ [\AA]")
	ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.98))
	#f.savefig("spectra_z-5_c3_carbon.pdf")
	#f.savefig("spectra_z-7_c5_1d3d.pdf")
	#f.savefig("spectra_z-5_c3_1d3d_2.pdf")
	#f.savefig("spectra_z-5_c3_carbon_corrected.pdf")

	f
end

# ╔═╡ ffc98088-a033-4020-97c0-0ffd53efbdbb
md"""
The corrections is (qualitatively):

$\rm [Fe/H]=-5, [C/Fe]=3.0: \Delta_{3D}\ = - 0.75\ dex$
$\rm [Fe/H]=-4, [C/Fe]=0.75: \Delta_{3D}\ = - 0.75\ dex$
"""

# ╔═╡ 92d4f31f-4d11-4605-85f0-843abe5c0bb9


# ╔═╡ 67d13344-01e3-4f82-a198-1d3ca43ca6bb
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(10, 5))

	ax.set_ylim(3000, 12500)
	ax.set_xlim(-4, 4)
	
	for k in plot_order_4297_4303
		!haskey(labels_4297_4303, k) && continue
		s = spectra_4297_4303[k]

		x = MUST.pyconvert(Array, s.run.ltau)
		y = MUST.pyconvert(Array, s.run.temp)
		if ndims(x)==1
			ax.plot(x, y, color=colors[k], ls=ls[k], lw=lw[k], label=labels_4297_4303[k])
		else
			!haskey(labels, k) && continue
			H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=250)

			xe = MUST.pyconvert(Array, xedges)
			ye = MUST.pyconvert(Array, yedges)
			H = MUST.pyconvert(Array, H.T)
			
			ax.imshow(
				H, interpolation="bicubic", origin="lower", aspect="auto",
				extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
			)
			
			ax.plot([], [], color="k", label=labels_4297_4303[k], marker="s", ls="")
			#xx, yy = MUST.numpy.meshgrid(x, y)
			#ax.pcolormesh(xx, yy, h.T, norm=matplotlib.colors.LogNorm())
		end
	end
	
	#ax.axhline(1.0, color="k", ls=":", alpha=0.5, lw=0.8)
	ax.set_xlabel(L"\rm \tau_{500}")
	ax.set_ylabel(L"\rm temperature\ [K]")
	ax.legend(loc="upper left", ncol=2)
	#f.savefig("temperature_MARCS_z-5_c3.0_2.pdf")
	
	f
end

# ╔═╡ a977f2d4-1089-4bea-b809-0caa6b750793


# ╔═╡ 73b20893-32b6-49f7-9eb3-d71482b57943
md"# Molecule formations"

# ╔═╡ 75564ab3-e9ed-4600-9c7b-5eb62ce2c405
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(10, 8))
	ax.set_xlim(-4, 4)
	
	plot_abund_3D(s, name; ax) = begin
		totn = MUST.pyconvert(Array, s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0])
		toth = MUST.pyconvert(Array, s.run.get_toth())
		C1 = log10.(totn ./ toth) .+ 12

		#= Plot =#
		x = MUST.pyconvert(Array, s.run.ltau)
		y = C1
		
		H, xedges, yedges = MUST.numpy.histogram2d(reshape(x, :), reshape(y, :), bins=250)

		xe = MUST.pyconvert(Array, xedges)
		ye = MUST.pyconvert(Array, yedges)
		H = MUST.pyconvert(Array, H.T)
		
		ax.imshow(
			H, interpolation="bicubic", origin="lower", aspect="auto",
			extent=[xe[1], xe[end], ye[1], ye[end]], cmap="Greys", norm=matplotlib.colors.LogNorm()
		)
	end

	plot_abund_av_3D(s, name; ax, kwargs...) = begin
		try
			totn = MUST.pyconvert(Array, s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0])
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
			totn = MUST.pyconvert(Array, s.run.read_patch_save(name, lazy=false, concat=true, fdim=0, zdim=2)[0][0][0])
			toth = MUST.pyconvert(Array, s.run.get_toth())
			C1 = log10.(totn ./ toth) .+ 12
	
			#= Plot =#
			x = MUST.pyconvert(Array, s.run.ltau)
			y = C1

			ax.plot(x, y; kwargs...)
	
			x, y
		end
	end
	
	for k in plot_order_molecules
		!haskey(labels_4297_4303, k) && continue
		s = spectra_4297_4303[k]
		name = replace(k, ' '=>'_')
	
		x, y = plot_abund_av_3D(s, "C_I", ax=ax, label=labels_4297_4303[k]*" C I", color="k", ls=ls[k], lw=lw[k])
		#x, y = MUST.pyconvert(Array, s.run.zz), MUST.pyconvert(Array, s.run.temp)
		#x = [MUST.mean(x[:, :, z]) for z in axes(x, 3)]
		#y = [MUST.mean(y[:, :, z]) for z in axes(y, 3)]
		#x, y = s.run.get_ltau_mean(s.run.temp)
		#x = MUST.pyconvert(Array, x)
		#y = MUST.pyconvert(Array, y)
		#ax.plot(x, log10.(y),  label=labels_4297_4303[k])
		MUST.writedlm(name*"_C_I.txt", [x y])
		
		x, y = plot_abund_av_3D(s, "C_II", ax=ax, label=labels_4297_4303[k]*" C II", color="cyan", ls=ls[k], lw=lw[k])
		MUST.writedlm(name*"_C_II.txt", [x y])
		
		x, y = plot_abund_av_3D(s, "CH", ax=ax, label=labels_4297_4303[k]*" CH", color="magenta", ls=ls[k], lw=lw[k])
		MUST.writedlm(name*"_CH.txt", [x y])
	end
	
	ax.set_xlabel(L"\rm \log{\tau_{500}}")
	ax.set_ylabel(L"\rm [X/H]")
	ax.legend(loc="lower center", ncol=1, framealpha=0)

	#f.savefig("molecules_CEMP-nonCEMP-MARCS.pdf")
	f
end

# ╔═╡ a260989a-93c3-4d39-8a8c-c5449320d9b5


# ╔═╡ 8b9d2ad3-3456-4cba-be40-a52a46dd8fcf
md"# Yoon-Beers 2019"

# ╔═╡ b7a63c29-b438-4b28-ae20-8ada77b87850
yb19 = MUST.readdlm("yoon_beers19.txt", ';')

# ╔═╡ 2bb06471-dca7-4baf-a1cb-f08e7d18eb81
feh = yb19[4:end, 6]

# ╔═╡ 5f8ef6ad-94a6-4c96-9098-5b7fb53c8553
c = yb19[4:end, 8]

# ╔═╡ d40a7e1d-8856-415a-88d4-3314689bee78
begin
	p1 = (-4.0, 0.7-4+8.560, -0.55)
	p2 = (-5.0, 3.5-5+8.560, -0.7+0.15)
	p3 = (-7.0, 5.5-7+8.560, -0.4+0.15)
	
	scatter_int(x, y) = MUST.pyconvert(typeof(x),
		MUST.scipy_interpolate.griddata(
		([p1[1], p2[1], p3[1]], 
		[p1[2], p2[2], p3[2]]), 
		[p1[3], p2[3], p3[3]], 
		(x, y), method="linear")
	)
	corrections = scatter_int(feh, c)
	mask = isnan.(corrections)
	scatter_int2(x, y) = MUST.pyconvert(typeof(x),
		MUST.scipy_interpolate.griddata(
		([p1[1], p2[1], p3[1]], 
		[p1[2], p2[2], p3[2]]), 
		[p1[3], p2[3], p3[3]], 
		(x[mask], y[mask]), method="nearest")
	)
	corrections[mask] = scatter_int2(feh, c)
end

# ╔═╡ df61824d-8045-40c8-8273-6f9fae669595
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
	ax.scatter(feh, c+corrections, s=10, marker="o", c="red", label="corrected")

	cemp_line(x) = 0.75 + x + 8.39
	x = range(-9, 0, length=100)|>collect
	ax.plot(x, cemp_line.(x), color="k", alpha=0.5, ls="--")

	ax.set_xlim(-8.5, -1)
	ax.set_ylim(4, 10)
	ax.set_xlabel("[Fe/H]")
	ax.set_ylabel("A(C)")
	ax.legend(loc="upper left")
	#f.savefig("yoon-beers_19-corrected_slope.pdf")

	f
end

# ╔═╡ Cell order:
# ╠═cc6eeb7e-5974-11ef-0436-fd682756392a
# ╟─a73ba74d-30ca-4716-9a0c-455aa24cf975
# ╟─39703f8e-40b2-476a-8ccb-cf5d208866fa
# ╠═5e8a2581-4cb0-4337-a913-4c4df3121455
# ╠═d0815aec-c693-4b7f-ab3d-2ba6b8f7176c
# ╠═cea7a484-48f5-4c38-9966-00e0f4409b50
# ╠═f17e97d8-3239-46d5-9883-8124cb0170a3
# ╠═2d84b97b-b835-421b-8b8a-2800da47e130
# ╠═b6f2f6e8-67e6-4bc8-936d-6d0f4a6a5111
# ╠═a7ee3a6e-5855-4dc1-a683-7f32f1bccb18
# ╟─307e6d9a-d7ac-435c-b4eb-c1c9e7825110
# ╟─05895238-a947-429d-9ecb-751ebe199ca9
# ╟─b737400b-0e76-40ae-98d4-01c00e8f2dd3
# ╠═4c85536c-87d8-4d0b-87d1-b2fc32ca5118
# ╟─4ad62be3-8d9f-446e-bf61-7d77d0475d25
# ╟─f00cc719-1090-46fc-ba84-31bb1351c726
# ╠═44429920-277b-4a09-a8a4-a6e143056628
# ╟─37bcf775-4ccc-44a7-af80-b10ffa8a26c7
# ╟─1f6071d5-ab91-4eef-9dd2-dea9237f1ec8
# ╠═21e58231-22b5-48e6-a8aa-c4bbe30a8c4c
# ╠═91bb18f2-3f87-4c10-a722-02232bc71ea5
# ╠═7f14c7c4-b39e-4fef-8e30-0381ac27cb8c
# ╟─6ca44b98-6cbd-44ab-9ca6-4356591b4fac
# ╟─2e6e1727-29c6-4608-852b-527a76b4cec6
# ╠═448e6a35-04da-4d51-bbb7-a60e5c3c00f1
# ╠═437ab9e3-1eea-4a2a-95db-0ad9317e654f
# ╠═59bab38c-5a3f-46a1-b006-2c3a466cf10f
# ╠═24f1cca3-59d3-4d01-823b-75c63dea045b
# ╠═fce9aeb5-8c56-4774-b686-1d93b5857276
# ╟─5bca5b7c-5261-458d-92de-958f907d9531
# ╟─218bdee1-929f-413b-bdab-bcc0886c4208
# ╟─150d3049-ab7f-4d59-8e3f-3931fc0bff90
# ╠═ee9da77b-cfdb-45c5-86de-8475b2c26bb4
# ╟─5291439a-7d01-4361-912d-1fe9f48eccfd
# ╠═e058beb7-be6e-4ce8-afef-93e2af62d55a
# ╟─c4982179-9b76-424e-8428-b0fbbdcf9196
# ╟─f75b8038-5a97-4173-b11f-4d0a339f9630
# ╟─c4763ca6-ee61-4cc0-a201-d6efb73a9c15
# ╟─b6aa69ca-36ad-4ff3-93fa-bce19c6d0dda
# ╟─cd57629b-9337-49b3-ba3c-74b2eccca8be
# ╟─2a885218-86b8-4286-9346-1ba4eefc5471
# ╟─d3059747-d6ab-4690-9751-0e50bc6869ce
# ╟─e80482ee-abb4-404d-95e8-2a9ff7d00c76
# ╟─8ee15ed7-4582-4967-ba7f-4f67f121a606
# ╟─cf8bfd7b-c571-444f-b57a-8bc35a549a6e
# ╟─1bb13217-808d-43f8-9705-0792d968f86d
# ╟─da40ed04-2c17-4057-9803-3a04175a80ec
# ╟─b99fc273-2918-4f4b-bacd-3b897291638e
# ╟─202a1cb6-b7b3-446c-8c44-da438b645180
# ╟─ffc98088-a033-4020-97c0-0ffd53efbdbb
# ╟─92d4f31f-4d11-4605-85f0-843abe5c0bb9
# ╟─67d13344-01e3-4f82-a198-1d3ca43ca6bb
# ╟─a977f2d4-1089-4bea-b809-0caa6b750793
# ╟─73b20893-32b6-49f7-9eb3-d71482b57943
# ╟─75564ab3-e9ed-4600-9c7b-5eb62ce2c405
# ╟─a260989a-93c3-4d39-8a8c-c5449320d9b5
# ╟─8b9d2ad3-3456-4cba-be40-a52a46dd8fcf
# ╠═b7a63c29-b438-4b28-ae20-8ada77b87850
# ╠═2bb06471-dca7-4baf-a1cb-f08e7d18eb81
# ╠═5f8ef6ad-94a6-4c96-9098-5b7fb53c8553
# ╠═d40a7e1d-8856-415a-88d4-3314689bee78
# ╟─df61824d-8045-40c8-8273-6f9fae669595
