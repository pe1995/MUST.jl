### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ b0b8ef88-d23d-11ee-1e5f-4b51fc39a67b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ d2f37232-47da-44a6-ae85-3d64b4acee9d
plt = matplotlib.pyplot

# ╔═╡ caf3e0e8-3792-4fbb-8c29-fcd32512a03a
md"# M3D Setup"

# ╔═╡ 50a53c7e-07e1-4a09-bf02-7f3521e33b75
MUST.@import_m3dis "../../../Multi3D"

# ╔═╡ 82c1bb90-8647-4091-9d19-f522d76dcd05
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ 8d43e13b-5033-4e23-b9a7-24bf42f34aa9
linelists = String[
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list",
	"/home/eitner/shared/StAt/LINE-LISTS/25500-200000_cut-4/atom_25500-200000.list",
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata",
	#"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/12CH_multi.list"
]

# ╔═╡ 61b5121d-670c-4baf-b021-c1ed259d8b4c
md"# Creating opacity table input"

# ╔═╡ 77792400-dabb-4dbe-a0cd-36b80ebbe25b
md"Creating opacity tables with the new Multi design is more straigt forward. The input atmosphere only needs to contain the minimum and maximum temperature and density."

# ╔═╡ ab618049-b942-4ae3-b2a7-703e5f36df1a
modelatmosfolder = "input_multi3d/test_opac_table/"

# ╔═╡ e9e182bd-e296-43ce-a309-2c89df055cbc
begin
	minT = 1000.
	maxT = 50000.
	minρ = 1e-20
	maxρ = 1e-2
	nT   = 150
	nρ   = 150

	λs = 1000
	λe = 200000
	nλ = 100000

	vmic = 2.0
end

# ╔═╡ b87b7824-bd8d-43f2-ae88-b454d293acaa
function eosTableInput(wheretosave; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3, vmic=1.0)
	path = MUST.@in_m3dis(wheretosave)

	if isdir(path)
		# clear the folder
		rm(path, recursive=true)
	end
	mkdir(path)

	z = [-1.0, 0.0, 1.0]
	TSO.saveAsText(
		joinpath(path, "TSO-M3D"), 
		z=z, 
		T=[minT, (maxT-minT)/2, maxT], 
		ρ=[minρ, (maxρ-minρ)/2, maxρ],
		vmic=vmic
	)

	"TSO-M3D"
end

# ╔═╡ 930978f2-f898-420d-8da8-d395be0748e9
model = eosTableInput(
	modelatmosfolder; 
	minT=minT, maxT=maxT, 
	minρ=minρ, maxρ=maxρ,
	vmic=vmic
)

# ╔═╡ 6e98979c-c748-4272-b4a9-dd7c2dbaeb96
md"# Run M3D"

# ╔═╡ 18f2d770-4803-4ee7-860b-aa87db9245f7
compute = true

# ╔═╡ 15a2ab2c-71e1-4778-8d22-759ad16a0bbe
eosTable(model; folder, linelist, λs, λe, δλ, δlnT, δlnρ, FeH=0.0, nν=10,
				in_log=true, slurm=false, m3dis_kwargs=Dict(), kwargs...) = begin
    MUST.whole_spectrum(
		model, 
		namelist_kwargs=(
			:model_folder=>folder,
			:linelist=>nothing,
			:absmet=>nothing,
			:linelist_params=>(:line_lists=>linelist,),
			:atom_params=>(:atom_file=>"", ),
			:spectrum_params=>(:daa=>δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>in_log),
			:atmos_params=>(
				:dims=>1, 
				:atmos_format=>"Text",
				:use_density=>true, 
				:use_ne=>false,
				:FeH=>FeH
			),
			:m3d_params=>(
				:n_nu=>nν, 
				:ilambd=>0,
				:quad_scheme=>"disk_center",
				:long_scheme=>"none",
				:make_eos=>true
			),
			:composition_params=>(
				:absdat_file=>"./input_multi3d/TS_absdat.dat",
                :abund_file=>"./input_multi3d/abund_magg_c3_a4",
				:ldtemp=>δlnT,
				:ldrho=>δlnρ,
				:tmolim=>10000.0
			),
            kwargs...
		),
		m3dis_kwargs=m3dis_kwargs,
		slurm=slurm
	)
end

# ╔═╡ d88497d9-d66a-492f-90f1-2581688ce5ec
if compute 
	eosTable(
		model; 
		folder=modelatmosfolder, 
		linelist=linelists,
		λs=log(λs),
		λe=log(λe),
		δλ=(log(λe)-log(λs))/nλ,
		in_log=true,
		δlnT=(log(maxT)-log(minT))/nT, 
		δlnρ=(log(maxρ)-log(minρ))/nρ,
		slurm=false,
		nν=20,
		FeH=-4.0
		#=m3dis_kwargs=Dict(
			:threads=>32,
			:memMB=>90000
		)=#
	)
end

# ╔═╡ cf6b3809-a161-4093-8285-2cb303c380f4
md"# Read EoS Table"

# ╔═╡ fe9f9763-0b3e-4a47-9a2c-e5a24d2bfd67
md"Decide where to save the table"

# ╔═╡ 270b37a1-d0c9-48b1-bd05-891a8c83cf72
extension = "magg_m4_a4_c3_vmic2"

# ╔═╡ e3f7e4bf-17c2-4c7d-b086-8ef846fdff32
md"EoS Versions:
-------------------------------------------------------
	v1.0: First test, no H lines, no Molecules
- v1.1: + H lines (possibly issues with Lalpha)
- v1.2: same as v1.1, but after Richard fixed H lines
- v1.3: same as v1.0, but with smaller range
- v1.4: small range without molecules in EoS + Hlines
- v1.5: small range without molecules in EoS 
-------------------------------------------------------
	v2.X: w/o scattering in absorptio + sep. table, no molecular lines yet
- v2.0:
- v2.0.1: Parallel rosseland_opacity 
- v2.1: low density table
-------------------------------------------------------
	v3.X: FeH = -1, still no molecular lines!
- v3.0:	Carbon enhanced <-> Carbon solar
-------------------------------------------------------
	v4.X: FeH = -4, still no molecular lines!
- v4.0:	Carbon enhanced <-> Carbon solar
"

# ╔═╡ 934be5d3-a7c5-46f2-870d-8ba7d8c134dc
eos_folder = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(extension)_v4.0"

# ╔═╡ d84c4140-1702-4fa4-8fc5-955a1e9c0d78
!isdir(eos_folder) && mkdir(eos_folder)

# ╔═╡ 4611660d-f561-457a-80d4-fd37630e5c3b
run = MUST.M3DISRun("data/$(model)")

# ╔═╡ e44b4fab-f33f-4d93-afea-c80b4a1849e1
eos, opa, scat = TSO.collect_opacity(
	run, 
	compute_ross=true
)

# ╔═╡ da9cbcec-5e3d-47e6-8b52-e9aad755de20
begin
	save(eos, joinpath(eos_folder, "combined_eos_$(extension).hdf5"))
	save(opa, joinpath(eos_folder, "combined_opacities_$(extension).hdf5"))
	if !isnothing(scat)
		save(scat, joinpath(eos_folder, "combined_sopacities_$(extension).hdf5"))
	end
end

# ╔═╡ 8c4c378c-0d0a-4d01-bcff-0b8201fdd402
size(eos)

# ╔═╡ f06e07d1-7989-4e45-bbe2-e4b25f77c36b
aos = @axed eos

# ╔═╡ 1f320b5e-30d5-4275-a99a-e2ac07eaa910
let
	tt, dd = TSO.meshgrid(aos)
	
	plt.close()
	
	im = plt.scatter(tt, dd, s=1.0, c=eos.lnRoss)

	plt.colorbar(im)

	gcf()
end

# ╔═╡ 0e818674-be73-4e6e-9de8-2ad7951c8ed0
let
	plt.close()

	plt.title("ρ=$(exp(eos.lnRho[10]))")
	im = plt.plot(eos.lnT, eos.lnEi[:, 10])

	plt.ylim(28, 35)

	gcf()
end

# ╔═╡ d03b0d8b-aa2d-49dc-a8bb-6064564c11ca
let
	tt, dd = TSO.meshgrid(aos)
	
	plt.close()

	i = 130
	j = 130

	plt.title("T=$(exp.(eos.lnT[i])), ρ=$(exp.(eos.lnRho[j]))")
	
	im = plt.plot(opa.λ, opa.κ[i, j, :])
	#plt.xlim(1000, 10000)
	#plt.ylim(-1, 100)
	

	gcf()
end

# ╔═╡ Cell order:
# ╠═b0b8ef88-d23d-11ee-1e5f-4b51fc39a67b
# ╠═d2f37232-47da-44a6-ae85-3d64b4acee9d
# ╟─caf3e0e8-3792-4fbb-8c29-fcd32512a03a
# ╠═50a53c7e-07e1-4a09-bf02-7f3521e33b75
# ╠═82c1bb90-8647-4091-9d19-f522d76dcd05
# ╠═8d43e13b-5033-4e23-b9a7-24bf42f34aa9
# ╟─61b5121d-670c-4baf-b021-c1ed259d8b4c
# ╟─77792400-dabb-4dbe-a0cd-36b80ebbe25b
# ╠═ab618049-b942-4ae3-b2a7-703e5f36df1a
# ╠═e9e182bd-e296-43ce-a309-2c89df055cbc
# ╠═b87b7824-bd8d-43f2-ae88-b454d293acaa
# ╠═930978f2-f898-420d-8da8-d395be0748e9
# ╟─6e98979c-c748-4272-b4a9-dd7c2dbaeb96
# ╠═18f2d770-4803-4ee7-860b-aa87db9245f7
# ╟─15a2ab2c-71e1-4778-8d22-759ad16a0bbe
# ╠═d88497d9-d66a-492f-90f1-2581688ce5ec
# ╟─cf6b3809-a161-4093-8285-2cb303c380f4
# ╟─fe9f9763-0b3e-4a47-9a2c-e5a24d2bfd67
# ╠═270b37a1-d0c9-48b1-bd05-891a8c83cf72
# ╟─e3f7e4bf-17c2-4c7d-b086-8ef846fdff32
# ╠═934be5d3-a7c5-46f2-870d-8ba7d8c134dc
# ╠═d84c4140-1702-4fa4-8fc5-955a1e9c0d78
# ╠═4611660d-f561-457a-80d4-fd37630e5c3b
# ╠═e44b4fab-f33f-4d93-afea-c80b4a1849e1
# ╠═da9cbcec-5e3d-47e6-8b52-e9aad755de20
# ╠═8c4c378c-0d0a-4d01-bcff-0b8201fdd402
# ╠═f06e07d1-7989-4e45-bbe2-e4b25f77c36b
# ╠═1f320b5e-30d5-4275-a99a-e2ac07eaa910
# ╠═0e818674-be73-4e6e-9de8-2ad7951c8ed0
# ╠═d03b0d8b-aa2d-49dc-a8bb-6064564c11ca
