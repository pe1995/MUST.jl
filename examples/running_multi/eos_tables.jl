### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ 61b5121d-670c-4baf-b021-c1ed259d8b4c
md"# Creating opacity table input"

# ╔═╡ 77792400-dabb-4dbe-a0cd-36b80ebbe25b
md"Creating opacity tables with the new Multi design is more straigt forward. The input atmosphere only needs to contain the minimum and maximum temperature and density."

# ╔═╡ ab618049-b942-4ae3-b2a7-703e5f36df1a
modelatmosfolder = "input_multi3d/test_opac_table/"

# ╔═╡ fb801c36-1de2-49a8-847e-99f05bc1e2cc
function absdat_abundances(;α=0.0, default="./input_multi3d/abund_magg", eles...)
	# read the default abundances
	abund_default = MUST.readdlm(@in_m3dis(default))
	abund_new = deepcopy(abund_default)
	new_name = default
	if α != 0.0
		for ele in TSO.α_elements
			iele = findfirst(ele .== abund_new[:, 1])
			if !isnothing(iele)
				abund_new[iele, 2] = abund_default[iele, 2] + α
			else
				@warn "element $(ele) not found in absdat $(default)."
			end
		end
		new_name *= "_a$(α)"
	end

	for (eleS, val) in eles
		ele = string(eleS)
		iele = findfirst(ele .== abund_new[:, 1])
		if !isnothing(iele)
			abund_new[iele, 2] = abund_default[iele, 2] + val
			new_name *= "_$(ele)$(val)"
		else
			@warn "element $(ele) not found in absdat $(default)."
		end
	end

	if new_name != default
		open(@in_m3dis(new_name), "w") do f
			for i in axes(abund_new, 1)
				line = MUST.@sprintf "%-4s %-.3f\n" abund_new[i, 1] abund_new[i, 2]
				write(f, line)
			end
		end
	else
		default
	end

	new_name
end

# ╔═╡ e9e182bd-e296-43ce-a309-2c89df055cbc
begin
	# Dimensions of the EoS Table
	minT = 1000.
	maxT = 100000.
	minρ = 1e-18
	maxρ = 1e-2
	nT = 150
	nρ = 150
	
	# Wavelength coverage and resolution
	λs = 1000
	λe = 200000
	nλ = 100000
	vmic = 1.0

	# pick linelists you want to use
	linelists = String[
		"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
		"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list",
		"/home/eitner/shared/StAt/LINE-LISTS/25500-200000_cut-4/atom_25500-200000.list",
		"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata",
		"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/12CH_multi.list",
		"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/13CH_multi.list",
		"/home/eitner/shared/StAt/LINE-LISTS/combined_molecules/H2O_multi.list",
		MUST.glob("*.list", "/home/eitner/shared/StAt/LINE-LISTS/combined_molecules/most_relevant/")...
	]

	# modify chemical composition by element name (beyond [Fe/H] with abund file, give as [X/Fe])
	FeH = -2.0
	abund_file = absdat_abundances(
		α=0.0, 
		#C=3.0,
		default="./input_multi3d/abund_magg"
	)

	# name of the new EoS table
	extension = "magg_m2_a0_vmic1"
	version = "v5.1"

	# computation setup
	compute = true  # set to false if you only want to collect the output from M3D
	nν = 34
	threads = 34
	slurm = false
end;

# ╔═╡ b87b7824-bd8d-43f2-ae88-b454d293acaa
function eosTableInput(wheretosave; minT=1000., maxT=5.5e5, minρ=1e-30, maxρ=1e-3, vmic=1.0)
	path = MUST.@in_m3dis(wheretosave)

	if isdir(path)
		# clear the folder
		rm(path, recursive=true)
	end
	mkdir(path)

	z = [-1.0, 1.0]
	TSO.saveAsText(
		joinpath(path, "TSO-M3D"), 
		z=z, 
		T=[minT, maxT], 
		ρ=[minρ, maxρ],
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
				:FeH=>FeH,
				:nz=>2,
				:amr=>false
			),
			:m3d_params=>(
				:n_nu=>nν, 
				:ilambd=>0,
				:short_scheme=>"disk_center",
				:long_scheme=>"none",
				:make_eos=>true
			),
			:composition_params=>(
				:absdat_file=>"./input_multi3d/TS_absdat.dat",
                :abund_file=>abund_file,
				:ldtemp=>δlnT,
				:ldrho=>δlnρ,
				:tmolim=>100000.0,
				:mhd_eos=>true
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
		slurm=slurm,
		nν=nν,
		FeH=FeH,
		m3dis_kwargs=Dict(
			:threads=>threads,
			#:memMB=>90000
		)
	)
end

# ╔═╡ cf6b3809-a161-4093-8285-2cb303c380f4
md"# Read EoS Table"

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
	v4.X: FeH = -X, still no molecular lines!
- v4.0:	Carbon enhanced <-> Carbon solar
-------------------------------------------------------
	v5.X: with molecules, OP formalism for H in EoS 
-v5.0: 50,000 K
-v5.1: 100,000 K
-------------------------------------------------------
"

# ╔═╡ 934be5d3-a7c5-46f2-870d-8ba7d8c134dc
eos_folder = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(extension)_$(version)"

# ╔═╡ d84c4140-1702-4fa4-8fc5-955a1e9c0d78
!isdir(eos_folder) && mkdir(eos_folder)

# ╔═╡ 4611660d-f561-457a-80d4-fd37630e5c3b
begin
	run = MUST.M3DISRun("data/$(model)", read_atmos=false)
end

# ╔═╡ e44b4fab-f33f-4d93-afea-c80b4a1849e1
eos, opa, scat = TSO.collect_opacity(
	run, 
	compute_ross=true
)

# ╔═╡ 05c98599-76e6-4747-8689-1cc4cc249994
begin
	eos_mono = deepcopy(eos)
	TSO.smoothAccumulate!(eos_mono)
end

# ╔═╡ da9cbcec-5e3d-47e6-8b52-e9aad755de20
begin
	save(eos_mono, joinpath(eos_folder, "combined_eos_$(extension).hdf5"))
	save(opa, joinpath(eos_folder, "combined_opacities_$(extension).hdf5"))
	if !isnothing(scat)
		save(scat, joinpath(eos_folder, "combined_sopacities_$(extension).hdf5"))
	end
end

# ╔═╡ 8c4c378c-0d0a-4d01-bcff-0b8201fdd402
size(eos)

# ╔═╡ f06e07d1-7989-4e45-bbe2-e4b25f77c36b
aos = @axed eos_mono

# ╔═╡ 1f320b5e-30d5-4275-a99a-e2ac07eaa910
let
	tt, dd = TSO.meshgrid(aos)
	
	plt.close()
	
	im = plt.pcolormesh(log10.(exp.(tt)), log10.(exp.(dd)), log10.(exp.(eos_mono.lnEi)))

	c = plt.colorbar(im)
	c.set_label("log E")

	plt.xlabel("log T")
	plt.ylabel("log ρ")
	
	plt.savefig("internalenergy.png")

	gcf()
end

# ╔═╡ 0e818674-be73-4e6e-9de8-2ad7951c8ed0
let
	plt.close()

	plt.title("ρ=$(exp(eos_mono.lnRho[100]))")
	im = plt.plot(eos_mono.lnT, eos_mono.lnEi[:, 100])

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
	
	plt.savefig("opacity.png")

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
# ╟─fb801c36-1de2-49a8-847e-99f05bc1e2cc
# ╠═e9e182bd-e296-43ce-a309-2c89df055cbc
# ╟─b87b7824-bd8d-43f2-ae88-b454d293acaa
# ╟─930978f2-f898-420d-8da8-d395be0748e9
# ╟─6e98979c-c748-4272-b4a9-dd7c2dbaeb96
# ╠═18f2d770-4803-4ee7-860b-aa87db9245f7
# ╟─15a2ab2c-71e1-4778-8d22-759ad16a0bbe
# ╠═d88497d9-d66a-492f-90f1-2581688ce5ec
# ╟─cf6b3809-a161-4093-8285-2cb303c380f4
# ╟─e3f7e4bf-17c2-4c7d-b086-8ef846fdff32
# ╠═934be5d3-a7c5-46f2-870d-8ba7d8c134dc
# ╠═d84c4140-1702-4fa4-8fc5-955a1e9c0d78
# ╠═4611660d-f561-457a-80d4-fd37630e5c3b
# ╠═e44b4fab-f33f-4d93-afea-c80b4a1849e1
# ╠═05c98599-76e6-4747-8689-1cc4cc249994
# ╠═da9cbcec-5e3d-47e6-8b52-e9aad755de20
# ╠═8c4c378c-0d0a-4d01-bcff-0b8201fdd402
# ╠═f06e07d1-7989-4e45-bbe2-e4b25f77c36b
# ╠═1f320b5e-30d5-4275-a99a-e2ac07eaa910
# ╠═0e818674-be73-4e6e-9de8-2ad7951c8ed0
# ╠═d03b0d8b-aa2d-49dc-a8bb-6064564c11ca
