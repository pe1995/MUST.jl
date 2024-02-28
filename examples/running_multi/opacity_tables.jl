### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ dc194f68-641c-11ee-3bdf-758771368780
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using MUST
	using TSO
	using PythonPlot
	using PythonCall 
end

# ╔═╡ cc74395c-b2ab-4575-af5a-b00071957ad9
np = pyimport("numpy")

# ╔═╡ 80c4f1c0-d4e2-445b-af7e-41353c8d5151
plt = matplotlib.pyplot

# ╔═╡ 399ad336-8066-404e-a6df-cf515f55236e
md"# M3D Setup"

# ╔═╡ e6c0702b-240a-4da8-a791-e360471bb990
# ╠═╡ show_logs = false
MUST.@import_m3dis "../../../Multi3D"

# ╔═╡ 473f9d7c-f0b6-4d44-bef2-207dfb42f2d4
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ c4a8a402-f700-48c4-b13f-61f6a84b4dd5
linelists = String[
	#"./input_multi3d/nlte_ges_linelist_jmg25jan2023_I_II",
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list",
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata"
]

# ╔═╡ e5f3464e-fe56-45c0-8e47-ca4f1321b1b5
md"# Creating opacity tables"

# ╔═╡ bc0d7f88-f34b-4323-9e1a-4809ceac9efd
md"We split up the desired density-temperature grid in 1D atmospheres consisting of one density and all temperatures. Those are saved so that M3D can read them."

# ╔═╡ 3e09464b-a09f-4cc8-8533-d2d8cdf40509
modelatmosfolder = "input_multi3d/test_opac_table/"

# ╔═╡ acd24cd4-3676-4a54-af7b-690564f97425
models = TSO.opacityTableInput(
	MUST.@in_m3dis(modelatmosfolder),
	lnT = range(log(1.1e3), log(5.5e5); length=50) |> collect, 
    lnρ = range(log(1e-30), log(1e-3); length=50) |> collect
)

# ╔═╡ d32c241f-97c8-46db-9512-1c7111e03e99
md"# Running M3D"

# ╔═╡ 9bce1738-531c-459e-95ac-ba5526df59be
md"Execute M3D in the given wavelength range"

# ╔═╡ 4c9becb2-07d4-4126-ab76-a741e97abb85
compute = false 

# ╔═╡ 078700b7-0b47-43a1-85d2-4567392e3479
opacityTable(models; folder, linelist, λs, λe, δλ, 
				in_log=true, slurm=false, m3dis_kwargs=Dict(), kwargs...) = begin
    MUST.whole_spectrum(
		models, 
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
				:use_ne=>false
			),
			:m3d_params=>(
				:n_nu=>1, 
				:ilambd=>0,
				:quad_scheme=>"disk_center",
				:long_scheme=>"disk_center",
				:make_opac_table=>true
			),
			:composition_params=>(
				:absdat_file=>"./input_multi3d/TS_absdat.dat",
                :abund_file=>"./input_multi3d/abund_magg"
			),
            kwargs...
		),
		m3dis_kwargs=m3dis_kwargs,
		slurm=slurm
	)
end

# ╔═╡ 9b48e5c5-ef7a-4bd2-9b8c-213d706d44c8
if compute 
	opacityTable(
		models; 
		folder=modelatmosfolder, 
		linelist=linelists,
		λs=log(1000), λe=log(100000), δλ=0.00003, in_log=true,
		slurm=true,
		m3dis_kwargs=Dict(
			:threads=>64,
			:memMB=>90000
		)
	)
end

# ╔═╡ 9e2be75a-7d0b-4c65-b750-845798e9223d
pyconvert

# ╔═╡ 711972c6-dfb3-4143-b619-0cba4e96a2e4
md"# Collecting the output"

# ╔═╡ bcff1c95-3d8e-40a1-87e0-3c7a36cd4206
m3dis_models = [MUST.M3DISRun("data/$(m)") for m in models]

# ╔═╡ 170c9ca0-bd1c-41d8-908e-4b02e58a081d
#=begin
	r1 = m3dis_models[1]

	f = joinpath(join(r1.sfolder), "chi_meta.txt")
	open(f, "r") do f
		for line in readlines(f)
			if !(first(strip(line))=='*')
				info = split(line)
				file = first(info)
				prec = parse(Int, info[2])
                dtype = TSO.@sprintf "<f%d" prec//8
                dims = Tuple([parse(Int, d) for d in info[3:end]])

				key = split(split(file, ".")[1], "_")[end]
				key = [parse(Int, k) for k in [key[1:4], key[4:7], key[7:9]]]

				 if occursin("_out", file)
                    key[2:end] = key[2:end] .*-1
				 end

				m = np.memmap(
					joinpath(join(r1.sfolder), file), dtype=dtype, mode="r", shape=dims, order="F"
				)

				@show m
			end
		end
	end
end=#

# ╔═╡ 2c3eef18-c8fb-42a0-ba56-607cc99d1b36
begin
	eos, opa = TSO.collect_opacity(m3dis_models)
	dirname = "m3dis_a0_m0_magg22"
	!isdir(dirname) && mkdir(dirname)
	save(eos, joinpath(dirname, "combined_eos.hdf5"))
	save(opa, joinpath(dirname, "combined_opacities.hdf5"))
end

# ╔═╡ 056c6791-ed98-4c49-a424-e1ac72e6b22a
aos = @axed eos

# ╔═╡ 2c9bb5de-0fc6-4f27-a5ce-0fca1f871baa
begin
	plt.close()

	tt, rr = meshgrid(@axed eos)
	plt.scatter(
		tt, rr, c=log10.(opa.κ[:, :, 10]), s=20
	)

	plt.colorbar()
	plt.savefig(joinpath(dirname,"opacity.png"))
	
	gcf()
end

# ╔═╡ 10642e82-6dbf-46a6-9914-fe25148365e6
begin
	plt.close()

	plt.scatter(
		tt, rr, c=eos.lnEi, s=20
	)

	plt.colorbar()
	plt.savefig(joinpath(dirname,"lnEi.png"))
	
	gcf()
end

# ╔═╡ 1d879108-07c7-4027-b300-e7855d1e5ea9
begin
	plt.close()

	plt.scatter(
		log10.(opa.λ), log10.(opa.κ[10, 10, :]), s=2
	)

	@show length(opa.λ)

	plt.savefig(joinpath(dirname,"lambda_opacity.png"))
	
	gcf()
end

# ╔═╡ d3464b51-4652-49ba-87e4-545cec44dafe
begin
	plt.close()

	plt.scatter(
		log10.(opa.λ), log10.(opa.src[10, 10, :]), s=2
	)

	plt.savefig(joinpath(dirname,"lambda_sourcefunction.png"))
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═dc194f68-641c-11ee-3bdf-758771368780
# ╠═cc74395c-b2ab-4575-af5a-b00071957ad9
# ╠═80c4f1c0-d4e2-445b-af7e-41353c8d5151
# ╟─399ad336-8066-404e-a6df-cf515f55236e
# ╠═e6c0702b-240a-4da8-a791-e360471bb990
# ╠═473f9d7c-f0b6-4d44-bef2-207dfb42f2d4
# ╠═c4a8a402-f700-48c4-b13f-61f6a84b4dd5
# ╟─e5f3464e-fe56-45c0-8e47-ca4f1321b1b5
# ╟─bc0d7f88-f34b-4323-9e1a-4809ceac9efd
# ╠═3e09464b-a09f-4cc8-8533-d2d8cdf40509
# ╠═acd24cd4-3676-4a54-af7b-690564f97425
# ╟─d32c241f-97c8-46db-9512-1c7111e03e99
# ╟─9bce1738-531c-459e-95ac-ba5526df59be
# ╠═4c9becb2-07d4-4126-ab76-a741e97abb85
# ╠═078700b7-0b47-43a1-85d2-4567392e3479
# ╠═9b48e5c5-ef7a-4bd2-9b8c-213d706d44c8
# ╠═9e2be75a-7d0b-4c65-b750-845798e9223d
# ╟─711972c6-dfb3-4143-b619-0cba4e96a2e4
# ╠═bcff1c95-3d8e-40a1-87e0-3c7a36cd4206
# ╟─170c9ca0-bd1c-41d8-908e-4b02e58a081d
# ╠═2c3eef18-c8fb-42a0-ba56-607cc99d1b36
# ╠═056c6791-ed98-4c49-a424-e1ac72e6b22a
# ╟─2c9bb5de-0fc6-4f27-a5ce-0fca1f871baa
# ╟─10642e82-6dbf-46a6-9914-fe25148365e6
# ╟─1d879108-07c7-4027-b300-e7855d1e5ea9
# ╟─d3464b51-4652-49ba-87e4-545cec44dafe
