### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 6f2ac3f4-d1bb-4564-b35e-20c2a98b386b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ f8bc39c3-c284-46ce-987f-c8fa810dc207
md"# MUST & M3D"

# ╔═╡ 44248cdb-a3e0-4a3a-862f-9384dcdc1a3d
PythonPlot.matplotlib.rcParams["font.size"] = 12

# ╔═╡ 1c1b28b5-8319-4cd1-a72a-750c7aed124a
plt = PythonPlot.matplotlib.pyplot

# ╔═╡ 24364b76-d62e-4da1-9c1f-e8b9fa53d60d
visual = MUST.ingredients("visual.jl")

# ╔═╡ 93d55469-8b09-48c1-9138-5783933a96ff
pc(x) = MUST.pyconvert(Any, x)

# ╔═╡ 78fdd0d7-2853-4f52-8a3b-39a62818779d
md"## Dispatch & M3D
We load the Python modules for handling Dispatch, M3D and loading the correct module locations."

# ╔═╡ 031de813-3a84-489b-9279-dcfb3f27a6bb
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 9ce0786f-fef0-4d1e-9a38-07ea96bd147f
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 2a19b30b-1f6c-4b50-8137-9a748d2dd028
md"## Model Atmospheres
Pick the model atmosphere. It should either be located in input_multi3d/MUST or the correct folder should be given to the namelist function."

# ╔═╡ 04f7c167-3b1d-4ad8-be9f-58fbac9a20a0
#modelatmos = "m3dis_sun_magg22_10x10x299_2" 
#modelatmos = "marcs_sun.mod"
modelatmos = "DIS_MARCS_E_t5777g44m00_v0.1.m3da"

# ╔═╡ b7790d90-97bf-4f89-8c15-0d70f855471d
#modelatmosfolder = "input_multi3d/magg2022_150x300/"
#modelatmosfolder = "input_multi3d/magg2022_150x300/"
modelatmosfolder = "input_multi3d/grd/"

# ╔═╡ 3fa3d345-62c7-4f8a-ab08-a570c0aeca82
model_type = "Text"

# ╔═╡ 4b3cb27a-b0c0-40e8-b7a9-81e0c06a2391
md"The initial condition to Dispatch and to M3D is slightly different, one can convert is by simply inverting (top to bottom) and change the minus sign. The electron pressure field can be left empty or filled from an EoS."

# ╔═╡ 567f1478-7691-4d33-b9e0-4f9af78572a8
function average3DBox(path)
	z, T, r = Float64[], Float64[], Float64[]
	zTprp = TSO.readdlm(path, skipstart=2)

	x = zeros(1, 1, size(zTprp, 1))
	y = zeros(1, 1, size(zTprp, 1))
	z = reshape(zTprp[:, 1], 1, 1, :)
	T = reshape(zTprp[:, 2], 1, 1, :)
	r = reshape(zTprp[:, 4], 1, 1, :)
	
	data = Dict{Symbol, typeof(z)}(
        :T=>T, :d=>r
    )
	
	b = Box(x, y, z, data, MUST.AtmosphericParameters())
	flip!(b)
end

# ╔═╡ 6634146d-1c22-4f65-8563-7f327878de7e
begin
	function average3DMulti(path, path_new; eos=nothing)
		# Load the the model 
		m = Average3D(path)
		k = length(m.z)
		z = - reverse(m.z)
		T = exp.(reverse(m.lnT))
		r = exp.(reverse(m.lnρ))
		p = zeros(eltype(T), k)

		if !isnothing(eos)
			p .= exp.(lookup(eos, :lnNe, r, T)) .* MUST.KBoltzmann .* T
		end
		
		open(path_new, "w") do f
			write(f, "$(path)\n")
			write(f, "$(k)\n")
			for i in eachindex(z)
				write(f, "$(z[i]) $(T[i]) $(p[i]) $(r[i]) $(0.0)")
				if i < length(z)
					write(f, "\n")
				end
			end
		end
	
		path_new
	end

	average3DMulti(path; kwargs...) = begin
		actual_path = first(split(path, ".m3da"))
		old_path = joinpath(actual_path, "inim.dat")

		average3DMulti(old_path, path; kwargs...)
	end
end

# ╔═╡ 6ccb7b73-2ea1-40d8-9ade-55948682cde9
md"## Line List
One can either use no line list to compute Teff, use a line list in the Turbospectrum format, or use an absmet file. The absmet file contains sampled backgournd opacities, and might be good enough for the spectrum computation in terms of Teff. To activate the absmet, one has to set the linelist to nothing, in order for the namelist to be created correctly."

# ╔═╡ 03ed4c44-da3b-45bc-946b-634ca23cf30a
linelist = "./input_multi3d/vald_2490-25540.list"

# ╔═╡ 484b7245-8e34-41b2-b9f8-2a9b3c1ba9e4
absmet = "./input_multi3d/absmet"

# ╔═╡ fe9351c2-c97b-4a66-a230-f2a0ab8067ce
md"## Opacity tables
We compute the heating from the unbinned and binned opacity table and then compare the results. Both is possible within M3D."

# ╔═╡ 5185341e-9299-40d5-9881-9e002c926bfc
eos = reload(
	SqEoS, 
	"/u/peitner/DISPATCH/binned/DIS_MARCS_v1.7.2/eos.hdf5"
)

# ╔═╡ fa7355e5-9072-4c4e-be9a-d5412f6cde5c
opa = reload(
	SqOpacity, 
	"/u/peitner/DISPATCH/binned/DIS_MARCS_v1.7.2/binned_opacities.hdf5"
)

# ╔═╡ 69ee39a1-ae7d-4ec9-bae2-b830f80836ed
md"These tables need to be combined and saved in the Multi format in order to be able to run them in M3D. We save them in common folder directly in Multi to make things easier."

# ╔═╡ e830f3c8-ce12-4e46-94fe-33b847ae9b61
folder_eos_rel = "input_multi3d/M3D_MARCS_v1.7.2"

# ╔═╡ fb609e98-ab06-487e-8bf1-91e232e5687b
folder_eos = MUST.@in_m3dis folder_eos_rel

# ╔═╡ 7d38779c-c2e3-4ff2-b30c-f4229f35545f
!isdir(folder_eos) && mkdir(folder_eos)

# ╔═╡ 3c76bd08-a807-4398-b277-a622e010d05e
m3d_table = save(eos, opa, joinpath(folder_eos, "eos_opa"))

# ╔═╡ 778ece71-8bb6-493c-bb4f-390f778a947c
if model_type == "Text"
	average3DMulti(MUST.@in_m3dis(joinpath(modelatmosfolder, modelatmos)), eos=eos)
end

# ╔═╡ 81ef72b6-eafa-49ef-9923-fb290cc4091d
md"## Running M3D
First we need to run M3D in the unbinned case, compute the heating, and the run it in the binned case for comparison. If you have computed this already, you can also read an existing output."

# ╔═╡ f9e7cd2c-7a09-4eb5-b0a1-1247ff4fbdf3
compute = true 

# ╔═╡ cbe117c0-f6c7-44bb-a86e-3d57362666d7
m3d_binned = if compute 
	MUST.heating(
		modelatmos, 
		joinpath(folder_eos_rel, "eos_opa"),
			namelist_kwargs=(
				:model_folder=>modelatmosfolder,
				:linelist=>nothing,
				:absmet=>nothing,
				:atom_params=>(:atom_file=>"", ),
				:spectrum_params=>(:daa=>1., :aa_blue=>1500, :aa_red=>4000),
				:atmos_params=>(
					:dims=>1, 
					:atmos_format=>model_type,
					:use_density=>true, 
					:use_ne=>false
				),
				:m3d_params=>(:n_nu=>1, )
			),
		slurm=false
	)
else
	MUST.M3DISRun("data/$(modelatmos)_binned")
end

# ╔═╡ 4668ca71-bd41-413a-8c67-d7e7dfe9b6ba
m3d_unbinned = if compute
	MUST.heating(
		modelatmos, 
		namelist_kwargs=(
			:model_folder=>modelatmosfolder,
			:linelist=>nothing,
			:absmet=>absmet,
			:atom_params=>(:atom_file=>"", ),
			:spectrum_params=>(:daa=>1., :aa_blue=>1500, :aa_red=>9000),
			:atmos_params=>(
				:dims=>1, 
				:atmos_format=>model_type,
				:use_density=>true, 
				:use_ne=>false
			),
			:m3d_params=>(:n_nu=>1, )
		),
		slurm=false
	)
else
	MUST.M3DISRun("data/$(modelatmos)_unbinned")
end

# ╔═╡ fb0a997b-f2dd-4a19-bdd6-c5293159eb2c
md"## Heating Rate"

# ╔═╡ 44f48b12-f884-49a6-94b9-728715daa2ab
function get_heating(run)
	qrad = MUST.pyconvert(
		Array, run.run.read_patch_save("Qrad")[0][0].compute()
	)[:, :, :, 1, 1]
	ltau = MUST.pyconvert(Array, run.run.tau)

	ltau, qrad
end

# ╔═╡ 4fe50a65-7729-4a13-ab39-d81c62fae60d
tau_unbinned, qrad_unbinned = get_heating(m3d_unbinned)

# ╔═╡ 12e97194-909a-493a-a242-9bf21e736713
tau_binned, qrad_binned = get_heating(m3d_binned)

# ╔═╡ 858b7fc4-9d95-4a18-aabb-7621ccf90b91
md"## API Integration
The heating can be integrated into the ```Box``` API by loading the corresponding Multi Box"

# ╔═╡ 0d6bc811-6e93-4f79-bbb2-1081c9314e7b
modelBox = flip!(
	if model_type=="Marcs"
		marcsBox(
			MUST.@in_m3dis(joinpath(modelatmosfolder, modelatmos))
		)
	elseif model_type=="Text"
		average3DBox(
			MUST.@in_m3dis(joinpath(modelatmosfolder, modelatmos))
		)
	else
		multiBox(
			MUST.@in_m3dis(joinpath(modelatmosfolder, modelatmos))
		)	
	end 
)

# ╔═╡ 2295240c-cca2-441a-b901-437e1d7f8924
keys(modelBox.data)

# ╔═╡ 5d155f12-b42c-43da-918e-6714bfbb0eac
md"For average comparison, we need to compare them on the optical depth scale. For this the cube has to be interpolated to this scale."

# ╔═╡ ae25d7f5-afad-4e41-8125-336c67047d33
function attach_heating(box, run)
	t, qr = get_heating(run)

	box_new = if length(run.zz) != size(box, 3)
		z_new = range(0, 1, length=size(box, 3)) |> collect
		
		ip = ginterpolate(MUST.Grid(run.zz), MUST.Grid(z_new))
		
		box_new = deepcopy(box)
		box_new[:qrad] = reshape(gevaluate!(ip, reshape(qr, :)), size(box)...)
		box_new[:tau] = reshape(gevaluate!(ip, reshape(t, :)), size(box)...)

		box_new
	else
		box_new = deepcopy(box)
		box_new[:qrad] = reshape(qr, size(box)...)
		box_new[:tau] = reshape(t, size(box)...)

		box_new
	end

	# Interpolate to new optical depth scale
	# for this use the logtau that was computed from Multi
	MUST.height_scale_fast(box_new, :tau)
end

# ╔═╡ e85054fc-fa63-48cc-959a-1fa4c4df6965
modelbinned = attach_heating(modelBox, m3d_binned)

# ╔═╡ 06168a77-8a58-41a6-8305-dcc32a71ab7f
modelunbinned = attach_heating(modelBox, m3d_unbinned)

# ╔═╡ 638dd690-d2f0-4381-b371-8fb444c18e5f
profile(MUST.std, modelbinned, :log10tau, :qrad)

# ╔═╡ 134144b6-6f21-42d5-a816-fd5bb57b9e55
begin
	plt.close()
	fH, axH = plt.subplots(2, 1, figsize=(7, 6), sharex=true)
	visual.basic_plot!.(axH)

	rms(x) = sqrt(MUST.mean(x .^2))
	stat = MUST.mean
	
	axH[0].plot(
		profile(stat, modelbinned, :log10tau, :qrad)..., 
		label="binned",
		color="r"
	)
	
	axH[1].plot(
		profile(stat, modelunbinned, :log10tau, :qrad)..., 
		label="unbinned",
		color="k"
	)

	axH[0].legend(framealpha=0)
	axH[1].legend(framealpha=0)

	gcf()
end

# ╔═╡ 6618b520-84f2-4fc6-a18a-dcec79cfff79
begin
	plt.close()
	fH2, axH2 = plt.subplots(2, 1, figsize=(7, 6), sharex=true)
	visual.basic_plot!.(axH2)

	i_specific = 1
	
	axH2[0].plot(
		log10.(modelbinned[:tau][i_specific, i_specific, :]),
		modelbinned[:qrad][i_specific, i_specific, :],
		label="binned",
		color="r"
	)
	
	axH2[1].plot(
		log10.(modelbinned[:tau][i_specific, i_specific, :]),
		modelunbinned[:qrad][i_specific, i_specific, :],
		label="unbinned",
		color="k"
	)

	axH2[0].legend(framealpha=0)
	axH2[1].legend(framealpha=0)

	gcf()
end

# ╔═╡ Cell order:
# ╟─f8bc39c3-c284-46ce-987f-c8fa810dc207
# ╠═6f2ac3f4-d1bb-4564-b35e-20c2a98b386b
# ╠═44248cdb-a3e0-4a3a-862f-9384dcdc1a3d
# ╠═1c1b28b5-8319-4cd1-a72a-750c7aed124a
# ╠═24364b76-d62e-4da1-9c1f-e8b9fa53d60d
# ╠═93d55469-8b09-48c1-9138-5783933a96ff
# ╟─78fdd0d7-2853-4f52-8a3b-39a62818779d
# ╠═031de813-3a84-489b-9279-dcfb3f27a6bb
# ╠═9ce0786f-fef0-4d1e-9a38-07ea96bd147f
# ╟─2a19b30b-1f6c-4b50-8137-9a748d2dd028
# ╠═04f7c167-3b1d-4ad8-be9f-58fbac9a20a0
# ╠═b7790d90-97bf-4f89-8c15-0d70f855471d
# ╠═3fa3d345-62c7-4f8a-ab08-a570c0aeca82
# ╟─4b3cb27a-b0c0-40e8-b7a9-81e0c06a2391
# ╟─567f1478-7691-4d33-b9e0-4f9af78572a8
# ╠═6634146d-1c22-4f65-8563-7f327878de7e
# ╟─6ccb7b73-2ea1-40d8-9ade-55948682cde9
# ╠═03ed4c44-da3b-45bc-946b-634ca23cf30a
# ╠═484b7245-8e34-41b2-b9f8-2a9b3c1ba9e4
# ╟─fe9351c2-c97b-4a66-a230-f2a0ab8067ce
# ╠═5185341e-9299-40d5-9881-9e002c926bfc
# ╠═fa7355e5-9072-4c4e-be9a-d5412f6cde5c
# ╟─69ee39a1-ae7d-4ec9-bae2-b830f80836ed
# ╠═e830f3c8-ce12-4e46-94fe-33b847ae9b61
# ╠═fb609e98-ab06-487e-8bf1-91e232e5687b
# ╠═7d38779c-c2e3-4ff2-b30c-f4229f35545f
# ╠═3c76bd08-a807-4398-b277-a622e010d05e
# ╠═778ece71-8bb6-493c-bb4f-390f778a947c
# ╟─81ef72b6-eafa-49ef-9923-fb290cc4091d
# ╠═f9e7cd2c-7a09-4eb5-b0a1-1247ff4fbdf3
# ╠═cbe117c0-f6c7-44bb-a86e-3d57362666d7
# ╠═4668ca71-bd41-413a-8c67-d7e7dfe9b6ba
# ╟─fb0a997b-f2dd-4a19-bdd6-c5293159eb2c
# ╟─44f48b12-f884-49a6-94b9-728715daa2ab
# ╠═4fe50a65-7729-4a13-ab39-d81c62fae60d
# ╠═12e97194-909a-493a-a242-9bf21e736713
# ╟─858b7fc4-9d95-4a18-aabb-7621ccf90b91
# ╠═0d6bc811-6e93-4f79-bbb2-1081c9314e7b
# ╠═2295240c-cca2-441a-b901-437e1d7f8924
# ╟─5d155f12-b42c-43da-918e-6714bfbb0eac
# ╟─ae25d7f5-afad-4e41-8125-336c67047d33
# ╠═e85054fc-fa63-48cc-959a-1fa4c4df6965
# ╠═06168a77-8a58-41a6-8305-dcc32a71ab7f
# ╠═638dd690-d2f0-4381-b371-8fb444c18e5f
# ╟─134144b6-6f21-42d5-a816-fd5bb57b9e55
# ╟─6618b520-84f2-4fc6-a18a-dcec79cfff79
