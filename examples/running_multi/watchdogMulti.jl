### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 16244768-3acc-11ef-2fbf-0d708f18e84d
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PythonPlot
end

# ╔═╡ 003095e1-f6e9-4bd8-bfdc-6008ffa505af
plt = matplotlib.pyplot

# ╔═╡ 6ab17505-ff25-495b-a546-852da7e5619c
MUST.@import_m3dis "../../../Multi3D"

# ╔═╡ a70742e9-c450-4ae6-a80a-b5a620acbdc1
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ a32be968-0714-47f8-8623-e559a0e3fc31
begin
	λs = 6100
	λe = 6790
	Δλ = 0.5
end

# ╔═╡ deedf225-cedf-409d-8833-46def77393c8
(λe - λs) / Δλ

# ╔═╡ 5ca326db-5903-4947-a6b5-af9f8792fbf9
name = "grid_t5777g44m00_h1.1"

# ╔═╡ 0b00f487-b3a5-44f6-af94-c0e0f1955823
get_teff(name) = begin
	names = split(name, "_", keepempty=false)
	i_name = 0
	for (i, split) in enumerate(names)
		if occursin("t", split) & occursin("g", split) & occursin("m", split)
			i_name = i
		end
	end

	relevantPart = names[i_name]
	t1 = findfirst(x->x=='t', relevantPart) + 1
	t2 = findfirst(x->x=='g', relevantPart) - 1
	t3 = findfirst(x->x=='m', relevantPart) - 1

	t = parse(Float64, relevantPart[t1:t2])
	g = parse(Float64, relevantPart[t2+2:t3]) ./10
	m = parse(Float64, relevantPart[t3+2:end])

	t, g, m = if i_name == 0
		@warn "could not guess atmospheric parameters from name."
		-1, -1, 0.0
	else
		t, g, m
	end

	t, g, m
end

# ╔═╡ d9b21458-2a4a-4759-b5ef-84c74647bc24
get_teff(name)

# ╔═╡ 1e53a540-17b7-463b-bc4e-51c4aebbdb72
run = @in_dispatch("data/$(name)")

# ╔═╡ f56644af-c01c-4d3a-8e78-d675560432b1
isnap = 302

# ╔═╡ e1aaad21-9340-4e31-9644-7f11b446e859
multi_name = joinpath(run, "m3dis_$(isnap)")

# ╔═╡ b3b1f71c-dfff-4c10-af63-967634b3d92e
b, bτ = pick_snapshot(run, isnap)

# ╔═╡ cc4090d0-b6bf-49ab-8b9f-dbcc27a3ec37


# ╔═╡ 46021094-c218-462d-9b4c-a14b06c618d5
begin
	bDown = MUST.gresample(b; nx=40, ny=40, nz=size(b, 3) * 2 - 1)
	bDown.data[:ne] = pop!(bDown.data, :Ne)
end;

# ╔═╡ 9ba6635c-d71e-4ca1-8ee9-2f2004be10cb
MUST.multiBox(bDown, joinpath(run, "m3dis_$(isnap)"))

# ╔═╡ bc8b742e-34e9-4148-917c-b054d9a0e952


# ╔═╡ 4d92338e-7b8d-47f0-8b07-d838ba4cad00
linelists = String[
	"/home/eitner/shared/StAt/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata",
	"input_multi3d/nlte_ges_linelist_jmg04sep2023_I_II"
]

# ╔═╡ 2941966e-ddf5-453e-bc89-63068ccc8eff
spectrum_namelist = Dict(
	:model_folder=>run,
	:linelist=>nothing,
	:absmet=>nothing,
	:linelist_params=>(:line_lists=>linelists,),
	:atom_params=>(:atom_file=>"",),
	:spectrum_params=>(:daa=>Δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>false),
	:atmos_params=>(
		:dims=>1, 
		:atmos_format=>"must",
		:use_density=>true, 
		:use_ne=>false,
		:FeH=>0.0
	),
	:m3d_params=>(
		:save_resolved=>true,
		:long_scheme=>"disk_center"
	),
	
)

# ╔═╡ b18dc1cc-153e-49e3-bc97-e03c306c6dc5
m3dis_kwargs = Dict(
	:threads=>2
)

# ╔═╡ eed8ef16-ee05-4245-9ae6-a57ec85f4eec
result = MUST.spectrum(
	"m3dis_$(isnap)"; name=name, NLTE=false, slurm=false, namelist_kwargs=spectrum_namelist, m3dis_kwargs=m3dis_kwargs
)

# ╔═╡ 2f718c66-6c8c-4fae-953d-83c8cc676e06
#result = MUST.M3DISRun("m3dis_300_$(name)", @in_m3dis("data"))

# ╔═╡ 309a3b42-0f48-4618-a1f1-83467df6e02d
i = MUST.pyconvert(Array, result.i3)

# ╔═╡ fa14ada8-02d3-4202-99a0-a9b7f4339b66
c = MUST.pyconvert(Array, result.c3)

# ╔═╡ 729c8ba1-2460-417c-a72a-f1a2ddc505c5
let
	plt.close()
	f, ax = plt.subplots()

	ax.plot(result.lam, i[:,10,10]./c[:,10,10])

	#ax.set_xlim(6500, 6600)
	f
end

# ╔═╡ 22ad2a7e-a358-4aac-8f6d-4b6ce97f4d0b
let
	plt.imshow(result.temp[:,:,end]')

	gcf()
end

# ╔═╡ 782d5099-c19d-44ff-8be3-bea4c340ea43
l = MUST.pyconvert(Array, result.xx)

# ╔═╡ Cell order:
# ╠═16244768-3acc-11ef-2fbf-0d708f18e84d
# ╠═003095e1-f6e9-4bd8-bfdc-6008ffa505af
# ╠═6ab17505-ff25-495b-a546-852da7e5619c
# ╠═a70742e9-c450-4ae6-a80a-b5a620acbdc1
# ╠═a32be968-0714-47f8-8623-e559a0e3fc31
# ╠═deedf225-cedf-409d-8833-46def77393c8
# ╠═5ca326db-5903-4947-a6b5-af9f8792fbf9
# ╠═0b00f487-b3a5-44f6-af94-c0e0f1955823
# ╠═d9b21458-2a4a-4759-b5ef-84c74647bc24
# ╠═1e53a540-17b7-463b-bc4e-51c4aebbdb72
# ╠═f56644af-c01c-4d3a-8e78-d675560432b1
# ╠═e1aaad21-9340-4e31-9644-7f11b446e859
# ╠═b3b1f71c-dfff-4c10-af63-967634b3d92e
# ╟─cc4090d0-b6bf-49ab-8b9f-dbcc27a3ec37
# ╠═46021094-c218-462d-9b4c-a14b06c618d5
# ╠═9ba6635c-d71e-4ca1-8ee9-2f2004be10cb
# ╟─bc8b742e-34e9-4148-917c-b054d9a0e952
# ╠═4d92338e-7b8d-47f0-8b07-d838ba4cad00
# ╠═2941966e-ddf5-453e-bc89-63068ccc8eff
# ╠═b18dc1cc-153e-49e3-bc97-e03c306c6dc5
# ╠═eed8ef16-ee05-4245-9ae6-a57ec85f4eec
# ╠═2f718c66-6c8c-4fae-953d-83c8cc676e06
# ╠═309a3b42-0f48-4618-a1f1-83467df6e02d
# ╠═fa14ada8-02d3-4202-99a0-a9b7f4339b66
# ╠═729c8ba1-2460-417c-a72a-f1a2ddc505c5
# ╠═22ad2a7e-a358-4aac-8f6d-4b6ce97f4d0b
# ╠═782d5099-c19d-44ff-8be3-bea4c340ea43
