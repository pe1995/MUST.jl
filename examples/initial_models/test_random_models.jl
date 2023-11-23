### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 15678f68-882c-11ee-05ce-fd760daebd8d
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using Plots
	using Printf
	using DataFrames
	using DelimitedFiles 
	using LazySets
	using Polyhedra
end

# ╔═╡ dc53b3e9-ab5f-4d14-828b-5fc961faf37d
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ bc1caffd-1787-40b3-8ba0-68a3eafe09b3
staggergrid = MUST.StaggerGrid("stagger_grid_full.mgrid")

# ╔═╡ bf868bbc-a51f-4d3a-9357-64d2eb05addd
randomgrid = MUST.StaggerGrid("random_models.mgrid")

# ╔═╡ 8a44b83c-eeea-4224-a7ef-52d98895fa9f
eos_folder = "/u/peitner/DISPATCH/TSO.jl/examples/binning_opacities/tables/TSO_MARCS_v1.6"

# ╔═╡ 9de39e0c-7525-42d3-86ac-2e5c06e37759
i_test = findfirst(randomgrid["name"] .== "t45.86g44.00m0.000")

# ╔═╡ 4ce7b47f-f8b8-400e-9c41-4c1b0b09046a
edges = [
	i for i in 1:nrow(staggergrid.info) 
	if staggergrid["name", i] ∈ ["t45g40m00", "t45g45m00", "t50g40m00", "t50g45m00", "t50g45p05"]
]

# ╔═╡ 01d636c2-db4d-45af-8a04-0af50443fe64


# ╔═╡ 3b36e263-fef2-4c16-8c68-b20dbca3c34d
eos = reload(SqEoS, joinpath(eos_folder, "combined_ross_eos_magg22.hdf5"))

# ╔═╡ c8578ef3-4eaa-4b51-8802-2557b78c6aa5
opa = reload(
	SqOpacity, 
	joinpath(eos_folder, "combined_opacities_magg22.hdf5"), 
	mmap=true
)

# ╔═╡ 6fe9c410-f5f4-459a-bfe7-7adde931838c


# ╔═╡ fc9e065d-c19a-4edb-95a7-322ac5a575b4
model = @optical(
	Average3D(eos, randomgrid["av_path", i_test], logg=randomgrid["logg", i_test]), eos, opa
)

# ╔═╡ a2da7fef-8afb-4f73-898f-944122dff5f3
models_stagger = [
	@optical(
		Average3D(eos, staggergrid["av_path", i], logg=staggergrid["logg", i]), eos, opa
	)
	for i in edges
]

# ╔═╡ 171e7f10-8d6d-4670-9c92-5027efb0ce29


# ╔═╡ e3f4417e-9aa2-401d-8b59-b9611c7e8bad
function interpolate_average_test(grid; teff, logg, feh, common_size=1000)
	logg_gr = grid.info[!, "logg"]
	teff_gr = grid.info[!, "teff"]
	feh_gr  = grid.info[!, "feh"]
	femask = trues(length(logg_gr)) #feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]
	
	# read all average models and interpolate them to the same number of points
	models = [modelgrids.initial_model(grid, String(name))
					for name in grid.info[!, "name"]]
	r_models = [minimum(abs.(diff(m.z))) for m in models]
	
	models = [TSO.upsample(m, common_size) for m in models]

	# now we interpolate all points to one common point, for every point
	scatter_int(v, x, y, z) = MUST.pyconvert(typeof(x),
		first(MUST.scipy_interpolate.griddata((teff_gr[femask], logg_gr[femask],  feh_gr[femask]), v, ([x], [y], [z]), method="linear"))
	)
	
	points = zeros(eltype(models[1].z), length(models), 3)
	z = zeros(eltype(models[1].z), common_size)
	t = zeros(eltype(models[1].z), common_size)
	d = zeros(eltype(models[1].z), common_size)
	
	for i in 1:common_size
		for j in eachindex(models)
			points[j, 1] = models[j].z[i]
			points[j, 2] = models[j].lnT[i]
			points[j, 3] = models[j].lnρ[i]
		end
		z[i] = scatter_int(points[femask, 1], teff, logg, feh)
		t[i] = scatter_int(points[femask, 2], teff, logg, feh)
		d[i] = scatter_int(points[femask, 3], teff, logg, feh)
	end

	r_target = scatter_int(r_models[femask], teff, logg, feh)
	npoints  = ceil(Int, abs(maximum(z) - minimum(z)) / r_target)

	TSO.upsample(Model1D(z=z, lnT=t, lnρ=d, logg=logg), npoints)
end

# ╔═╡ dded0e6d-d0f4-4eff-8a91-8ab32791c61c
m_ip = interpolate_average_test(
	staggergrid,
	teff=randomgrid["teff", i_test], 
	logg=randomgrid["logg", i_test], 
	feh=randomgrid["feh", i_test]
)

# ╔═╡ 86d785dc-befc-47f4-aac5-d7e23652a714


# ╔═╡ 0ec025ac-a545-4dcb-94b4-8c366e604a17
md"An adiabat for this model would look like this"

# ╔═╡ 9d8528ed-de0f-4f74-99fe-bd864442a9ba
ad = modelgrids.adiabat(
	eos_folder, 
	randomgrid["av_path", i_test], 
	randomgrid["logg", i_test],
	1000,
	ee_min=randomgrid["ee_min", i_test],
	eospath=joinpath(randomgrid["binned_E_tables", i_test], "eos.hdf5")
)

# ╔═╡ bdf65b83-2b93-471b-922f-21464b4dcea1


# ╔═╡ 602803ec-1801-4b4e-a18d-4c41de0451ab
model_stagger_tau = """
-5.0000    3.52048E+03    3.24575E+10    0.00000E+00    6.95897E-01
     -4.9002    3.53195E+03    3.66174E+10    0.00000E+00    6.76037E-01
     -4.7999    3.54419E+03    4.13680E+10    0.00000E+00    6.57935E-01
     -4.7001    3.55721E+03    4.67822E+10    0.00000E+00    6.41492E-01
     -4.5998    3.57105E+03    5.29323E+10    0.00000E+00    6.26633E-01
     -4.5000    3.58572E+03    5.99309E+10    0.00000E+00    6.13354E-01
     -4.4001    3.60124E+03    6.79273E+10    0.00000E+00    6.01640E-01
     -4.2999    3.61760E+03    7.70710E+10    0.00000E+00    5.91493E-01
     -4.2001    3.63479E+03    8.75081E+10    0.00000E+00    5.82891E-01
     -4.0999    3.65281E+03    9.94247E+10    0.00000E+00    5.75761E-01
     -4.0000    3.67168E+03    1.13078E+11    0.00000E+00    5.69995E-01
     -3.9001    3.69135E+03    1.28768E+11    0.00000E+00    5.65460E-01
     -3.7999    3.71176E+03    1.46774E+11    0.00000E+00    5.62029E-01
     -3.7001    3.73283E+03    1.67363E+11    0.00000E+00    5.59638E-01
     -3.5999    3.75444E+03    1.90893E+11    0.00000E+00    5.58306E-01
     -3.5000    3.77645E+03    2.17818E+11    0.00000E+00    5.58117E-01
     -3.4001    3.79868E+03    2.48576E+11    0.00000E+00    5.59198E-01
     -3.2999    3.82100E+03    2.83497E+11    0.00000E+00    5.61655E-01
     -3.2001    3.84323E+03    3.22937E+11    0.00000E+00    5.65558E-01
     -3.0999    3.86528E+03    3.67461E+11    0.00000E+00    5.70957E-01
     -3.0000    3.88706E+03    4.17736E+11    0.00000E+00    5.77920E-01
     -2.9000    3.90849E+03    4.74346E+11    0.00000E+00    5.86492E-01
     -2.7999    3.92955E+03    5.37766E+11    0.00000E+00    5.96757E-01
     -2.7001    3.95021E+03    6.08656E+11    0.00000E+00    6.08779E-01
     -2.6000    3.97054E+03    6.88100E+11    0.00000E+00    6.22544E-01
     -2.5000    3.99066E+03    7.77438E+11    0.00000E+00    6.38016E-01
     -2.4000    4.01073E+03    8.78055E+11    0.00000E+00    6.55131E-01
     -2.2999    4.03100E+03    9.91291E+11    0.00000E+00    6.73759E-01
     -2.2001    4.05174E+03    1.11914E+12    0.00000E+00    6.93718E-01
     -2.1001    4.07333E+03    1.26473E+12    0.00000E+00    7.14746E-01
     -2.0000    4.09615E+03    1.43215E+12    0.00000E+00    7.36549E-01
     -1.8999    4.12064E+03    1.62566E+12    0.00000E+00    7.58815E-01
     -1.8000    4.14724E+03    1.84958E+12    0.00000E+00    7.81247E-01
     -1.7000    4.17638E+03    2.11027E+12    0.00000E+00    8.03549E-01
     -1.6001    4.20849E+03    2.41712E+12    0.00000E+00    8.25460E-01
     -1.5000    4.24400E+03    2.78203E+12    0.00000E+00    8.46802E-01
     -1.3999    4.28330E+03    3.21703E+12    0.00000E+00    8.67539E-01
     -1.3000    4.32685E+03    3.73639E+12    0.00000E+00    8.87665E-01
     -1.2000    4.37516E+03    4.36004E+12    0.00000E+00    9.07161E-01
     -1.1000    4.42891E+03    5.11445E+12    0.00000E+00    9.25993E-01
     -1.0000    4.48891E+03    6.03152E+12    0.00000E+00    9.44164E-01
     -0.9000    4.55610E+03    7.13516E+12    0.00000E+00    9.61737E-01
     -0.8000    4.63155E+03    8.44901E+12    0.00000E+00    9.78770E-01
     -0.7000    4.71640E+03    1.00417E+13    0.00000E+00    9.95258E-01
     -0.6000    4.81182E+03    1.20245E+13    0.00000E+00    1.01117E+00
     -0.5000    4.91916E+03    1.44861E+13    0.00000E+00    1.02649E+00
     -0.4000    5.04013E+03    1.73973E+13    0.00000E+00    1.04123E+00
     -0.3000    5.17687E+03    2.08332E+13    0.00000E+00    1.05536E+00
     -0.2000    5.33169E+03    2.51234E+13    0.00000E+00    1.06889E+00
     -0.1000    5.50683E+03    3.06881E+13    0.00000E+00    1.08181E+00
      0.0000    5.70407E+03    3.82205E+13    0.00000E+00    1.09407E+00
      0.1000    5.92383E+03    4.97700E+13    0.00000E+00    1.10616E+00
      0.2000    6.16502E+03    6.87113E+13    0.00000E+00    1.11903E+00
      0.3000    6.42576E+03    1.01397E+14    0.00000E+00    1.13352E+00
      0.4000    6.70395E+03    1.57348E+14    0.00000E+00    1.14972E+00
      0.5000    6.99708E+03    2.50833E+14    0.00000E+00    1.16723E+00
      0.6000    7.30120E+03    4.00832E+14    0.00000E+00    1.18491E+00
      0.7000    7.61019E+03    6.27517E+14    0.00000E+00    1.20156E+00
      0.8000    7.91553E+03    9.45449E+14    0.00000E+00    1.21612E+00
      0.9000    8.20706E+03    1.35463E+15    0.00000E+00    1.22789E+00
      1.0000    8.47527E+03    1.83660E+15    0.00000E+00    1.23624E+00
      1.1000    8.71380E+03    2.36141E+15    0.00000E+00    1.24088E+00
      1.2000    8.92098E+03    2.90109E+15    0.00000E+00    1.24186E+00
      1.3000    9.09954E+03    3.44038E+15    0.00000E+00    1.23954E+00
      1.3999    9.25489E+03    3.97931E+15    0.00000E+00    1.23441E+00
      1.5000    9.39308E+03    4.52831E+15    0.00000E+00    1.22689E+00
      1.6001    9.51944E+03    5.10247E+15    0.00000E+00    1.21728E+00
      1.7000    9.63805E+03    5.71715E+15    0.00000E+00    1.20599E+00
      1.8000    9.75192E+03    6.38727E+15    0.00000E+00    1.19341E+00
      1.8999    9.86319E+03    7.12750E+15    0.00000E+00    1.17982E+00
      2.0000    9.97332E+03    7.95177E+15    0.00000E+00    1.16519E+00
      2.1001    1.00835E+04    8.87601E+15    0.00000E+00    1.14922E+00
      2.2001    1.01943E+04    9.91547E+15    0.00000E+00    1.13238E+00
      2.2999    1.03064E+04    1.10886E+16    0.00000E+00    1.11466E+00
      2.4000    1.04202E+04    1.24142E+16    0.00000E+00    1.09615E+00
      2.5000    1.05362E+04    1.39159E+16    0.00000E+00    1.07710E+00
      2.6000    1.06543E+04    1.56164E+16    0.00000E+00    1.05760E+00
      2.7001    1.07747E+04    1.75426E+16    0.00000E+00    1.03741E+00
      2.7999    1.08973E+04    1.97253E+16    0.00000E+00    1.01592E+00
      2.9000    1.10225E+04    2.21996E+16    0.00000E+00    9.93776E-01
      3.0000    1.11501E+04    2.50034E+16    0.00000E+00    9.71223E-01
      3.0999    1.12804E+04    2.81869E+16    0.00000E+00    9.48223E-01
      3.2001    1.14133E+04    3.17983E+16    0.00000E+00    9.24106E-01
      3.2999    1.15486E+04    3.58890E+16    0.00000E+00    8.99319E-01
      3.4001    1.16863E+04    4.05252E+16    0.00000E+00    8.74476E-01
      3.5000    1.18269E+04    4.57809E+16    0.00000E+00    8.49563E-01
      3.5999    1.19701E+04    5.17316E+16    0.00000E+00    8.25113E-01
      3.7001    1.21159E+04    5.84654E+16    0.00000E+00    8.01058E-01
      3.7999    1.22645E+04    6.60893E+16    0.00000E+00    7.77811E-01
      3.9001    1.24158E+04    7.47227E+16    0.00000E+00    7.55199E-01
      4.0000    1.25701E+04    8.44967E+16    0.00000E+00    7.32756E-01
      4.0999    1.27271E+04    9.55497E+16    0.00000E+00    7.10665E-01
      4.2001    1.28868E+04    1.08036E+17    0.00000E+00    6.89607E-01
      4.2999    1.30495E+04    1.22146E+17    0.00000E+00    6.69654E-01
      4.4001    1.32155E+04    1.38090E+17    0.00000E+00    6.50448E-01
      4.5000    1.33843E+04    1.56098E+17    0.00000E+00    6.31527E-01
      4.5998    1.35563E+04    1.76449E+17    0.00000E+00    6.13193E-01
      4.7001    1.37317E+04    1.99460E+17    0.00000E+00    5.94336E-01
      4.7999    1.39107E+04    2.25485E+17    0.00000E+00    5.75421E-01
      4.9002    1.40928E+04    2.54865E+17    0.00000E+00    5.56634E-01
      5.0000    1.42786E+04    2.87999E+17    0.00000E+00    5.38068E-01
"""

# ╔═╡ a593e6fd-6702-4aaf-8917-64ae7290aeb5
open("test_model.txt", "w") do f
	write(f, model_stagger_tau)
end

# ╔═╡ 123b165f-0a6c-4894-9157-a538452cf139
data = readdlm("test_model.txt")

# ╔═╡ 20df94f0-1d93-4604-9c2c-193856d1f4f9


# ╔═╡ 250c95e4-e4c5-4272-b8d8-a1b4fd3d72d1
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	scatter!(staggergrid["teff"], staggergrid["logg"])
	scatter!([randomgrid["teff", i_test]], [randomgrid["logg", i_test]])
	
	xflip!(true)
	yflip!(true)
end

# ╔═╡ e5c4f17b-da5a-4b72-8de1-6a215804e7b0
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	#=for i in eachindex(edges)
		plot!(models_stagger[i].z, models_stagger[i].lnT, label=staggergrid["name", edges[i]])
	end=#

	plot!(model.z, model.lnT, label=randomgrid["name", i_test])
	plot!(ad.z, ad.lnT, label="adiabat")
	
	plot!()
end

# ╔═╡ 124f13ff-e0cc-4d85-9854-0b68c1061776
begin
	plot()
	
	plot!(data[:, 1], data[:, 2])
	
	for i in eachindex(edges)
		(i != length(edges)) && continue
		plot!(log10.(models_stagger[i].τ), exp.(models_stagger[i].lnT), 
				label=staggergrid["name", edges[i]])
	end

	plot!(xlim=(-8, 5), ylim=(2000, 15000))
end

# ╔═╡ Cell order:
# ╠═15678f68-882c-11ee-05ce-fd760daebd8d
# ╠═dc53b3e9-ab5f-4d14-828b-5fc961faf37d
# ╠═bc1caffd-1787-40b3-8ba0-68a3eafe09b3
# ╠═bf868bbc-a51f-4d3a-9357-64d2eb05addd
# ╠═8a44b83c-eeea-4224-a7ef-52d98895fa9f
# ╠═9de39e0c-7525-42d3-86ac-2e5c06e37759
# ╠═4ce7b47f-f8b8-400e-9c41-4c1b0b09046a
# ╟─01d636c2-db4d-45af-8a04-0af50443fe64
# ╠═3b36e263-fef2-4c16-8c68-b20dbca3c34d
# ╠═c8578ef3-4eaa-4b51-8802-2557b78c6aa5
# ╟─6fe9c410-f5f4-459a-bfe7-7adde931838c
# ╠═fc9e065d-c19a-4edb-95a7-322ac5a575b4
# ╠═a2da7fef-8afb-4f73-898f-944122dff5f3
# ╟─171e7f10-8d6d-4670-9c92-5027efb0ce29
# ╟─e3f4417e-9aa2-401d-8b59-b9611c7e8bad
# ╠═dded0e6d-d0f4-4eff-8a91-8ab32791c61c
# ╟─86d785dc-befc-47f4-aac5-d7e23652a714
# ╟─0ec025ac-a545-4dcb-94b4-8c366e604a17
# ╠═9d8528ed-de0f-4f74-99fe-bd864442a9ba
# ╟─bdf65b83-2b93-471b-922f-21464b4dcea1
# ╟─602803ec-1801-4b4e-a18d-4c41de0451ab
# ╟─a593e6fd-6702-4aaf-8917-64ae7290aeb5
# ╟─123b165f-0a6c-4894-9157-a538452cf139
# ╟─20df94f0-1d93-4604-9c2c-193856d1f4f9
# ╟─250c95e4-e4c5-4272-b8d8-a1b4fd3d72d1
# ╠═e5c4f17b-da5a-4b72-8de1-6a215804e7b0
# ╟─124f13ff-e0cc-4d85-9854-0b68c1061776
