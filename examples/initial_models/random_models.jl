### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ f89573bc-8b9b-4897-8ad6-7b10fbdf9b4d
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

# ╔═╡ d3b3d5a1-82ac-494f-ad15-684be62555d3
using LaTeXStrings

# ╔═╡ a0516377-218a-4260-ae15-acf6ac36f2c1
md"# Random Models
The goal is to read an existing grid of average 3D models and interpolate within this grid to obtain new initial conditions."

# ╔═╡ a0e60cf6-7268-11ee-2e25-17dd31433d8f
md"## Setup"

# ╔═╡ fb5cf860-1fbd-472a-acf1-ad66ea43cddb
begin
	gr()
	default(
		framestyle=:box, grid=false, minorticks=true, 
		legendforegroundcolor=nothing,legendbackgroundcolor=nothing
	)
end

# ╔═╡ 4e200464-64d6-48d2-9c80-b4e91e5b2d3b
md"## Interpolation Grid"

# ╔═╡ 4ce94434-db5a-4323-81e8-c0c5340bab18
grid = MUST.StaggerGrid("stagger_grid_full.mgrid")

# ╔═╡ 5165efd2-ecbd-486d-acba-a11e555184d8
#grid = MUST.StaggerGrid("MARCS_z0/marcs_grid_z0.mgrid")

# ╔═╡ b3b31521-9967-4ab8-9bf0-3bde45e2db6b
deleteat!(grid.info, .!isfile.(grid["av_path"]));

# ╔═╡ d3dd31b5-237e-4b51-9706-8b33269ce5dd
deleteat!(grid.info, grid["logg"] .<= 2.0);

# ╔═╡ f5ff3e9c-44ca-4442-b587-064272205dcc
#deleteat!(grid.info, grid["feh"] .<= -0.1);

# ╔═╡ 5701903e-59f3-4495-9dcc-4e730ed2e15f
md"## Get any model within
The simplest way to get any model within this grid, is to interpolate in Teff, logg and FeH. The corresponding 1D initial models can also be interpolated point-wise to the new point in the grid. The more average 3D models available the better. Alternatively, also adiabts can be used for this.
Because this includes MUST and TSO, it is only available as ingredient, to keep both repos separate"

# ╔═╡ c72e1a5a-d13a-481a-b732-3b6be3326326
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ 6dffcedb-ae0b-4bb4-9490-94c22ec3e953
begin
	function random_paramters(grid, N; 
		teff=[minimum(grid["teff"]), maximum(grid["teff"])], 
		logg=[minimum(grid["logg"]), maximum(grid["logg"])], 
		feh=[minimum(grid["feh"]), maximum(grid["feh"])])
		
		points = [[t, l, f] for (t, l, f) in zip(grid["teff"], grid["logg"], grid["feh"])]
		hull = convex_hull(points)
		P = VPolytope(hull)
		
		rand_points = zeros(N, 3)
		found = 0
		while found < N
			t = MUST.randrange(teff...)
			l = MUST.randrange(logg...) 
			f = MUST.randrange(feh...)
		
			if [t, l, f] ∈ P
				found += 1
				rand_points[found, 1] = t
				rand_points[found, 2] = l
				rand_points[found, 3] = f
			end
		end
	
		rand_points
	end
	
	function random_paramters(grid, N, lowerlimit, upperlimit; 
		teff=[minimum(grid["teff"]), maximum(grid["teff"])], 
		logg=[minimum(grid["logg"]), maximum(grid["logg"])], 
		feh=[minimum(grid["feh"]), maximum(grid["feh"])])
		
		points = [[t, l, f] for (t, l, f) in zip(grid["teff"], grid["logg"], grid["feh"])]
		hull = convex_hull(points)
		P = VPolytope(hull)
		
		rand_points = zeros(N, 3)
		found = 0
		while found < N
			t = MUST.randrange(teff...)
			l = MUST.randrange(logg...) 
			f = MUST.randrange(feh...)
		
			if [t, l, f] ∈ P
				# check if it is below the lower or above the upper limit
				lgLowLim = min(lowerlimit(t), upperlimit(t))
				lgUpLim = max(lowerlimit(t), upperlimit(t))
	
				if (l>=lgLowLim) & (l<=lgUpLim)
					found += 1
					rand_points[found, 1] = t
					rand_points[found, 2] = l
					rand_points[found, 3] = f
				end
			end
		end
	
		rand_points
	end
end

# ╔═╡ 77758e3f-5d6a-4b4d-9239-0b575d78ef76
begin
	lineLimit(teff, x1, x2) = begin
			m = (x2[2] - x1[2]) / (x2[1] - x1[1]) 
			m* teff + x2[2] - m*x2[1]
		end
	
	lowerLimit(x) = lineLimit(x, (5000.0, 4.5), (7000.0, 4.3))
	upperLimit(x) = lineLimit(x, (5000.0, 3.6), (7000.0, 3.4))
end

# ╔═╡ d1415b10-403e-4736-9930-14c451f4f366
begin
	
	#=paras = random_paramters(
		grid, 
		60, 
		lowerLimit,
		upperLimit,
		teff=[4500, 7000], 
		logg=[3.0, 4.7], 
		feh=[0.0, 0.0]
	)=#
	

	# or reload parameters already sampled
	
	#=parasList = MUST.StaggerGrid("CEMP/random_CEMP_2.mgrid")
	paras = zeros(nrow(parasList.info), 3)
	paras[:, 1] = parasList["teff"]
	paras[:, 2] = parasList["logg"]
	paras[:, 3] = parasList["feh"]=#

	#=parasList = MUST.StaggerGrid("random_MS_5.mgrid")
	paras = zeros(nrow(parasList.info), 3)
	paras[:, 1] = parasList["teff"]
	paras[:, 2] = parasList["logg"]
	paras[:, 3] = parasList["feh"]=#

	
	# select a few sub-giants for CEMP tests
	
	#=paras = random_paramters(
		grid, 
		10, 
		teff=[4500, 5500], 
		logg=[3.0, 4.0], 
		feh=[-4.0, -4.0]
	)=#


	# select a few special stars
	
	#=paras = zeros(2, 3)
	paras[:, 1] = [6250.0, 5250.0]
	paras[:, 2] = [4.0, 4.0]
	paras[:, 3] = [0.0, 0.0]=#


	# test case for MARCS
	
	paras = random_paramters(
		grid, 
		5, 
		teff=[4500, 7000], 
		logg=[3.0, 4.5], 
		feh=[-4.0, -4.0]
	)
	
end

# ╔═╡ 0bffe69d-50fa-4b1c-adba-de962e8f38cc
begin
	#= add hot Jupyter testcase
	paras_extended = zeros(eltype(paras), size(paras, 1)+1, 3)
	paras_extended[1:end-1, :] .= paras
	paras_extended[end, :] = [4500.0, 4.8, 0.0]=#

	# add solar test case
	paras_extended = zeros(eltype(paras), size(paras, 1)+1, 3)
	paras_extended[1:end-1, :] .= paras
	paras_extended[end, :] = [5777.0, 4.44, -4.0]
	
	#paras_extended = paras
end

# ╔═╡ 5fbd389f-9ea5-4b09-839b-59aad0c72a0d
begin
	plot(xflip=true, yflip=true, legend=:topleft)
	
	scatter!(
		paras_extended[:, 1],
		paras_extended[:, 2],
		color="red",
		label="New models",
		marker=:square,
		markersize=6
	)
	scatter!(
		grid["teff"][grid["feh"].==0], 
		grid["logg"][grid["feh"].==0], 
		label="Grid",
		color="black",
		markersize=8
	)

	x = range(3000, 8000, length=100)
	#plot!(x, lowerLimit.(x), color=:blue, lw=3, label=nothing)
	#plot!(x, upperLimit.(x), color=:blue, lw=3, label=nothing)
	plot!(ylabel="log g")
	plot!(xlabel="T"*L"\textrm{\textbf{_{eff}}}")
end

# ╔═╡ ed48ac99-eb46-4386-b1a5-31f3e78c2574


# ╔═╡ 8ab3ea9f-8c5e-4ff2-8b5b-631beb9f6438
md"## EoS
We can add information about the EoS. This will make the interpolation better because it will use the rosseland optical depth information to interpolate at constant optical depth between models. If not, just position is used for the interpolation. Depending on the chemical composition of the model you need to pick a different EoS here! So it might be best to include the mother table of the respective model in the grid, so that for each metallicity the correct one can be used."

# ╔═╡ 78e843fc-42c1-49ba-8d8b-0a8ada2b2a3a
#extension="magg_m4_a4_c3_vmic2"

# ╔═╡ 6dd1a861-7654-4f5a-b871-3f2526154d26
extension="magg_m0_a0"

# ╔═╡ cb7d1e6e-5462-44e9-88d0-96a27a85990c
#extension="magg_m10_vmic2"

# ╔═╡ e2598dbe-f097-428f-987c-5c36f62845e7
extension="magg_m10_a4_c3_vmic2"

# ╔═╡ 1ebff512-8f72-48e0-9cb6-c97f23073893


# ╔═╡ e17316a1-a29e-4ce7-b651-65a75881180f
#mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_$(extension)_v1.8"

# ╔═╡ 1191c34a-c27e-4eff-875a-57f7fc043ee9
mother_table_path = "../../../opacity_tables/TSO_MARCS_$(extension)_v1.8"

# ╔═╡ 5e777d3b-2db8-4e82-a125-564c7db64281
#=mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m1_a0_c1_v3.0"=#

# ╔═╡ ad7d6660-fea7-4808-bff3-7288c20966b2
#mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(extension)_v4.0"

# ╔═╡ 453fc753-4298-4909-97fa-5624783d9ef6
mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(extension)_v5.0"

# ╔═╡ 775aff22-9bb6-4171-aa9d-3b8daff4c248


# ╔═╡ 1ca259f3-54b6-4c96-8930-eab6454a2cd2
eos_mother_path = "ross_combined_eos_$(extension).hdf5"

# ╔═╡ b40c3cee-303b-4a24-9692-fc92a69a9b50
#eos_mother_path = "combined_eos_magg_m1_a0_c1.hdf5"

# ╔═╡ 03cde4ab-bd80-4476-a27a-750a13935a96
#eos_mother_path = "combined_eos_$(extension).hdf5"

# ╔═╡ af910237-0a7d-4b58-a870-1fbcfd5e883e


# ╔═╡ 2f101c34-f3f7-46d8-b1e3-fb855281c6ac
eos = [
	reload(SqEoS, joinpath(mother_table_path, eos_mother_path))
	for _ in 1:size(paras_extended, 1)
]

# ╔═╡ 2b9e746b-94c9-4a5a-9838-9cccaf96e406


# ╔═╡ 5bc6df18-c39d-4cc3-9010-7de7b986935d
let
	plot()
	for i in 1:nrow(grid.info)
		m = TSO.Average3D(grid.info[i, "av_path"]; logg=grid.info[i, "logg"]) 
		m = @optical m eos[1]
		znew = TSO.rosseland_depth(eos[1], m)
		m.z .= znew
		TSO.optical_surface!(m)
		TSO.flip!(m)
		
		plot!(m.z, m.lnT, label=nothing)
	end

	plot!()
end

# ╔═╡ 52bdf098-3209-4703-b74d-3de5886ad8fe


# ╔═╡ ee45e7db-2fc3-43f3-83cf-bd2dca6d52e7
md"## Interpolation"

# ╔═╡ 0364f91e-e768-4126-b347-ec6684b9c033
md"We can add a matching eos for each of the grid nodes. In the case of non-existing optical depth averaged models, this will be used to determine the optical dept."

# ╔═╡ 2dee306a-531f-4eb2-bac7-e4a5c8a219a7
if !("matching_eos" in names(grid.info))
	@info "Add matching EoS."
	
	grid.info[!, "matching_eos"] = [
		joinpath(mother_table_path, eos_mother_path) 
		for _ in 1:nrow(grid.info)
	]
else
	@info "Matching EoS already present."
end

# ╔═╡ 4aedcbe9-ee05-48d7-8919-51514b09f018
begin
	newModelCollection = "StaggerExtrapolated"
	newdir = "StaggerExtrapolated/av_$(extension)/"
	newSet = "SE2"
end

# ╔═╡ 268b93aa-b610-467e-be1c-0306ca8ac31e
begin
	if !isdir(newModelCollection)
		mkdir(newModelCollection)
	end
	if !isdir(newdir)
		mkdir(newdir)
	end
end;

# ╔═╡ 27f3dc38-c32c-4e35-b5a6-cce2443e64d3
ig = modelgrids.interpolate_from_grid(
	grid, 
	teff=round.(paras_extended[:, 1], sigdigits=4), 
	logg=round.(paras_extended[:, 2], sigdigits=4), 
	feh=paras_extended[:, 3],
	eos=eos,
	folder=newdir,
	adiabatic_extrapolation=true,
	τbottom=2
)

# ╔═╡ a47fa3e9-bd97-45e4-bf55-8de2e0d75c64
ig.info[!, "eos_root"] = [mother_table_path for _ in 1:size(paras_extended, 1)]

# ╔═╡ 47f9397e-a5fc-4c7a-a2f6-cf2eb45653e0
md"## Validate interpolation"

# ╔═╡ 2536778d-66e7-4c08-a5a4-541f4ad46795
i_select = 1

# ╔═╡ 213cb3c7-7e08-4929-9c7c-fee75f004ded
begin
	m1 = @optical Average3D(ig["av_path", i_select]) eos[i_select]
	other_models = [Average3D(grid["av_path", i]) for i in 1:nrow(grid.info)]
end

# ╔═╡ 6618d5a5-3c28-4b73-8fd9-a386c4e422d6
begin
	plot(framestyle=:box, grid=false)

	plot!(
		-m1.z, m1.lnT, linewidth=3, color=:red, 
		label="T:$(ig["teff", i_select]) G:$(ig["logg", i_select]) M:$(ig["feh", i_select])"
	)

	for i in eachindex(other_models)
		#=if i == 1
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label="grid nodes"
			)
		else
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label=nothing
			)
		end=#
		g = grid
		t = paras[i_select, 1]
		mask = (g["teff", i] == t+250) && (g["logg", i] == 4.0) && (g["feh", i] == 0.0)
		mask2 = (g["teff", i] == t-250) && (g["logg", i] == 4.0) && (g["feh", i] == 0.0)
		
		mask = mask | mask2
		if mask
			m = @optical other_models[i] eos[i_select]
			znew = TSO.rosseland_depth(eos[i_select], m)
			m = deepcopy(m)
			m.z .= znew
			TSO.optical_surface!(m)
			plot!(
				-m.z, m.lnT, linewidth=1, color=:black,
				label="T:$(g["teff", i]) G:$(g["logg", i]) M:$(g["feh", i])"
			)
		end
	end

	

	#plot!(xlim=[minimum(-m1.z), maximum(-m1.z)])
	#plot!(ylim=[minimum(m1.lnT), maximum(m1.lnT)])
	plot!(xlabel="z", ylabel="log T")
	
	
	plot!()
end

# ╔═╡ dbee67e3-ee0a-4954-9ead-033141c8bf2f
begin
	plot(framestyle=:box, grid=false)

	

	#plot!(xlim=[minimum(-m1.z), maximum(-m1.z)])
	#plot!(ylim=[minimum(m1.lnT), maximum(m1.lnT)])
	plot!(xlabel="z", ylabel="log T")
	
	
	plot!()
end

# ╔═╡ dbee67e3-ee0a-4954-9ead-033141c8bf2f
begin
	plot(framestyle=:box, grid=false)

	plot!(
		log10.(m1.τ), m1.lnT, linewidth=3, color=:red, 
		label="T:$(ig["teff", i_select]) G:$(ig["logg", i_select]) M:$(ig["feh", i_select])"
	)

	for i in eachindex(other_models)
		#=if i == 1
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label="grid nodes"
			)
		else
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label=nothing
			)
		end=#
		g = grid
		t = paras[i_select, 1]
		mask = (g["teff", i] == t+250) && (g["logg", i] == 4.0) && (g["feh", i] == 0.0)
		mask2 = (g["teff", i] == t-250) && (g["logg", i] == 4.0) && (g["feh", i] == 0.0)
		
		mask = mask | mask2
		if mask
			ml = @optical other_models[i] eos[i_select]
			plot!(
				log10.(ml.τ), ml.lnT, linewidth=1, color=:black,
				label="T:$(g["teff", i]) G:$(g["logg", i]) M:$(g["feh", i])"
			)
		end
	end



	#plot!(xlim=[minimum(-m1.z), maximum(-m1.z)])
	#plot!(ylim=[minimum(m1.lnT), maximum(m1.lnT)])
	plot!(xlabel="z", ylabel="log T")
	
	plot!()
end

# ╔═╡ ba92c5a2-fc64-431f-bb0f-0e8af8437b0f
begin
	m1t = @optical Average3D(eos[1], ig["av_path", 1]) eos[1]
	m2t = @optical Average3D(eos[2], ig["av_path", 2]) eos[2]
end

# ╔═╡ 1d80590e-fca6-4200-a9e2-c6e2299330a2
let
	plot(framestyle=:box, grid=false)

	plot!(log10.(m1t.τ), exp.(m1t.lnT), label="interpolated")

	Tselect = [6500.0, 6000.0, 5777.]
	loggselect = [4.0, 4.0, 4.44]
	for (t, l) in zip(Tselect, loggselect)
		mask = (grid["teff", !] .== t) .& (grid["logg", !] .== l) .& (grid["feh", !] .== 0.0) 
		i = findfirst(mask)
		mi = @optical other_models[i] eos[1]
		plot!(log10.(mi.τ), exp.(mi.lnT), label="T:$(t), logg:$(l)")
	end
	
	plot!(xlim=(-6, 6.5), ylim=(1000,40000))
end

# ╔═╡ d9e75d15-d74c-4d7a-8f8b-c71380d3bd44
let 
	mask = (grid["teff", !] .== 5777.0) .& (grid["logg", !] .== 4.44) .& (grid["feh", !] .== 0.0) 
	i = findfirst(mask)
	mi = @optical other_models[i] eos[1]
	writedlm("test_sun.csv", [log10.(mi.τ) exp.(mi.lnT) mi.lnρ mi.z])
end

# ╔═╡ db086ed6-641b-47df-a5df-bcc6df2cbd84
MUST.save(ig, joinpath(newModelCollection, "$(newSet)_$(extension).mgrid"))

# ╔═╡ Cell order:
# ╟─a0516377-218a-4260-ae15-acf6ac36f2c1
# ╟─a0e60cf6-7268-11ee-2e25-17dd31433d8f
# ╠═f89573bc-8b9b-4897-8ad6-7b10fbdf9b4d
# ╠═d3b3d5a1-82ac-494f-ad15-684be62555d3
# ╠═fb5cf860-1fbd-472a-acf1-ad66ea43cddb
# ╟─4e200464-64d6-48d2-9c80-b4e91e5b2d3b
# ╠═4ce94434-db5a-4323-81e8-c0c5340bab18
# ╠═5165efd2-ecbd-486d-acba-a11e555184d8
# ╠═b3b31521-9967-4ab8-9bf0-3bde45e2db6b
# ╠═d3dd31b5-237e-4b51-9706-8b33269ce5dd
# ╠═f5ff3e9c-44ca-4442-b587-064272205dcc
# ╟─5701903e-59f3-4495-9dcc-4e730ed2e15f
# ╠═c72e1a5a-d13a-481a-b732-3b6be3326326
# ╟─6dffcedb-ae0b-4bb4-9490-94c22ec3e953
# ╠═77758e3f-5d6a-4b4d-9239-0b575d78ef76
# ╠═d1415b10-403e-4736-9930-14c451f4f366
# ╠═0bffe69d-50fa-4b1c-adba-de962e8f38cc
# ╟─5fbd389f-9ea5-4b09-839b-59aad0c72a0d
# ╟─ed48ac99-eb46-4386-b1a5-31f3e78c2574
# ╟─8ab3ea9f-8c5e-4ff2-8b5b-631beb9f6438
# ╠═78e843fc-42c1-49ba-8d8b-0a8ada2b2a3a
# ╠═6dd1a861-7654-4f5a-b871-3f2526154d26
# ╠═cb7d1e6e-5462-44e9-88d0-96a27a85990c
# ╠═e2598dbe-f097-428f-987c-5c36f62845e7
# ╟─1ebff512-8f72-48e0-9cb6-c97f23073893
# ╠═e17316a1-a29e-4ce7-b651-65a75881180f
# ╠═1191c34a-c27e-4eff-875a-57f7fc043ee9
# ╠═5e777d3b-2db8-4e82-a125-564c7db64281
# ╠═ad7d6660-fea7-4808-bff3-7288c20966b2
# ╠═453fc753-4298-4909-97fa-5624783d9ef6
# ╟─775aff22-9bb6-4171-aa9d-3b8daff4c248
# ╠═1ca259f3-54b6-4c96-8930-eab6454a2cd2
# ╠═b40c3cee-303b-4a24-9692-fc92a69a9b50
# ╠═03cde4ab-bd80-4476-a27a-750a13935a96
# ╟─af910237-0a7d-4b58-a870-1fbcfd5e883e
# ╟─2f101c34-f3f7-46d8-b1e3-fb855281c6ac
# ╟─2b9e746b-94c9-4a5a-9838-9cccaf96e406
# ╟─5bc6df18-c39d-4cc3-9010-7de7b986935d
# ╟─52bdf098-3209-4703-b74d-3de5886ad8fe
# ╟─ee45e7db-2fc3-43f3-83cf-bd2dca6d52e7
# ╟─0364f91e-e768-4126-b347-ec6684b9c033
# ╟─2dee306a-531f-4eb2-bac7-e4a5c8a219a7
# ╠═4aedcbe9-ee05-48d7-8919-51514b09f018
# ╟─268b93aa-b610-467e-be1c-0306ca8ac31e
# ╠═27f3dc38-c32c-4e35-b5a6-cce2443e64d3
# ╠═a47fa3e9-bd97-45e4-bf55-8de2e0d75c64
# ╟─47f9397e-a5fc-4c7a-a2f6-cf2eb45653e0
# ╠═2536778d-66e7-4c08-a5a4-541f4ad46795
# ╟─213cb3c7-7e08-4929-9c7c-fee75f004ded
# ╟─6618d5a5-3c28-4b73-8fd9-a386c4e422d6
# ╟─dbee67e3-ee0a-4954-9ead-033141c8bf2f
# ╠═db086ed6-641b-47df-a5df-bcc6df2cbd84
