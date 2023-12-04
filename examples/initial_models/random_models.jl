### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ a0516377-218a-4260-ae15-acf6ac36f2c1
md"# Random Models
The goal is to read an existing grid of average 3D models and interpolate within this grid to obtain new initial conditions."

# ╔═╡ a0e60cf6-7268-11ee-2e25-17dd31433d8f
md"## Setup"

# ╔═╡ 4e200464-64d6-48d2-9c80-b4e91e5b2d3b
md"## Interpolation Grid"

# ╔═╡ 4ce94434-db5a-4323-81e8-c0c5340bab18
grid = MUST.StaggerGrid("stagger_grid_full.mgrid")

# ╔═╡ b3b31521-9967-4ab8-9bf0-3bde45e2db6b
deleteat!(grid.info, .!isfile.(grid["av_path"]))

# ╔═╡ 5701903e-59f3-4495-9dcc-4e730ed2e15f
md"## Get any model within
The simplest way to get any model within this grid, is to interpolate in Teff, logg and FeH. The corresponding 1D initial models can also be interpolated point-wise to the new point in the grid. The more average 3D models available the better. Alternatively, also adiabts can be used for this.
Because this includes MUST and TSO, it is only available as ingredient, to keep both repos separate"

# ╔═╡ c72e1a5a-d13a-481a-b732-3b6be3326326
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ 6e0c166b-58a5-4a25-a726-a66cc357731c
begin
	plot(framestyle=:box, grid=false)
	scatter!(
		grid["logg"][grid["feh"].==0], grid["teff"][grid["feh"].==0], 
		label="FeH=0"
	)
end

# ╔═╡ 6dffcedb-ae0b-4bb4-9490-94c22ec3e953
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

# ╔═╡ d1415b10-403e-4736-9930-14c451f4f366
begin
	paras = random_paramters(grid, 10, teff=[4500, 6500], logg=[4.0, 4.5], feh=[0.0, 0.0])
	#=paras = zeros(1, 3)
	paras[:, 1] .= 5000.0
	paras[:, 2] .= 4.0
	paras[:, 3] .= 0.0=#
end

# ╔═╡ 27f3dc38-c32c-4e35-b5a6-cce2443e64d3
ig = modelgrids.interpolate_from_grid(
	grid, 
	teff=round.(paras[:, 1], sigdigits=4), 
	logg=round.(paras[:, 2], sigdigits=4), 
	feh=paras[:, 3]
)

# ╔═╡ 47f9397e-a5fc-4c7a-a2f6-cf2eb45653e0
md"## Validate interpolation"

# ╔═╡ 213cb3c7-7e08-4929-9c7c-fee75f004ded
begin
	m1 = Average3D(ig["av_path", 1])
	other_models = [Average3D(grid["av_path", i]) for i in 1:nrow(grid.info)]
end

# ╔═╡ 6618d5a5-3c28-4b73-8fd9-a386c4e422d6
begin
	plot(framestyle=:box, grid=false)

	for i in eachindex(other_models)
		if i == 1
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label="grid nodes"
			)
		else
			plot!(
				-other_models[i].z, other_models[i].lnT, linewidth=1, color=:black,
				label=nothing
			)
		end
	end

	plot!(
		-m1.z, m1.lnT, linewidth=7, color=:red, 
		label="T:$(ig["teff", 1]) G:$(ig["logg", 1]) M:$(ig["feh", 1])"
	)

	plot!(xlim=[minimum(-m1.z), maximum(-m1.z)])
	plot!(ylim=[minimum(m1.lnT), maximum(m1.lnT)])
	plot!(xlabel="z", ylabel="log T")
	

	plot!()
end

# ╔═╡ db086ed6-641b-47df-a5df-bcc6df2cbd84
MUST.save(ig, "random_grid.mgrid")

# ╔═╡ Cell order:
# ╟─a0516377-218a-4260-ae15-acf6ac36f2c1
# ╟─a0e60cf6-7268-11ee-2e25-17dd31433d8f
# ╠═f89573bc-8b9b-4897-8ad6-7b10fbdf9b4d
# ╟─4e200464-64d6-48d2-9c80-b4e91e5b2d3b
# ╠═4ce94434-db5a-4323-81e8-c0c5340bab18
# ╠═b3b31521-9967-4ab8-9bf0-3bde45e2db6b
# ╟─5701903e-59f3-4495-9dcc-4e730ed2e15f
# ╠═c72e1a5a-d13a-481a-b732-3b6be3326326
# ╟─6e0c166b-58a5-4a25-a726-a66cc357731c
# ╟─6dffcedb-ae0b-4bb4-9490-94c22ec3e953
# ╠═d1415b10-403e-4736-9930-14c451f4f366
# ╠═27f3dc38-c32c-4e35-b5a6-cce2443e64d3
# ╟─47f9397e-a5fc-4c7a-a2f6-cf2eb45653e0
# ╠═213cb3c7-7e08-4929-9c7c-fee75f004ded
# ╟─6618d5a5-3c28-4b73-8fd9-a386c4e422d6
# ╠═db086ed6-641b-47df-a5df-bcc6df2cbd84
