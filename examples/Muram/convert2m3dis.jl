### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ f617333a-53be-11ee-38b0-8dc017f5e8b8
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	using MUST
	using TSO
	using Glob
end

# ╔═╡ f2ecb12c-e939-4bbd-9e57-37761de41d78
eos = reload(
	SqEoS, "/u/peitner/DISPATCH/TSO.jl/examples/creating_table/eos_asplund07_v5.0.hdf5"
)

# ╔═╡ c6388487-c6c9-4283-9ec6-8475e27c9a22
must2multi = MUST.ingredients("convert2multi.jl")

# ╔═╡ a21d5e49-c1b2-40c5-a0d1-8b7aa1bceaa9
paths = glob("*", "feh2")

# ╔═╡ 9cebddeb-42bd-4cda-946f-fafbc92fb1f8
models = Box.(paths)

# ╔═╡ 53537391-1fd9-4ffb-be4f-46fdee9d7110
labels = last.(split.(first.(split.(last.(split.(paths, "/")), ".hdf5")), "_"))

# ╔═╡ 64b1fbb2-730b-40b3-b95e-4bebd21a4b22
for model in models
	if (model[:d][1, 1, 1] < model[:d][1, 1, end])
		@info "cube has lower density at the bottom than at the top. Reversing..."
		for k in keys(model.data)
			model.data[k] .= reverse(model.data[k], dims=3)
		end
		model.x .= reverse(model.x, dims=3)
		model.y .= reverse(model.y, dims=3)
		model.z .= reverse(model.z, dims=3)
	end
	
	if (model[:d][1, 1, 1] > model[:d][1, 1, end]) & 
		(model.z[1, 1, 1] > model.z[1, 1, end])
		@info "Switching sign of zscale to be negative in the interior..."
		model.z = -model.z
	end
end

# ╔═╡ d1ad1725-2819-44ab-8825-8c23bae1623c
begin
	example_model = models[1]
	example_model.z = - example_model.z
	n = gresample(example_model; nz=200)
end

# ╔═╡ 27aa01be-cdf9-496a-8ecf-3266f558c0c8
#xres = 50

# ╔═╡ 6ec71eaf-5d10-4dee-add8-ad633a5f78ff
#=must2multi.snaps2multi(
	models..., 
	eos=eos, 
	label=["$(l)_$(xres)x$(xres)" for (i, l) in enumerate(labels)],
	name="muram_m2",
	n_horizontal=xres,
	n_vertical=200,
	outfolder="muram_m2"
)=#

# ╔═╡ Cell order:
# ╠═f617333a-53be-11ee-38b0-8dc017f5e8b8
# ╠═f2ecb12c-e939-4bbd-9e57-37761de41d78
# ╠═c6388487-c6c9-4283-9ec6-8475e27c9a22
# ╠═a21d5e49-c1b2-40c5-a0d1-8b7aa1bceaa9
# ╠═9cebddeb-42bd-4cda-946f-fafbc92fb1f8
# ╠═53537391-1fd9-4ffb-be4f-46fdee9d7110
# ╠═64b1fbb2-730b-40b3-b95e-4bebd21a4b22
# ╠═d1ad1725-2819-44ab-8825-8c23bae1623c
# ╠═27aa01be-cdf9-496a-8ecf-3266f558c0c8
# ╠═6ec71eaf-5d10-4dee-add8-ad633a5f78ff
