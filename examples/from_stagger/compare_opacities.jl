### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ fe22481c-f644-11ed-1b92-dd3e526e7845
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using Plots
end

# ╔═╡ 201b667f-a7f1-457d-bf91-3f41e1943589
MUST.@import_dispatch "/home/eitner/shared/model_grid/dispatch2"

# ╔═╡ ef5e0076-3c49-4dc9-9ecf-56ce63b011b9
tables = [
	MUST.@in_dispatch("input_data/grd/DIS_MARCS_E_t5777g44m00_v0.1"),
	MUST.@in_dispatch("input_data/DIS_MARCS_E_v0.5.1")
]

# ╔═╡ 47692886-50b1-44d3-99d5-eda5afae147d
eos = [reload(SqEoS, joinpath(tables[i], "eos.hdf5")) for i in eachindex(tables)]

# ╔═╡ 19f00ac1-5290-47b0-afd4-818a2d27ab3e
opa = [
	reload(SqOpacity, joinpath(tables[i], "binned_opacities.hdf5")) 
			for i in eachindex(tables)
]

# ╔═╡ 8893ea31-175c-4a39-b3a0-151debf17058
# ╠═╡ show_logs = false
models = [
	Average3D(eos[i], joinpath(tables[1], "inim.dat")) 
		for i in eachindex(tables)
]

# ╔═╡ Cell order:
# ╠═fe22481c-f644-11ed-1b92-dd3e526e7845
# ╠═201b667f-a7f1-457d-bf91-3f41e1943589
# ╠═ef5e0076-3c49-4dc9-9ecf-56ce63b011b9
# ╠═47692886-50b1-44d3-99d5-eda5afae147d
# ╠═19f00ac1-5290-47b0-afd4-818a2d27ab3e
# ╠═8893ea31-175c-4a39-b3a0-151debf17058
