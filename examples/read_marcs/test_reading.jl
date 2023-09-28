### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 61c9ba44-5d1f-11ee-032a-df5d2ad76a14
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using Revise
	import MUST
end

# ╔═╡ 06a0a5f1-7e90-444c-b83d-918743ff6b5b
model = MUST.MARCSModel("test_model.txt")

# ╔═╡ d9decb50-6341-45e2-8ade-10bb2a776d54
for q in keys(model.structure)
	@info "$(q) $(length(model.structure[q]))"
end

# ╔═╡ a80c13ca-ed98-4eee-9d3d-4338a6f68177
for q in keys(model.pressures)
	@info "$(q) $(length(model.pressures[q]))"
end

# ╔═╡ f4bf9544-0017-4bdb-85b9-558085f41064
model.structure["Density"]

# ╔═╡ 7fdd1ac9-ee66-481d-9946-629f82e2f1b3
model.structure["T"]

# ╔═╡ Cell order:
# ╠═61c9ba44-5d1f-11ee-032a-df5d2ad76a14
# ╠═06a0a5f1-7e90-444c-b83d-918743ff6b5b
# ╠═d9decb50-6341-45e2-8ade-10bb2a776d54
# ╠═a80c13ca-ed98-4eee-9d3d-4338a6f68177
# ╠═f4bf9544-0017-4bdb-85b9-558085f41064
# ╠═7fdd1ac9-ee66-481d-9946-629f82e2f1b3
