### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ ff1e868e-5941-11ef-1cda-77db6c33abba
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
end

# ╔═╡ 2f3d8a85-a7b5-4f74-9f94-95bbc3d4b03b
iniCond = MUST.ingredients("initial_condition.jl")

# ╔═╡ 33544257-c51a-47a3-9e2e-041a3ff7d42a


# ╔═╡ e493ee91-73cc-4062-bd0d-ed2b88843ee5
md"# Load Grid"

# ╔═╡ 3ed4fb6f-3766-4891-b9a8-c37948bced88
grid = iniCond.staggergrid

# ╔═╡ 7b87b39f-fcdd-473f-9bed-8ddcc5e355db
md"Next, we expand the relative paths within this grid to absolute paths. They are assumed to be replative to the `MUST.jl/src` directory by default."

# ╔═╡ 9988a070-5951-4f7f-bcbf-3cdd53d63d89
MUST.absolute_path!(grid)

# ╔═╡ 07619a18-d647-4f75-8627-5fc8451823f8
grid

# ╔═╡ 976b86b7-cc84-442c-885c-5d65065b674a


# ╔═╡ 3903fea2-3c5c-4b73-a1b6-4a503d9e65d3
md"# New EoS
You can now define any rule you like to pick for each row of the grid a EoS at any path. You need to use the absolute paths here!"

# ╔═╡ 3377dc64-141a-40ce-a9b1-dff24b238507
new_name = "stagger_grid_v5.1"

# ╔═╡ ff93980e-d013-4d14-bc9f-e858ae824cf5
oproot = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/"

# ╔═╡ e6929cbd-0908-469b-ae31-c55b95394dcf
eos_name(ext, version) = joinpath(oproot, "TSO_M3D_$(ext)_$(version)/combined_eos_$(ext).hdf5")

# ╔═╡ fb3b92d2-c689-45b4-a7c5-6ef2a702fbb7
allEoS = Dict(
	0=>eos_name("magg_m0_a0_vmic1", "v5.1"),
	-1=>eos_name("magg_m1_a0_vmic1", "v5.1"),
	-2=>eos_name("magg_m2_a0_vmic1", "v5.1"),
	-3=>eos_name("magg_m3_a0_vmic1", "v5.1"),
	-4=>eos_name("magg_m4_a0_vmic1", "v5.1"),
	-5=>eos_name("magg_m5_a0_vmic1", "v5.1")
)

# ╔═╡ def67aaa-d41d-4a47-a77a-b679ca3d8c65
closest_eos(feh) = begin
	k = sort(keys(allEoS) |> collect)
	ik = k[argmin(abs.(k .- feh))]
	allEoS[ik]
end

# ╔═╡ 0f651973-cc43-4a0f-8f1b-75f7762f1624
new_eos = closest_eos.(grid["feh"])

# ╔═╡ 7d35fd29-caae-45b9-a823-0d4a98010b33
if any(.!isfile.(new_eos))
	@info unique(new_eos[.!isfile.(new_eos)])
end

# ╔═╡ 07d04544-1c25-43a3-841e-6d17b922e373


# ╔═╡ 98fdc9ca-3b95-46a6-98ea-851c4c7bfadd
md"# Replace EoS
Replace the `matching_eos` and `eos_root` columns of the table with the new EoS."

# ╔═╡ 317f87f8-9862-48cb-8924-cb75bb118f02
grid2 = deepcopy(grid)

# ╔═╡ 9a3a65fa-c866-4d44-a613-748f8c5e58be
grid2.info[!, "eos_root"] = dirname.(new_eos)

# ╔═╡ 4d09e9ed-885c-4e59-9bd3-461c58932872
grid2.info[!, "matching_eos"] = new_eos

# ╔═╡ 7462b440-a580-4a58-adee-a170e1788d44
#grid2.info[!, "av_path"] = joinpath.(abspath(pwd(), "MARCS/av_models"), basename.(grid2.info[!, "av_path"]))

# ╔═╡ 21242730-5fd7-4c73-8951-deaee78fd287
#grid2.info[!, "avo_path"] = joinpath.(abspath(pwd(), "MARCS/av_models"), basename.(grid2.info[!, "avo_path"]))

# ╔═╡ 3c03580b-d97f-459f-ae72-694f7d2c0453
#grid2.info[!, "abs_av_path"] = joinpath.(abspath(pwd(), "MARCS/av_models"), basename.(grid2.info[!, "av_path"]))

# ╔═╡ 9c6f4603-1a89-4614-ba23-07a52f4996bc
MUST.relative_path!(grid2)

# ╔═╡ 88489ab8-44d8-4726-baa6-a52fcf10a3ae
path_new = joinpath(dirname(grid2.path), new_name*".mgrid")

# ╔═╡ 15cc88df-ed46-47a5-901f-fa9cdade8072
grid3 = MUST.Atmos1DGrid(new_name, path_new, grid2.info)

# ╔═╡ a5f3ba9b-7b19-4dbc-bfd1-6a4b983da19e
MUST.save(grid3)

# ╔═╡ Cell order:
# ╠═ff1e868e-5941-11ef-1cda-77db6c33abba
# ╠═2f3d8a85-a7b5-4f74-9f94-95bbc3d4b03b
# ╟─33544257-c51a-47a3-9e2e-041a3ff7d42a
# ╟─e493ee91-73cc-4062-bd0d-ed2b88843ee5
# ╠═3ed4fb6f-3766-4891-b9a8-c37948bced88
# ╟─7b87b39f-fcdd-473f-9bed-8ddcc5e355db
# ╠═9988a070-5951-4f7f-bcbf-3cdd53d63d89
# ╠═07619a18-d647-4f75-8627-5fc8451823f8
# ╟─976b86b7-cc84-442c-885c-5d65065b674a
# ╟─3903fea2-3c5c-4b73-a1b6-4a503d9e65d3
# ╠═3377dc64-141a-40ce-a9b1-dff24b238507
# ╠═ff93980e-d013-4d14-bc9f-e858ae824cf5
# ╠═e6929cbd-0908-469b-ae31-c55b95394dcf
# ╠═fb3b92d2-c689-45b4-a7c5-6ef2a702fbb7
# ╠═def67aaa-d41d-4a47-a77a-b679ca3d8c65
# ╠═0f651973-cc43-4a0f-8f1b-75f7762f1624
# ╠═7d35fd29-caae-45b9-a823-0d4a98010b33
# ╟─07d04544-1c25-43a3-841e-6d17b922e373
# ╟─98fdc9ca-3b95-46a6-98ea-851c4c7bfadd
# ╠═317f87f8-9862-48cb-8924-cb75bb118f02
# ╠═9a3a65fa-c866-4d44-a613-748f8c5e58be
# ╠═4d09e9ed-885c-4e59-9bd3-461c58932872
# ╠═7462b440-a580-4a58-adee-a170e1788d44
# ╠═21242730-5fd7-4c73-8951-deaee78fd287
# ╠═3c03580b-d97f-459f-ae72-694f7d2c0453
# ╠═9c6f4603-1a89-4614-ba23-07a52f4996bc
# ╠═88489ab8-44d8-4726-baa6-a52fcf10a3ae
# ╠═15cc88df-ed46-47a5-901f-fa9cdade8072
# ╠═a5f3ba9b-7b19-4dbc-bfd1-6a4b983da19e
