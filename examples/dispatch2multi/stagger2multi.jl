### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6d7bb008-84da-11ee-0801-9fe63eee8429
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PlutoUI
end

# ╔═╡ 90b5f877-f02b-43cb-94a0-9e833218e138
md"# Setup"

# ╔═╡ 5755ef9a-4478-4ab7-96b0-d057092b1142
TableOfContents()

# ╔═╡ cd117e11-2412-42c3-8cf3-d4f8a7e079fb
all_models_at(path) = split.(MUST.glob("t*/", path), "/", keepempty=false) .|> last

# ╔═╡ f2e2f930-ff80-485f-b24a-c306c718962e
all_snapshots(path) = split.(basename.(MUST.glob("*.aux", path)), ".aux") .|> first

# ╔═╡ f76fa48c-c41e-40fe-bc4b-d92ae5b1fb6e
md"# Stagger cubes
## Pick Snapshot"

# ╔═╡ a043761b-5033-4a49-a549-40edb3ac0fcb
md"Pick location of Stagger grid: $(@bind datafolder TextField(default=\"/mnt/beegfs/gemini/groups/bergemann/users/shared-storage/bergemann-data/Stagger_zaz/grid/\"))"

# ╔═╡ f3e0f3e6-fab5-4fe8-9085-263b1f251818
md"Select from available models: $(@bind stagger_model Select(all_models_at(datafolder)))"

# ╔═╡ dfa4df7d-136a-451f-9c81-943ad5ca4b16
md"Select from available snapshots: $(@bind snapshot_name Select(all_snapshots(joinpath(datafolder, stagger_model))))"

# ╔═╡ 7f706c09-e8ed-4e7b-a2d8-9ad7d1d169ab


# ╔═╡ ebf5816b-60c1-424c-b6da-830b605cceb1
md"## Convert to Box
By defaul, we do not interpolate it to a uniform grid in the vertical with a low resulution. To pick the current number of points instead we have to pick `max_resolution=false`"

# ╔═╡ f35b1eaa-9876-4e84-8bbb-63154860b380
md"Tick box to start converting snapshots: $(@bind start_conversion CheckBox(default=false))"

# ╔═╡ d7300bbe-51eb-462a-9b30-d20b85e73452
if start_conversion
	box = Box(
		MUST.StaggerSnap(snapshot_name, joinpath(datafolder, stagger_model)),
		skip_interpolation=false, 
		max_resolution=false
	)
end

# ╔═╡ 81ac358d-ca8c-46ca-b34b-ddfa58b54103


# ╔═╡ 277a1a92-df1d-49f4-9c72-47e231ca52ee
md"## Convert to M3D
We can now convert this to M3D format. It can be downsampled directly."

# ╔═╡ 4e450df8-b244-4bd9-958d-757668eaa9e6
must2multi = MUST.ingredients("convert2multi.jl")

# ╔═╡ 86a90fa6-efd6-4ecb-b361-0cf3988c318d
if start_conversion
	must2multi.snaps2multi(
		box,
	    label=snapshot_name, 
		n_horizontal=10,
		outfolder="", 
		name="m3dis_stagger"
	)
end

# ╔═╡ Cell order:
# ╟─90b5f877-f02b-43cb-94a0-9e833218e138
# ╠═6d7bb008-84da-11ee-0801-9fe63eee8429
# ╟─5755ef9a-4478-4ab7-96b0-d057092b1142
# ╟─cd117e11-2412-42c3-8cf3-d4f8a7e079fb
# ╟─f2e2f930-ff80-485f-b24a-c306c718962e
# ╟─f76fa48c-c41e-40fe-bc4b-d92ae5b1fb6e
# ╟─a043761b-5033-4a49-a549-40edb3ac0fcb
# ╟─f3e0f3e6-fab5-4fe8-9085-263b1f251818
# ╟─dfa4df7d-136a-451f-9c81-943ad5ca4b16
# ╟─7f706c09-e8ed-4e7b-a2d8-9ad7d1d169ab
# ╟─ebf5816b-60c1-424c-b6da-830b605cceb1
# ╟─f35b1eaa-9876-4e84-8bbb-63154860b380
# ╠═d7300bbe-51eb-462a-9b30-d20b85e73452
# ╟─81ac358d-ca8c-46ca-b34b-ddfa58b54103
# ╟─277a1a92-df1d-49f4-9c72-47e231ca52ee
# ╟─4e450df8-b244-4bd9-958d-757668eaa9e6
# ╠═86a90fa6-efd6-4ecb-b361-0cf3988c318d
