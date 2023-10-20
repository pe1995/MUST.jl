### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 136196ea-6f25-11ee-0b3c-ffd40a627be2
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
end

# ╔═╡ 73be8ecf-bb0a-4361-8839-cc0b21df09be
# ╠═╡ show_logs = false
@import_dispatch "../../../dispatch2"

# ╔═╡ 7cf60e1b-ba47-4e56-b346-2cb255496676
folder = @in_dispatch "data/grid_t5777g44m00"

# ╔═╡ 57fbb4cd-2478-4cfa-8ce5-3df0b83aaf79
md"What snapshots are available"

# ╔═╡ b1cd3611-f1e4-4e7d-814a-e4c6656ea6ce
snapshots = list_of_snapshots(folder)

# ╔═╡ cce53fe6-8d5e-49cb-bcce-52877a7d898e
md"Pick a couple (or one specific) to convert them. This takes a while and can, depending on the size of the model, require some memory. This is because at the moment every Patch hat to be loaded into RAM first. When reloading, this is not the case anymore and mmaps are used. This may be optimized in the future."

# ╔═╡ 85ba291b-7456-45b2-a926-09865a77b3fa
b, bτ = snapshotBox(snapshots[end-1], folder=folder)

# ╔═╡ Cell order:
# ╠═136196ea-6f25-11ee-0b3c-ffd40a627be2
# ╠═73be8ecf-bb0a-4361-8839-cc0b21df09be
# ╠═7cf60e1b-ba47-4e56-b346-2cb255496676
# ╟─57fbb4cd-2478-4cfa-8ce5-3df0b83aaf79
# ╠═b1cd3611-f1e4-4e7d-814a-e4c6656ea6ce
# ╟─cce53fe6-8d5e-49cb-bcce-52877a7d898e
# ╠═85ba291b-7456-45b2-a926-09865a77b3fa
