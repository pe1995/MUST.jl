### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ ab2ddb62-8737-11ee-1799-dd19d7d06046
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PythonPlot

	plt = matplotlib.pyplot
end

# ╔═╡ d962459b-3797-4452-bd50-e2eff8a72064
md"# Modules"

# ╔═╡ 2d624dd2-b15f-4b0b-9879-9f3b01dbf655
@import_dispatch "../../../dispatch2"

# ╔═╡ df5c119e-7841-4e95-b84d-7cb1f06392bb
md"# Test Model"

# ╔═╡ 336c5497-58ce-4b85-9c82-bb67648dd30f
folder = @in_dispatch("data/grid_t5777g44m00")

# ╔═╡ 8d59704b-05ee-496e-9588-24fc3a9ecd0f
md"# Converting function"

# ╔═╡ a3ec9236-0e38-45ca-b163-acfd3ec5e2d2
add_if_not_exists!(arr, data) = begin
	for entry in data
		if !(entry in arr)
			append!(arr, entry)
		end
	end
end

# ╔═╡ d780dcfa-9d72-4a1f-803e-f1bcde4b4262
function build_cube(patch_data)
	# (x, y, z), (li, ui), (patches...)
	patch_range = zeros(Int, 3, 2, length(patch_data))

	global_x = []
	global_y = []
	global_z = []
	for i in eachindex(patch_data)
		add_if_not_exists!(global_x, patch_data[i][:x])
		add_if_not_exists!(global_y, patch_data[i][:y])
		add_if_not_exists!(global_z, patch_data[i][:z])
	end

	global_x .= sort(global_x)
	global_y .= sort(global_y)
	global_z .= sort(global_z)

	p = zeros(3)
	for i in eachindex(patch_data)
		p[1] = first(patch_data[i][:x])
		p[2] = first(patch_data[i][:y])
		p[3] = first(patch_data[i][:z])
		patch_range[:, 1, i] .= MUST._find_in_meshgrid(p, global_x, global_y, global_z)

		p[1] = last(patch_data[i][:x])
		p[2] = last(patch_data[i][:y])
		p[3] = last(patch_data[i][:z])
		patch_range[:, 2, i] .= MUST._find_in_meshgrid(p, global_x, global_y, global_z)
	end
	
	global_x, global_y, global_z, patch_range
end

# ╔═╡ 2fa1d233-a7df-42cf-99f5-14ccd042d7cf
"""
	Box(snap::Py)

Convert a snapshot to a `Box`. Assumes that the patch arangement is cubic, i.e. there is no local mesh refinement. If this is not true, convert to `Space` first, and the interpolate to `Box`.
"""
function Box_fast(snap, quantities...; density=:d, use_mmap=true)
	# first we loop through and get the sizes of all patches
	patch_data = []
	for (i, patch) in enumerate(snap.patches)
		x = MUST.pyconvert.(Float32, patch.xi)
        y = MUST.pyconvert.(Float32, patch.yi) 
        z = MUST.pyconvert.(Float32, patch.zi)

        fname = MUST.pyconvert(Any, patch.filename)
		shp = Tuple(MUST.pyconvert.(Any, patch.ncell))
		off = MUST.pyconvert(Vector{Int}, patch.offset)
		li = MUST.pyconvert(Vector{Int}, patch.li)
		ui = MUST.pyconvert(Vector{Int}, patch.ui)
        idxd = patch.idx.__dict__

		meta = Dict(
			:x=>x,
			:y=>y,
			:z=>z,
			:fname=>fname,
			:shp=>shp,
            :off=>off,
            :li=>li,
            :ui=>ui,
            :idxd=>idxd
		)
		append!(patch_data, [meta])
	end

	x, y, z, patch_range = build_cube(patch_data)
	
	# now we create the data arrays
	data = if use_mmap
		Dict{Symbol,Array{Float32,3}}(
			q=>Array{Float32, 3}(undef, length(x), length(y), length(z)) 
			for q in quantities
		)
	else
		pname = "test.bin" #tempname(pwd())
		io = open(pname, "w+")
		d = MUST.mmap(io, Array{Float32, 4}, (length(x), length(y), length(z), length(quantities)))
		Dict{Symbol,Array{Float32,3}}(
			q=>@view d[:, :, :, i] 
			for (i, q) in enumerate(quantities)
		)
	end

	# and we loop through the patches, read the mmaps, and fill them in the data arrays
	r = zeros(Int, 3, 2)
	for (i, pd) in enumerate(patch_data)
		r .= patch_range[:, :, i]
		for (j, q) in enumerate(quantities)
			data[q][r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] .= MUST._var_from_patch(
				q, 
				pd[:fname], 
				pd[:shp], 
				pd[:off], 
				pd[:li], 
				pd[:ui], 
				pd[:idxd], 
				density=density
			)
		end
	end

	time = MUST.pyconvert(Any, snap.nml_list["snapshot_nml"]["time"])
	xx, yy, zz = MUST.meshgrid(x, y, z)
	MUST.Box(
		xx, yy, zz, data,
		MUST.AtmosphericParameters(
			time, 
			Base.convert(typeof(time), -99.0), 
			Base.convert(typeof(time), -99.0), 
			Dict{Symbol, typeof(time)}()
		)
	)
end

# ╔═╡ 1b517457-f7db-411f-8765-4a600a26c75e
b = Box_fast(dispatch.snapshot(1, data=folder), :d, :e, :ee, use_mmap=false)

# ╔═╡ a4fe72a9-9a5b-42be-8186-c3a418f8d1a9
b2 = Box(Space(dispatch.snapshot(1, data=folder), :d, :e, :ee))

# ╔═╡ ada89a62-d14f-4e0e-bfa7-747785fbf94e
let
	plt.close()
	
	f, ax = plt.subplots(1, 1)

	ax.plot(profile(MUST.mean, b, :z, :log10d)..., label="fast")
	ax.plot(profile(MUST.mean, b2, :z, :log10d)..., label="slow")

	gcf()
end

# ╔═╡ fa6555c6-6a43-470a-b5ad-b443b727190f
let 
	plt.close()
	
	f, ax = plt.subplots(1, 2)

	iclose = MUST.closest(MUST.axis(b, :z), 0.0)
	iclose2 = MUST.closest(MUST.axis(b2, :z), 0.0)
	
	ax[0].imshow(b[:ee][:, :, iclose], label="fast")
	ax[1].imshow(b[:ee][:, :, iclose2], label="slow")

	gcf()
end

# ╔═╡ Cell order:
# ╟─d962459b-3797-4452-bd50-e2eff8a72064
# ╠═ab2ddb62-8737-11ee-1799-dd19d7d06046
# ╠═2d624dd2-b15f-4b0b-9879-9f3b01dbf655
# ╟─df5c119e-7841-4e95-b84d-7cb1f06392bb
# ╠═336c5497-58ce-4b85-9c82-bb67648dd30f
# ╟─8d59704b-05ee-496e-9588-24fc3a9ecd0f
# ╠═a3ec9236-0e38-45ca-b163-acfd3ec5e2d2
# ╠═d780dcfa-9d72-4a1f-803e-f1bcde4b4262
# ╠═2fa1d233-a7df-42cf-99f5-14ccd042d7cf
# ╠═1b517457-f7db-411f-8765-4a600a26c75e
# ╠═a4fe72a9-9a5b-42be-8186-c3a418f8d1a9
# ╟─ada89a62-d14f-4e0e-bfa7-747785fbf94e
# ╟─fa6555c6-6a43-470a-b5ad-b443b727190f
