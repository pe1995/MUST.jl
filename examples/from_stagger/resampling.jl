### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 13e36110-0f52-11ee-0142-bd1afe3d9a36
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using Revise
	import MUST
	using Plots
end

# ╔═╡ 82119670-88be-4e8e-88e8-8c8c2e706b4f
md"# 3D-Cube Interpolation"

# ╔═╡ 72bf89a2-6ef5-4112-abe1-e7a8c26897f7
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 12883f0f-d6f6-40c8-b389-8a18bdc1606b
model_location = MUST.@in_dispatch "data/sun_hllc"

# ╔═╡ ce272b96-8e58-41d9-852c-44971ce370b7
snap, snap_τ = MUST.pick_snapshot(model_location, :recent)

# ╔═╡ 04fb3d30-9bd5-47f9-b58f-f0a8124b7593
md"## Current Grid
We can construct a `RegularGrid` object from a `Box` object by converting its x, y and z axis into `RegularAxis` objects."

# ╔═╡ 6c4f3c3b-302c-480d-ab00-c257dc4ca9f5
grid = MUST.Grid(snap)

# ╔═╡ 6154b79e-405f-47ab-a23a-5805b5fad950
md"We then can decide onto what grid we want to interpolate the whole cube. We can e.g. increase the resolution by a factor of 2."

# ╔═╡ ceb266ff-7069-48b4-b395-d89f8e1d4158
function scale_axis(axis, factor)
	axis_new = range(first(axis), last(axis), length=Int(ceil(length(axis).*factor)))
	Base.convert.(eltype(axis), axis_new)
end

# ╔═╡ 58d7df98-b3a9-4846-8d0d-d2f1506953ec
begin
	factor = 0.5
	x_new = scale_axis(MUST.axis(snap, :x), factor)
	y_new = scale_axis(MUST.axis(snap, :y), factor)
	z_new = scale_axis(MUST.axis(snap, :z), factor)
end

# ╔═╡ 3cad6318-3fb1-4061-bf98-e347f8cda6d5
begin
@info "Old size vs. new size (x): $(length(MUST.axis(snap, :x)))-$(length(x_new))"
@info "Old size vs. new size (y): $(length(MUST.axis(snap, :y)))-$(length(y_new))" 
@info "Old size vs. new size (z): $(length(MUST.axis(snap, :z)))-$(length(z_new))" 
end

# ╔═╡ ca9a1b49-26e0-4925-a186-c27903ac3c75
target_grid = MUST.Grid(x_new, y_new, z_new)

# ╔═╡ f558c152-9e91-4b2f-8738-128804cd9d0e
md"Now interpolate the old cube to this new cube. First build the interpolator and then evaluate it."

# ╔═╡ b48a0e3e-5358-4028-b536-3009223ef7f4
ip = MUST.ginterpolate(grid, target_grid)

# ╔═╡ f1601c61-a499-4e59-ae84-74278b3d85e6
T_old = log.(deepcopy(snap[:T]));

# ╔═╡ 92b25162-7264-4cf9-bda4-37847fb2bde3
T_new = MUST.gevaluate!(ip, T_old);

# ╔═╡ 05807cbb-e21d-4635-b88d-e7d10ead9ce3
function haver(cube)
	mv = similar(cube, size(cube, 3))
	for k in axes(cube, 3)
		mv[k] = MUST.mean(@view(cube[:,:,k]))
	end

	mv
end

# ╔═╡ 502fb3e2-9036-4637-ac17-759aaeea675a
begin
	plot(framestyle=:box, grid=false)
	
	scatter!(z_new, haver(T_new), 
				markershape=:circle, label="resample", markersize=4)
	plot!(MUST.axis(snap, :z), haver(T_old), 
				label="original", lw=2)
	
end

# ╔═╡ 36f3d298-2aec-48dc-a305-18b19fb6d3b3
begin
	i0_new = argmin(abs.(z_new))
	i0_old = argmin(abs.(MUST.axis(snap, :z)))
		
	h1 = heatmap(T_new[:,:,i0_new])
	h2 = heatmap(T_old[:,:,i0_old])

	plot(h1, h2, size=(1300,600))
	
end

# ╔═╡ Cell order:
# ╟─82119670-88be-4e8e-88e8-8c8c2e706b4f
# ╠═13e36110-0f52-11ee-0142-bd1afe3d9a36
# ╟─72bf89a2-6ef5-4112-abe1-e7a8c26897f7
# ╠═12883f0f-d6f6-40c8-b389-8a18bdc1606b
# ╠═ce272b96-8e58-41d9-852c-44971ce370b7
# ╟─04fb3d30-9bd5-47f9-b58f-f0a8124b7593
# ╠═6c4f3c3b-302c-480d-ab00-c257dc4ca9f5
# ╟─6154b79e-405f-47ab-a23a-5805b5fad950
# ╠═ceb266ff-7069-48b4-b395-d89f8e1d4158
# ╠═58d7df98-b3a9-4846-8d0d-d2f1506953ec
# ╟─3cad6318-3fb1-4061-bf98-e347f8cda6d5
# ╠═ca9a1b49-26e0-4925-a186-c27903ac3c75
# ╟─f558c152-9e91-4b2f-8738-128804cd9d0e
# ╠═b48a0e3e-5358-4028-b536-3009223ef7f4
# ╠═f1601c61-a499-4e59-ae84-74278b3d85e6
# ╠═92b25162-7264-4cf9-bda4-37847fb2bde3
# ╠═05807cbb-e21d-4635-b88d-e7d10ead9ce3
# ╟─502fb3e2-9036-4637-ac17-759aaeea675a
# ╟─36f3d298-2aec-48dc-a305-18b19fb6d3b3
