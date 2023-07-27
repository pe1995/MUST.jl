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
model_location = MUST.@in_dispatch "data/pretty_good_sun_new_magg3"

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
"""
    scale_axis(axis, factor)

Scale the number of points on an axis by the given factor or to a specific number.
"""
function scale_axis(axis; factor=nothing, N=nothing)
    if isnothing(factor) & isnothing(N)
        error("N or factor required.")
    end

    len = if isnothing(N)
        Int(ceil(length(axis).*factor))
    else
        N
    end

	axis_new = range(first(axis), last(axis), length=len)
	Base.convert.(eltype(axis), axis_new)
end

# ╔═╡ 58d7df98-b3a9-4846-8d0d-d2f1506953ec
begin
	factor = 2.0
	x_new = scale_axis(MUST.axis(snap, :x), factor=1.0)
	y_new = scale_axis(MUST.axis(snap, :y), factor=1.0)
	z_new = scale_axis(MUST.axis(snap, :z), factor=factor)
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
ip = MUST.ginterpolate(grid, target_grid, method=:linear)

# ╔═╡ 2f9c0cc7-85ed-4804-8ef0-d4399dd34f71
MUST.method(ip)

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

	plot!(legendforegroundcolor=nothing)
	
end

# ╔═╡ 36f3d298-2aec-48dc-a305-18b19fb6d3b3
begin
	i0_new = argmin(abs.(z_new))
	i0_old = argmin(abs.(MUST.axis(snap, :z)))
		
	h1 = heatmap(T_new[:,:,i0_new])
	h2 = heatmap(T_old[:,:,i0_old])

	plot(h1, h2, size=(1800,600))
	
end

# ╔═╡ b9391c54-6b34-488c-8b92-e58625c36cb9
md"## Resampling the Box
In line with those tools, the entire box can be re-sampled in the following manner"

# ╔═╡ 720b6b6f-b7e0-4e21-9ef8-6068490d095b
function resample(b::MUST.Box; nx=size(b, 1), ny=size(b, 2), nz=size(b, 3), 
					method=:linear)
	if (nx==size(b, 1)) & (ny==size(b, 2)) & (nz==size(b, 3))
		@warn "Size of new box = size of old box."
		deepcopy(b)
	end

	# The coordinate grid of the input box
	grid = MUST.Grid(b)

	# build the new axis
	x_new = scale_axis(MUST.axis(snap, :x), N=nx)
	y_new = scale_axis(MUST.axis(snap, :y), N=ny)
	z_new = scale_axis(MUST.axis(snap, :z), N=nz)
	target_grid = MUST.Grid(x_new, y_new, z_new)

	# interpolate 
	ip = MUST.ginterpolate(grid, target_grid, method=method)

	# compute the interpolate quantities
	data_new = Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}()

	for f in keys(snap.data)
		d = if all(snap.data[f] .> 0.0) & (eltype(snap.data[f]) <: AbstractFloat)
			log.(snap.data[f])
		else
			snap.data[f]
		end
		
		data_new[f] = MUST.gevaluate!(ip, d)

		data_new[f] .= if all(snap.data[f] .> 0.0) & 
							(eltype(snap.data[f]) <: AbstractFloat)
			exp.(data_new[f])
		else
			data_new[f]
		end
	end

	xx, yy, zz = MUST.meshgrid(x_new, y_new, z_new)
	MUST.Box(xx, yy, zz, data_new, deepcopy(snap.parameter))
end

# ╔═╡ d51b3b07-382c-4a21-af9e-d38be320d5b1
keys(snap.data)

# ╔═╡ 4ce857f4-d78e-451c-a1b9-d7fba54d0343
is_log(x) = begin
	sx = String(x)
	
	xnew, f = if occursin("log10", sx)
		Symbol(sx[findfirst("log10", sx)[end]+1:end]), log10
	elseif occursin("log", sx)
		Symbol(sx[findfirst("log", sx)[end]+1:end]), log
	else
		x, identity
	end

	xnew, f
end

# ╔═╡ 19f1e6e8-d574-47b1-bcbd-bb61e7e2804a
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if xs == :τ_ross
		logx.(MUST.axis(model, xs, 3)), logy.(MUST.plane_statistic(f, model, ys)) 
	else
		logx.(MUST.axis(model, xs)), logy.(MUST.plane_statistic(f, model, ys))
	end
end

# ╔═╡ ee4db621-4d20-4173-b708-a9018db1bc47
md"Is there a difference when computing the optical depth of the resampled model?"

# ╔═╡ 48e2d238-03b6-443a-b32f-cc880dc1132a
box_resample0 = MUST.gresample(
	snap, nz=150, method=:linear
)

# ╔═╡ 27c0fb5c-15b1-423d-88ea-e58bcea4e12e
keys(box_resample0.data)

# ╔═╡ 3391cd19-8747-433e-be1a-f0fe570bf97b
box_resample = MUST.gresample(
	snap, nz=280, method=:linear
)

# ╔═╡ 12988bcb-3c0c-471a-9afb-3ba199f49361
begin
	plot(framestyle=:box, grid=false)
	
	plot!(profile(MUST.mean, snap, :z, :T)...,
				label="original", lw=2)
	plot!(profile(MUST.mean, box_resample, :z, :T)...,
				label="resample", lw=2, ls=:dash)

	plot!(legendforegroundcolor=nothing)
end

# ╔═╡ 485a35b0-ff12-4804-bf8c-82a62d9eba93
box_resample2 = MUST.gresample(
	snap, nz=299, method=:linear
)

# ╔═╡ 8fb3bbb2-3baa-4eff-a378-b5fb059e0eb5
b_s0 = MUST.height_scale(box_resample0, :τ_ross)

# ╔═╡ bb44545f-75c5-432d-b65a-5ba45b283ab9
b_s = MUST.height_scale(box_resample, :τ_ross)

# ╔═╡ 75170e02-209d-42f7-a159-2a8ee480f148
b_s2 = MUST.height_scale(box_resample2, :τ_ross)

# ╔═╡ 0f0b091c-f796-4d82-b687-74459bc85750
begin
	plot(framestyle=:box, grid=false)

	common_r = range(-4, 4, length=200) |> collect

	x, y = profile(MUST.mean, b_s0, :log10τ_ross, :T)
	m = sortperm(x)
	ip_res0   = MUST.linear_interpolation(x[m], y[m]).(common_r)

	x, y = profile(MUST.mean, b_s, :log10τ_ross, :T)
	m = sortperm(x)
	ip_res   = MUST.linear_interpolation(x[m], y[m]).(common_r)
	
	x, y = profile(MUST.mean, b_s2, :log10τ_ross, :T)
	m = sortperm(x)
	ip_res2   = MUST.linear_interpolation(x[m], y[m]).(common_r)

	x, y = profile(MUST.mean, snap_τ, :log10τ_ross, :T)
	m = sortperm(x)
	ip_ori   = MUST.linear_interpolation(x[m], y[m]).(common_r)

	
	plot!(common_r, (ip_res0 - ip_ori), #./ ip_ori,
				lw=2, ls=:solid, 
				label="resample (150) - original (150)",)
	
	plot!(common_r, (ip_res - ip_ori), #./ ip_ori,
				label="resample (280) - original (150)", lw=2, ls=:solid)

	plot!(common_r, (ip_res2 - ip_ori), # ./ ip_ori,
				label="resample (299) - original (150)", lw=2, ls=:solid)


	plot!(legendforegroundcolor=nothing, legendposition=:bottomleft)
	plot!(ylabel="ΔT [K]", xlabel="τ-ross")
end

# ╔═╡ a2c88970-91bd-4a50-a1f2-af713c4d380e
begin
	plot(framestyle=:box, grid=false)

	common_r2 = range(-.4e8, .4e8, length=1000) |> collect

	x2, y2 = profile(MUST.mean, box_resample0, :z, :e)
	m2 = sortperm(x2)
	ip_res02   = MUST.linear_interpolation(x2[m2], y2[m2]).(common_r2)

	x2, y2 = profile(MUST.mean, box_resample, :z, :e)
	m2 = sortperm(x2)
	ip_res12   = MUST.linear_interpolation(x2[m2], y2[m2]).(common_r2)
	
	x2, y2 = profile(MUST.mean, box_resample2, :z, :e)
	m2 = sortperm(x2)
	ip_res22   = MUST.linear_interpolation(x2[m2], y2[m2]).(common_r2)

	x2, y2 = profile(MUST.mean, snap, :z, :e)
	m2 = sortperm(x2)
	ip_ori2   = MUST.linear_interpolation(x2[m2], y2[m2]).(common_r2)


	plot!(common_r2, (ip_res02 - ip_ori2) ./ ip_ori2 .*100,
				label="resample (150) - original (150)", lw=2, ls=:solid)
	
	plot!(common_r2, (ip_res12 - ip_ori2) ./ ip_ori2 .*100,
				label="resample (280) - original (150)", lw=2, ls=:solid)

	plot!(common_r2, (ip_res22 - ip_ori2) ./ ip_ori2 .*100,
				label="resample (299) - original (150)", lw=2, ls=:solid)


	plot!(legendforegroundcolor=nothing, legendposition=:bottomright)
	plot!(ylabel="ΔT [%]", xlabel="z [cm]")
	
end

# ╔═╡ ebd8b320-3568-415c-b4bc-ad1706f09c39
md"The answer seems to be: No, it is not important for the post processing what the resilution is."

# ╔═╡ 72d86e29-c236-48ab-b5ca-383f8cd5c2b0
md"## Test: pchip interpolation instead of linear
It turns out that it might be relevant to change the interpolation technique to an alternative procedure, especially in the vertical direction."

# ╔═╡ 2f7bb7d0-f90e-4166-bd45-9a9e6b784736
ippchip = MUST.ginterpolate(grid, target_grid, method=:pchip)

# ╔═╡ 566f5a7c-d065-45ee-bc04-8328bf2716c9
T_new_pchip = MUST.gevaluate!(ippchip, T_old);

# ╔═╡ 68a08685-9d40-4ee1-a56c-4fb9af7f8e07
begin
	plot(framestyle=:box, grid=false)
	
	scatter!(z_new, haver(T_new_pchip), 
				markershape=:circle, label="resample", markersize=4)
	plot!(MUST.axis(snap, :z), haver(T_old), 
				label="original", lw=2)

	plot!(legendforegroundcolor=nothing)
end

# ╔═╡ 344a6492-f86a-4756-9395-8e217c5829ec
begin
	i0_new2 = argmin(abs.(z_new))
	i0_old2 = argmin(abs.(MUST.axis(snap, :z)))
		
	h12 = heatmap(T_new_pchip[:,:,i0_new2])
	h22 = heatmap(T_old[:,:,i0_old2])

	plot(h12, h22, size=(1800,600))
end

# ╔═╡ 6f749fc7-26a5-4168-9e1f-b8bf57b91ef2
md"## Test: What does resolution do?"

# ╔═╡ e784defd-6a29-42c6-9837-23e648fb7c36
box_resample_res = MUST.gresample(
	snap, nx=10, ny=10, nz=280, method=:linear
)

# ╔═╡ 75228348-9c13-440d-918a-5fb1d92c96f8
begin
	plot(framestyle=:box, grid=false)
	
	scatter!(
		profile(MUST.mean, box_resample_res, :z, :T)...,
		markershape=:circle, 
		label="resample", 
		markersize=4
	)
	
	plot!(
		profile(MUST.mean, snap, :z, :T)..., 
		label="original", 
		lw=2
	)

	p_line = plot!(legendforegroundcolor=nothing, xlabel="z", ylabel="<T>")

	i0_new3 = argmin(abs.(MUST.axis(box_resample_res, :z)))
	i0_old3 = argmin(abs.(MUST.axis(snap, :z)))
		
	h13 = heatmap(box_resample_res[:T][:,:,i0_new3], title="resampled")
	h23 = heatmap(snap[:T][:,:,i0_old3], title="original")


	plot(h13, h23, p_line, layout=(2, 2), size=(900,600))
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
# ╟─ceb266ff-7069-48b4-b395-d89f8e1d4158
# ╠═58d7df98-b3a9-4846-8d0d-d2f1506953ec
# ╟─3cad6318-3fb1-4061-bf98-e347f8cda6d5
# ╠═ca9a1b49-26e0-4925-a186-c27903ac3c75
# ╟─f558c152-9e91-4b2f-8738-128804cd9d0e
# ╠═b48a0e3e-5358-4028-b536-3009223ef7f4
# ╠═2f9c0cc7-85ed-4804-8ef0-d4399dd34f71
# ╠═f1601c61-a499-4e59-ae84-74278b3d85e6
# ╠═92b25162-7264-4cf9-bda4-37847fb2bde3
# ╟─05807cbb-e21d-4635-b88d-e7d10ead9ce3
# ╟─502fb3e2-9036-4637-ac17-759aaeea675a
# ╟─36f3d298-2aec-48dc-a305-18b19fb6d3b3
# ╟─b9391c54-6b34-488c-8b92-e58625c36cb9
# ╟─720b6b6f-b7e0-4e21-9ef8-6068490d095b
# ╟─2322b69e-57e5-4e0e-b7fa-f34240f19d42
# ╠═d51b3b07-382c-4a21-af9e-d38be320d5b1
# ╟─4ce857f4-d78e-451c-a1b9-d7fba54d0343
# ╟─19f1e6e8-d574-47b1-bcbd-bb61e7e2804a
# ╟─12988bcb-3c0c-471a-9afb-3ba199f49361
# ╟─ee4db621-4d20-4173-b708-a9018db1bc47
# ╠═27c0fb5c-15b1-423d-88ea-e58bcea4e12e
# ╠═48e2d238-03b6-443a-b32f-cc880dc1132a
# ╠═3391cd19-8747-433e-be1a-f0fe570bf97b
# ╠═485a35b0-ff12-4804-bf8c-82a62d9eba93
# ╠═8fb3bbb2-3baa-4eff-a378-b5fb059e0eb5
# ╠═bb44545f-75c5-432d-b65a-5ba45b283ab9
# ╠═75170e02-209d-42f7-a159-2a8ee480f148
# ╠═0f0b091c-f796-4d82-b687-74459bc85750
# ╠═a2c88970-91bd-4a50-a1f2-af713c4d380e
# ╟─ebd8b320-3568-415c-b4bc-ad1706f09c39
# ╟─72d86e29-c236-48ab-b5ca-383f8cd5c2b0
# ╠═2f7bb7d0-f90e-4166-bd45-9a9e6b784736
# ╠═566f5a7c-d065-45ee-bc04-8328bf2716c9
# ╟─68a08685-9d40-4ee1-a56c-4fb9af7f8e07
# ╟─344a6492-f86a-4756-9395-8e217c5829ec
# ╟─6f749fc7-26a5-4168-9e1f-b8bf57b91ef2
# ╠═e784defd-6a29-42c6-9837-23e648fb7c36
# ╠═75228348-9c13-440d-918a-5fb1d92c96f8
