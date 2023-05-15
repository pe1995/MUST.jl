### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ ccfef480-e409-11ed-050d-4fec63105069
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using CSV
	using MUST  # For the models
	using TSO   # For the EoS
	using PyPlot
	using DelimitedFiles
	using DataFrames
	using Glob
end

# ╔═╡ 781a1f12-16ba-4dcd-9fa4-55692ca2a954
grid = CSV.read("dispatch_grid.mgrid", DataFrame)

# ╔═╡ ee9af5e2-5d1b-4d0a-9b26-662fb8dfb840
giant = MUST.Box(MUST.StaggerSnap(grid[1, "snapshot"], grid[1, "folder"]));

# ╔═╡ 86b35722-669c-45e2-bccf-75d4e4db1f44
begin
	vmin_giant = min(abs(minimum(giant[:ux])), 
						abs(minimum(giant[:uy])), 
						abs(minimum(giant[:uz])))
	vmax_giant = max(abs(maximum(giant[:ux])), 
						abs(maximum(giant[:uy])), 
						abs(maximum(giant[:uz])));

	@info "vmin giant: $(vmin_giant)"
	@info "vmax giant: $(vmax_giant)"
end

# ╔═╡ 06c37b81-0aae-44f0-bae4-78cf0ba6bb91
dwarf = MUST.Box(MUST.StaggerSnap(grid[4, "snapshot"], grid[4, "folder"]));

# ╔═╡ 84d422ab-3807-46c1-8744-d9cf367cfa87
md"First the global velocity field boundary"

# ╔═╡ c9509809-8d76-417a-a6d4-b97d567d58a9
begin
	vmin_dwarf = min(abs(minimum(dwarf[:ux])), 
						abs(minimum(dwarf[:uy])), 
						abs(minimum(dwarf[:uz])))
	vmax_dwarf = max(abs(maximum(dwarf[:ux])), 
						abs(maximum(dwarf[:uy])), 
						abs(maximum(dwarf[:uz])));

	@info "vmin dwarf: $(vmin_dwarf)"
	@info "vmax dwarf: $(vmax_dwarf)"
end

# ╔═╡ cf61fc4b-41ae-49da-97c4-ac25912bad7f
begin
	resolution(model) = minimum(diff(MUST.axis(model, :z)))
	Δt(R, u, c=0.3) = c * R / u
end

# ╔═╡ edff72fa-d804-427e-9e7c-3dff1ded6734
velocity(model) = max(  abs(MUST.mean(model[:ux])), 
						abs(MUST.mean(model[:uy])), 
						abs(MUST.mean(model[:uz])) )

# ╔═╡ 81eb5c00-49eb-4485-8a6a-cc02ed2102eb
begin
	@info "Δt (giant, vmin) = $(Δt(resolution(giant), vmin_giant)) s"
	@info "Δt (giant, vmax) = $(Δt(resolution(giant), vmax_giant)) s"
	@info "Δt (dwarf, vmin) = $(Δt(resolution(dwarf), vmin_dwarf)) s"
	@info "Δt (dwarf, vmax) = $(Δt(resolution(dwarf), vmax_dwarf)) s"
end

# ╔═╡ bb38e83b-c943-40b2-82fb-9be37d8b46f0
velocity(giant)

# ╔═╡ bb643f3d-d96a-4bad-8f93-1f33bd2ac46e
md"Assuming that we want the timestep to be of the order of 1e-1 in order to be on the same scale as the Newton cooling and the damping time scale, we get the following sclaing factors:"

# ╔═╡ c3fb20f9-5143-42a2-b48f-1a63fac62f17
scaling(δt, goal_scale=1e-3) = round(exp10.(round(log10(δt / goal_scale), sigdigits=3)), sigdigits=1)

# ╔═╡ 93e3b3af-d410-40d1-928c-48f7970646c9
begin
	@info "Scaling (giant, vmin): $(scaling(Δt(resolution(giant), 
								vmin_giant, 0.1)))"
	@info "Scaling (giant, vmax): $(scaling(Δt(resolution(giant), 
								vmax_giant, 0.1)))"
	@info "Scaling (giant, mean): $(scaling(Δt(resolution(giant), 
								velocity(giant), 0.1)))"
	
	@info "Scaling (dwarf, vmin): $(scaling(Δt(resolution(dwarf), 
								vmin_dwarf, 0.3)))"
	@info "Scaling (dwarf, vmax): $(scaling(Δt(resolution(dwarf), 
								vmax_dwarf, 0.3)))"
	@info "Scaling (dwarf, mean): $(scaling(Δt(resolution(dwarf), 
								velocity(dwarf), 0.3)))"
end

# ╔═╡ Cell order:
# ╠═ccfef480-e409-11ed-050d-4fec63105069
# ╟─781a1f12-16ba-4dcd-9fa4-55692ca2a954
# ╠═ee9af5e2-5d1b-4d0a-9b26-662fb8dfb840
# ╠═86b35722-669c-45e2-bccf-75d4e4db1f44
# ╠═06c37b81-0aae-44f0-bae4-78cf0ba6bb91
# ╟─84d422ab-3807-46c1-8744-d9cf367cfa87
# ╟─c9509809-8d76-417a-a6d4-b97d567d58a9
# ╟─cf61fc4b-41ae-49da-97c4-ac25912bad7f
# ╠═edff72fa-d804-427e-9e7c-3dff1ded6734
# ╟─81eb5c00-49eb-4485-8a6a-cc02ed2102eb
# ╠═bb38e83b-c943-40b2-82fb-9be37d8b46f0
# ╟─bb643f3d-d96a-4bad-8f93-1f33bd2ac46e
# ╠═c3fb20f9-5143-42a2-b48f-1a63fac62f17
# ╟─93e3b3af-d410-40d1-928c-48f7970646c9
