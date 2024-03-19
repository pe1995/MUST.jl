### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ b91f494a-ccc7-11ee-0013-cd40032fda95
begin
	using Pkg; Pkg.activate(".")
	using Plots
end

# ╔═╡ 87affc67-4fb6-4b51-b1d8-7d8648c00d65
begin
	x = range(1, 10)
	y = x .^2
	
	plot(framestyle=:box, grid=false, minorticks=true)
	plot!(x, y)
end

# ╔═╡ Cell order:
# ╠═b91f494a-ccc7-11ee-0013-cd40032fda95
# ╠═87affc67-4fb6-4b51-b1d8-7d8648c00d65
