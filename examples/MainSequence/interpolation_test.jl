### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 1a26f95a-ee9a-11ef-06d1-71f3dcdd45a2
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
	using NaturalNeighbours
end

# ╔═╡ 54b63978-54af-49a2-8bbe-199dc7055d75
plt = matplotlib.pyplot

# ╔═╡ 63c8fc1b-36bd-46ad-a340-4a9b4e3ea155
@import_dispatch "../../../dispatch2"

# ╔═╡ a0a877af-e797-45a4-a869-fcb1206562bc
model1 = TSO.Average3D(@in_dispatch("interpolatedModels/IPTEST_E_t57.50g45.00m-5.000_v1.0/inim.dat"), logg=4.5)

# ╔═╡ 2ac891b4-195a-4d0b-9ae3-c1621fbdfcab
model2 = TSO.Average3D(@in_dispatch("interpolatedModels/VIDM5_E_t57.50g45.00m-5.000_v1.0/inim.dat"), logg=4.5)

# ╔═╡ 928a9f11-d20d-465d-9d60-a286ebede6dd


# ╔═╡ 0f50c982-8eac-4a75-b07a-564535ba84e7
let
	x = rand(200) .*2
	y = rand(200)

	z = x.^2 .+ y.^2

	itp = interpolate(x, y, z, method=Triangle())
	x_test = LinRange(0, 1, 50) .*2
	y_test = LinRange(0, 1, 50)

	xx, yy = meshgrid(x_test, y_test)
	zz = itp.(xx, yy)

	plt.close()
	f, ax = plt.subplots(1, 1)

	
	im = ax.scatter(xx, yy, c=log10.(abs.((zz .- (xx.^2 .+ yy.^2)) ./ zz)))
	ax.scatter(x, y, color="red")
	
	f.colorbar(im, ax=ax)

	f
end

# ╔═╡ 49f82e98-82cc-42e2-bee7-ae76a350e9d3


# ╔═╡ 185314b6-8c70-4f0c-a5e3-18a1711554fe
let
	f, ax = plt.subplots(1, 1)

	ax.plot(model1.z, model1.lnT, label="model 1")
	ax.plot(model2.z, model2.lnT, label="model 2")

	ax.legend()
	f
end

# ╔═╡ 7431b3ea-1058-4b80-af82-2b6d5cc87989
let
	f, ax = plt.subplots(1, 1)

	ax.plot(model1.z, model1.lnρ, label="model 1")
	ax.plot(model2.z, model2.lnρ, label="model 2")

	ax.legend()
	f
end

# ╔═╡ Cell order:
# ╠═1a26f95a-ee9a-11ef-06d1-71f3dcdd45a2
# ╠═54b63978-54af-49a2-8bbe-199dc7055d75
# ╠═63c8fc1b-36bd-46ad-a340-4a9b4e3ea155
# ╠═a0a877af-e797-45a4-a869-fcb1206562bc
# ╠═2ac891b4-195a-4d0b-9ae3-c1621fbdfcab
# ╟─928a9f11-d20d-465d-9d60-a286ebede6dd
# ╠═0f50c982-8eac-4a75-b07a-564535ba84e7
# ╟─49f82e98-82cc-42e2-bee7-ae76a350e9d3
# ╟─185314b6-8c70-4f0c-a5e3-18a1711554fe
# ╟─7431b3ea-1058-4b80-af82-2b6d5cc87989
