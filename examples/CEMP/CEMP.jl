### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ cc6eeb7e-5974-11ef-0436-fd682756392a
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using PythonPlot
	using LaTeXStrings
end

# ╔═╡ a73ba74d-30ca-4716-9a0c-455aa24cf975
begin
	@import_dispatch "../../../dispatch2"
	@import_m3dis "../../../Multi3D"
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
end;

# ╔═╡ 39703f8e-40b2-476a-8ccb-cf5d208866fa
md"# DISPATCH Models
Load a couple of Dispatch models that can be analysed."

# ╔═╡ 5e8a2581-4cb0-4337-a913-4c4df3121455
models = Dict(
	"metal-poor sun, canbon enhanced" => "M1_E_t57.77g44.40m-4.000_v1.0"
)

# ╔═╡ f17e97d8-3239-46d5-9883-8124cb0170a3
snapshots = Dict(
	"metal-poor sun, canbon enhanced" => 390
)

# ╔═╡ 05895238-a947-429d-9ecb-751ebe199ca9


# ╔═╡ b737400b-0e76-40ae-98d4-01c00e8f2dd3
spectra_4297_4303 = Dict(
	k=>M3DISRun(
		joinpath(
			@in_dispatch(joinpath("data", models[k])), 
			"spectra_sn$(isnap)_lam_4297-4303"
		)
	) for (k, isnap) in snapshots
)

# ╔═╡ 4c85536c-87d8-4d0b-87d1-b2fc32ca5118
cubes = Dict(
	k=>pick_snapshot(@in_dispatch(joinpath("data", m)), snapshots[k]) 
	for (k, m) in models
)

# ╔═╡ 37bcf775-4ccc-44a7-af80-b10ffa8a26c7


# ╔═╡ 1f6071d5-ab91-4eef-9dd2-dea9237f1ec8
md"## Style"

# ╔═╡ 21e58231-22b5-48e6-a8aa-c4bbe30a8c4c
colors = Dict(
	"metal-poor sun, canbon enhanced" => "k"
)

# ╔═╡ 91bb18f2-3f87-4c10-a722-02232bc71ea5
ls = Dict(
	"metal-poor sun, canbon enhanced" => "-"
)

# ╔═╡ 7f14c7c4-b39e-4fef-8e30-0381ac27cb8c
lw = Dict(
	"metal-poor sun, canbon enhanced" => 1.5
)

# ╔═╡ 448e6a35-04da-4d51-bbb7-a60e5c3c00f1
labels = Dict(
	"metal-poor sun, canbon enhanced" => L"\rm [C/Fe]_{model} = +3.0"
)

# ╔═╡ 5bca5b7c-5261-458d-92de-958f907d9531


# ╔═╡ 8ee15ed7-4582-4967-ba7f-4f67f121a606
md"# Average profiles"

# ╔═╡ cf8bfd7b-c571-444f-b57a-8bc35a549a6e
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	for (name, boxes) in cubes
		ax.plot(profile(mean, boxes[2], :log10τ_ross, :T)..., color=colors[name], ls=ls[name], lw=lw[name], label=labels[name])
	end

	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm\ temperature\ [K]")

	ax.set_xlim(-3, 4)
	ax.set_ylim(3500, 12000)
	ax.legend()
	
	f
end

# ╔═╡ 1bb13217-808d-43f8-9705-0792d968f86d
let
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(6, 6))

	for (name, boxes) in cubes
		ax.plot(profile(mean, boxes[2], :log10τ_ross, :log10d)..., color=colors[name], ls=ls[name], lw=lw[name], label=labels[name])
	end

	ax.set_xlabel(L"\rm \log \tau_{ross}")
	ax.set_ylabel(L"\rm\ \log density\ [g\ cm^{-3}]")

	ax.set_xlim(-3, 4)
	ax.set_ylim(-6.7, -5.6)
	ax.legend()

	f
end

# ╔═╡ da40ed04-2c17-4057-9803-3a04175a80ec


# ╔═╡ b99fc273-2918-4f4b-bacd-3b897291638e
md"# Spectra"

# ╔═╡ 202a1cb6-b7b3-446c-8c44-da438b645180
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(10, 5))

	for (k, s) in spectra_4297_4303
		λ, F = flux(s)

		ax.plot(λ, F, color=colors[k], ls=ls[k], lw=lw[k], label=labels[k])
	end

	ax.set_ylim(-0.1, 1.1)
	ax.axhline(1.0, color="k", ls=":")
	ax.set_ylabel(L"\rm flux\ [normalized]")
	ax.set_xlabel(L"\rm wavelength\ [\AA]")
	ax.legend(loc="lower center")
	
	f
end

# ╔═╡ Cell order:
# ╠═cc6eeb7e-5974-11ef-0436-fd682756392a
# ╟─a73ba74d-30ca-4716-9a0c-455aa24cf975
# ╟─39703f8e-40b2-476a-8ccb-cf5d208866fa
# ╠═5e8a2581-4cb0-4337-a913-4c4df3121455
# ╠═f17e97d8-3239-46d5-9883-8124cb0170a3
# ╟─05895238-a947-429d-9ecb-751ebe199ca9
# ╟─b737400b-0e76-40ae-98d4-01c00e8f2dd3
# ╟─4c85536c-87d8-4d0b-87d1-b2fc32ca5118
# ╟─37bcf775-4ccc-44a7-af80-b10ffa8a26c7
# ╟─1f6071d5-ab91-4eef-9dd2-dea9237f1ec8
# ╠═21e58231-22b5-48e6-a8aa-c4bbe30a8c4c
# ╠═91bb18f2-3f87-4c10-a722-02232bc71ea5
# ╠═7f14c7c4-b39e-4fef-8e30-0381ac27cb8c
# ╠═448e6a35-04da-4d51-bbb7-a60e5c3c00f1
# ╟─5bca5b7c-5261-458d-92de-958f907d9531
# ╟─8ee15ed7-4582-4967-ba7f-4f67f121a606
# ╟─cf8bfd7b-c571-444f-b57a-8bc35a549a6e
# ╟─1bb13217-808d-43f8-9705-0792d968f86d
# ╟─da40ed04-2c17-4057-9803-3a04175a80ec
# ╟─b99fc273-2918-4f4b-bacd-3b897291638e
# ╟─202a1cb6-b7b3-446c-8c44-da438b645180
