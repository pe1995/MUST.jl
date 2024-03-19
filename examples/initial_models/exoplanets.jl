### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 2e580624-98f8-11ee-1a97-2be453e8f5e0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using DataFrames
	using PythonPlot
	using PlutoUI
	using LaTeXStrings
	
	plt = matplotlib.pyplot
end;

# ╔═╡ 1d15a3bb-1f6e-4b3a-a448-741da9119908
TableOfContents()

# ╔═╡ b28233a2-a806-449c-ab57-604a9741a949
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ 277a9baf-4c46-4f35-ae13-42177e176e5d
matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"));

# ╔═╡ f3f8e1bf-f901-438e-a2ae-b338dd369f58
md"# PLATO Catalog"

# ╔═╡ e74b67ed-c86b-407c-af46-bd11d46b8825
md"## Targets"

# ╔═╡ 6734f522-8956-4cf6-a0bc-33f032dbdd94
catalogPath = "PLATO/exoplanet.eu_catalog_19-10-2022.csv"

# ╔═╡ 4615e75e-02bd-4366-a16b-bb0168bdea5e
catalog = MUST.CSV.read(catalogPath, DataFrame)

# ╔═╡ defaf6d7-eb72-4cdb-8ad3-2b9ea923746b
md"## Masks"

# ╔═╡ a5643fc2-59bb-400e-a8a1-f8a0b77c3c1d
massMask = .!ismissing.(catalog[!, "star_mass"]);

# ╔═╡ f387975d-fb8a-49c0-a2d6-a8cd108b9e09
radiusMask = .!ismissing.(catalog[!, "star_radius"]);

# ╔═╡ f15bc1c9-a00d-4601-a6ca-e263ed2abd3b
planetStatusMask = catalog[!, "planet_status"] .== "Confirmed";

# ╔═╡ 0389b1b9-86fa-4d5e-a5b9-a0d805a7603d
pmassMask = .!ismissing.(catalog[!, "mass"]);

# ╔═╡ 96213042-7b81-4582-abef-878e51426929
pradiusMask = .!ismissing.(catalog[!, "radius"]);

# ╔═╡ b6512166-eff6-4d47-866a-b9a5c887e1a1
orbitalPeriodMask = .!ismissing.(catalog[!, "orbital_period"]);

# ╔═╡ 82378d19-71a9-4e6d-9438-354684d651b9
semiMajorAxisMask = .!ismissing.(catalog[!, "semi_major_axis"]);

# ╔═╡ f36553ff-3ba1-4ad5-a114-9e144f70c562
starMetallicityMask = .!ismissing.(catalog[!, "star_metallicity"]);

# ╔═╡ d728c980-7794-4c38-8a79-cbc589992781
starTeffMask = .!ismissing.(catalog[!, "star_teff"]);

# ╔═╡ 75c87521-4fcd-4081-858a-da7d0c3a4b60
starDistanceMask = .!ismissing.(catalog[!, "star_distance"]);

# ╔═╡ 72e766fe-afeb-4034-96f7-16462ff1f106
catalogClean = catalog[massMask .& radiusMask .& planetStatusMask .& pmassMask .& pradiusMask .& orbitalPeriodMask .& semiMajorAxisMask .& starMetallicityMask .& starTeffMask .& starDistanceMask, :]

# ╔═╡ 8ae46a65-b690-4d3e-827f-cc79711cd2bf
md"## Additional Parameters"

# ╔═╡ 09b3eae9-d1cc-4c93-8d32-fea0fc4c9319
logg = log10.(6.67e-8 .* (catalogClean[!, "star_mass"] *3.955e33) ./ (catalogClean[!, "star_radius"] .* 6.9634e10) .^2)

# ╔═╡ 853e5542-b4b8-4758-b394-cef13f8c93b1
catalogClean[!, "logg"] = logg

# ╔═╡ 1680701c-3d1c-4c6b-ac28-357823983470
md"## Sub-selection"

# ╔═╡ 8b4d8ba4-39aa-4b86-a4f8-b9cb2ea99f1c
loggSelect = (logg .< 4.5) .& (logg .> 4.0)

# ╔═╡ c3410b35-adb7-43dc-9e79-7c0bca4379b3
fehSelect = (catalogClean[!, "star_metallicity"] .< 0.05) .& (catalogClean[!, "star_metallicity"] .> -0.05)

# ╔═╡ 6897b883-298b-4f3a-9c45-33ab23364ec8
teffSelect = (catalogClean[!, "star_teff"] .< 6000.0) .& (catalogClean[!, "star_teff"] .> 5500.0)

# ╔═╡ 88ee579a-4632-4089-8dd0-55e8cbb27a96
catalogSelect = catalogClean[loggSelect .& fehSelect .& teffSelect, :]

# ╔═╡ d6d0d25e-86fd-4281-a271-55a85277e570
md"# Stagger Grid
As initial conditions, we rely on the Stagger grid."

# ╔═╡ cc585845-d526-41b8-a4f0-16999bb585a4
grid = MUST.StaggerGrid("stagger_grid_full.mgrid")

# ╔═╡ 273cf3bf-804e-4987-a7fd-106f2142b800
deleteat!(grid.info, .!isfile.(grid["av_path"]));

# ╔═╡ 95ef9861-643d-40e9-a7e6-fe999cb6eba0
md"# Final Sample"

# ╔═╡ 6ab4b465-f8f7-4fb5-b510-e905368cdf23
md"## HRD location"

# ╔═╡ 0b5991a3-e7d9-4256-9ebc-4b9550e2f935
let
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	ax.scatter(
		grid["teff"],
		grid["logg"], 
		color="k", label="Stagger grid", marker="s", s=50
	)

	ax.scatter(
		catalogSelect[!, "star_teff"],
		catalogSelect[!, "logg"],
		color="r", marker="x", label="PLATO targets", s=70
	)

	ax.set_xlabel(L"\rm T_{eff}\ [K]")
	ax.set_ylabel(L"\rm \log(g)")
	ax.invert_xaxis()
	ax.invert_yaxis()
	ax.legend(framealpha=0)

	gcf()
end

# ╔═╡ 8ced1b68-f0a4-43b7-810a-3a121e5ebb3e
md"# EoS
We need to add an EoS to each model at this point, because opacities are required in order to create the height scale after interpolating."

# ╔═╡ c6ea5cfe-c64c-4887-8d14-52c253620b2f
mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"

# ╔═╡ bc3054de-c877-467b-87de-2b31e972f178
eos = [
	reload(SqEoS, joinpath(mother_table_path, "ross_combined_eos_magg_m0_a0.hdf5"))
	for _ in 1:nrow(catalogSelect)
]

# ╔═╡ 32a0ed24-e9c9-4489-b3b5-1d2d64770085
md"Also make a note of the given EoS in the grid we interpolate in. NOTE: This is not required if the table already contains information about the matching EoS. This is only needed to recompute the rosseland optical depth scale consistently. You can comment this out if you made sure the column is present, or if this does not concern you. Note that ideally you list the CORRECT EoS for every composition. The following uses the same table for all models! This is not relevant if you attempt to interpolate only one metallicity, but otherwise it is crutial!"

# ╔═╡ d5264ffe-ff35-484e-be53-fb1c1f39c345
grid.info[!, "matching_eos"] = [
	joinpath(mother_table_path, "ross_combined_eos_magg_m0_a0.hdf5") 
	for _ in 1:nrow(grid.info)
]

# ╔═╡ d961226b-fba6-4641-8deb-d304c38cea18
md"## Interpolated Initial Condition
We can guess the initial condition by looking at the corresponding Stagger models and interpolated between them. The initial model will always be an average model, however it can be replaced by an adiabat later, when the binning is done."

# ╔═╡ 0fe8fc46-0fca-4c38-a068-2e4669c61e1b
begin
	teffFinal = zeros(nrow(catalogSelect))
	loggFinal = zeros(nrow(catalogSelect))
	fehFinal = zeros(nrow(catalogSelect))

	teffFinal .= catalogSelect[!, "star_teff"]
	loggFinal .= catalogSelect[!, "logg"]
end

# ╔═╡ 369477df-f078-4458-95f9-df0f0597f61e
begin
	ig = modelgrids.interpolate_from_grid(
		grid, 
		teff=teffFinal,
		logg=loggFinal, 
		feh=fehFinal,
		folder="PLATO",
		eos=eos
	)
end

# ╔═╡ 45c4f9f6-6216-40fc-9fed-bd0158421a71


# ╔═╡ 488b10f3-5ed7-441b-bbd6-ef3c05c1dfa1
md"Also make a note of the used EoS table for the opacity binning later."

# ╔═╡ c4e647a7-2086-4475-a8e8-9e57e8bf8f57
ig.info[!, "eos_root"] = [mother_table_path for _ in 1:nrow(ig.info)]

# ╔═╡ 9f6804f0-0248-4fe4-80b4-a272b003d7b5
md"Save the grid for later."

# ╔═╡ 8670e802-7835-4d5f-96f0-1241ae102af1
MUST.save(ig, "PLATO/plato_initial_121223.mgrid")

# ╔═╡ Cell order:
# ╠═2e580624-98f8-11ee-1a97-2be453e8f5e0
# ╟─1d15a3bb-1f6e-4b3a-a448-741da9119908
# ╟─b28233a2-a806-449c-ab57-604a9741a949
# ╟─277a9baf-4c46-4f35-ae13-42177e176e5d
# ╟─f3f8e1bf-f901-438e-a2ae-b338dd369f58
# ╟─e74b67ed-c86b-407c-af46-bd11d46b8825
# ╠═6734f522-8956-4cf6-a0bc-33f032dbdd94
# ╠═4615e75e-02bd-4366-a16b-bb0168bdea5e
# ╟─defaf6d7-eb72-4cdb-8ad3-2b9ea923746b
# ╠═a5643fc2-59bb-400e-a8a1-f8a0b77c3c1d
# ╠═f387975d-fb8a-49c0-a2d6-a8cd108b9e09
# ╠═f15bc1c9-a00d-4601-a6ca-e263ed2abd3b
# ╠═0389b1b9-86fa-4d5e-a5b9-a0d805a7603d
# ╠═96213042-7b81-4582-abef-878e51426929
# ╠═b6512166-eff6-4d47-866a-b9a5c887e1a1
# ╠═82378d19-71a9-4e6d-9438-354684d651b9
# ╠═f36553ff-3ba1-4ad5-a114-9e144f70c562
# ╠═d728c980-7794-4c38-8a79-cbc589992781
# ╠═75c87521-4fcd-4081-858a-da7d0c3a4b60
# ╠═72e766fe-afeb-4034-96f7-16462ff1f106
# ╟─8ae46a65-b690-4d3e-827f-cc79711cd2bf
# ╠═09b3eae9-d1cc-4c93-8d32-fea0fc4c9319
# ╠═853e5542-b4b8-4758-b394-cef13f8c93b1
# ╟─1680701c-3d1c-4c6b-ac28-357823983470
# ╠═8b4d8ba4-39aa-4b86-a4f8-b9cb2ea99f1c
# ╠═c3410b35-adb7-43dc-9e79-7c0bca4379b3
# ╠═6897b883-298b-4f3a-9c45-33ab23364ec8
# ╠═88ee579a-4632-4089-8dd0-55e8cbb27a96
# ╟─d6d0d25e-86fd-4281-a271-55a85277e570
# ╠═cc585845-d526-41b8-a4f0-16999bb585a4
# ╠═273cf3bf-804e-4987-a7fd-106f2142b800
# ╟─95ef9861-643d-40e9-a7e6-fe999cb6eba0
# ╟─6ab4b465-f8f7-4fb5-b510-e905368cdf23
# ╟─0b5991a3-e7d9-4256-9ebc-4b9550e2f935
# ╟─8ced1b68-f0a4-43b7-810a-3a121e5ebb3e
# ╠═c6ea5cfe-c64c-4887-8d14-52c253620b2f
# ╠═bc3054de-c877-467b-87de-2b31e972f178
# ╟─32a0ed24-e9c9-4489-b3b5-1d2d64770085
# ╠═d5264ffe-ff35-484e-be53-fb1c1f39c345
# ╟─d961226b-fba6-4641-8deb-d304c38cea18
# ╠═0fe8fc46-0fca-4c38-a068-2e4669c61e1b
# ╠═369477df-f078-4458-95f9-df0f0597f61e
# ╟─45c4f9f6-6216-40fc-9fed-bd0158421a71
# ╟─488b10f3-5ed7-441b-bbd6-ef3c05c1dfa1
# ╠═c4e647a7-2086-4475-a8e8-9e57e8bf8f57
# ╟─9f6804f0-0248-4fe4-80b4-a272b003d7b5
# ╠═8670e802-7835-4d5f-96f0-1241ae102af1
