### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 4ab3d223-fbd7-4324-b324-2f907cefd5f5
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using TSO
	using Plots
	using LaTeXStrings
end

# ╔═╡ 88df19df-fe1b-4fe0-b07d-b0007c2e107b
md"# Investigate Grid Opacities
We analyze formation opacities and binned opacities created by the automatic grid computation procedures."

# ╔═╡ 6f4a50b4-0463-4879-a756-9e959b0007ba
md"## Formation opacity"

# ╔═╡ e7ae7228-9f1e-4b94-8bdd-8477a1153495
table = "/home/eitner/shared/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_v1.4"

# ╔═╡ 547dd066-f6c8-46b0-943e-809baaf4ecd1
begin
	star = ARGS[1] 
end

# ╔═╡ 9cab4a03-4ca6-4ea8-9c6c-efd73954cb67
formOpa = reload(
	SqOpacity, 
	joinpath(table, "combined_formation_opacities_$(star).hdf5"),
	mmap=true
)

# ╔═╡ 9e2a0907-bea3-4fb3-a34c-cc1eec637f59
any(isnan.(formOpa.κ))

# ╔═╡ 46273bb4-c4fc-4292-b4a0-be3f1cdc5494
all(isfinite.(formOpa.κ))

# ╔═╡ 513c51d7-8e85-409c-8e48-809726ef260a
size(formOpa.κ)

# ╔═╡ f2d0b2b5-d671-4dfd-828c-cde7cc721abb
begin
	plot(ylabel=L"\rm -\tau_{ross}\ (\tau_{\lambda}=1)\ [log]", 
		 	xlabel=L"\rm \lambda\ [log\ \AA]", 
			grid=false, framestyle=:box, title=star)
	
	plot!(log10.(formOpa.λ), -log10.(formOpa.κ_ross), 
		label=nothing, color=:black, ticklength=20)

	savefig("form_opacity_$(star).png")
end

# ╔═╡ fc83f1c9-18cf-469a-98cd-9d67134575de
md"## Binned Opacity
We compute the heating rates from the unbinned solution and compare the results"

# ╔═╡ 564c2e52-5bb2-4203-a0ce-d135b0e55210
eos_raw = reload(SqEoS, 
		joinpath(table, "combined_ross_eos.hdf5"))

# ╔═╡ 735888bb-c117-477c-9d50-ac351be2dc9a
opa_raw = reload(SqOpacity, 
		joinpath(table, "combined_opacities.hdf5"), mmap=true)

# ╔═╡ fedfd32c-a3af-4f05-8233-ec0afb399dff
star_folder = "DIS_MARCS_$(star)_v0.1"

# ╔═╡ a2be5613-51c3-4605-b920-aee1aa3b3716
star_folder_E = "DIS_MARCS_E_$(star)_v0.1"

# ╔═╡ 1033c4a3-9755-48a4-b9be-f2fcaaf06e96
eos = reload(SqEoS, joinpath(star_folder, "eos.hdf5"))

# ╔═╡ 7587fc40-91f3-4669-a38c-c4813a51a0e1
opa = reload(SqOpacity, joinpath(star_folder, "binned_opacities.hdf5"))

# ╔═╡ 7b7d0ae7-bf3e-4278-855d-b8b41dee8b9d
weights = ω_midpoint(opa_raw)

# ╔═╡ 9fd3dcc2-5696-46e7-9399-a7c4c5dc8739
opacities = @binned(opa)

# ╔═╡ 5136a956-f1d4-4c77-a49c-78b24cee0c6c
opacities_raw = @binned opa_raw eos_raw

# ╔═╡ 0849c0be-e7b9-4a14-84b4-1671628161bb
md"The model atmosphere we want to use for comparison. We also compute the optical depth based on the unbinned table."

# ╔═╡ 0c1121c7-58ed-4937-a4d2-121731bf58d4
# ╠═╡ show_logs = false
solar_model = @optical Average3D(eos_raw, joinpath(star_folder_E, "inim.dat")) eos_raw opa_raw

# ╔═╡ 63f7c3a0-4dae-460b-85ca-92589306bfa0
md"Next we construct the solvers for these tables."

# ╔═╡ bfa9457f-be16-4c2c-a779-ed9370af152c
solver = Solver(solar_model, @axed(eos), opacities=opacities)

# ╔═╡ 091eb332-fb0a-40a8-842b-82d15731e505
solver_raw = Solver(solar_model, @axed(eos_raw), opacities=opacities_raw)

# ╔═╡ 8d53ccdb-331a-4552-b47c-cadcd1f92ce8
md"Now solve"

# ╔═╡ 09d7d91e-ad1a-4c94-9844-ff17ea884371
q = Qr(solver)

# ╔═╡ 3fe242c9-ca58-46a8-b66a-422ef47e693d
md"In the unbinned case, we have to additionally provide λ integration weights, because the λ integration in the binned case happens during the binning already."

# ╔═╡ 00a0ab81-0a74-449f-a6f1-c0c7ff2cf676
q_raw = Qr(solver_raw, weights)

# ╔═╡ 995d8f32-e4d2-433f-93b5-38fdade719b1
z, lnT, τ = solver_raw.model[:, 1], solver_raw.model[:, 2], reverse(solar_model.τ)

# ╔═╡ 36beeffb-cbaf-4ec3-8965-4a22ded9b024
begin
	mask = log10.(τ) .< 5
	plot(log10.(τ)[mask], q_raw[mask], label="unbinned", lw=3, 
		ylabel="heating rate", xlabel="log τ ross", framestyle=:box, grid=false,
		color=:black)


	plot!(log10.(τ), q, label="binned", lw=3, ls=:dash, color=:red)

	plot!(xlim=(-8, 5.5))
	plot!(legend=:bottomleft)

	savefig(joinpath(star_folder, "heating_rate_$(star).png"))
end

# ╔═╡ b9a4e5c2-ac32-46c7-9ab0-9cf4357d0ea6
begin
	plot(ylabel="heating rate (rel. difference)", xlabel="log τ ross", 
		framestyle=:box, grid=false)
	

	plot!(log10.(τ[mask]), (q[mask] .- q_raw[mask]) ./ q_raw[mask] , lw=3, 
			label=nothing, color=:black)

	plot!(xlim=(-8, 5.5), ylim=(-.5,0.5))
	plot!(legend=:topright)

	savefig(joinpath(star_folder, "heating_rate_rel_$(star).png"))
end

# ╔═╡ Cell order:
# ╟─88df19df-fe1b-4fe0-b07d-b0007c2e107b
# ╠═4ab3d223-fbd7-4324-b324-2f907cefd5f5
# ╟─6f4a50b4-0463-4879-a756-9e959b0007ba
# ╠═e7ae7228-9f1e-4b94-8bdd-8477a1153495
# ╠═547dd066-f6c8-46b0-943e-809baaf4ecd1
# ╠═9cab4a03-4ca6-4ea8-9c6c-efd73954cb67
# ╠═9e2a0907-bea3-4fb3-a34c-cc1eec637f59
# ╠═46273bb4-c4fc-4292-b4a0-be3f1cdc5494
# ╠═513c51d7-8e85-409c-8e48-809726ef260a
# ╠═f2d0b2b5-d671-4dfd-828c-cde7cc721abb
# ╟─fc83f1c9-18cf-469a-98cd-9d67134575de
# ╟─564c2e52-5bb2-4203-a0ce-d135b0e55210
# ╟─735888bb-c117-477c-9d50-ac351be2dc9a
# ╠═fedfd32c-a3af-4f05-8233-ec0afb399dff
# ╠═a2be5613-51c3-4605-b920-aee1aa3b3716
# ╟─1033c4a3-9755-48a4-b9be-f2fcaaf06e96
# ╟─7587fc40-91f3-4669-a38c-c4813a51a0e1
# ╟─7b7d0ae7-bf3e-4278-855d-b8b41dee8b9d
# ╟─9fd3dcc2-5696-46e7-9399-a7c4c5dc8739
# ╟─5136a956-f1d4-4c77-a49c-78b24cee0c6c
# ╟─0849c0be-e7b9-4a14-84b4-1671628161bb
# ╠═0c1121c7-58ed-4937-a4d2-121731bf58d4
# ╟─63f7c3a0-4dae-460b-85ca-92589306bfa0
# ╠═bfa9457f-be16-4c2c-a779-ed9370af152c
# ╠═091eb332-fb0a-40a8-842b-82d15731e505
# ╟─8d53ccdb-331a-4552-b47c-cadcd1f92ce8
# ╠═09d7d91e-ad1a-4c94-9844-ff17ea884371
# ╟─3fe242c9-ca58-46a8-b66a-422ef47e693d
# ╠═00a0ab81-0a74-449f-a6f1-c0c7ff2cf676
# ╠═995d8f32-e4d2-433f-93b5-38fdade719b1
# ╟─36beeffb-cbaf-4ec3-8965-4a22ded9b024
# ╟─b9a4e5c2-ac32-46c7-9ab0-9cf4357d0ea6
