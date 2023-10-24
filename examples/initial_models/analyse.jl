### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 4a837d4c-724c-11ee-0b6b-21a88a56df5b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using PythonPlot
	using TSO 
	using LaTeXStrings
	using DelimitedFiles
end;

# ╔═╡ 9e9c5903-762d-43d7-adf2-0eccee134974
begin
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	MUST.@import_dispatch "../../../dispatch2"
end;

# ╔═╡ 4abf1ef3-e08c-46db-a4b5-a986434c2938
names = [
	"DIS_MARCS_E_t5777g44m00_v0.1",
]

# ╔═╡ 40e718a0-bed0-46bd-a77f-4c7039290b45
out_folder = [
	@in_dispatch("data/grid_t50g40m00")
]

# ╔═╡ ce9f7f32-078b-44b7-93d9-d10330c9d26a
eos_folder = [
	@in_dispatch("input_data/grd/DIS_MARCS_E_t50g40m00_v0.1"),
]

# ╔═╡ b3df2bd4-8ed3-4ecb-97ff-dcd1e86c053d
labels = [
	"test subgiant"
]

# ╔═╡ 17ca4c8c-9563-4be4-8ab7-a7446ba28dec
colors = ["black"]

# ╔═╡ cfc0ab3e-5989-414a-8631-97acee2ae886
list_of_snapshots(out_folder[1])

# ╔═╡ ab79e004-4abb-449b-a50c-d8f093efa8e6
begin
	eos = [reload(SqEoS, joinpath(eos_folder[i], "eos.hdf5")) 
				for i in eachindex(eos_folder)]
	opa = [reload(SqOpacity, joinpath(eos_folder[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder)]
end

# ╔═╡ b633c06f-afaf-429d-8f41-4aea94100853
begin
	snapshots = []
	snapshots_τ = []
	
	for i in eachindex(names)
		snapshot, snapshot_τ = pick_snapshot(
			converted_snapshots(out_folder[i]), :recent
		)
		
		append!(snapshots, [snapshot])
		append!(snapshots_τ, [snapshot_τ])
	end
end

# ╔═╡ 3c9cecd9-8c59-4311-93fa-5137655fdfef
initial_model = [
	@optical(Average3D(eos[i], joinpath(eos_folder[i], "inim.dat")), eos[i], opa[i])
	for i in eachindex(eos)
]

# ╔═╡ f3ffa698-bcea-4ab3-9b63-84d518c14068
md"# Average figures"

# ╔═╡ ef9b0e78-00f0-41d0-8429-2d36f54c0a02
begin
	plt.close()

	fA, axA = plt.subplots(1, 1, figsize=(5, 6))

	for i in eachindex(names)
		axA.plot(
			-initial_model[i].z, exp.(initial_model[i].lnT), 
			lw=1., 
			color=colors[i],
			alpha=0.7,
			ls="--"
		)
		
		axA.plot(
			profile(mean, snapshots[i], :z, :T)..., 
			lw=2., 
			color=colors[i], 
			label=labels[i]
		)
	end

	axA.legend()
	axA.set_ylabel("temperature [K]")
	axA.set_xlabel("z [cm]")
	
	
	gcf()
end

# ╔═╡ 36eabb43-1660-4e11-8bca-e7f09568d695
begin
	plt.close()

	fB, axB = plt.subplots(1, 1, figsize=(5, 6))

	for i in eachindex(names)
		axB.plot(
			-initial_model[i].z, log10.(exp.(initial_model[i].lnρ)), 
			lw=1., 
			color=colors[i],
			alpha=0.7,
			ls="--"
		)
		
		axB.plot(
			profile(mean, snapshots[i], :z, :log10d)..., 
			lw=2., 
			color=colors[i], 
			label=labels[i]
		)
	end

	axB.legend()
	axB.set_ylabel(L"\rm density\ [g\ \times\ cm^{-3}]")
	axB.set_xlabel("z [cm]")
	
	
	gcf()
end

# ╔═╡ e99b4bed-f093-4fe6-ad1d-af0f2235470b
begin
	plt.close()

	fC, axC = plt.subplots(1, 1, figsize=(5, 6))

	for i in eachindex(names)
		axC.plot(
			log10.(exp.(initial_model[i].lnρ)), 
			exp.(initial_model[i].lnT),
			lw=1., 
			color=colors[i],
			alpha=0.7,
			ls="--"
		)

		_, T = profile(mean, snapshots[i], :z, :T)
		_, d = profile(mean, snapshots[i], :z, :log10d)
		
		axC.plot(
			d, T,
			lw=2., 
			color=colors[i], 
			label=labels[i]
		)
	end

	axC.legend()
	axC.set_xlabel(L"\rm density\ [g\ \times\ cm^{-3}]")
	axC.set_ylabel(L"\rm temperature\ [K]")	
	
	gcf()
end

# ╔═╡ 6c4f5ac6-a03b-4237-aa8f-fe50f77bde6f
begin
	plt.close()

	fE, axE = plt.subplots(1, 1, figsize=(5, 6))

	for i in eachindex(names)
		axE.plot(
			log10.(initial_model[i].τ), exp.(initial_model[i].lnT), 
			lw=1., 
			color=colors[i],
			alpha=0.7,
			ls="--"
		)
		
		axE.plot(
			profile(mean, snapshots_τ[i], :log10τ_ross, :T)..., 
			lw=2., 
			color=colors[i], 
			label=labels[i]
		)
	end

	axE.legend()
	axE.set_ylabel("temperature [K]")
	axE.set_xlabel(L"\rm \tau_{ross}\ [cm]")
	
	
	gcf()
end

# ╔═╡ 292e9c4b-10d9-4cf9-b5e7-2594eec4bcfd
md"# Optical surface"

# ╔═╡ 12446032-6cc4-4eac-94cc-6ccb8de46c5d
extent(snap) = begin
	x = MUST.axis(snap, :x) ./1e8
	y = MUST.axis(snap, :x) ./1e8

	[minimum(x), maximum(x), minimum(y), maximum(y)]
end

# ╔═╡ 488b2eec-daf8-4920-9140-7d67c6ca3de1
begin
	uz = []
	for (i,snapshot_τ) in enumerate(snapshots_τ)
		uz_τ_surf = if !isnothing(snapshot_τ)
			isurf = MUST.closest(log10.(MUST.axis(snapshot_τ, :τ_ross, 3)), 0)
			snapshot_τ[:uz][:, :, isurf]
		else
			@info "Interpolating Uz..."
			uz_τ_surf = MUST.interpolate_to(snapshots[i], :uz, τ_ross=1)
			MUST.axis(uz_τ_surf, :uz, 3)
		end

		append!(uz, [uz_τ_surf])
	end

	uz
end

# ╔═╡ 1ffcf11f-58dd-41dd-921f-995d0a84f0d0
begin

	fDs, axDs = [], []

	vmin = minimum([minimum(u) for u in uz]) ./1e5
	vmax = maximum([maximum(u) for u in uz]) ./1e5

	vmax = min(abs.([vmin, vmax])...)
	vmin = -vmax

	for i in eachindex(names)
		plt.close()
		
		fD, axD = plt.subplots(1, 1, figsize=(5, 6))
		
		i = axD.imshow(
			uz[i] ./1e5,
			origin="lower",
			vmin=vmin, vmax=vmax,
			extent=extent(snapshots[i]),
			cmap="coolwarm"
		)

		cb = fD.colorbar(i, ax=axD, fraction=0.046, pad=0.04)
		cb.set_label(L"\rm U_z\ [km\ \times\ s^{-1}]")

		axD.set_xlabel("x [cm]")
		axD.set_ylabel("y [cm]")	

		append!(fDs, [fD])
		append!(axDs, [axD])
	end
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═4a837d4c-724c-11ee-0b6b-21a88a56df5b
# ╟─9e9c5903-762d-43d7-adf2-0eccee134974
# ╠═4abf1ef3-e08c-46db-a4b5-a986434c2938
# ╠═40e718a0-bed0-46bd-a77f-4c7039290b45
# ╠═ce9f7f32-078b-44b7-93d9-d10330c9d26a
# ╠═b3df2bd4-8ed3-4ecb-97ff-dcd1e86c053d
# ╠═17ca4c8c-9563-4be4-8ab7-a7446ba28dec
# ╠═cfc0ab3e-5989-414a-8631-97acee2ae886
# ╟─ab79e004-4abb-449b-a50c-d8f093efa8e6
# ╠═b633c06f-afaf-429d-8f41-4aea94100853
# ╟─3c9cecd9-8c59-4311-93fa-5137655fdfef
# ╟─f3ffa698-bcea-4ab3-9b63-84d518c14068
# ╟─ef9b0e78-00f0-41d0-8429-2d36f54c0a02
# ╟─36eabb43-1660-4e11-8bca-e7f09568d695
# ╟─e99b4bed-f093-4fe6-ad1d-af0f2235470b
# ╟─6c4f5ac6-a03b-4237-aa8f-fe50f77bde6f
# ╟─292e9c4b-10d9-4cf9-b5e7-2594eec4bcfd
# ╟─12446032-6cc4-4eac-94cc-6ccb8de46c5d
# ╟─488b2eec-daf8-4920-9140-7d67c6ca3de1
# ╟─1ffcf11f-58dd-41dd-921f-995d0a84f0d0
