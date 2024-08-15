### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2f1edd2a-b56e-11ee-29e7-c353938e7088
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PythonPlot
	using Plots
	using Images
	using PlutoUI
	using ProgressLogging
end

# ╔═╡ c7dc3b15-6555-4824-872a-d487fe5145ea
md"""
# Dispatch Monitoring Board
Monitoring Board to visualize DISPATCH stellar atmosphere simulations while they are running or after completion. All visualizations rely on the computations from a `MUST.WatchDog`, that converts, computes statistics and deletes converted snapshots afterwards to avoid wasting disk space. 
"""

# ╔═╡ c1ad6c59-fae6-49e7-bfc0-8426c553aa2d
md"# Setup"

# ╔═╡ 7cd7d6f0-8498-44ff-b59c-d298365d6416
TableOfContents()

# ╔═╡ 78c88a26-1e84-4ba2-a8f2-d4c7f6468dd3
@import_dispatch "../../../dispatch2"

# ╔═╡ e2e0b39b-9c60-4630-817b-f180c2631a08
datafolder = @in_dispatch "data" 

# ╔═╡ 409ff57f-8d9d-419b-b448-fdf40c0843b4
begin
	mean = MUST.mean
	plt = matplotlib.pyplot
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	#matplotlib.style.use("dark_background")
	scipy_fft = MUST.pyimport("scipy.fft")

	# create a tempdir for all the movie content
	# this is needed if multiple people try to run the same code at the same time
	v_opt_folder_name = "v_opt" #mktempdir()
end;

# ╔═╡ 6754b2c3-d205-4a12-88b3-53fe62c5637f
md"## Select Run"

# ╔═╡ 41f0864e-26ee-46a6-b4ab-c401a4712941
md"Press the button to check for new monitored runs: $(@bind reload_data Button(\"Refresh\"))"

# ╔═╡ c596d1b3-32c5-4651-a8c1-3100fcd6cd59
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	mask = isdir.(joinpath.(allruns, "monitoring"))
	split.(allruns[mask], "/", keepempty=false) .|> last
end;

# ╔═╡ 8af0c339-237c-42ca-bda5-0b0e44f11c30
begin
	reload_data
	md"""Select one of the available monitored runs:\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""
end

# ╔═╡ 498f998c-9996-43dc-a647-3b08abf9786d
md"Click 'Submit' again to look for new snapshots."

# ╔═╡ 3c05fe1e-1c19-4c30-b34e-a874f63b91bc


# ╔═╡ ed9cc79f-0161-4178-b6f0-81a2bccbf188
md"## Load Monitoring"

# ╔═╡ 28d708db-45a7-4fad-a226-e907ebf88b43
wd = MUST.WatchDog(selectedRun, folder=datafolder)

# ╔═╡ 6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
monitoring = MUST.reload!(wd; mmap=true);

# ╔═╡ 1fcefd1e-1c50-43b5-b203-62920944344a
begin
	"""
		timeevolution(m, group, field) 
	
	Return an array containing the time evolution from all available snapshots in `m` for the given topic `group` and variable `field` in that topic.
	"""
	timeevolution(m, group, field) = [
		moni[group][field] 
		for moni in m
	]
	
	"""
		timeevolution(m, group) 
	
	Return an Dictionary containing arrays of the time evolution from all available snapshots in `m` for the given topic `group` with all variables in that topic as dictionary entried.
	"""
	timeevolution(m, group) = begin
		mg = m[1][group]
		Dict(
			k => [moni[group][k] for moni in m]
			for k in keys(mg)
		)
	end
end

# ╔═╡ 60199001-8f0e-44a6-ae50-e829687c045c


# ╔═╡ 0e4df1b3-52ef-4fce-8f32-bddc966b0516
md"## Deleting Snapshots"

# ╔═╡ 3aadeb1c-3585-46af-8b9d-daf33bfcb7f3
snapshotnames(moni) = [nothing, wd.snapshotsCompleted...];

# ╔═╡ 2d580186-96a6-41b8-91aa-da527691ea1d
md"""You can delete a snapshot by picking its number:\
$(@bind selectDeleted confirm(Select(snapshotnames(monitoring))))"""

# ╔═╡ e3f00981-eb96-49ba-8451-9af7560d3556
md"A deleted snapshot will be removed from the monitoring only. If the watchdog is still running in the background it will be added again at a later time automatically."

# ╔═╡ 73df16bf-b033-494d-b321-608ffff33467
if !isnothing(selectDeleted)
	f = joinpath(datafolder, selectedRun, "monitoring", "snap_$(selectDeleted).hdf5")
	if isfile(f)
		rm(f)
	@info "snap $(selectDeleted) removed from monitoring." 
	end
end

# ╔═╡ bd936d7d-e79f-4f9b-ba54-e0694c6a83f0
md"# General Properties"

# ╔═╡ 2c64fcf2-1a0b-49cf-a3f1-890f152d0650
time = timeevolution(monitoring, "atmosphericParameters", "time")

# ╔═╡ 65c9d5ea-45a2-4ab8-99e8-5eee29935589
snapshots = snapshotnames(monitoring)[2:end]

# ╔═╡ 0f9a1524-9c47-49d6-a554-66db88663093


# ╔═╡ cf85b913-bfd0-4a52-b363-7a70b8675134
begin
	logg = last(monitoring)["atmosphericParameters"]["logg"]
	effective_temperature = last(monitoring)["atmosphericParameters"]["teff"]
	@info "Current Teff, logg" effective_temperature logg
end

# ╔═╡ 835470b0-73d6-4c15-8d13-46abe85e2ed3


# ╔═╡ 2ebb80e6-1966-4673-bc79-50d17add3969
md"## Available Statistics"

# ╔═╡ 9265061d-eaa5-4fc1-a7b9-51392a357c91
md"Available Groups that contain a number of fields with respect to the same topic:"

# ╔═╡ 63b0d9d6-f27a-492e-9018-876db8091914
let
	a = sort(keys(monitoring[1]) |> collect)
	for name in a
		fields = keys(monitoring[1][name]) |> collect
		@info name fields
	end
end

# ╔═╡ c7bb05cd-0798-4989-94de-f7238fea8d89


# ╔═╡ b76f0efd-6b2f-47d3-84e7-ca9fb866292e
md"## Resolution"

# ╔═╡ 4eef1012-c8cd-4579-9a7e-ac2b4d4c7af6
let
	z = timeevolution(monitoring, "geometricalAverages", "z") |> first
	HD_dims = (length(z)*2, length(z)*2, length(z))
	HD_res = maximum(diff(z)) / 1e5
	
	@info "[$selectedRun] HD Resolution (Nx, Ny, Nz), dz [km]" HD_dims HD_res
end

# ╔═╡ 33c1da97-0760-495e-abd5-65531d5e1170


# ╔═╡ 8d6d674b-153b-4357-9f2d-c4e3cb05d059
md"# Individual Snapshots"

# ╔═╡ b9a721cf-46ef-4e3c-a37c-8b35653e31cb
md"Pick the time for which you want to see the status: $(@bind timeSurface confirm(Slider(snapshots, show_value=true, default=last(snapshots))))
"

# ╔═╡ 725f7500-c31b-4efa-9df9-00c51640914a
md"Pick a second time for comparison: $(@bind timeSurface2 confirm(Slider(snapshots, show_value=true)))
"

# ╔═╡ a34793ae-db12-4dac-b2f8-348a88092815
itimeSurface = findfirst(snapshots .== timeSurface);

# ╔═╡ 63f58f6c-954c-44e9-84b1-1aa799323586
itimeSurface2 = findfirst(snapshots .== timeSurface2);

# ╔═╡ a54433ec-14b4-4e5e-a0a7-2e6e50a5fa49
md"Snapshot Info:"

# ╔═╡ 695e28be-31c8-4f44-85e3-72ce387d9da5
@info "You are looking at snapshots (A) $(timeevolution(monitoring, "general", "snapshot")[itimeSurface]) and (B) $(timeevolution(monitoring, "general", "snapshot")[itimeSurface2])"

# ╔═╡ 76a1714e-a5bb-488f-ad93-b8552e4531fd


# ╔═╡ f5ccc550-6bf3-4f0f-baf3-d8de3d216737
md"## Optical Surface"

# ╔═╡ b757013d-aee5-41a4-ab0a-7715ba47bd97
topticalsurfaces = timeevolution(monitoring, "opticalSurfaces")

# ╔═╡ 0639ce7d-955d-448f-84a0-353dfd4e93a3
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[0].imshow(
		topticalsurfaces["uzplane"][itimeSurface] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)


	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[1].imshow(
		topticalsurfaces["uzplane"][itimeSurface2] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")


	ax[0].set_xlabel("x [Mm]")
	ax[1].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")
	
	gcf()
end

# ╔═╡ c24a3b12-c7c0-449b-9a76-6e9c5d475344
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		log10.(exp.(topticalsurfaces["lnDplane"][itimeSurface])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		log10.(exp.(topticalsurfaces["lnDplane"][itimeSurface2])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \log_{10} \rho\ [g\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 3d371088-2322-462b-ab93-7cb49fcdf75f
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		topticalsurfaces["Tplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		topticalsurfaces["Tplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 9758eaf6-9d57-4847-a608-2ba81e1192d0
haskey(topticalsurfaces, "fluxplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		topticalsurfaces["fluxplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		topticalsurfaces["fluxplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")
	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 80cf7a65-87c6-48d1-89d9-90f07ec25d51
haskey(topticalsurfaces, "qrplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		topticalsurfaces["qrplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		topticalsurfaces["qrplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")
	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 5817821d-67f5-4e41-a0d3-7ca12961b0c7


# ╔═╡ 7f77f259-505d-4344-8ee4-8628387f2401
md"## Upper Boundary Surface"

# ╔═╡ cc5fbd5a-c8a0-471a-a56b-0512e4c3989b
tuppersurfaces = timeevolution(monitoring, "upperBoundarySurface")

# ╔═╡ 495e3733-d290-40ab-af63-0eb10a033b53
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tuppersurfaces["x"][itimeSurface] ./1e8
	y = tuppersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[0].imshow(
		tuppersurfaces["uzplane"][itimeSurface] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)


	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[1].imshow(
		tuppersurfaces["uzplane"][itimeSurface2] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")


	ax[0].set_xlabel("x [Mm]")
	ax[1].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")
	
	gcf()
end

# ╔═╡ ed29d53f-00bc-4295-93f6-864a44f92ccb
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tuppersurfaces["x"][itimeSurface] ./1e8
	y = tuppersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		log10.(exp.(tuppersurfaces["lnDplane"][itimeSurface])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		log10.(exp.(tuppersurfaces["lnDplane"][itimeSurface2])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \log_{10} \rho\ [g\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ fa8161aa-0e1c-404f-8c7e-0e3914917df4
let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tuppersurfaces["x"][itimeSurface] ./1e8
	y = tuppersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tuppersurfaces["Tplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tuppersurfaces["Tplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ ef776fe2-5f0f-49a8-b2d5-2e1235d41ec1
haskey(topticalsurfaces, "fluxplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tuppersurfaces["x"][itimeSurface] ./1e8
	y = tuppersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tuppersurfaces["fluxplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tuppersurfaces["fluxplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")
	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ ff3d827b-a0dd-4d86-bfd4-d25fc5f4ab3e
haskey(topticalsurfaces, "qrplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tuppersurfaces["x"][itimeSurface] ./1e8
	y = tuppersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tuppersurfaces["qrplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tuppersurfaces["qrplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]")
	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 0e3d2723-1ecb-4e5a-8125-78f6f27407e2
if haskey(tuppersurfaces, "dtplane") 
	let 
		plt.close()
		f, ax = plt.subplots(1, 2, figsize=(10, 6))
	
		x = tuppersurfaces["x"][itimeSurface] ./1e8
		y = tuppersurfaces["y"][itimeSurface] ./1e8
	
		extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

		#=i = ax[0].imshow(
			tuppersurfaces["Tplane"][itimeSurface],
			origin="lower",
			extent=extent,
			cmap="gist_heat"
		)=#
		i = ax[0].imshow(
			tuppersurfaces["dtplane"][itimeSurface],
			origin="lower",
			extent=extent,
			cmap="rainbow_r", alpha=1
		)
		
		cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)
	
		ax[0].set_xlabel("x [Mm]")
		ax[0].set_ylabel("y [Mm]")
	
	
	
		x = tuppersurfaces["x"][itimeSurface2] ./1e8
		y = tuppersurfaces["y"][itimeSurface2] ./1e8
	
		extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		
		i = ax[1].imshow(
			tuppersurfaces["dtplane"][itimeSurface2],
			origin="lower",
			extent=extent,
			cmap="rainbow_r"
		)
		
		cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
		cb.set_label(L"\rm \Delta t\ [s]")
	
		ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
		ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")
	
		ax[1].set_xlabel("x [Mm]")
		
		gcf()
	end
end

# ╔═╡ 3403a014-2441-4ad2-95f6-2e686ae99ba8


# ╔═╡ 77f5cb10-ce91-44c2-91c6-c01432655121
md"## Lower Boundary Surface"

# ╔═╡ 7532d800-c183-42d5-8bd0-9970ca507cfd
if "lowerBoundarySurface" in keys(monitoring[1])
	tlowersurfaces = timeevolution(monitoring, "lowerBoundarySurface")
end

# ╔═╡ 9504d2ae-a060-4b02-a873-27f4c01d3364
("lowerBoundarySurface" in keys(monitoring[1])) && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[0].imshow(
		tlowersurfaces["uzplane"][itimeSurface] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)


	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[1].imshow(
		tlowersurfaces["uzplane"][itimeSurface2] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm_r"
	)
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")


	ax[0].set_xlabel("x [Mm]")
	ax[1].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")
	
	gcf()
end

# ╔═╡ ab5e5713-7dde-4592-9371-88e2a9c00400
("lowerBoundarySurface" in keys(monitoring[1])) && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		log10.(exp.(tlowersurfaces["lnDplane"][itimeSurface])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		log10.(exp.(tlowersurfaces["lnDplane"][itimeSurface2])),
		origin="lower",
		extent=extent,
		cmap="gist_heat_r"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \log_{10} \rho\ [g\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 2b627607-4d4f-447f-a29a-daacf7977c28
("lowerBoundarySurface" in keys(monitoring[1])) && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tlowersurfaces["Tplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tlowersurfaces["Tplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ b5a51746-f8e3-4760-b73a-a760b2e27c3c
("lowerBoundarySurface" in keys(monitoring[1])) && haskey(tlowersurfaces, "fluxplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tlowersurfaces["fluxplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tlowersurfaces["fluxplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")
	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 40b66681-4ee6-485c-a94f-6d09cf4d56e6
("lowerBoundarySurface" in keys(monitoring[1])) && haskey(tlowersurfaces, "qrplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[0].imshow(
		tlowersurfaces["qrplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tlowersurfaces["qrplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="jet"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]")
	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 97e5e57c-3828-4781-b93e-172a1c3c2cce
("lowerBoundarySurface" in keys(monitoring[1])) && haskey(tlowersurfaces, "qrplane") && let 
	plt.close()
	f, ax = plt.subplots(1, 2, figsize=(10, 6))

	x = tlowersurfaces["x"][itimeSurface] ./1e8
	y = tlowersurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	#=i = ax[0].imshow(
		tuppersurfaces["Tplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="gist_heat"
	)=#
	i = ax[0].imshow(
		tlowersurfaces["dtplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="rainbow_r", alpha=1
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [Mm]")
	ax[0].set_ylabel("y [Mm]")



	x = tlowersurfaces["x"][itimeSurface2] ./1e8
	y = tlowersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tlowersurfaces["dtplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="rainbow_r"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \Delta t\ [s]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [Mm]")
	
	gcf()
end

# ╔═╡ 82d82155-981c-4fdb-b0ae-d8d40ce016f5


# ╔═╡ b0c40c50-3361-4b01-ae87-45ae30387526
md"## Geometrical Profiles"

# ╔═╡ 159019c1-b183-4d25-83a9-0ae27d9eb0ee
begin
	tGeoAv = timeevolution(monitoring, "geometricalAverages")
	tGeoMin = timeevolution(monitoring, "geometricalMinimum")
	tGeoMax = timeevolution(monitoring, "geometricalMaximum")
	tGeoRms = timeevolution(monitoring, "geometricalRMS")
	tGeoQuantiles = Dict(
		k=>timeevolution(
			monitoring,
			"geometrical$(k)thQuantile"
		) 
		for k in [15, 30, 45, 60, 75, 90]
		if "geometrical$(k)thQuantile" in keys(monitoring[1])
	)
end

# ╔═╡ 29803b42-3f55-4fb1-8853-8848bbdd5bd9
md"### Averages"

# ╔═╡ 68477423-a6f7-424f-8fd4-849d63648b57
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoMax["z"][itimeSurface] ./1e8, tGeoMax["T"][itimeSurface]
	ax.plot(
		x, y,
		color="cyan", marker="", ls="-"
	) 
	x, y = tGeoMin["z"][itimeSurface] ./1e8, tGeoMin["T"][itimeSurface]
	ax.plot(
		x, y,
		color="magenta", marker="", ls="-"
	) 
	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["T"][itimeSurface]
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 



	
	x, y = tGeoMax["z"][itimeSurface2] ./1e8, tGeoMax["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="cyan", marker="", ls="--"
	) 
	x, y = tGeoMin["z"][itimeSurface2] ./1e8, tGeoMin["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="magenta", marker="", ls="--"
	) 
	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
	) 

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm T\ [K]")
	ax.legend()

	gcf()
end

# ╔═╡ 8ef763ae-d1bd-40de-a26a-f91f529c03bf
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoMax["z"][itimeSurface] ./1e8, tGeoMax["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="cyan", marker="", ls="-"
	)
	x, y = tGeoMin["z"][itimeSurface] ./1e8, tGeoMin["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="magenta", marker="", ls="-"
	)
	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 



	x, y = tGeoMax["z"][itimeSurface2] ./1e8, tGeoMax["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="cyan", marker="", ls="--",
	) 
	x, y = tGeoMin["z"][itimeSurface2] ./1e8, tGeoMin["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="magenta", marker="", ls="--",
	) 
	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
	) 

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm \log \rho\ [g\ cm^{-3}]")
	ax.legend()

	gcf()
end

# ╔═╡ d9546636-2dbb-4b14-9954-76872b95fd06
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoMax["d"][itimeSurface], tGeoMax["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="cyan", marker="", ls="-"
	) 
	x, y = tGeoMin["d"][itimeSurface], tGeoMin["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="magenta", marker="", ls="-"
	) 
	x, y = tGeoAv["d"][itimeSurface], tGeoAv["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 

	

	x, y = tGeoMax["d"][itimeSurface2], tGeoMax["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="cyan", marker="", ls="--"
	) 
	x, y = tGeoMin["d"][itimeSurface2], tGeoMin["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="magenta", marker="", ls="--"
	) 
	x, y = tGeoAv["d"][itimeSurface2], tGeoAv["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2
	) 

	ax.set_xlabel(L"\rm \log \rho\ [g\ cm^{-3}]")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()

	gcf()
end

# ╔═╡ 7fd50fcd-8a9d-417d-a28a-febd9f0f68f3
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoMax["z"][itimeSurface] ./1e8, tGeoMax["uz"][itimeSurface]./1e5
	ax.plot(
		x, y,
		color="cyan", marker="", ls="-"
	)
	x, y = tGeoMin["z"][itimeSurface] ./1e8, tGeoMin["uz"][itimeSurface]./1e5
	ax.plot(
		x, y,
		color="magenta", marker="", ls="-"
	)
	x, y = tGeoRms["z"][itimeSurface] ./1e8, tGeoRms["uz"][itimeSurface] ./1e5
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 



	x, y = tGeoMax["z"][itimeSurface2] ./1e8, tGeoMax["uz"][itimeSurface2]./1e5
	ax.plot(
		x, y,
		color="cyan", marker="", ls="--",
	) 
	x, y = tGeoMin["z"][itimeSurface2] ./1e8, tGeoMin["uz"][itimeSurface2]./1e5
	ax.plot(
		x, y,
		color="magenta", marker="", ls="--",
	) 
	x, y = tGeoRms["z"][itimeSurface2] ./1e8, tGeoRms["uz"][itimeSurface2] ./1e5
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
	) 

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm rms\ v_z\ [km\ s^{-1}]")
	ax.legend()

	gcf()
end

# ╔═╡ cdfb5217-e298-4ee1-ac5d-f294d04f5abc


# ╔═╡ 274d42e7-20eb-4560-a0f0-e4a631aedcec
md"### Quantiles"

# ╔═╡ 46db420f-2f46-4989-846e-145a70001ddf
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	colors = plt.cm.coolwarm.(Vector(0.0:0.01:1.0))
	sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=0, vmax=100))

	getcolor(q) = begin
		colors[ceil(Int, q)]
	end

	for (q, data) in tGeoQuantiles
		x, y = data["z"][itimeSurface] ./1e8, data["T"][itimeSurface]
		ax.plot(
			x, y,
			color=getcolor(q), 
			marker="", ls="-", lw=2.
		) 
	
		x, y = data["z"][itimeSurface2] ./1e8, data["T"][itimeSurface2]
		ax.plot(
			x, y,
			color=getcolor(q), 
			marker="", ls="--", lw=2.,
		)
	end

	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["T"][itimeSurface]
	ax.plot(
		x, y,
		color="k", 
		marker="", ls="-", lw=2.5,
		label="t = $(time[itimeSurface]) s"
	) 

	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="k", 
		marker="", ls="--", lw=2.,
		label="t = $(time[itimeSurface2]) s"
	)

	f.colorbar(sm, ax=ax, location="right", label="quantile [%]")

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm T\ [K]")
	ax.legend()

	gcf()
end

# ╔═╡ b046b3a0-05d4-4b6c-be80-7396f475be3d
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))
	colors = plt.cm.coolwarm.(Vector(0.0:0.01:1.0))
	sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=0, vmax=100))

	getcolor(q) = begin
		colors[ceil(Int, q)]
	end

	for (q, data) in tGeoQuantiles
		x, y = data["z"][itimeSurface] ./1e8, data["d"][itimeSurface]
		ax.plot(
			x, log10.(y),
			color=getcolor(q), 
			marker="", ls="-", lw=2.
		) 
	
		x, y = data["z"][itimeSurface2] ./1e8, data["d"][itimeSurface2]
		ax.plot(
			x, log10.(y),
			color=getcolor(q), 
			marker="", ls="--", lw=2.,
		)
	end

	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="k", 
		marker="", ls="-", lw=2.5,
		label="t = $(time[itimeSurface]) s"
	) 

	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="k", 
		marker="", ls="--", lw=2.,
		label="t = $(time[itimeSurface2]) s"
	)

	f.colorbar(sm, ax=ax, location="right", label="quantile [%]")

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm \log \rho\ [g\ cm^{-3}]")
	ax.legend()

	gcf()
end

# ╔═╡ b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5


# ╔═╡ db1a9405-e93c-476d-a43b-f11f3138b57a
md"## Optical Profiles"

# ╔═╡ a575dee8-6283-469b-a6ad-365ca9d1b5e9
begin
    tOptAv = timeevolution(monitoring, "opticalAverages")
	tOptMin = timeevolution(monitoring, "opticalMinimum")
	tOptMax = timeevolution(monitoring, "opticalMaximum")
	tOptRms = timeevolution(monitoring, "opticalRMS")
end

# ╔═╡ f5c4d377-09f2-4ecf-9cd6-458c10be2943
md"### Averages"

# ╔═╡ 0fdf3055-4d11-4dea-8a50-e595ef1c112d
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	# First one
	x, y = tOptMin["log10τ_ross"][itimeSurface], tOptMin["T"][itimeSurface]
	ax.plot(
		x, y,
		color="magenta", marker="", ls="-"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface], tOptMax["T"][itimeSurface]
	ax.plot(
		x, y,
		color="cyan", marker="", ls="-"
	) 
	x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["T"][itimeSurface]
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 

	# Second one
	x, y = tOptMin["log10τ_ross"][itimeSurface2], tOptMin["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="magenta", marker="", ls="--"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface2], tOptMax["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="cyan", marker="", ls="--"
	)
	x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
	) 

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	ax.set_ylabel(L"\rm T\ [K]")
	ax.legend()

	gcf()
end

# ╔═╡ eddc02bf-d7ca-41e1-878e-ef1103bf1b0f
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	# first one
	x, y = tOptMin["log10τ_ross"][itimeSurface], tOptMin["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="magenta", marker="", ls="-"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface], tOptMax["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="cyan", marker="", ls="-"
	) 
	x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 

	# second one
	x, y = tOptMin["log10τ_ross"][itimeSurface2], tOptMin["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="magenta", marker="", ls="--"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface2], tOptMax["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="cyan", marker="", ls="--"
	) 
	x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.5
	) 

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	ax.set_ylabel(L"\rm \log \rho\ [g\ cm^{-3}]")
	ax.legend()

	gcf()
end

# ╔═╡ 25c3d608-9440-4ac9-8277-3855ba3b6a7b
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tOptMax["d"][itimeSurface], tOptMax["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="cyan", marker="", ls="-"
	) 
	x, y = tOptMin["d"][itimeSurface], tOptMin["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="magenta", marker="", ls="-"
	) 
	x, y = tOptAv["d"][itimeSurface], tOptAv["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 


	
	x, y = tOptMax["d"][itimeSurface2], tOptMax["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="cyan", marker="", ls="--"
	) 
	x, y = tOptMin["d"][itimeSurface2], tOptMin["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="magenta", marker="", ls="--"
	) 
	x, y = tOptAv["d"][itimeSurface2], tOptAv["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.5
	) 

	ax.set_xlabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()

	gcf()
end

# ╔═╡ e59ab649-3b30-4af1-ab28-6f6d6aacee1c
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	# First one
	x, y = tOptMin["log10τ_ross"][itimeSurface], tOptMin["uz"][itimeSurface]./1e5
	ax.plot(
		x, y,
		color="magenta", marker="", ls="-"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface], tOptMax["uz"][itimeSurface]./1e5
	ax.plot(
		x, y,
		color="cyan", marker="", ls="-"
	) 
	x, y = tOptRms["log10τ_ross"][itimeSurface], tOptRms["uz"][itimeSurface]./1e5
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
	) 

	# Second one
	x, y = tOptMin["log10τ_ross"][itimeSurface2], tOptMin["uz"][itimeSurface2]./1e5
	ax.plot(
		x, y,
		color="magenta", marker="", ls="--"
	) 
	x, y = tOptMax["log10τ_ross"][itimeSurface2], tOptMax["uz"][itimeSurface2]./1e5
	ax.plot(
		x, y,
		color="cyan", marker="", ls="--"
	)
	x, y = tOptRms["log10τ_ross"][itimeSurface2], tOptRms["uz"][itimeSurface2]./1e5
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
	) 

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	ax.set_ylabel(L"\rm rms\ v_z\ [km\ s^{-1}]")
	ax.legend()

	gcf()
end

# ╔═╡ 23c0607b-d70e-4bbd-85ba-7f2b322d93dc


# ╔═╡ 0da614f4-9021-44f3-af89-bd6dab57dc1b
md"## Time Steps
Optionally the radiative timestep can be saved in dispatch. If it is, it should be included in the monitoring by default."

# ╔═╡ d55d7d42-c78c-447c-9959-3689f5341655
if ("dt_rt" in keys(tOptAv)) 
	let
		plt.close()
	
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
		# First one
		x, y = tOptMin["log10τ_ross"][itimeSurface], tOptMin["dt_rt"][itimeSurface]
		ax.plot(
			x, y,
			color="magenta", marker="", ls="-"
		) 
		x, y = tOptMax["log10τ_ross"][itimeSurface], tOptMax["dt_rt"][itimeSurface]
		ax.plot(
			x, y,
			color="cyan", marker="", ls="-"
		) 
		x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["dt_rt"][itimeSurface]
		ax.plot(
			x, y,
			color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s", lw=2.5
		) 
	
		# Second one
		x, y = tOptMin["log10τ_ross"][itimeSurface2], tOptMin["dt_rt"][itimeSurface2]
		ax.plot(
			x, y,
			color="magenta", marker="", ls="--"
		) 
		x, y = tOptMax["log10τ_ross"][itimeSurface2], tOptMax["dt_rt"][itimeSurface2]
		ax.plot(
			x, y,
			color="cyan", marker="", ls="--"
		)
		x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["dt_rt"][itimeSurface2]
		ax.plot(
			x, y,
			color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s", lw=2.
		) 
	
		ax.set_yscale("log")
		ax.set_xlabel(L"\rm \log\ \tau_{ross}")
		ax.set_ylabel(L"\rm \Delta t_{RT}\ [s]")
		ax.legend()
	
		gcf()
	end
end

# ╔═╡ c3fb528f-4bd8-4ae7-bb4c-ac90899fbf21


# ╔═╡ 10e5e8cb-0881-40d0-b392-2f66c0cbfc7c
md"""
## Bolometric Flux
It is also possible to monitor the evolution of the bolometric flux. The effective temperature then is computed as 

$\rm T_{eff} = (F_{bol}/\sigma_{SB})^{1/4}$

where we divide the angle weights by 2 to only account for the outgoing flux in DISPATCH. This is needed because $\rm q_r$ has been integrated over all angles, including downwards radiation, within the RT solver of DISPATCH. 
"""

# ╔═╡ 4c3b9b8a-03a5-494d-a321-7f5c400ed054
teff(f) = (abs.(f) /MUST.σ_S) ^0.25

# ╔═╡ dc0315d8-0616-433b-8182-33840afb0b0f
teff_str(f) = L"\rm T_{eff} =\ "*"$(round(teff(f), sigdigits=5))"*L"\rm \ K"

# ╔═╡ aef86067-68b0-48b8-8bb2-0d410a7521c2
if haskey(tOptAv, "flux")
	let
		plt.close()
	
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
		x, y = tOptMax["log10τ_ross"][itimeSurface], tOptMax["flux"][itimeSurface]
		ax.plot(
			x, y,
			color="cyan", marker="", ls="-"
		) 
		x, y = tOptMin["log10τ_ross"][itimeSurface], tOptMin["flux"][itimeSurface]
		ax.plot(
			x, y,
			color="magenta", marker="", ls="-"
		) 
		x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["flux"][itimeSurface]
		ax.plot(
			x, y,
			color="k", marker="", ls="-", 
			label="t = $(time[itimeSurface]) s, "*teff_str(tGeoAv["flux"][itimeSurface][end]), lw=2.5
		) 
	
	
	
		
		x, y = tOptMax["log10τ_ross"][itimeSurface2], tOptMax["flux"][itimeSurface2]
		ax.plot(
			x, y,
			color="cyan", marker="", ls="--"
		) 
		x, y = tOptMin["log10τ_ross"][itimeSurface2], tOptMin["flux"][itimeSurface2]
		ax.plot(
			x, y,
			color="magenta", marker="", ls="--"
		) 
		x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["flux"][itimeSurface2]
		ax.plot(
			x, y,
			color="k", marker="", ls="--", 
			label="t = $(time[itimeSurface2]) s, "*teff_str(tGeoAv["flux"][itimeSurface2][end]), lw=2.
		) 

		ax.set_yscale("log")
		ax.set_xlabel(L"\rm \log \tau_{ross}")
		ax.set_ylabel(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")
		ax.legend()

		gcf()
	end
end

# ╔═╡ 321e3dda-cd15-4787-95e6-f928125535d5
md"""
# Time evolution
In the following the time evolution of different quantities is shown by computing statistics for every snapshot.
"""

# ╔═╡ 643eaee3-a1f8-4be2-b032-08fcb325fbe8


# ╔═╡ aae74666-7fdd-4db9-8869-68e4597f6940
md"__Snapshot selector__\
You can select a snapshot by moving the vertical line in the plots:"

# ╔═╡ 15f1108e-6f72-4700-89d9-fc2dfaf7ea47
md"$(@bind snapshotSelector confirm(Slider(snapshots, show_value=true, default=first(snapshots))))"

# ╔═╡ e452c94a-3b4b-4dd1-a286-3c32a5c8e995
md"Click to show the selector in the plots: $(@bind activateSnapshotSelector CheckBox(default=false))"

# ╔═╡ c370ebd5-fc25-4e81-b3f9-37bf99b35598
snapshotSelected = findfirst(i->i==snapshotSelector, snapshots);

# ╔═╡ fc550522-d506-4f0e-b582-e7c54cd23be0
plot_time_selector!(ax) = begin
	# Selector
	if activateSnapshotSelector
		#xlim = ax.get_xlim()
		xloc = time[snapshotSelected]./(60*60) 
		ylim = ax.get_ylim()
		yloc = abs(ylim[1] -ylim[0])*0.01 + ylim[1]
		ax.axvline(time[snapshotSelected]./(60*60), color="0.5", ls="--", lw=1)
		
		ax.text(
			xloc, yloc, 
			"$snapshotSelector", 
			color="0.5",
			ha="center", va="bottom"
		)
	end
end

# ╔═╡ 35f64e1d-2273-4178-879e-187b86b24043
md"## Optical Surface"

# ╔═╡ cb3e21b3-0197-4685-a9ae-fe4040106222
begin
	ttempsurface = timeevolution(monitoring, "opticalSurfaces", "Tplane")
	tdsurface = timeevolution(monitoring, "opticalSurfaces", "lnDplane")
	tuzsurface = timeevolution(monitoring, "opticalSurfaces", "uzplane")
end

# ╔═╡ eb335a4d-e662-499c-bb80-8bf38c84329f
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceTemp = [mean(t) for t in ttempsurface]
	surfaceMax = [maximum(t) for t in ttempsurface]
	surfaceMin = [minimum(t) for t in ttempsurface]
	
	ax.plot(
		time./(60*60), surfaceMax, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceMin, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceTemp, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm <T_{\tau=1}>\ [K]")

	plot_time_selector!(ax)

	gcf()
end

# ╔═╡ 35b4aa7e-e41b-4444-af3b-74684c723d5c
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceRho = [log10.(mean(exp.(t))) for t in tdsurface]
	surfaceMax = [log10.(maximum(exp.(t))) for t in tdsurface]
	surfaceMin = [log10.(minimum(exp.(t))) for t in tdsurface]

	ax.plot(
		time./(60*60), surfaceMax, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceMin, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceRho, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 	

	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm \log < \rho_{\tau=1}>\ [g\ cm^{-3}]")

	plot_time_selector!(ax)
	
	gcf()
end

# ╔═╡ 87d43815-9733-40ac-b4e7-81a2c5dbd0d1
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceUz  = [sqrt.(mean(t .^2)) for t in tuzsurface]
	surfaceMax = [maximum(t) for t in tuzsurface]
	surfaceMin = [minimum(t)  for t in tuzsurface]
	
	ax.plot(
		time./(60*60), surfaceMax ./1e5, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceMin ./1e5, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceUz ./1e5, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm rms\ v_z^{\tau=1}\ [km\ s^{-1}]")

	plot_time_selector!(ax)

	gcf()
end

# ╔═╡ c8492ca5-893f-4939-a256-f0872f351c5d


# ╔═╡ df9b58cf-fe4f-4858-a296-879bb0385ba7
md"## Upper Boundary"

# ╔═╡ ccf5f4e6-adc6-419d-a717-4b3b597c2233
ttopgeo = timeevolution(monitoring, "geometricalAverages")

# ╔═╡ 913c1cbe-fedc-4e7d-a8a7-c2ad416a21e6
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceT = (ttopgeo["T"] .|> last)
	surfaceMa = (tGeoMax["T"] .|> last)
	surfaceMi = (tGeoMin["T"] .|> last)
	
	ax.plot(
		time./(60*60), surfaceMi, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceMa, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), surfaceT, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	#ax.set_title("Upper Boundary")
	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm <T_{top}>\ [K]")

	plot_time_selector!(ax)
	
	gcf()
end

# ╔═╡ 2f76eaac-0793-4fa8-9ee1-67dc81e9a1ac
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceT = (ttopgeo["d"] .|> last)
	surfaceMa = (tGeoMax["d"] .|> last)
	surfaceMi = (tGeoMin["d"] .|> last)

	ax.plot(
		time./(60*60), log10.(surfaceMi), 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), log10.(surfaceMa), 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time./(60*60), log10.(surfaceT), 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	#ax.set_title("Upper Boundary")
	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm  <\rho_{top}>\ [g\ cm^{-3}]")

	plot_time_selector!(ax)
	
	gcf()
end

# ╔═╡ aa0376dc-5982-4274-9b4d-4aef1d1de896
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceT = (ttopgeo["uz"] .|> last) ./1e5
	surfaceMa = (tGeoMax["uz"] .|> last) ./1e5
	surfaceMi = (tGeoMin["uz"] .|> last) ./1e5

	ax.plot(
		time./(60*60), surfaceMi, 
		color="magenta", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 
	ax.plot(
		time./(60*60), surfaceMa, 
		color="cyan", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 
	ax.plot(
		time./(60*60), surfaceT, 
		color="k", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	#ax.set_title("Upper Boundary")
	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm  <v_z^{top}>\ [km\ s^{-1}]")

	plot_time_selector!(ax)
	
	gcf()
end

# ╔═╡ 0378dcdf-0291-4e9c-a30b-5dacd346707c
if haskey(ttopgeo, "flux")
	let
		plt.close()
	
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
	
		surfaceT = teff.((ttopgeo["flux"] .|> last))
		surfaceMa = teff.((tGeoMax["flux"] .|> last))
		surfaceMi = teff.((tGeoMin["flux"] .|> last))
		
		ax.plot(
			time./(60*60), surfaceMi, 
			color="magenta", marker="s", markerfacecolor="w", ls="-"
		) 
		ax.plot(
			time./(60*60), surfaceMa, 
			color="cyan", marker="s", markerfacecolor="w", ls="-"
		) 
		ax.plot(
			time./(60*60), surfaceT, 
			color="k", marker="s", markerfacecolor="w", ls="-"
		) 
	
		#ax.set_title("Upper Boundary")
		ax.set_xlabel("time [h]")
		ax.set_ylabel(L"\rm T_{eff}\ [K]")

		plot_time_selector!(ax)
		
	
		gcf()
	end
end

# ╔═╡ b4fd0320-74b2-4ed4-b7b0-beca4eccd553


# ╔═╡ 9e883b44-f225-4566-9d76-ecd3e34d3b5b
md"## Mass Flux"

# ╔═╡ def11de6-b7af-446b-ada7-e61f19b6503b
begin
	tmassfluxgeo = timeevolution(monitoring, "geoMassFlux")
	tmassfluxopt = timeevolution(monitoring, "optMassFlux")
end

# ╔═╡ 01026bad-d93a-4246-a940-e9eb894ccf67
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceMF =[
		tmassfluxopt["massFlux"][i][end] ./1e5
		for i in eachindex(tmassfluxopt["massFlux"])
	]
	
	ax.plot(
		time./(60*60), surfaceMF, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm <\rho\ v_z>\ /\ <\rho>\ [km\ s^{-1}]")

	plot_time_selector!(ax)

	gcf()
end

# ╔═╡ dd1297cd-62d0-4b39-b0a3-e41fab79bc8b


# ╔═╡ 43590370-0a2e-4d6e-9531-f93eb1630acd
md"You can also pick a optical depth: $(@bind opplotmass confirm(TextField(default=\"-3.0\")))"

# ╔═╡ cd9ff780-66d2-4bcb-8b89-145d9d857306
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	masks = [
		sortperm(tmassfluxopt["log10τ_ross"][i]) 
		for i in eachindex(tmassfluxopt["log10τ_ross"])
	]
	
	ip = [
		MUST.linear_interpolation(
			tmassfluxopt["log10τ_ross"][i][masks[i]],
			tmassfluxopt["massFlux"][i][masks[i]]./1e5,
			extrapolation_bc=NaN
		) 
		for i in eachindex(tmassfluxopt["log10τ_ross"])]
	
	surfaceMF = [ip[i].(parse(Float64, opplotmass)) for i in eachindex(ip)]
	
	ax.plot(
		time./(60*60), surfaceMF, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	ax.set_xlabel("time [h]")
	ax.set_ylabel(L"\rm <\rho\ v_z>\ /\ <\rho>\ [km\ s^{-1}]")

	plot_time_selector!(ax)

	gcf()
end

# ╔═╡ 82daeb7f-1a3b-420a-8193-148d7314f557


# ╔═╡ b59908ac-50a9-49fd-bcc2-94c59b205543
md"## Time Step (minimum)"

# ╔═╡ 787b1cc4-9e1f-4421-874d-bdff7c231bf9
if ("dt_rt" in keys(tGeoAv)) 
	let
		plt.close()
	
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
		
		surfaceT = (tGeoAv["dt_rt"] .|> minimum)
		surfaceMa = (tGeoMax["dt_rt"] .|> minimum)
		surfaceMi = (tGeoMin["dt_rt"] .|> minimum)
		
		ax.plot(
			time./(60*60), surfaceMi, 
			color="magenta", marker="s", markerfacecolor="w", ls="-"
		) 
		ax.plot(
			time./(60*60), surfaceMa, 
			color="cyan", marker="s", markerfacecolor="w", ls="-"
		) 
		ax.plot(
			time./(60*60), surfaceT, 
			color="k", marker="s", markerfacecolor="w", ls="-"
		) 

		ax.set_yscale("log")
		ax.set_xlabel("time [h]")
		ax.set_ylabel(L"\rm min\left( \Delta t_{RT}\right)\ [s]")

		plot_time_selector!(ax)
		
		gcf()
	end
end

# ╔═╡ 8ae48e0a-cddf-4775-9a1a-8eada0c07ca7


# ╔═╡ 43b11d24-2496-44f2-9e9f-fc30ff8d285d
md"# Fourier transformation"

# ╔═╡ dbeee710-c9f3-4561-969a-321f256f99f5
md"You can select any number of snapshots you want to include in the Fourier transformation. Make sure that you include enough snapshots so that the time sampling is dense enough, but also make sure you include a long enough time sequence so that you resolve long-wavelength waves."

# ╔═╡ 408f2044-cac0-48b2-8b33-f7ddf9b66849
md"__First snapshot to include:__ $(@bind start_fft confirm(Slider(snapshots, default=snapshots[max(length(snapshots)-200, 1)], show_value=true)))"

# ╔═╡ 7a1d3d74-fea2-402b-aee1-728ff89eb548
md"__Frequency limit [μHz]:__ $(@bind fft_xlim confirm(TextField(default=\"5000\")))"

# ╔═╡ e9b7dfd5-1034-4b95-941c-df003436e26e
md"__Frequency selector [μHz]:__ $(@bind freqSelect confirm(TextField(default=\"3000\")))"

# ╔═╡ 39d3636d-1161-4cbd-a11e-b900e75e77ce


# ╔═╡ c5550bd5-725d-4d3e-9990-324860af7f67
fft(args...; kwargs...) = MUST.pyconvert(Array, scipy_fft.fftn(args...; kwargs...))

# ╔═╡ 43ab1632-3df6-4b76-9385-627452ce97aa
fftfreq(axis) = begin
	s = length(axis)
	d = abs(diff(axis) |> first)
	MUST.pyconvert(Array, scipy_fft.fftfreq(s, d=d))
end

# ╔═╡ 587cd5a6-2117-4a65-b50e-c8d2107adfd9
rms(data) = sqrt(mean(data .^2))

# ╔═╡ 47f967fb-e063-4b23-97e9-02ca76985fa9
let 
	iend_fft = length(time)
	istart_fft = findfirst(snapshots.==start_fft)
	if istart_fft>=1
		data, ylabel = if haskey(topticalsurfaces, "fluxplane")
			dat = topticalsurfaces["fluxplane"][istart_fft:iend_fft]
			[sum(d)/prod(size(d)) for d in dat],  L"\rm fft(<F_{bol}^{\tau=1}>)\ [normalized]"
		else
			dat = topticalsurfaces["uzplane"][istart_fft:iend_fft] ./ 1e5
			[mean(d) for d in dat], L"\rm fft(<v_z^{\tau=1}>)\ [normalized]"
		end
		
		t = time[istart_fft:iend_fft]
		surfaces = zeros(length(data))
		for (i, d) in enumerate(data)
			surfaces[i] = d
		end

		q = fft(surfaces)
		surfaceFFT = sqrt.(real(q).^2 .+ imag(q).^2)
		timeFFT = fftfreq(t)
		meanSurfaceFFT = mean(surfaceFFT, dims=1)
		sortmask = sortperm(timeFFT)
	
		plt.close()
		f, ax = plt.subplots(1, 1, figsize=(5, 6))
		x = timeFFT[sortmask] .* 1e6
		y = surfaceFFT[sortmask] ./ maximum(surfaceFFT[sortmask])

		first_half = 1#floor(Int, length(x) /2) 
		ax.plot(x[first_half:end], y[first_half:end], color="k")

		if haskey(tGeoAv, "flux")
			ax.axvline(exp10(logg) / exp10(4.44) * (5777.0 / teff(tGeoAv["flux"][itimeSurface2][end]))^0.5 * 3166, ls="--", alpha=0.5, color="k")
		end
		
		ax.set_xlabel(L"\rm frequency\ [\mu Hz]")
		ax.set_ylabel(ylabel)
		#ax.set_xlim(1, last(x))
		#ax.set_xlim(-parse(Float64, fft_xlim), parse(Float64, fft_xlim))
		ax.set_xlim(1, parse(Float64, fft_xlim))
		#ax.set_ylim(0, 0.01)
		

		#ax.set_xscale("log")
		ax.set_yscale("log")
		xpos = x[x.>100.0]
		ypos = y[x.>100.0]

		Δt = first(diff(t)) ./(60)
		t_max = (last(t) - first(t)) ./(60*60)
		@info "Fourier transformation including $(length(t)) snapshots."
		@info "Time spacing [min], time sequence length [h]" Δt t_max
		@info "Peak around $(xpos[argmax(ypos)]) μHz = $(1.0/(1e-6*xpos[argmax(ypos)])) s"

		ax.axvline(parse(Float64, freqSelect), ls=":", alpha=0.5, color="r")
		

		f
	end
end

# ╔═╡ 3d530773-bd52-41d6-a3dd-509f237ae439
if "flux" in keys(tGeoAv)
	νmax_expected = exp10(logg) / exp10(4.44) * (5777.0 / teff(tGeoAv["flux"][itimeSurface2][end]))^0.5 * 3166
	νtimescale = 1.0/(νmax_expected*1e-6)
	@info "Expected νmax from parameters (Themeßl et al. 2018) [μHz], [s]:" νmax_expected νtimescale
end

# ╔═╡ 61b3e522-2af4-4413-9c69-ee30b0d4673f
md"# Movies"

# ╔═╡ 713164e6-4a7a-4475-beda-bda33d16e206
md"Pick the snapshot from which you want to start: $(@bind startTimeMovie confirm(Slider(snapshots, show_value=true, default=first(snapshots))))
"

# ╔═╡ fb94a874-a34b-4721-be63-e443f92330af
md"Pick the snapshot on which you want to end: $(@bind endTimeMovie confirm(Slider(snapshots, show_value=true, default=last(snapshots))))
"

# ╔═╡ 2e99b941-935b-4379-89bf-55a4635df686
begin
	i_start_fps = findfirst(i->i==startTimeMovie, snapshots)
	i_end_fps = findfirst(i->i==endTimeMovie, snapshots)
	
	md"Pick the frames-per-second (FPS) of the movies: $(@bind fps confirm(Slider(1:(i_end_fps-i_start_fps+1), default=4, show_value=true)))"
end

# ╔═╡ 48b4180b-6891-4c2a-a22e-7761b8e1f439


# ╔═╡ 9d0348a6-24a1-447c-9af4-5b648a381370
md"""
__Click to create GIFs__: $(@bind createGifImages CheckBox(default=false))
"""

# ╔═╡ c30611f4-1d2e-41c9-b8db-f85b3c99990d
md"From snapshot $(startTimeMovie) to $(endTimeMovie) with $(fps) FPS."

# ╔═╡ 4564cbd9-34bb-4b21-9fca-1027d5c95518


# ╔═╡ 40bf37a4-318b-4263-8705-a566a3321f87
md"
## 2D Surface Movie
You can select any of the 2D plane monitoring results you wish, as well as the image resolution and FPS of the final movie. There are a couple of movie setting already available. If you want to plot something else, please have a look at the cell below (hidden, click the eye symbol on the left to show cell content)"

# ╔═╡ c1faa953-4b76-43fa-a958-a5476e325afc
begin
	struct SurfaceMovieContent
		movie
		label
		xAxis
		yAxis
	end
	surfacesMovieSelection = Dict()
	surfacesMovieSelection["upper boundary - temperature"] = SurfaceMovieContent(
		tuppersurfaces["Tplane"], 
		L"\rm temperature\ [K]", 
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	surfacesMovieSelection["upper boundary - density"] = SurfaceMovieContent(
		[exp.(x) for x in tuppersurfaces["lnDplane"]], 
		L"\rm density\ [g\ cm^{-3}]",
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	if haskey(tuppersurfaces, "fluxplane")
		surfacesMovieSelection["upper boundary - flux"] = SurfaceMovieContent(
			tuppersurfaces["fluxplane"], 
			L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
	surfacesMovieSelection["upper boundary - vertical velocity"] = SurfaceMovieContent(
		tuppersurfaces["uzplane"] ./1e5, 
		L"\rm v_z\ [km\ s^{-1}]", 
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	if haskey(tuppersurfaces, "dtplane")
		surfacesMovieSelection["upper boundary - timestep"] = SurfaceMovieContent(
			tuppersurfaces["dtplane"], 
			L"\rm timestep\ [s]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
	if haskey(tuppersurfaces, "qrplane")
		surfacesMovieSelection["upper boundary - heating"] = SurfaceMovieContent(
			tuppersurfaces["qrplane"], 
			L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
		
		
	surfacesMovieSelection["optical surface - temperature"] = SurfaceMovieContent(
		topticalsurfaces["Tplane"], 
		L"\rm temperature\ [K]", 
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	surfacesMovieSelection["optical surface - density"] = SurfaceMovieContent(
		[exp.(x) for x in topticalsurfaces["lnDplane"]], 
		L"\rm density\ [g\ cm^{-3}]",
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	if haskey(tuppersurfaces, "fluxplane")
		surfacesMovieSelection["optical surface - flux"] = SurfaceMovieContent(
			topticalsurfaces["fluxplane"], 
			L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
	surfacesMovieSelection["optical surface - vertical velocity"] = SurfaceMovieContent(
		topticalsurfaces["uzplane"] ./1e5, 
		L"\rm v_z\ [km\ s^{-1}]", 
		tuppersurfaces["x"][1] ./1e8, 
		tuppersurfaces["y"][1] ./1e8
	)
	if haskey(tuppersurfaces, "dtplane")
		surfacesMovieSelection["optical surface - timestep"] = SurfaceMovieContent(
			topticalsurfaces["dtplane"], 
			L"\rm timestep\ [s]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
	if haskey(tuppersurfaces, "qrplane")
		surfacesMovieSelection["optical surface - heating"] = SurfaceMovieContent(
			topticalsurfaces["qrplane"], 
			L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end
end;

# ╔═╡ b901bb2a-1f3c-46f5-b3f7-1d530975d767
md"""
	Select 2D surface movie topic
$(@bind surfaceMovie_choice Select(sort(keys(surfacesMovieSelection)|>collect), default="optical surface - temperature"))
"""

# ╔═╡ 3c3c8847-ae3e-40d5-97ef-d69e3e5ccd14
begin
	surfaceMovieSelected = surfacesMovieSelection[surfaceMovie_choice]
	surfacesMovie = surfaceMovieSelected.movie
	labelsurfaceMovie = surfaceMovieSelected.label
	xAxis = surfaceMovieSelected.xAxis
	yAxis = surfaceMovieSelected.xAxis
end;

# ╔═╡ 6575c66a-fe2a-4427-bd4c-f23eef344caf
md"""
	Modify movie appearence
"""

# ╔═╡ 5646dc57-693f-4605-9edc-46b0648f3ab5
md"""
figure resolution (DPI): $(@bind dpi confirm(Slider(10:600, default=75, show_value=true)))\

colormap (matplotlib): $(@bind cmap confirm(TextField(default=\"gist_heat\")))\

log scale for color axis: $(@bind logy CheckBox(default=true))
"""

# ╔═╡ 8945b7b7-f3e5-4846-b76e-49d9c864391e
begin		
	i_start = findfirst(i->i==startTimeMovie, snapshots)
	i_end = findfirst(i->i==endTimeMovie, snapshots)
	v_min_movie = minimum(minimum, surfacesMovie[i_start:i_end])
	v_max_movie = maximum(maximum, surfacesMovie[i_start:i_end])	
	
	if !isdir(v_opt_folder_name)
		mkdir(v_opt_folder_name)
	end
	
	f_movie, ax_movie, fnames_movie = [], [], []
end;

# ╔═╡ ba5049f0-b015-41ad-907e-b86edbcbf7fb
begin
	redoMovie = true
	if createGifImages
		rm.(fnames_movie)
		for i in eachindex(fnames_movie)
			pop!(fnames_movie)
			pop!(f_movie)
			pop!(ax_movie)
		end
		let 
			@progress for i in eachindex(snapshots)
				if (snapshots[i] > endTimeMovie) | (snapshots[i] < startTimeMovie)
					continue
				end
				
				plt.close()
				f, ax = plt.subplots(1, 1, figsize=(5, 6))
				
				x = xAxis
				y = yAxis
			
				extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
			
				im = if !logy
					ax.imshow(
						surfacesMovie[i]',
						origin="lower",
						extent=extent,
						cmap=cmap,
						vmin=v_min_movie,
						vmax=v_max_movie,
						aspect="equal"
					)
				else
					ax.imshow(
						surfacesMovie[i]',
						origin="lower",
						extent=extent,
						cmap=cmap,
						aspect="equal",
						norm=matplotlib.colors.LogNorm(vmin=v_min_movie, vmax=v_max_movie)
					)
				end
				cb = f.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
				cb.set_label(labelsurfaceMovie)
			
				
				ax.set_title("t = $(MUST.@sprintf("%.2f", time[i]/(60*60))) h")
				
				ax.set_xlabel("x [Mm]")
				ax.set_ylabel("y [Mm]")
				
				append!(f_movie, [f])
				append!(ax_movie, [ax])
				append!(fnames_movie, [joinpath(v_opt_folder_name, "im_$i.png")])
				f.savefig(joinpath(v_opt_folder_name, "im_$i.png"), dpi=dpi)
			end
		end
	end;
end

# ╔═╡ 7308fc51-5236-4b07-b773-df5451852fbb
(length(fnames_movie) > 0) && begin
	redoMovie
	
	v_images = Images.load.(fnames_movie)
	anim = @animate for i ∈ eachindex(v_images)
		Plots.plot(v_images[i], axis=([], false), background_color=:transparent)
	end every 1
	gif(anim, joinpath(v_opt_folder_name, "v_opt.gif"), fps=fps)
end

# ╔═╡ fe34b982-e39b-4c54-827f-53b27b2d6db2


# ╔═╡ 5d5886f6-565f-47d5-be94-7e7bf5b0e89d


# ╔═╡ bdf0fb65-7a95-471a-99dc-b2ae8ce9a97e
md"## Vertical yz, x-center slice"

# ╔═╡ a2718d44-abdc-47fd-8ad1-6ad2aece1718
vertMovieSelection = if "centerVerticalCut" in keys(monitoring[1])
	vertMovieSelection = Dict()
	tvc = timeevolution(monitoring, "centerVerticalCut")
	vertMovieSelection["temperature"] = SurfaceMovieContent(
		tvc["Tplane"], 
		L"\rm temperature\ [K]", 
		tvc["y"][1] ./1e8, 
		tvc["z"][1] ./1e8
	)
	vertMovieSelection["density"] = SurfaceMovieContent(
		[exp.(x) for x in tvc["lnDplane"]], 
		L"\rm density\ [g\ cm^{-3}]",
		tvc["y"][1] ./1e8, 
		tvc["z"][1] ./1e8
	)
	if haskey(tvc, "fluxplane")
		vertMovieSelection["flux"] = SurfaceMovieContent(
			tvc["fluxplane"], 
			L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]",
			tvc["y"][1] ./1e8, 
			tvc["z"][1] ./1e8
		)
	end
	vertMovieSelection["vertical velocity"] = SurfaceMovieContent(
		tvc["uzplane"] ./1e5, 
		L"\rm v_z\ [km\ s^{-1}]", 
		tvc["y"][1] ./1e8, 
		tvc["z"][1] ./1e8
	)
	if haskey(tvc, "dtplane")
		vertMovieSelection["timestep"] = SurfaceMovieContent(
			tvc["dtplane"], 
			L"\rm timestep\ [s]",
			tvc["y"][1] ./1e8, 
			tvc["z"][1] ./1e8
		)
	end
	if haskey(tvc, "qrplane")
		vertMovieSelection["heating"] = SurfaceMovieContent(
			tvc["qrplane"], 
			L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]",
			tvc["y"][1] ./1e8, 
			tvc["z"][1] ./1e8
		)
	end
	vertMovieSelection
else
	Dict("Temperature"=>nothing)
end;

# ╔═╡ 97a4f8d2-6c3b-48c0-a5ff-ac4f23cfcaf9
("centerVerticalCut" in keys(monitoring[1])) && md"""
	Select 2D vertical movie topic
$(@bind vertMovie_choice Select(sort(keys(vertMovieSelection)|>collect), default="temperature"))
"""

# ╔═╡ 7ccca88c-23f4-4440-8346-0d3cc25ebae9
if "centerVerticalCut" in keys(monitoring[1])
	vertMovieSelected = vertMovieSelection[vertMovie_choice]
	surfacesMovie_vert = vertMovieSelected.movie
	labelsurfaceMovie_vert = vertMovieSelected.label
	xAxis_vert = vertMovieSelected.xAxis
	yAxis_vert = vertMovieSelected.yAxis
end;

# ╔═╡ 5540fbd5-67b4-4d18-8d3b-b12f1249e2ec
md"""
	Modify movie appearence
"""

# ╔═╡ d0e41ab4-6d23-4fa2-89dd-5253013af2d2
md"""
figure resolution (DPI): $(@bind dpi_vert confirm(Slider(10:600, default=75, show_value=true)))\

colormap (matplotlib): $(@bind cmap_vert confirm(TextField(default=\"rainbow\")))\

log scale for color axis: $(@bind logy_vert CheckBox(default=true))
"""

# ╔═╡ 3f5f399f-e012-4a74-8e72-a6ed45c66b7a
if ("centerVerticalCut" in keys(monitoring[1]))
	v_min_movie_vert = minimum(minimum, surfacesMovie_vert[i_start:i_end])
	v_max_movie_vert = maximum(maximum, surfacesMovie_vert[i_start:i_end])
	f_movie_vert, ax_movie_vert, fnames_movie_vert = [], [], []
else
	f_movie_vert, ax_movie_vert, fnames_movie_vert = [], [], []
end;

# ╔═╡ fe3a6c8f-1135-479e-bf5f-00a7671abd3e
("centerVerticalCut" in keys(monitoring[1])) && begin
	redoMovie_vert = true
	if createGifImages
		for (i, fp) in enumerate(fnames_movie_vert)
			isfile(fp) && rm(fp)
			pop!(fnames_movie_vert)
			pop!(f_movie_vert)
			pop!(ax_movie_vert)
		end
		let 
			@progress for i in eachindex(snapshots)
				if (snapshots[i] > endTimeMovie) | (snapshots[i] < startTimeMovie)
					continue
				end
				
				plt.close()
				f, ax = plt.subplots(1, 1, figsize=(6, 4))
				
				x = xAxis_vert
				y = yAxis_vert
			
				extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

				im = if !logy_vert
					ax.imshow(
						surfacesMovie_vert[i]',
						origin="lower",
						extent=extent,
						cmap=cmap_vert,
						vmin=v_min_movie_vert,
						vmax=v_max_movie_vert,
						aspect="auto"
					)
				else
					ax.imshow(
						surfacesMovie_vert[i]',
						origin="lower",
						extent=extent,
						cmap=cmap_vert,
						aspect="auto",
						norm=matplotlib.colors.LogNorm(vmin=v_min_movie_vert, vmax=v_max_movie_vert)
					)
				end
				cb = f.colorbar(im, ax=ax)
				cb.set_label(labelsurfaceMovie_vert)
				
				ax.set_title("t = $(MUST.@sprintf("%.2f", time[i]/(60*60))) h")
				
				ax.set_xlabel("y [Mm]")
				ax.set_ylabel("z [Mm]")
				
				append!(f_movie_vert, [f])
				append!(ax_movie_vert, [ax])
				append!(fnames_movie_vert, [joinpath(v_opt_folder_name, "im_vert_$i.png")])
				f.savefig(joinpath(v_opt_folder_name, "im_vert_$i.png"), dpi=dpi_vert)
			end
		end
	end
end;

# ╔═╡ 0ddaddb6-784f-45ae-adc2-9217f0f52996
(length(fnames_movie_vert) > 0) && let
	redoMovie_vert
	
	v_images = Images.load.(fnames_movie_vert)
	anim = @animate for i ∈ eachindex(v_images)
		Plots.plot(v_images[i], axis=([], false), background_color=:transparent)
	end every 1
	gif(anim, joinpath(v_opt_folder_name, "v_opt_vert.gif"), fps=fps)
end

# ╔═╡ 401758b0-9d90-4e90-bb44-36337afe4425


# ╔═╡ 0580f869-e419-4ad3-890e-0844b10d7887


# ╔═╡ f2923fdf-d717-46e5-bd64-87cf599d9a11
md"## 1D stratification movie"

# ╔═╡ d07278c5-f95a-4e3a-9709-9a1aa52c1500
begin
	struct AverageMovieContent
		xmovie
		ymovie
		xlabel
		ylabel
	end
	averageMovieSelection = Dict()
	averageMovieSelection["geometrical - temperature"] = AverageMovieContent(
		[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
		[tGeoAv["T"], tGeoMin["T"], tGeoMax["T"]],
		L"\rm z\ [km]", 
		L"\rm temperature\ [K]", 
	)
	averageMovieSelection["geometrical - density"] = AverageMovieContent(
		[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
		[tGeoAv["d"], tGeoMin["d"], tGeoMax["d"]],
		L"\rm z\ [km]", 
		L"\rm density\ [g\ cm^{-3}]",
	)
	if haskey(tGeoMin, "flux")
		averageMovieSelection["geometrical - flux"] = AverageMovieContent(
			[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
			[tGeoAv["flux"], tGeoMin["flux"], tGeoMax["flux"]],
			L"\rm z\ [km]", 
			L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]",
		)
	end
	averageMovieSelection["geometrical - vertical velocity"] = AverageMovieContent(
		[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
		[tGeoAv["uz"], tGeoMin["uz"], tGeoMax["uz"]] ./1e5,
		L"\rm z\ [km]", 
		L"\rm v_z\ [km\ s^{-1}]", 
	)
	if haskey(tGeoMin, "dt_rt")
		averageMovieSelection["geometrical - timestep"] = AverageMovieContent(
			[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
			[tGeoAv["dt_rt"], tGeoMin["dt_rt"], tGeoMax["dt_rt"]],
			L"\rm z\ [km]", 
			L"\rm timestep\ [s]",
		)
	end
	if haskey(tGeoMin, "qr")
		averageMovieSelection["geometrical - heating"] = AverageMovieContent(
			[tGeoAv["z"], tGeoMin["z"], tGeoMax["z"]] ./1e8,
			[tGeoAv["qr"], tGeoMin["qr"], tGeoMax["qr"]],
			L"\rm z\ [km]", 
			L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]",
		)
	end


	averageMovieSelection["optical - temperature"] = AverageMovieContent(
		[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
		[tOptAv["T"], tOptMin["T"], tOptMax["T"]],
		L"\rm \log \tau_{ross}", 
		L"\rm temperature\ [K]", 
	)
	averageMovieSelection["optical - density"] = AverageMovieContent(
		[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
		[tOptAv["d"], tOptMin["d"], tOptMax["d"]],
		L"\rm \log \tau_{ross}", 
		L"\rm density\ [g\ cm^{-3}]",
	)
	if haskey(tGeoMin, "dt_rt")
		averageMovieSelection["optical - flux"] = AverageMovieContent(
			[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
			[tOptAv["flux"], tOptMin["flux"], tOptMax["flux"]],
			L"\rm \log \tau_{ross}", 
			L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]",
		)
	end
	averageMovieSelection["optical - vertical velocity"] = AverageMovieContent(
		[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
		[tOptAv["uz"], tOptMin["uz"], tOptMax["uz"]] ./1e5,
		L"\rm \log \tau_{ross}", 
		L"\rm v_z\ [km\ s^{-1}]", 
	)
	if haskey(tGeoMin, "dt_rt")
		averageMovieSelection["optical - timestep"] = AverageMovieContent(
			[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
			[tOptAv["dt_rt"], tOptMin["dt_rt"], tOptMax["dt_rt"]],
			L"\rm \log \tau_{ross}", 
			L"\rm timestep\ [s]",
		)
	end
	if haskey(tGeoMin, "qr")
		averageMovieSelection["optical - heating"] = AverageMovieContent(
			[tOptAv["log10τ_ross"], tOptMin["log10τ_ross"], tOptMax["log10τ_ross"]],
			[tOptAv["qr"], tOptMin["qr"], tOptMax["qr"]],
			L"\rm \log \tau_{ross}", 
			L"\rm Q_r\ [erg\ s^{-1}\ cm^{-3}]",
		)
	end
end;

# ╔═╡ 12767600-b8d3-4a46-bea2-8325f91da347
md"""
	Select 1D average movie topic
$(@bind averageMovie_choice Select(sort(keys(averageMovieSelection)|>collect), default="optical - temperature"))
"""

# ╔═╡ 6f82fd7d-8be1-4065-a0d8-2100cada740c
begin
	averageMovieSelected = averageMovieSelection[averageMovie_choice]
	xAxis_averageMovie = averageMovieSelected.xmovie
	yAxis_averageMovie = averageMovieSelected.ymovie
	xlabelaverageMovie_vert = averageMovieSelected.xlabel
	ylabelaverageMovie_vert = averageMovieSelected.ylabel
end;

# ╔═╡ 0d9f3118-7029-4c89-97a1-6f1989584a50
md"""
	Modify movie appearence
"""

# ╔═╡ ff6d80d9-aadb-41a4-9cd8-0c601bfea714
md"""
figure resolution (DPI): $(@bind dpi_averageMovie confirm(Slider(10:800, default=200, show_value=true)))\
log scale for vertical axis: $(@bind logy_averageMovie CheckBox(default=true))
"""

# ╔═╡ a926b689-8009-42f5-acf7-b35c5d5ed53f
begin
	color_averageMovie = ["black", "magenta", "cyan"]
	kwargs_averageMovie = [
		Dict(:ls=>"-", :lw=>2), Dict(:ls=>"", :lw=>2, :marker=>"."), Dict(:ls=>"", :lw=>2, :marker=>".")
	]
end;

# ╔═╡ be0fc7d6-0145-4f7a-97fb-ee43fa61c099
begin
	x_min_movie_averageMovie = Inf
	x_max_movie_averageMovie = -Inf
	y_min_movie_averageMovie = Inf
	y_max_movie_averageMovie = -Inf
	for (x, y) in zip(xAxis_averageMovie, yAxis_averageMovie)
		x_min_loc = minimum(minimum, x[i_start:i_end])
		x_max_loc = maximum(maximum, x[i_start:i_end])
		y_min_loc = minimum(minimum, y[i_start:i_end])
		y_max_loc = maximum(maximum, y[i_start:i_end])

		x_min_movie_averageMovie = min(x_min_loc, x_min_movie_averageMovie)
		x_max_movie_averageMovie = max(x_max_loc, x_max_movie_averageMovie)
		y_min_movie_averageMovie = min(y_min_loc, y_min_movie_averageMovie)
		y_max_movie_averageMovie = max(y_max_loc, y_max_movie_averageMovie)
	end
	
	f_movie_averageMovie, ax_movie_averageMovie, fnames_movie_averageMovie = [], [], []
end;

# ╔═╡ ea8df11b-7355-4564-bcfc-6271e8211ef1
begin
	redoMovie_averageMovie = true
	if createGifImages
		for (i, fp) in enumerate(fnames_movie_averageMovie)
			isfile(fp) && rm(fp)
			pop!(fnames_movie_averageMovie)
			pop!(f_movie_averageMovie)
			pop!(ax_movie_averageMovie)
		end
		let 
			@progress for i in eachindex(snapshots)
				if (snapshots[i] > endTimeMovie) | (snapshots[i] < startTimeMovie)
					continue
				end
				
				plt.close()
				f, ax = plt.subplots(1, 1, figsize=(6, 4))

				for j in eachindex(xAxis_averageMovie)
					x = xAxis_averageMovie[j][i]
					y = yAxis_averageMovie[j][i]

					im = ax.plot(
						x, y,
						color=color_averageMovie[j]; 
						kwargs_averageMovie[j]...
					)
				end
				
				ax.set_title("t = $(MUST.@sprintf("%.2f", time[i]/(60*60))) h")
				
				ax.set_xlabel(xlabelaverageMovie_vert)
				ax.set_ylabel(ylabelaverageMovie_vert)

				ax.set_xlim(x_min_movie_averageMovie, x_max_movie_averageMovie)
				ax.set_ylim(y_min_movie_averageMovie, y_max_movie_averageMovie)

				if logy_averageMovie
					ax.set_yscale("log")
				end
				
				append!(f_movie_averageMovie, [f])
				append!(ax_movie_averageMovie, [ax])
				append!(fnames_movie_averageMovie, [joinpath(v_opt_folder_name, "im_averageMovie_$i.png")])
				f.savefig(joinpath(v_opt_folder_name, "im_averageMovie_$i.png"), dpi=dpi_averageMovie)
			end
		end
	end
end;

# ╔═╡ 260616da-d48e-4a7e-9a2b-1923c16408be
(length(fnames_movie_averageMovie) > 0) && let
	redoMovie_averageMovie
	
	v_images = Images.load.(fnames_movie_averageMovie)
	anim = @animate for i ∈ eachindex(v_images)
		Plots.plot(v_images[i], axis=([], false), background_color=:transparent)
	end every 1
	gif(anim, joinpath(v_opt_folder_name, "v_opt_averageMovie.gif"), fps=fps)
end

# ╔═╡ Cell order:
# ╟─c7dc3b15-6555-4824-872a-d487fe5145ea
# ╟─c1ad6c59-fae6-49e7-bfc0-8426c553aa2d
# ╠═2f1edd2a-b56e-11ee-29e7-c353938e7088
# ╟─7cd7d6f0-8498-44ff-b59c-d298365d6416
# ╠═78c88a26-1e84-4ba2-a8f2-d4c7f6468dd3
# ╠═e2e0b39b-9c60-4630-817b-f180c2631a08
# ╟─409ff57f-8d9d-419b-b448-fdf40c0843b4
# ╟─6754b2c3-d205-4a12-88b3-53fe62c5637f
# ╟─41f0864e-26ee-46a6-b4ab-c401a4712941
# ╟─c596d1b3-32c5-4651-a8c1-3100fcd6cd59
# ╟─8af0c339-237c-42ca-bda5-0b0e44f11c30
# ╟─498f998c-9996-43dc-a647-3b08abf9786d
# ╟─3c05fe1e-1c19-4c30-b34e-a874f63b91bc
# ╟─ed9cc79f-0161-4178-b6f0-81a2bccbf188
# ╟─28d708db-45a7-4fad-a226-e907ebf88b43
# ╟─6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
# ╟─1fcefd1e-1c50-43b5-b203-62920944344a
# ╟─60199001-8f0e-44a6-ae50-e829687c045c
# ╟─0e4df1b3-52ef-4fce-8f32-bddc966b0516
# ╟─3aadeb1c-3585-46af-8b9d-daf33bfcb7f3
# ╟─2d580186-96a6-41b8-91aa-da527691ea1d
# ╟─e3f00981-eb96-49ba-8451-9af7560d3556
# ╟─73df16bf-b033-494d-b321-608ffff33467
# ╟─bd936d7d-e79f-4f9b-ba54-e0694c6a83f0
# ╟─2c64fcf2-1a0b-49cf-a3f1-890f152d0650
# ╟─65c9d5ea-45a2-4ab8-99e8-5eee29935589
# ╟─0f9a1524-9c47-49d6-a554-66db88663093
# ╟─cf85b913-bfd0-4a52-b363-7a70b8675134
# ╟─835470b0-73d6-4c15-8d13-46abe85e2ed3
# ╟─2ebb80e6-1966-4673-bc79-50d17add3969
# ╟─9265061d-eaa5-4fc1-a7b9-51392a357c91
# ╟─63b0d9d6-f27a-492e-9018-876db8091914
# ╟─c7bb05cd-0798-4989-94de-f7238fea8d89
# ╟─b76f0efd-6b2f-47d3-84e7-ca9fb866292e
# ╟─4eef1012-c8cd-4579-9a7e-ac2b4d4c7af6
# ╟─33c1da97-0760-495e-abd5-65531d5e1170
# ╟─8d6d674b-153b-4357-9f2d-c4e3cb05d059
# ╟─b9a721cf-46ef-4e3c-a37c-8b35653e31cb
# ╟─725f7500-c31b-4efa-9df9-00c51640914a
# ╟─a34793ae-db12-4dac-b2f8-348a88092815
# ╟─63f58f6c-954c-44e9-84b1-1aa799323586
# ╟─a54433ec-14b4-4e5e-a0a7-2e6e50a5fa49
# ╟─695e28be-31c8-4f44-85e3-72ce387d9da5
# ╟─76a1714e-a5bb-488f-ad93-b8552e4531fd
# ╟─f5ccc550-6bf3-4f0f-baf3-d8de3d216737
# ╟─b757013d-aee5-41a4-ab0a-7715ba47bd97
# ╟─0639ce7d-955d-448f-84a0-353dfd4e93a3
# ╟─c24a3b12-c7c0-449b-9a76-6e9c5d475344
# ╟─3d371088-2322-462b-ab93-7cb49fcdf75f
# ╟─9758eaf6-9d57-4847-a608-2ba81e1192d0
# ╟─80cf7a65-87c6-48d1-89d9-90f07ec25d51
# ╟─5817821d-67f5-4e41-a0d3-7ca12961b0c7
# ╟─7f77f259-505d-4344-8ee4-8628387f2401
# ╟─cc5fbd5a-c8a0-471a-a56b-0512e4c3989b
# ╟─495e3733-d290-40ab-af63-0eb10a033b53
# ╟─ed29d53f-00bc-4295-93f6-864a44f92ccb
# ╟─fa8161aa-0e1c-404f-8c7e-0e3914917df4
# ╟─ef776fe2-5f0f-49a8-b2d5-2e1235d41ec1
# ╟─ff3d827b-a0dd-4d86-bfd4-d25fc5f4ab3e
# ╟─0e3d2723-1ecb-4e5a-8125-78f6f27407e2
# ╟─3403a014-2441-4ad2-95f6-2e686ae99ba8
# ╟─77f5cb10-ce91-44c2-91c6-c01432655121
# ╟─7532d800-c183-42d5-8bd0-9970ca507cfd
# ╟─9504d2ae-a060-4b02-a873-27f4c01d3364
# ╟─ab5e5713-7dde-4592-9371-88e2a9c00400
# ╟─2b627607-4d4f-447f-a29a-daacf7977c28
# ╟─b5a51746-f8e3-4760-b73a-a760b2e27c3c
# ╟─40b66681-4ee6-485c-a94f-6d09cf4d56e6
# ╟─97e5e57c-3828-4781-b93e-172a1c3c2cce
# ╟─82d82155-981c-4fdb-b0ae-d8d40ce016f5
# ╟─b0c40c50-3361-4b01-ae87-45ae30387526
# ╟─159019c1-b183-4d25-83a9-0ae27d9eb0ee
# ╟─29803b42-3f55-4fb1-8853-8848bbdd5bd9
# ╟─68477423-a6f7-424f-8fd4-849d63648b57
# ╟─8ef763ae-d1bd-40de-a26a-f91f529c03bf
# ╟─d9546636-2dbb-4b14-9954-76872b95fd06
# ╟─7fd50fcd-8a9d-417d-a28a-febd9f0f68f3
# ╟─cdfb5217-e298-4ee1-ac5d-f294d04f5abc
# ╟─274d42e7-20eb-4560-a0f0-e4a631aedcec
# ╟─46db420f-2f46-4989-846e-145a70001ddf
# ╟─b046b3a0-05d4-4b6c-be80-7396f475be3d
# ╟─b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5
# ╟─db1a9405-e93c-476d-a43b-f11f3138b57a
# ╟─a575dee8-6283-469b-a6ad-365ca9d1b5e9
# ╟─f5c4d377-09f2-4ecf-9cd6-458c10be2943
# ╟─0fdf3055-4d11-4dea-8a50-e595ef1c112d
# ╟─eddc02bf-d7ca-41e1-878e-ef1103bf1b0f
# ╟─25c3d608-9440-4ac9-8277-3855ba3b6a7b
# ╟─e59ab649-3b30-4af1-ab28-6f6d6aacee1c
# ╟─23c0607b-d70e-4bbd-85ba-7f2b322d93dc
# ╟─0da614f4-9021-44f3-af89-bd6dab57dc1b
# ╟─d55d7d42-c78c-447c-9959-3689f5341655
# ╟─c3fb528f-4bd8-4ae7-bb4c-ac90899fbf21
# ╟─10e5e8cb-0881-40d0-b392-2f66c0cbfc7c
# ╟─4c3b9b8a-03a5-494d-a321-7f5c400ed054
# ╟─dc0315d8-0616-433b-8182-33840afb0b0f
# ╟─aef86067-68b0-48b8-8bb2-0d410a7521c2
# ╟─321e3dda-cd15-4787-95e6-f928125535d5
# ╟─643eaee3-a1f8-4be2-b032-08fcb325fbe8
# ╟─aae74666-7fdd-4db9-8869-68e4597f6940
# ╟─15f1108e-6f72-4700-89d9-fc2dfaf7ea47
# ╟─e452c94a-3b4b-4dd1-a286-3c32a5c8e995
# ╟─c370ebd5-fc25-4e81-b3f9-37bf99b35598
# ╟─fc550522-d506-4f0e-b582-e7c54cd23be0
# ╟─35f64e1d-2273-4178-879e-187b86b24043
# ╟─cb3e21b3-0197-4685-a9ae-fe4040106222
# ╟─eb335a4d-e662-499c-bb80-8bf38c84329f
# ╟─35b4aa7e-e41b-4444-af3b-74684c723d5c
# ╟─87d43815-9733-40ac-b4e7-81a2c5dbd0d1
# ╟─c8492ca5-893f-4939-a256-f0872f351c5d
# ╟─df9b58cf-fe4f-4858-a296-879bb0385ba7
# ╟─ccf5f4e6-adc6-419d-a717-4b3b597c2233
# ╟─913c1cbe-fedc-4e7d-a8a7-c2ad416a21e6
# ╟─2f76eaac-0793-4fa8-9ee1-67dc81e9a1ac
# ╟─aa0376dc-5982-4274-9b4d-4aef1d1de896
# ╟─0378dcdf-0291-4e9c-a30b-5dacd346707c
# ╟─b4fd0320-74b2-4ed4-b7b0-beca4eccd553
# ╟─9e883b44-f225-4566-9d76-ecd3e34d3b5b
# ╟─def11de6-b7af-446b-ada7-e61f19b6503b
# ╟─01026bad-d93a-4246-a940-e9eb894ccf67
# ╟─dd1297cd-62d0-4b39-b0a3-e41fab79bc8b
# ╟─43590370-0a2e-4d6e-9531-f93eb1630acd
# ╟─cd9ff780-66d2-4bcb-8b89-145d9d857306
# ╟─82daeb7f-1a3b-420a-8193-148d7314f557
# ╟─b59908ac-50a9-49fd-bcc2-94c59b205543
# ╟─787b1cc4-9e1f-4421-874d-bdff7c231bf9
# ╟─8ae48e0a-cddf-4775-9a1a-8eada0c07ca7
# ╟─43b11d24-2496-44f2-9e9f-fc30ff8d285d
# ╟─dbeee710-c9f3-4561-969a-321f256f99f5
# ╟─408f2044-cac0-48b2-8b33-f7ddf9b66849
# ╟─7a1d3d74-fea2-402b-aee1-728ff89eb548
# ╟─e9b7dfd5-1034-4b95-941c-df003436e26e
# ╟─39d3636d-1161-4cbd-a11e-b900e75e77ce
# ╟─c5550bd5-725d-4d3e-9990-324860af7f67
# ╟─43ab1632-3df6-4b76-9385-627452ce97aa
# ╟─587cd5a6-2117-4a65-b50e-c8d2107adfd9
# ╟─47f967fb-e063-4b23-97e9-02ca76985fa9
# ╟─3d530773-bd52-41d6-a3dd-509f237ae439
# ╟─61b3e522-2af4-4413-9c69-ee30b0d4673f
# ╟─713164e6-4a7a-4475-beda-bda33d16e206
# ╟─fb94a874-a34b-4721-be63-e443f92330af
# ╟─2e99b941-935b-4379-89bf-55a4635df686
# ╟─48b4180b-6891-4c2a-a22e-7761b8e1f439
# ╟─9d0348a6-24a1-447c-9af4-5b648a381370
# ╟─c30611f4-1d2e-41c9-b8db-f85b3c99990d
# ╟─4564cbd9-34bb-4b21-9fca-1027d5c95518
# ╟─40bf37a4-318b-4263-8705-a566a3321f87
# ╟─c1faa953-4b76-43fa-a958-a5476e325afc
# ╟─b901bb2a-1f3c-46f5-b3f7-1d530975d767
# ╟─3c3c8847-ae3e-40d5-97ef-d69e3e5ccd14
# ╟─6575c66a-fe2a-4427-bd4c-f23eef344caf
# ╟─5646dc57-693f-4605-9edc-46b0648f3ab5
# ╟─8945b7b7-f3e5-4846-b76e-49d9c864391e
# ╟─ba5049f0-b015-41ad-907e-b86edbcbf7fb
# ╟─7308fc51-5236-4b07-b773-df5451852fbb
# ╟─fe34b982-e39b-4c54-827f-53b27b2d6db2
# ╟─5d5886f6-565f-47d5-be94-7e7bf5b0e89d
# ╟─bdf0fb65-7a95-471a-99dc-b2ae8ce9a97e
# ╟─a2718d44-abdc-47fd-8ad1-6ad2aece1718
# ╟─97a4f8d2-6c3b-48c0-a5ff-ac4f23cfcaf9
# ╟─7ccca88c-23f4-4440-8346-0d3cc25ebae9
# ╟─5540fbd5-67b4-4d18-8d3b-b12f1249e2ec
# ╟─d0e41ab4-6d23-4fa2-89dd-5253013af2d2
# ╟─3f5f399f-e012-4a74-8e72-a6ed45c66b7a
# ╟─fe3a6c8f-1135-479e-bf5f-00a7671abd3e
# ╟─0ddaddb6-784f-45ae-adc2-9217f0f52996
# ╟─401758b0-9d90-4e90-bb44-36337afe4425
# ╟─0580f869-e419-4ad3-890e-0844b10d7887
# ╟─f2923fdf-d717-46e5-bd64-87cf599d9a11
# ╟─d07278c5-f95a-4e3a-9709-9a1aa52c1500
# ╟─12767600-b8d3-4a46-bea2-8325f91da347
# ╟─6f82fd7d-8be1-4065-a0d8-2100cada740c
# ╟─0d9f3118-7029-4c89-97a1-6f1989584a50
# ╟─ff6d80d9-aadb-41a4-9cd8-0c601bfea714
# ╠═a926b689-8009-42f5-acf7-b35c5d5ed53f
# ╟─be0fc7d6-0145-4f7a-97fb-ee43fa61c099
# ╟─ea8df11b-7355-4564-bcfc-6271e8211ef1
# ╟─260616da-d48e-4a7e-9a2b-1923c16408be
