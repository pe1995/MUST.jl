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
md"## Load monitoring"

# ╔═╡ 28d708db-45a7-4fad-a226-e907ebf88b43
wd = MUST.WatchDog(selectedRun, folder=datafolder)

# ╔═╡ 6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
monitoring = MUST.reload!(wd)

# ╔═╡ 62e93fea-6600-4619-bb66-5908be8c26e3


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
md"Pick the time for which you want to see the status: $(@bind timeSurface Slider(snapshots, show_value=true, default=last(snapshots)))
"

# ╔═╡ 725f7500-c31b-4efa-9df9-00c51640914a
md"Pick a second time for comparison: $(@bind timeSurface2 Slider(snapshots, show_value=true))
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
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")

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
		cmap="gist_heat"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]")
	
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


# ╔═╡ b0c40c50-3361-4b01-ae87-45ae30387526
md"## Geometrical Profiles"

# ╔═╡ 96a3ddc0-3d70-4b65-a7c2-7482c8817186
tGeoAv = timeevolution(monitoring, "geometricalAverages")

# ╔═╡ 8c03ccae-6a25-4f2b-a99d-bf537221d007
tGeoMin = timeevolution(monitoring, "geometricalMinimum")

# ╔═╡ 2f86b88f-b53a-489b-8a1e-7d5116492e34
tGeoMax = timeevolution(monitoring, "geometricalMaximum")

# ╔═╡ f3002fa3-b796-4479-8008-fc8149e797bf
tGeoRms = timeevolution(monitoring, "geometricalRMS")

# ╔═╡ 32603241-ed2c-4579-9319-40e50eb2172c
tGeoQuantiles = Dict(
	k=>timeevolution(monitoring, "geometrical$(k)thQuantile") 
	for k in [15, 30, 45, 60, 75, 90] 
	if "geometrical$(k)thQuantile" in keys(monitoring[1])
)

# ╔═╡ 53bed7de-b834-465e-8d13-2909240286b5


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

# ╔═╡ 3da231de-ff58-4161-a6a4-58162483825a
tOptAv = timeevolution(monitoring, "opticalAverages")

# ╔═╡ efc11ddf-99b9-4344-ab4f-1e5d7b7e2805
tOptMin = timeevolution(monitoring, "opticalMinimum")

# ╔═╡ e0d8aa2e-85cb-43fe-a5f4-65a7a757f19c
tOptMax = timeevolution(monitoring, "opticalMaximum")

# ╔═╡ 699ef294-a1b0-4ea5-8618-5cfd0b143c9a
tOptRms = timeevolution(monitoring, "opticalRMS")

# ╔═╡ 46944b1e-1590-4cba-9b9a-078f34d987ad


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
md"$(@bind snapshotSelector Slider(snapshots, show_value=true, default=first(snapshots)))"

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

# ╔═╡ c2b64d82-09a7-4d67-97d0-b51d96d30d25
ttempsurface = timeevolution(monitoring, "opticalSurfaces", "Tplane")

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

# ╔═╡ 71c4c81a-5567-45f0-925a-1dcad0f97082
tdsurface = timeevolution(monitoring, "opticalSurfaces", "lnDplane")

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

# ╔═╡ 1521d2b2-56a7-41eb-ad52-b661af7f30c6
tuzsurface = timeevolution(monitoring, "opticalSurfaces", "uzplane")

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

# ╔═╡ 3bfef187-d899-473d-9a92-acf5c65fa50d
tmassfluxgeo = timeevolution(monitoring, "geoMassFlux")

# ╔═╡ 786e489d-42fb-4ef0-83b6-0e6a7fce2d93
tmassfluxopt = timeevolution(monitoring, "optMassFlux")

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

# ╔═╡ 58b167bf-d79d-48e1-8118-ffc1a13ba913
begin
	iend_fft = length(time)
	istart_fft = iend_fft - 100
end

# ╔═╡ 47f967fb-e063-4b23-97e9-02ca76985fa9
let 
	if istart_fft>1
		data = topticalsurfaces["uzplane"][istart_fft:iend_fft] ./ 1e5
		data = [mean(d) for d in data]
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

		first_half = floor(Int, length(x) /2) 
		ax.plot(x[first_half:end], y[first_half:end], color="k")
		ax.set_xlabel(L"\rm frequency\ [\mu Hz]")
		ax.set_ylabel(L"\rm fft(<v_z^{\tau=1}>)\ [normalized]")

		ax.set_xlim(0, last(x))
		f
	end
end

# ╔═╡ 3d530773-bd52-41d6-a3dd-509f237ae439
begin
	νmax_expected = exp10(logg) / exp10(4.44) * (5777.0 / teff(tGeoAv["flux"][itimeSurface2][end]))^0.5 * 3090
	@info "Expected νmax from parameters [μHz]:" νmax_expected
end

# ╔═╡ f6916a12-cb30-4fa1-8d24-e75d1729ede2


# ╔═╡ 61b3e522-2af4-4413-9c69-ee30b0d4673f
md"# Movies"

# ╔═╡ 40bf37a4-318b-4263-8705-a566a3321f87
md"
## 2D Surface Movie
You can select any of the 2D plane monitoring results you wish, as well as the image resolution and FPS of the final movie."

# ╔═╡ 362e608f-4765-4f03-85f9-0d4673f6bd2b


# ╔═╡ 00876598-ab41-4741-9cb8-1538b0bd01a0
fps = 3

# ╔═╡ 1483a522-8d36-45f3-a7bc-8f8f8076085f
begin
	#=
	surfacesMovie = topticalsurfaces["Tplane"]
	surfacesMovie = topticalsurfaces["fluxplane"]
	surfacesMovie = topticalsurfaces["uzplane"] ./1e5
	surfacesMovie = tuppersurfaces["lnDplane"]
	surfacesMovie = tuppersurfaces["uzplane"] ./1e5
	surfacesMovie = tuppersurfaces["Tplane"] 
	surfacesMovie = tuppersurfaces["dtplane"]
	surfacesMovie = tuppersurfaces["fluxplane"] 
	surfacesMovie = timeevolution(monitoring, "minimumTempSurface")["uzplane"]
	=#
	surfacesMovie = topticalsurfaces["fluxplane"]

	#=
	labelsurfaceMovie = L"\rm v_z\ [km\ s^{-1}]"
	labelsurfaceMovie = L"\rm temperature\ [K]"
	labelsurfaceMovie = L"\rm density\ [g\ cm^{-3}]"
	labelsurfaceMovie = L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]"
	=#
	labelsurfaceMovie = L"\rm F_{bol}\ [erg\ s^{-1}\ cm^{-2}]"
	
	dpi = 72
	cmap = "jet"
end;

# ╔═╡ a84007ca-49d0-4558-ace3-cc42da9cb8eb


# ╔═╡ 98cfdafa-9078-4d0e-84d5-fcf97ab53a55
md"Pick the snapshot from which you want to start: $(@bind startTimeMovie Slider(snapshots, show_value=true, default=first(snapshots)))
"

# ╔═╡ 20feed40-879d-4786-93d3-e6d05b72fefc
md"Pick the snapshot on which you want to end: $(@bind endTimeMovie Slider(snapshots, show_value=true, default=last(snapshots)))
"

# ╔═╡ 8945b7b7-f3e5-4846-b76e-49d9c864391e
begin		
	i_start = findfirst(i->i==startTimeMovie, snapshots)
	i_end = findfirst(i->i==endTimeMovie, snapshots)
	v_min_movie = minimum(minimum, surfacesMovie[i_start:i_end])
	v_max_movie = maximum(maximum, surfacesMovie[i_start:i_end])	
	v_opt_folder_name = "v_opt"
	
	if isdir(v_opt_folder_name)
		rm(v_opt_folder_name, recursive=true)
	end
	mkdir(v_opt_folder_name)
	
	f_movie, ax_movie, fnames_movie = [], [], []
end;

# ╔═╡ 2e865c16-f7ae-4cc8-bcb8-5e4e13aabb9a
md"""
__Click to create GIF__: $(@bind createGifImages CheckBox(default=false))
"""

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
				
				x = topticalsurfaces["x"][i] ./1e8
				y = topticalsurfaces["y"][i] ./1e8
			
				extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
			
				im = ax.imshow(
					surfacesMovie[i],
					origin="lower",
					extent=extent,
					cmap=cmap,
					vmin=v_min_movie,
					vmax=v_max_movie
				)
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
# ╟─62e93fea-6600-4619-bb66-5908be8c26e3
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
# ╟─5817821d-67f5-4e41-a0d3-7ca12961b0c7
# ╟─7f77f259-505d-4344-8ee4-8628387f2401
# ╟─cc5fbd5a-c8a0-471a-a56b-0512e4c3989b
# ╟─495e3733-d290-40ab-af63-0eb10a033b53
# ╟─ed29d53f-00bc-4295-93f6-864a44f92ccb
# ╟─fa8161aa-0e1c-404f-8c7e-0e3914917df4
# ╟─ef776fe2-5f0f-49a8-b2d5-2e1235d41ec1
# ╟─0e3d2723-1ecb-4e5a-8125-78f6f27407e2
# ╟─3403a014-2441-4ad2-95f6-2e686ae99ba8
# ╟─b0c40c50-3361-4b01-ae87-45ae30387526
# ╟─96a3ddc0-3d70-4b65-a7c2-7482c8817186
# ╟─8c03ccae-6a25-4f2b-a99d-bf537221d007
# ╟─2f86b88f-b53a-489b-8a1e-7d5116492e34
# ╟─f3002fa3-b796-4479-8008-fc8149e797bf
# ╟─32603241-ed2c-4579-9319-40e50eb2172c
# ╟─53bed7de-b834-465e-8d13-2909240286b5
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
# ╟─3da231de-ff58-4161-a6a4-58162483825a
# ╟─efc11ddf-99b9-4344-ab4f-1e5d7b7e2805
# ╟─e0d8aa2e-85cb-43fe-a5f4-65a7a757f19c
# ╟─699ef294-a1b0-4ea5-8618-5cfd0b143c9a
# ╟─46944b1e-1590-4cba-9b9a-078f34d987ad
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
# ╠═4c3b9b8a-03a5-494d-a321-7f5c400ed054
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
# ╟─c2b64d82-09a7-4d67-97d0-b51d96d30d25
# ╟─eb335a4d-e662-499c-bb80-8bf38c84329f
# ╟─71c4c81a-5567-45f0-925a-1dcad0f97082
# ╟─35b4aa7e-e41b-4444-af3b-74684c723d5c
# ╟─1521d2b2-56a7-41eb-ad52-b661af7f30c6
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
# ╟─3bfef187-d899-473d-9a92-acf5c65fa50d
# ╟─786e489d-42fb-4ef0-83b6-0e6a7fce2d93
# ╟─01026bad-d93a-4246-a940-e9eb894ccf67
# ╟─dd1297cd-62d0-4b39-b0a3-e41fab79bc8b
# ╟─43590370-0a2e-4d6e-9531-f93eb1630acd
# ╟─cd9ff780-66d2-4bcb-8b89-145d9d857306
# ╟─82daeb7f-1a3b-420a-8193-148d7314f557
# ╟─b59908ac-50a9-49fd-bcc2-94c59b205543
# ╟─787b1cc4-9e1f-4421-874d-bdff7c231bf9
# ╟─8ae48e0a-cddf-4775-9a1a-8eada0c07ca7
# ╟─43b11d24-2496-44f2-9e9f-fc30ff8d285d
# ╟─c5550bd5-725d-4d3e-9990-324860af7f67
# ╟─43ab1632-3df6-4b76-9385-627452ce97aa
# ╟─587cd5a6-2117-4a65-b50e-c8d2107adfd9
# ╠═58b167bf-d79d-48e1-8118-ffc1a13ba913
# ╟─47f967fb-e063-4b23-97e9-02ca76985fa9
# ╟─3d530773-bd52-41d6-a3dd-509f237ae439
# ╟─f6916a12-cb30-4fa1-8d24-e75d1729ede2
# ╟─61b3e522-2af4-4413-9c69-ee30b0d4673f
# ╟─40bf37a4-318b-4263-8705-a566a3321f87
# ╟─362e608f-4765-4f03-85f9-0d4673f6bd2b
# ╠═00876598-ab41-4741-9cb8-1538b0bd01a0
# ╠═1483a522-8d36-45f3-a7bc-8f8f8076085f
# ╟─a84007ca-49d0-4558-ace3-cc42da9cb8eb
# ╟─98cfdafa-9078-4d0e-84d5-fcf97ab53a55
# ╟─20feed40-879d-4786-93d3-e6d05b72fefc
# ╟─8945b7b7-f3e5-4846-b76e-49d9c864391e
# ╟─2e865c16-f7ae-4cc8-bcb8-5e4e13aabb9a
# ╟─ba5049f0-b015-41ad-907e-b86edbcbf7fb
# ╟─7308fc51-5236-4b07-b773-df5451852fbb
