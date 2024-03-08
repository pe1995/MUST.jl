### A Pluto.jl notebook ###
# v0.19.37

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
	using PlutoUI
end

# ╔═╡ c7dc3b15-6555-4824-872a-d487fe5145ea
md"# Dispatch Monitoring Board"

# ╔═╡ c1ad6c59-fae6-49e7-bfc0-8426c553aa2d
md"## Setup"

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

# ╔═╡ 6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
monitoring = MUST.reload(MUST.WatchDog, selectedRun, folder=datafolder)

# ╔═╡ 1fcefd1e-1c50-43b5-b203-62920944344a
timeevolution(m, group, field) = [
	moni[group][field] 
	for moni in m
]

# ╔═╡ e752753c-b982-40c6-b45d-ad9bfcd3bbfe
timeevolution(m, group) = begin
	mg = m[1][group]
	Dict(
		k => [moni[group][k] for moni in m]
		for k in keys(mg)
	)
end

# ╔═╡ bcb4fd58-2b60-4b24-ad3e-ebc68392f320


# ╔═╡ e5b373f5-f565-4910-88b0-f5580880fac4
md"Available Groups:"

# ╔═╡ 53f4fdc7-5509-44dd-bf71-ceacb78e3e54
keys(monitoring[1]) |> collect

# ╔═╡ 60199001-8f0e-44a6-ae50-e829687c045c


# ╔═╡ 0e4df1b3-52ef-4fce-8f32-bddc966b0516
md"## Deleting Snapshots"

# ╔═╡ 3aadeb1c-3585-46af-8b9d-daf33bfcb7f3
snapshotnames(moni) = [nothing, timeevolution(moni, "general", "snapshot")...];

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
md"## General Properties"

# ╔═╡ 2c64fcf2-1a0b-49cf-a3f1-890f152d0650
time = timeevolution(monitoring, "atmosphericParameters", "time")

# ╔═╡ ef5ee4f6-96be-4669-ba68-3c1117605a4c


# ╔═╡ 8d6d674b-153b-4357-9f2d-c4e3cb05d059
md"## Individual Snapshots"

# ╔═╡ b9a721cf-46ef-4e3c-a37c-8b35653e31cb
md"Pick the time for which you want to see the status: $(@bind timeSurface Slider(time, show_value=true, default=last(time)))
"

# ╔═╡ 725f7500-c31b-4efa-9df9-00c51640914a
md"Pick a second time for comparison: $(@bind timeSurface2 Slider(time, show_value=true))
"

# ╔═╡ a34793ae-db12-4dac-b2f8-348a88092815
itimeSurface = findfirst(time .== timeSurface);

# ╔═╡ 63f58f6c-954c-44e9-84b1-1aa799323586
itimeSurface2 = findfirst(time .== timeSurface2);

# ╔═╡ a54433ec-14b4-4e5e-a0a7-2e6e50a5fa49


# ╔═╡ 695e28be-31c8-4f44-85e3-72ce387d9da5
@info "You are looking at snapshots (A) $(timeevolution(monitoring, "general", "snapshot")[itimeSurface]) and (B) $(timeevolution(monitoring, "general", "snapshot")[itimeSurface2])"

# ╔═╡ 76a1714e-a5bb-488f-ad93-b8552e4531fd


# ╔═╡ f5ccc550-6bf3-4f0f-baf3-d8de3d216737
md"### Optical Surface"

# ╔═╡ b757013d-aee5-41a4-ab0a-7715ba47bd97
topticalsurfaces = timeevolution(monitoring, "opticalSurfaces")

# ╔═╡ 09b4edcd-d333-42a2-b8e7-1814e85b446a


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
		cmap="coolwarm"
	)
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)


	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[1].imshow(
		topticalsurfaces["uzplane"][itimeSurface2] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	
	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")


	ax[0].set_xlabel("x [cm]")
	ax[1].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")
	
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
		exp.(topticalsurfaces["lnDplane"][itimeSurface]),
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		exp.(topticalsurfaces["lnDplane"][itimeSurface2]),
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \rho\ [g\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [cm]")
	
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
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")



	x = topticalsurfaces["x"][itimeSurface2] ./1e8
	y = topticalsurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		topticalsurfaces["Tplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [cm]")
	
	gcf()
end

# ╔═╡ 5817821d-67f5-4e41-a0d3-7ca12961b0c7


# ╔═╡ 7f77f259-505d-4344-8ee4-8628387f2401
md"### Upper Boundary Surface"

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
		cmap="coolwarm"
	)
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)


	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]

	i = ax[1].imshow(
		tuppersurfaces["uzplane"][itimeSurface2] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")


	ax[0].set_xlabel("x [cm]")
	ax[1].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")
	
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
		tuppersurfaces["lnDplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tuppersurfaces["lnDplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \rho\ [g\ cm^{-3}]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [cm]")
	
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
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[0], fraction=0.046, pad=0.04)

	ax[0].set_xlabel("x [cm]")
	ax[0].set_ylabel("y [cm]")



	x = tuppersurfaces["x"][itimeSurface2] ./1e8
	y = tuppersurfaces["y"][itimeSurface2] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax[1].imshow(
		tuppersurfaces["Tplane"][itimeSurface2],
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax[1], fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax[0].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface])) s")
	ax[1].set_title("t = $(MUST.@sprintf("%i", time[itimeSurface2])) s")

	ax[1].set_xlabel("x [cm]")
	
	gcf()
end

# ╔═╡ 3403a014-2441-4ad2-95f6-2e686ae99ba8


# ╔═╡ b0c40c50-3361-4b01-ae87-45ae30387526
md"### Geometrical Profiles"

# ╔═╡ 96a3ddc0-3d70-4b65-a7c2-7482c8817186
tGeoAv = timeevolution(monitoring, "geometricalAverages")

# ╔═╡ 8c03ccae-6a25-4f2b-a99d-bf537221d007
tGeoMin = timeevolution(monitoring, "geometricalMinimum")

# ╔═╡ 2f86b88f-b53a-489b-8a1e-7d5116492e34
tGeoMax = timeevolution(monitoring, "geometricalMaximum")

# ╔═╡ f3002fa3-b796-4479-8008-fc8149e797bf
tGeoRms = timeevolution(monitoring, "geometricalRMS")

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

	ax.set_xlabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
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

# ╔═╡ b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5


# ╔═╡ db1a9405-e93c-476d-a43b-f11f3138b57a
md"### Optical Profiles"

# ╔═╡ 3da231de-ff58-4161-a6a4-58162483825a
tOptAv = timeevolution(monitoring, "opticalAverages")

# ╔═╡ efc11ddf-99b9-4344-ab4f-1e5d7b7e2805
tOptMin = timeevolution(monitoring, "opticalMinimum")

# ╔═╡ e0d8aa2e-85cb-43fe-a5f4-65a7a757f19c
tOptMax = timeevolution(monitoring, "opticalMaximum")

# ╔═╡ 699ef294-a1b0-4ea5-8618-5cfd0b143c9a
tOptRms = timeevolution(monitoring, "opticalRMS")

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

# ╔═╡ d55d7d42-c78c-447c-9959-3689f5341655


# ╔═╡ 321e3dda-cd15-4787-95e6-f928125535d5
md"## Time evolution"

# ╔═╡ 35f64e1d-2273-4178-879e-187b86b24043
md"### Optical Surace"

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
		time, surfaceMax, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceMin, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceTemp, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <T_{\tau=1}>\ [K]")

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
		time, surfaceMax, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceMin, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceRho, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 	

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm \log < \rho_{\tau=1}>\ [g\ cm^{-3}]")

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
		time, surfaceMax ./1e5, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceMin ./1e5, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceUz ./1e5, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm rms\ v_z^{\tau=1}\ [km\ s^{-1}]")

	gcf()
end

# ╔═╡ c8492ca5-893f-4939-a256-f0872f351c5d


# ╔═╡ df9b58cf-fe4f-4858-a296-879bb0385ba7
md"### Upper Boundary"

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
		time, surfaceMi, 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceMa, 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, surfaceT, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_title("Upper Boundary")
	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <T_{top}>\ [K]")

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
		time, log10.(surfaceMi), 
		color="magenta", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, log10.(surfaceMa), 
		color="cyan", marker="s", markerfacecolor="w", ls="-"
	) 
	ax.plot(
		time, log10.(surfaceT), 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.set_title("Upper Boundary")
	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm  <\rho_{top}>\ [g\ cm^{-3}]")
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
		time, surfaceMi, 
		color="magenta", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 
	ax.plot(
		time, surfaceMa, 
		color="cyan", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 
	ax.plot(
		time, surfaceT, 
		color="k", marker="s", markerfacecolor="w", ls="-", markersize=5.5
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	

	ax.set_title("Upper Boundary")
	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm  <v_z^{top}>\ [km\ s^{-1}]")
	gcf()
end

# ╔═╡ b4fd0320-74b2-4ed4-b7b0-beca4eccd553


# ╔═╡ 9e883b44-f225-4566-9d76-ecd3e34d3b5b
md"### Mass flux"

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
		time, surfaceMF, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <\rho\ v_z>\ /\ <\rho>\ [km\ s^{-1}]")


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
		time, surfaceMF, 
		color="k", marker="s", markerfacecolor="w", ls="-"
	) 

	ax.axhline(0.0, color="k", alpha=0.2, ls="--")

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <\rho\ v_z>\ /\ <\rho>\ [km\ s^{-1}]")


	gcf()
end

# ╔═╡ Cell order:
# ╟─c7dc3b15-6555-4824-872a-d487fe5145ea
# ╟─c1ad6c59-fae6-49e7-bfc0-8426c553aa2d
# ╠═2f1edd2a-b56e-11ee-29e7-c353938e7088
# ╟─7cd7d6f0-8498-44ff-b59c-d298365d6416
# ╟─78c88a26-1e84-4ba2-a8f2-d4c7f6468dd3
# ╟─e2e0b39b-9c60-4630-817b-f180c2631a08
# ╟─409ff57f-8d9d-419b-b448-fdf40c0843b4
# ╟─6754b2c3-d205-4a12-88b3-53fe62c5637f
# ╟─41f0864e-26ee-46a6-b4ab-c401a4712941
# ╟─c596d1b3-32c5-4651-a8c1-3100fcd6cd59
# ╟─8af0c339-237c-42ca-bda5-0b0e44f11c30
# ╟─498f998c-9996-43dc-a647-3b08abf9786d
# ╟─3c05fe1e-1c19-4c30-b34e-a874f63b91bc
# ╟─ed9cc79f-0161-4178-b6f0-81a2bccbf188
# ╟─6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
# ╟─1fcefd1e-1c50-43b5-b203-62920944344a
# ╟─e752753c-b982-40c6-b45d-ad9bfcd3bbfe
# ╟─bcb4fd58-2b60-4b24-ad3e-ebc68392f320
# ╟─e5b373f5-f565-4910-88b0-f5580880fac4
# ╟─53f4fdc7-5509-44dd-bf71-ceacb78e3e54
# ╟─60199001-8f0e-44a6-ae50-e829687c045c
# ╟─0e4df1b3-52ef-4fce-8f32-bddc966b0516
# ╟─3aadeb1c-3585-46af-8b9d-daf33bfcb7f3
# ╟─2d580186-96a6-41b8-91aa-da527691ea1d
# ╟─e3f00981-eb96-49ba-8451-9af7560d3556
# ╟─73df16bf-b033-494d-b321-608ffff33467
# ╟─bd936d7d-e79f-4f9b-ba54-e0694c6a83f0
# ╟─2c64fcf2-1a0b-49cf-a3f1-890f152d0650
# ╟─ef5ee4f6-96be-4669-ba68-3c1117605a4c
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
# ╠═09b4edcd-d333-42a2-b8e7-1814e85b446a
# ╟─0639ce7d-955d-448f-84a0-353dfd4e93a3
# ╟─c24a3b12-c7c0-449b-9a76-6e9c5d475344
# ╟─3d371088-2322-462b-ab93-7cb49fcdf75f
# ╟─5817821d-67f5-4e41-a0d3-7ca12961b0c7
# ╟─7f77f259-505d-4344-8ee4-8628387f2401
# ╟─cc5fbd5a-c8a0-471a-a56b-0512e4c3989b
# ╟─495e3733-d290-40ab-af63-0eb10a033b53
# ╟─ed29d53f-00bc-4295-93f6-864a44f92ccb
# ╟─fa8161aa-0e1c-404f-8c7e-0e3914917df4
# ╟─3403a014-2441-4ad2-95f6-2e686ae99ba8
# ╟─b0c40c50-3361-4b01-ae87-45ae30387526
# ╟─96a3ddc0-3d70-4b65-a7c2-7482c8817186
# ╟─8c03ccae-6a25-4f2b-a99d-bf537221d007
# ╟─2f86b88f-b53a-489b-8a1e-7d5116492e34
# ╟─f3002fa3-b796-4479-8008-fc8149e797bf
# ╟─68477423-a6f7-424f-8fd4-849d63648b57
# ╟─8ef763ae-d1bd-40de-a26a-f91f529c03bf
# ╟─d9546636-2dbb-4b14-9954-76872b95fd06
# ╟─7fd50fcd-8a9d-417d-a28a-febd9f0f68f3
# ╟─b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5
# ╟─db1a9405-e93c-476d-a43b-f11f3138b57a
# ╟─3da231de-ff58-4161-a6a4-58162483825a
# ╟─efc11ddf-99b9-4344-ab4f-1e5d7b7e2805
# ╟─e0d8aa2e-85cb-43fe-a5f4-65a7a757f19c
# ╟─699ef294-a1b0-4ea5-8618-5cfd0b143c9a
# ╟─0fdf3055-4d11-4dea-8a50-e595ef1c112d
# ╟─eddc02bf-d7ca-41e1-878e-ef1103bf1b0f
# ╟─25c3d608-9440-4ac9-8277-3855ba3b6a7b
# ╟─e59ab649-3b30-4af1-ab28-6f6d6aacee1c
# ╟─d55d7d42-c78c-447c-9959-3689f5341655
# ╟─321e3dda-cd15-4787-95e6-f928125535d5
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
# ╟─b4fd0320-74b2-4ed4-b7b0-beca4eccd553
# ╟─9e883b44-f225-4566-9d76-ecd3e34d3b5b
# ╟─3bfef187-d899-473d-9a92-acf5c65fa50d
# ╟─786e489d-42fb-4ef0-83b6-0e6a7fce2d93
# ╟─01026bad-d93a-4246-a940-e9eb894ccf67
# ╟─dd1297cd-62d0-4b39-b0a3-e41fab79bc8b
# ╟─43590370-0a2e-4d6e-9531-f93eb1630acd
# ╟─cd9ff780-66d2-4bcb-8b89-145d9d857306
