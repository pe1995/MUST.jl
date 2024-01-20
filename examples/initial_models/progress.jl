### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ c596d1b3-32c5-4651-a8c1-3100fcd6cd59
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	mask = isdir.(joinpath.(allruns, "monitoring"))
	split.(allruns[mask], "/", keepempty=false) .|> last
end

# ╔═╡ f15835aa-8e63-441f-a951-3b606850346f


# ╔═╡ 8af0c339-237c-42ca-bda5-0b0e44f11c30
md"""Select one of the available monitored runs:\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""

# ╔═╡ 3c05fe1e-1c19-4c30-b34e-a874f63b91bc


# ╔═╡ ed9cc79f-0161-4178-b6f0-81a2bccbf188
md"## Load monitoring"

# ╔═╡ 6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
monitoring = MUST.reload(MUST.WatchDog, selectedRun, folder=datafolder)

# ╔═╡ 1fcefd1e-1c50-43b5-b203-62920944344a
timeevolution(m, group, field) = [moni[group][field] for moni in m]

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


# ╔═╡ bd936d7d-e79f-4f9b-ba54-e0694c6a83f0
md"## General Properties"

# ╔═╡ 2c64fcf2-1a0b-49cf-a3f1-890f152d0650
time = timeevolution(monitoring, "atmosphericParameters", "time")

# ╔═╡ ef5ee4f6-96be-4669-ba68-3c1117605a4c


# ╔═╡ 8d6d674b-153b-4357-9f2d-c4e3cb05d059
md"## Individual Snapshots"

# ╔═╡ b9a721cf-46ef-4e3c-a37c-8b35653e31cb
md"Pick the time for which you want to see the status: $(@bind timeSurface Slider(time, show_value=true))
"

# ╔═╡ 725f7500-c31b-4efa-9df9-00c51640914a
md"Pick a second time for comparison: $(@bind timeSurface2 Slider(time, show_value=true))
"

# ╔═╡ a34793ae-db12-4dac-b2f8-348a88092815
itimeSurface = findfirst(time .== timeSurface);

# ╔═╡ 63f58f6c-954c-44e9-84b1-1aa799323586
itimeSurface2 = findfirst(time .== timeSurface2);

# ╔═╡ f5ccc550-6bf3-4f0f-baf3-d8de3d216737
md"### Optical Surface"

# ╔═╡ b757013d-aee5-41a4-ab0a-7715ba47bd97
topticalsurfaces = timeevolution(monitoring, "opticalSurfaces")

# ╔═╡ 0639ce7d-955d-448f-84a0-353dfd4e93a3
let 
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax.imshow(
		topticalsurfaces["uzplane"][itimeSurface] ./1e5,
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax, fraction=0.046, pad=0.04)
	cb.set_label(L"\rm v_z\ [km\ s^{-1}]")

	ax.set_xlabel("x [cm]")
	ax.set_ylabel("y [cm]")
	
	gcf()
end

# ╔═╡ c24a3b12-c7c0-449b-9a76-6e9c5d475344
let 
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax.imshow(
		exp.(topticalsurfaces["lnDplane"][itimeSurface]),
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax, fraction=0.046, pad=0.04)
	cb.set_label(L"\rm \rho\ [g\ cm^{-3}]")

	ax.set_xlabel("x [cm]")
	ax.set_ylabel("y [cm]")
	
	gcf()
end

# ╔═╡ 3d371088-2322-462b-ab93-7cb49fcdf75f
let 
	plt.close()
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x = topticalsurfaces["x"][itimeSurface] ./1e8
	y = topticalsurfaces["y"][itimeSurface] ./1e8

	extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
	
	i = ax.imshow(
		topticalsurfaces["Tplane"][itimeSurface],
		origin="lower",
		extent=extent,
		cmap="coolwarm"
	)
	
	cb = f.colorbar(i, ax=ax, fraction=0.046, pad=0.04)
	cb.set_label(L"\rm T\ [K]")

	ax.set_xlabel("x [cm]")
	ax.set_ylabel("y [cm]")
	
	gcf()
end

# ╔═╡ 3403a014-2441-4ad2-95f6-2e686ae99ba8


# ╔═╡ b0c40c50-3361-4b01-ae87-45ae30387526
md"### Geometrical Profiles"

# ╔═╡ 96a3ddc0-3d70-4b65-a7c2-7482c8817186
tGeoAv = timeevolution(monitoring, "geometricalAverages")

# ╔═╡ 68477423-a6f7-424f-8fd4-849d63648b57
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["T"][itimeSurface]
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 

	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
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

	x, y = tGeoAv["z"][itimeSurface] ./1e8, tGeoAv["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 

	x, y = tGeoAv["z"][itimeSurface2] ./1e8, tGeoAv["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
	) 

	ax.set_xlabel("z [Mm]")
	ax.set_ylabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
	ax.legend()

	gcf()
end

# ╔═╡ d9546636-2dbb-4b14-9954-76872b95fd06
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tGeoAv["d"][itimeSurface], tGeoAv["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 

	x, y = tGeoAv["d"][itimeSurface2], tGeoAv["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
	) 

	ax.set_xlabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()

	gcf()
end

# ╔═╡ b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5


# ╔═╡ db1a9405-e93c-476d-a43b-f11f3138b57a
md"### Optical Profiles"

# ╔═╡ 3da231de-ff58-4161-a6a4-58162483825a
tOptAv = timeevolution(monitoring, "opticalAverages")

# ╔═╡ 0fdf3055-4d11-4dea-8a50-e595ef1c112d
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["T"][itimeSurface]
	ax.plot(
		x, y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 

	x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["T"][itimeSurface2]
	ax.plot(
		x, y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
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

	x, y = tOptAv["log10τ_ross"][itimeSurface], tOptAv["d"][itimeSurface]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 

	x, y = tOptAv["log10τ_ross"][itimeSurface2], tOptAv["d"][itimeSurface2]
	ax.plot(
		x, log10.(y),
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
	) 

	ax.set_xlabel(L"\rm \log\ \tau_{ross}")
	ax.set_ylabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
	ax.legend()

	gcf()
end

# ╔═╡ 25c3d608-9440-4ac9-8277-3855ba3b6a7b
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	x, y = tOptAv["d"][itimeSurface], tOptAv["T"][itimeSurface]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="-", label="t = $(time[itimeSurface]) s"
	) 
	
	x, y = tOptAv["d"][itimeSurface2], tOptAv["T"][itimeSurface2]
	ax.plot(
		log10.(x), y,
		color="k", marker="", ls="--", label="t = $(time[itimeSurface2]) s"
	) 

	ax.set_xlabel(L"\rm \log \rho\ [g\ cm^{-3}]]")
	ax.set_ylabel(L"\rm T\ [K]")
	
	ax.legend()

	gcf()
end

# ╔═╡ d55d7d42-c78c-447c-9959-3689f5341655


# ╔═╡ 321e3dda-cd15-4787-95e6-f928125535d5
md"## Time evolution"

# ╔═╡ 51858291-e78b-4d05-8441-64f33694d133
md"### Mass Flux"

# ╔═╡ 145f084a-a24b-4bcd-b5a2-42b114fa8df6
tmassfluxgeo = timeevolution(monitoring, "geoMassFlux")

# ╔═╡ fb9216dd-c811-46fe-96f6-316369a00a1e
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceMF = (tmassfluxgeo["massFlux"] .|> last) ./1e5
	
	ax.plot(
		time, surfaceMF, 
		color="k", marker="s", markerfacecolor="None", ls="-"
	) 

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <\rho\ v_z>\ /\ <\rho>\ [km\ s^{-1}]")


	gcf()
end

# ╔═╡ dc0cd97f-f4e9-4a91-8c4a-ed27571a7f48
md"### Temperature Optical Surface"

# ╔═╡ c2b64d82-09a7-4d67-97d0-b51d96d30d25
ttempsurface = timeevolution(monitoring, "opticalSurfaces", "Tplane")

# ╔═╡ eb335a4d-e662-499c-bb80-8bf38c84329f
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceTemp = [mean(t) for t in ttempsurface]
	
	ax.plot(
		time, surfaceTemp, 
		color="k", marker="s", markerfacecolor="None", ls="-"
	) 

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm <T_{\tau=1}>\ [K]")

	gcf()
end

# ╔═╡ f328cd09-36f4-4910-a2e5-bc839df18e1e
md"### Density Optical Surface"

# ╔═╡ 71c4c81a-5567-45f0-925a-1dcad0f97082
tdsurface = timeevolution(monitoring, "opticalSurfaces", "lnDplane")

# ╔═╡ 35b4aa7e-e41b-4444-af3b-74684c723d5c
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceRho = [log10.(mean(exp.(t))) for t in tdsurface]
	
	ax.plot(
		time, surfaceRho, 
		color="k", marker="s", markerfacecolor="None", ls="-"
	) 	

	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm \log < \rho_{\tau=1}>\ [g\ cm^{-3}]")

	gcf()
end

# ╔═╡ 3b3e1234-8efc-4bd1-9807-1b63548c08c2
md"### RMS Velocity Optical Surface"

# ╔═╡ 1521d2b2-56a7-41eb-ad52-b661af7f30c6
tuzsurface = timeevolution(monitoring, "opticalSurfaces", "uzplane")

# ╔═╡ 87d43815-9733-40ac-b4e7-81a2c5dbd0d1
let
	plt.close()

	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	surfaceUz = [sqrt.(mean(t .^2)) for t in tuzsurface]

	ax.plot(
		time, surfaceUz ./1e5, 
		color="k", marker="s", markerfacecolor="None", ls="-"
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
	
	ax.plot(
		time, surfaceT, 
		color="k", marker="s", markerfacecolor="None", ls="-"
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
	
	ax.plot(
		time, log10.(surfaceT), 
		color="k", marker="s", markerfacecolor="None", ls="-"
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
	
	ax.plot(
		time, surfaceT, 
		color="k", marker="s", markerfacecolor="None", ls="-"
	) 

	ax.set_title("Upper Boundary")
	ax.set_xlabel("time [s]")
	ax.set_ylabel(L"\rm  <v_z^{top}>\ [g\ cm^{-3}]")
	gcf()
end

# ╔═╡ d50b1713-c0a9-4da7-8397-95b1cd3059a7


# ╔═╡ Cell order:
# ╟─c7dc3b15-6555-4824-872a-d487fe5145ea
# ╠═2f1edd2a-b56e-11ee-29e7-c353938e7088
# ╟─7cd7d6f0-8498-44ff-b59c-d298365d6416
# ╟─78c88a26-1e84-4ba2-a8f2-d4c7f6468dd3
# ╠═e2e0b39b-9c60-4630-817b-f180c2631a08
# ╟─409ff57f-8d9d-419b-b448-fdf40c0843b4
# ╟─6754b2c3-d205-4a12-88b3-53fe62c5637f
# ╟─c596d1b3-32c5-4651-a8c1-3100fcd6cd59
# ╟─f15835aa-8e63-441f-a951-3b606850346f
# ╟─8af0c339-237c-42ca-bda5-0b0e44f11c30
# ╟─3c05fe1e-1c19-4c30-b34e-a874f63b91bc
# ╟─ed9cc79f-0161-4178-b6f0-81a2bccbf188
# ╟─6b68e1e5-fe97-4485-9c7b-c3e821f23a7c
# ╟─1fcefd1e-1c50-43b5-b203-62920944344a
# ╟─e752753c-b982-40c6-b45d-ad9bfcd3bbfe
# ╟─bcb4fd58-2b60-4b24-ad3e-ebc68392f320
# ╟─e5b373f5-f565-4910-88b0-f5580880fac4
# ╟─53f4fdc7-5509-44dd-bf71-ceacb78e3e54
# ╟─60199001-8f0e-44a6-ae50-e829687c045c
# ╟─bd936d7d-e79f-4f9b-ba54-e0694c6a83f0
# ╟─2c64fcf2-1a0b-49cf-a3f1-890f152d0650
# ╟─ef5ee4f6-96be-4669-ba68-3c1117605a4c
# ╟─8d6d674b-153b-4357-9f2d-c4e3cb05d059
# ╟─b9a721cf-46ef-4e3c-a37c-8b35653e31cb
# ╟─725f7500-c31b-4efa-9df9-00c51640914a
# ╟─a34793ae-db12-4dac-b2f8-348a88092815
# ╟─63f58f6c-954c-44e9-84b1-1aa799323586
# ╟─f5ccc550-6bf3-4f0f-baf3-d8de3d216737
# ╟─b757013d-aee5-41a4-ab0a-7715ba47bd97
# ╟─0639ce7d-955d-448f-84a0-353dfd4e93a3
# ╟─c24a3b12-c7c0-449b-9a76-6e9c5d475344
# ╟─3d371088-2322-462b-ab93-7cb49fcdf75f
# ╟─3403a014-2441-4ad2-95f6-2e686ae99ba8
# ╟─b0c40c50-3361-4b01-ae87-45ae30387526
# ╟─96a3ddc0-3d70-4b65-a7c2-7482c8817186
# ╟─68477423-a6f7-424f-8fd4-849d63648b57
# ╟─8ef763ae-d1bd-40de-a26a-f91f529c03bf
# ╟─d9546636-2dbb-4b14-9954-76872b95fd06
# ╟─b09b0c9e-3d76-4ff0-ae82-205cbc3b83b5
# ╟─db1a9405-e93c-476d-a43b-f11f3138b57a
# ╟─3da231de-ff58-4161-a6a4-58162483825a
# ╟─0fdf3055-4d11-4dea-8a50-e595ef1c112d
# ╟─eddc02bf-d7ca-41e1-878e-ef1103bf1b0f
# ╟─25c3d608-9440-4ac9-8277-3855ba3b6a7b
# ╟─d55d7d42-c78c-447c-9959-3689f5341655
# ╟─321e3dda-cd15-4787-95e6-f928125535d5
# ╟─51858291-e78b-4d05-8441-64f33694d133
# ╟─145f084a-a24b-4bcd-b5a2-42b114fa8df6
# ╟─fb9216dd-c811-46fe-96f6-316369a00a1e
# ╟─dc0cd97f-f4e9-4a91-8c4a-ed27571a7f48
# ╟─c2b64d82-09a7-4d67-97d0-b51d96d30d25
# ╟─eb335a4d-e662-499c-bb80-8bf38c84329f
# ╟─f328cd09-36f4-4910-a2e5-bc839df18e1e
# ╟─71c4c81a-5567-45f0-925a-1dcad0f97082
# ╟─35b4aa7e-e41b-4444-af3b-74684c723d5c
# ╟─3b3e1234-8efc-4bd1-9807-1b63548c08c2
# ╟─1521d2b2-56a7-41eb-ad52-b661af7f30c6
# ╟─87d43815-9733-40ac-b4e7-81a2c5dbd0d1
# ╟─c8492ca5-893f-4939-a256-f0872f351c5d
# ╟─df9b58cf-fe4f-4858-a296-879bb0385ba7
# ╠═ccf5f4e6-adc6-419d-a717-4b3b597c2233
# ╟─913c1cbe-fedc-4e7d-a8a7-c2ad416a21e6
# ╟─2f76eaac-0793-4fa8-9ee1-67dc81e9a1ac
# ╟─aa0376dc-5982-4274-9b4d-4aef1d1de896
# ╟─d50b1713-c0a9-4da7-8397-95b1cd3059a7
