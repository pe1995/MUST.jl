### A Pluto.jl notebook ###
# v0.19.43

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
	using Pkg; Pkg.activate("."); Pkg.instantiate()
	using MUST
	using PlutoPlotly
	using Images
	using PlutoUI
	using PlutoUI: combine
	using ProgressLogging
	using LaTeXStrings
	using Printf
	using ImageFiltering
	using DelimitedFiles
	using JSON

	#plotly()
	#PlutoPlotly.default(grid=false, framestyle=:box, minorticks=true)
	PlutoPlotly.default_plotly_template("plotly_dark")

	#datafolder = "/data/web2/public/atmos3D"
	datafolder = "../data/"
	if !(isdir(datafolder))
		mkdir(datafolder)
	end
	
	mean = MUST.mean
	#plt = matplotlib.pyplot
	#matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	#matplotlib.style.use("dark_background")
	
	scipy_fft = MUST.pyimport("scipy.fft")

	# create a tempdir for all the movie content
	# this is needed if multiple people try to run the same code at the same time
	v_opt_folder_name = mktempdir("../")
	#v_opt_folder_name = mktempdir(parent="/data/web2/public/results/")

	selectedPointsFile = joinpath(v_opt_folder_name, "selectedPoints.txt")
	
	TableOfContents(title="DISPATCH Spectroscopy 🌟")
end

# ╔═╡ c7dc3b15-6555-4824-872a-d487fe5145ea
md"""
# DISPATCH 3D Stellar Atmospheres
This is a monitoring board to visualize DISPATCH stellar atmosphere simulations while they are running or after completion. Within this interactive notebook you can visualize key physical aspects for selected 3D simulation of stars created within the DISPATCH framework. [More details can be found in Eitner et al. 2024.](https://ui.adsabs.harvard.edu/abs/2024arXiv240506338E/abstract)
"""

# ╔═╡ 6754b2c3-d205-4a12-88b3-53fe62c5637f
md"## Select Simulation
Simulations of multiple stars available within this interface. You can select a simulation from the dropdown menu below."

# ╔═╡ 41f0864e-26ee-46a6-b4ab-c401a4712941
md"Press this button to check for new simulations: $(@bind reload_data Button(\"Refresh\"))"

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

# ╔═╡ 3c05fe1e-1c19-4c30-b34e-a874f63b91bc
begin
	wd = MUST.WatchDog(selectedRun, folder=datafolder)
	monitoring = MUST.reload!(wd; mmap=true)

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

	snapshotnames(moni) = [nothing, wd.snapshotsCompleted...]
	time = timeevolution(monitoring, "atmosphericParameters", "time")
	snapshots = snapshotnames(monitoring)[2:end]
end;

# ╔═╡ 66e8e497-c3d4-48e7-9c2a-9de4d52b3dd4
begin
	dwnldFolder = joinpath("atmos3D", selectedRun, "snapshots")

	md"""
	You can download selected snapshots of your chosen simulation direclty by clicking the download button below. Note that the download button only appears if snapshots are available to download for the simulation that you chose. Note also that this may close the current window, so consider clicking it with ⌘+click.
	"""
end

# ╔═╡ 431bee5c-bfc1-4f53-8170-3ebee1277a2d
if isfile(dwnldFolder)
	dwnld_str = """<a href="$(dwnldFolder)", download="$(selectedRun).zip" >Download the snapshots by clicking here.</a>"""
	Docs.HTML(dwnld_str)
else
	@info "No snapshots available to download for the given simulation."
end

# ╔═╡ ae5ab0a3-89db-4f1b-8e36-027e84d65c45
md"Note that you can also download a static version of this notebook (including the plots you made) by clicking on the arrow to the top right of this page and then select 'static HTML'."

# ╔═╡ 8252972b-e496-4f13-b5d6-39e441fb9b84


# ╔═╡ b76f0efd-6b2f-47d3-84e7-ca9fb866292e
md"## Resolution
Your chosen model has the following resolution in the Hydrodynamics part. Please note that the resolution in the radiative transfer is independent from this resolution, and is usually twice as high."

# ╔═╡ 4eef1012-c8cd-4579-9a7e-ac2b4d4c7af6
let
	z = timeevolution(monitoring, "geometricalAverages", "z") |> first
	HD_dims = (length(z)*2, length(z)*2, length(z))
	HD_res = maximum(diff(z)) / 1e5
	
	@info "[$selectedRun] HD Resolution (Nx, Ny, Nz), dz [km]" HD_dims HD_res
end

# ╔═╡ 33c1da97-0760-495e-abd5-65531d5e1170


# ╔═╡ 8d6d674b-153b-4357-9f2d-c4e3cb05d059
md"# Investigating Snapshots
In the following you can investigate diiferent snapshots of your chosen simulation."

# ╔═╡ e938c3cf-e473-44ee-978c-4c0135856cab
md"Only show snapshots that have spectra available: $(@bind showOnlySpectra CheckBox(default=false))"

# ╔═╡ e61334ad-c74c-4334-87a3-46c15f622ce6
snapshots_pick = if showOnlySpectra
	isspectrum(m) = [occursin("resolvedSpectra", k) && haskey(m[k], "intensity") for k in keys(m)]
	mask = [any(isspectrum(m)) for m in monitoring]
	snapshots[mask]
else
	snapshots
end;

# ╔═╡ b9a721cf-46ef-4e3c-a37c-8b35653e31cb
md"Pick a snapshot: $(@bind timeSurface Slider(snapshots_pick, show_value=true, default=last(snapshots_pick)))
"

# ╔═╡ a34793ae-db12-4dac-b2f8-348a88092815
begin
	timeSurfacel = if length(timeSurface) > 0
		last(timeSurface)
	else
		timeSurface
	end
	itimeSurface = findfirst(snapshots .== timeSurfacel)
	topticalsurfaces = timeevolution(monitoring, "opticalSurfaces")
	tuppersurfaces = timeevolution(monitoring, "upperBoundarySurface")
	tlowersurfaces = if haskey(monitoring[1], "lowerBoundarySurface") 
		timeevolution(monitoring, "lowerBoundarySurface")
	else
		nothing
	end

	
	default_axis_template = attr(
		 showgrid = false, zeroline = false, showline=true, linecolor="white",
		 showticklabels = true
	)
	
	defaultLayout(;axis_template=default_axis_template, layout_kwargs...) = Layout(;
		xaxis=axis_template,
		yaxis=axis_template,
		layout_kwargs...
	)
	
	defaultPlot(f, args...; layout=defaultLayout(), kwargs...) = PlutoPlot(
		Plot(
			f(args...; kwargs...),
			layout
		)
	)
	defaultSubPlot(f, args...; layout=defaultLayout(), kwargs...) = f(args...; kwargs...)
end;

# ╔═╡ 27885149-1f0f-4c75-ae98-96abd9509ddc
begin
	struct SurfaceSlice
		zAxis
		label
		xAxis
		yAxis
	end
	surfacesSliceSelection = Dict()

	#==== Upper boundary ====#
	surfacesSliceSelection["upper boundary - temperature"] = SurfaceSlice(
		tuppersurfaces["Tplane"][itimeSurface], 
		"temperature [K]", 
		tuppersurfaces["x"][itimeSurface] ./1e8, 
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	surfacesSliceSelection["upper boundary - density"] = SurfaceSlice(
		exp.(tuppersurfaces["lnDplane"][itimeSurface]), 
		"density [g /cm3]",
		tuppersurfaces["x"][itimeSurface] ./1e8, 
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	if haskey(tuppersurfaces, "fluxplane")
		surfacesSliceSelection["upper boundary - flux"] = SurfaceSlice(
			tuppersurfaces["fluxplane"][itimeSurface], 
			"F [erg /s /cm2]",
			tuppersurfaces["x"][itimeSurface] ./1e8, 
			tuppersurfaces["y"][itimeSurface] ./1e8
		)
	end
	surfacesSliceSelection["upper boundary - vertical velocity"] = SurfaceSlice(
		tuppersurfaces["uzplane"][itimeSurface] ./1e5,
		"vz [km /s]", 
		tuppersurfaces["x"][itimeSurface] ./1e8,
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	#=if haskey(tuppersurfaces, "dtplane")
		surfacesSliceSelection["upper boundary - timestep"] = SurfaceSlice(
			tuppersurfaces["dtplane"][itimeSurface], 
			L"\rm timestep\ [s]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end=#
	if haskey(tuppersurfaces, "qrplane")
		surfacesSliceSelection["upper boundary - heating"] = SurfaceSlice(
			tuppersurfaces["qrplane"][itimeSurface], 
			"Qr [erg /s /cm3]",
			tuppersurfaces["x"][itimeSurface] ./1e8, 
			tuppersurfaces["y"][itimeSurface] ./1e8
		)
	end


	#==== Lower boundary ====#
	if !isnothing(tlowersurfaces)
		surfacesSliceSelection["lower boundary - temperature"] = SurfaceSlice(
			tlowersurfaces["Tplane"][itimeSurface], 
			"temperature [K]", 
			tlowersurfaces["x"][itimeSurface] ./1e8, 
			tlowersurfaces["y"][itimeSurface] ./1e8
		)
		surfacesSliceSelection["lower boundary - density"] = SurfaceSlice(
			exp.(tlowersurfaces["lnDplane"][itimeSurface]), 
			"density [g /cm3]",
			tlowersurfaces["x"][itimeSurface] ./1e8, 
			tlowersurfaces["y"][itimeSurface] ./1e8
		)
		if haskey(tlowersurfaces, "fluxplane")
			surfacesMovieSelection["lower boundary - flux"] = SurfaceSlice(
				tlowersurfaces["fluxplane"][itimeSurface], 
				"F [erg /s /cm2]",
				tlowersurfaces["x"][itimeSurface] ./1e8, 
				tlowersurfaces["y"][itimeSurface] ./1e8
			)
		end
		surfacesSliceSelection["lower boundary - vertical velocity"] = SurfaceSlice(
			tlowersurfaces["uzplane"][itimeSurface] ./1e5,
			"vz [km /s]", 
			tlowersurfaces["x"][itimeSurface] ./1e8,
			tlowersurfaces["y"][itimeSurface] ./1e8
		)
		#=if haskey(tlowersurfaces, "dtplane")
			surfacesSliceSelection["lower boundary - timestep"] = SurfaceSlice(
				tlowersurfaces["dtplane"][itimeSurface], 
				L"\rm timestep\ [s]",
				tlowersurfaces["x"][1] ./1e8, 
				tlowersurfaces["y"][1] ./1e8
			)
		end=#
		if haskey(tlowersurfaces, "qrplane")
			surfacesSliceSelection["lower boundary - heating"] = SurfaceSlice(
				tlowersurfaces["qrplane"][itimeSurface], 
				"Qr [erg /s /cm3]",
				tlowersurfaces["x"][itimeSurface] ./1e8, 
				tlowersurfaces["y"][itimeSurface] ./1e8
			)
		end
	end
		
	
	
	#==== Optical surface ====#
	surfacesSliceSelection["optical surface - temperature"] = SurfaceSlice(
		topticalsurfaces["Tplane"][itimeSurface], 
		"temperature [K]", 
		tuppersurfaces["x"][itimeSurface] ./1e8, 
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	surfacesSliceSelection["optical surface - density"] = SurfaceSlice(
		exp.(topticalsurfaces["lnDplane"][itimeSurface]), 
		"density [g /cm3]",
		tuppersurfaces["x"][itimeSurface] ./1e8, 
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	if haskey(topticalsurfaces, "fluxplane")
		surfacesSliceSelection["optical surface - flux"] = SurfaceSlice(
			topticalsurfaces["fluxplane"][itimeSurface], 
			"F [erg /s /cm2]",
			tuppersurfaces["x"][itimeSurface] ./1e8, 
			tuppersurfaces["y"][itimeSurface] ./1e8
		)
	end
	surfacesSliceSelection["optical surface - vertical velocity"] = SurfaceSlice(
		topticalsurfaces["uzplane"][itimeSurface] ./1e5, 
		"vz [km /s]", 
		tuppersurfaces["x"][itimeSurface] ./1e8, 
		tuppersurfaces["y"][itimeSurface] ./1e8
	)
	#=if haskey(topticalsurfaces, "dtplane")
		surfacesSliceSelection["optical surface - timestep"] = SurfaceSlice(
			topticalsurfaces["dtplane"][itimeSurface], 
			L"\rm timestep\ [s]",
			tuppersurfaces["x"][1] ./1e8, 
			tuppersurfaces["y"][1] ./1e8
		)
	end=#
	if haskey(topticalsurfaces, "qrplane")
		surfacesSliceSelection["optical surface - heating"] = SurfaceSlice(
			topticalsurfaces["qrplane"][itimeSurface], 
			"Qr [erg /s /cm3]",
			tuppersurfaces["x"][itimeSurface] ./1e8, 
			tuppersurfaces["y"][itimeSurface] ./1e8
		)
	end
end;

# ╔═╡ 333cbf0a-5a20-4612-8aa0-299e1ff40359
md"## Horizontal Slices"

# ╔═╡ 1ef75139-4ee0-44fc-a731-a341ccd48e74
md"Below horizintal slices of the chosen snapshot. Slices on the optical depth scale are created by interpolating snapshots from the geometrical to the rosseland optical depth scale."

# ╔═╡ aafbb147-e6ed-448f-bf4c-92f6c88d5199
md"""
	Select quantity to show
$(@bind surface_choice Select(sort(keys(surfacesSliceSelection)|>collect), default="optical surface - temperature"))
"""

# ╔═╡ 8a764271-7957-4a47-a979-0340a2ae7408
begin
	x = surfacesSliceSelection[surface_choice].xAxis[:, 1]
	y = surfacesSliceSelection[surface_choice].yAxis[1, :]
end;

# ╔═╡ ab72c112-5f17-4e71-a49e-5b2ffe2c6c20
md"""
	Adjust figure settings

__colormap settings__\


colormap: $(@bind sliceCmap confirm(TextField(32, default=\"linear_kryw_5_100_c64_n256\")))\

reverse colormap: $(@bind sliceCmapReverse CheckBox(default=false))\

log scale for color axis: $(@bind sliceCmapLog CheckBox(default=true))\

make colormap symmetric $(@bind sliceCmapSym CheckBox(default=false))\
"""

# ╔═╡ ada8d5a9-9a30-43a6-af66-050bcd000952


# ╔═╡ 2182de0b-bb9a-4fb2-b9df-45ef4ecd084e
md"""
Below you can find a list of colormaps that are available. Just expand the dictionary and have a look! A good choice for the velocity field is e.g.
`diverging_bwr_40_95_c42_n256` or `bwr`, and for the temperature `linear_kryw_5_100_c64_n256`.
"""

# ╔═╡ de173b50-fc57-4ac1-bd82-079921a29da5
begin
	cmaps = PlutoPlotly.PlotlyBase.colors.sequential
end

# ╔═╡ 0639ce7d-955d-448f-84a0-353dfd4e93a3
begin 
	@bind pickOptVSurf let 
		x = surfacesSliceSelection[surface_choice].xAxis
		y = surfacesSliceSelection[surface_choice].yAxis
		z = surfacesSliceSelection[surface_choice].zAxis'
		z = sliceCmapLog ? log10.(z) : z
		vmin, vmax = if sliceCmapSym
			vmax = maximum(abs.(z))
			-vmax, vmax
		else
			minimum(z), maximum(z)
		end
	
		extent = [minimum(x), maximum(x), minimum(y), maximum(y)]
		cmap = cmaps[Symbol(sliceCmap)]
		p1 = defaultPlot(
			heatmap, 
			z=z,
			x=x[:, 1], y=y[1, :],
			zmin=vmin, zmax=vmax,
			colorscale=cmap,
			colorbar_titleside="right",
			colorbar_title=surfacesSliceSelection[surface_choice].label,
			reversescale=sliceCmapReverse,
			layout = defaultLayout(
				xaxis_title="X [Mm]",
				yaxis_title="Y [Mm]",
				yaxis_scaleanchor="x",
				title=surface_choice
			)
		)

		add_plotly_listener!(
			p1, "plotly_click","
        (e) => {
            let dt = e.points[0];
            PLOT.value = [dt.x, dt.y];
            PLOT.dispatchEvent(new CustomEvent('input'));
        }
        "
    )
	
		p1
	end
end

# ╔═╡ f1e4cbf4-0adb-4156-b696-5ea3801c6ae2
function set_selected_names(selectedPoints, x, y) 
	return combine() do Child
		names_loc = ["(x=$(x[xlab]), y=$(y[ylab]))" for (xlab, ylab) in selectedPoints]
		namesChild = [
			Child(names_loc[ij], confirm(TextField(32, default=names_loc[ij])))
			for ij in eachindex(names_loc)
		]
		inputs = [
			md""" $(name): $(
				namesChild[ik]
			)"""
			for (ik, name) in enumerate(names_loc)
		]
		
		md"""
		$(inputs)
		"""
	end
end;

# ╔═╡ 585a336e-9c90-4b1b-93ac-8d81ff35d49d
md"""
You can __click__ in the figure to take a closer look at spatially-dependent quantities!
"""

# ╔═╡ c82a3030-2bba-49ee-9eb6-cf07e901438d
md"""
	Reset selection of points
__Click__ to remove selection, __click again__ to continue plotting: $(@bind removeSpectra2 CheckBox(default=false))
"""

# ╔═╡ 0584c2c0-de6e-43a9-a7cc-3dfc0742de21
begin
	if !removeSpectra2
		# check temp file
		newFile = if !isfile(selectedPointsFile)
			open(selectedPointsFile, "w") do f
			end
			true
		else
			open(selectedPointsFile, "a+") do f
				if !isnothing(pickOptVSurf) & !ismissing(pickOptVSurf)
					write(f, join(pickOptVSurf, ' ')*"\n")
				end
			end
			false
		end
	end
end;

# ╔═╡ c1ca0584-92cd-497b-a8bd-8ca50ac84631
begin
	if removeSpectra2
		isfile(selectedPointsFile) && rm(selectedPointsFile)
	end
end;

# ╔═╡ 3b070fdb-12c0-4ee8-b2dd-6d69fc972b81
begin
	pickOptVSurf
	removeSpectra2

	ix_pickOptVSurf_arr = []
	iy_pickOptVSurf_arr = []

	redo=true

	# read the file and see if there are points 
	isfile(selectedPointsFile) && !newFile && open(selectedPointsFile, "r") do f
		lines = readlines(f)
		
		for line in lines
			l = split(line, " ", keepempty=false)

			ix_pickOptVSurf = argmin(
				abs.(surfacesSliceSelection[surface_choice].xAxis[:,1] .- parse(Float64, l[1]))
			)
			iy_pickOptVSurf = argmin(
				abs.(surfacesSliceSelection[surface_choice].yAxis[1,:] .- parse(Float64, l[2]))
			)

			if !(ix_pickOptVSurf in ix_pickOptVSurf_arr) && !(iy_pickOptVSurf in iy_pickOptVSurf_arr)
				append!(ix_pickOptVSurf_arr, [ix_pickOptVSurf])
				append!(iy_pickOptVSurf_arr, [iy_pickOptVSurf])
			end
		end
	end

	ix_pickOptVSurf_arr, iy_pickOptVSurf_arr
end;

# ╔═╡ 4eb0810c-b724-4284-8e97-5eb0b3d357ed
begin
	#refresh
	#for i in eachindex(ix_pickOptVSurf_arr)
	#	pop!(ix_pickOptVSurf_arr)
	#end
	#for i in eachindex(iy_pickOptVSurf_arr)
	#	pop!(iy_pickOptVSurf_arr)
	#end

	if length(ix_pickOptVSurf_arr) > 0
		md"## Spatially Resolved Quantities"
	end
end

# ╔═╡ 0eb83755-4f86-4771-bcab-c7bbba096fb2
begin
	pickOptVSurf
	redo
	removeSpectra2
	
	struct ResolvedStatistik
		yAxis
		label
		xAxis
		xlabel
	end
	resolvedSelection = Dict()
	
	interpolate_spectrum(x_spec, y_spec, int, xlab, ylab) = begin
		ingrid = MUST.Grid(x_spec[:, 1], y_spec[1, :])
		outgrid = MUST.Grid([xlab], [ylab])
		
		ip = MUST.ginterpolate(ingrid, outgrid)
		intensity = [
			MUST.gevaluate!(ip, view(int, i, :, :))[1, 1] for i in axes(int, 1)
		]
	end

	time_point(arr, f=identity) = () -> [
		[f(arr[it][ix, iy]) for it in eachindex(time)] for (ix, iy) in zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr)
	]


	
	#=== Spectra ===#
	add_spectrum_key!(resolvedSelection, mon, name, output_name) = begin
		if haskey(mon, name)
			spec = mon
			x_spec = spec[name]["x"] ./1e8
			y_spec = spec[name]["y"] ./1e8
			int = spec[name]["intensity"]
			
			#=resolvedSelection["spectrum 6200-6790Å (4MOST)"] = ResolvedStatistik(
				[
					interpolate_spectrum(x_spec, y_spec, int, x[ix], y[iy]) for (ix, iy) in zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr)
				], 
				"intensity",
				spec["resolvedSpectra4MOST3"]["wavelength"],
				"λ [Å]"
			)=#
	
			ic = if haskey(spec[name], "continuum")
				cont = spec[name]["continuum"]
				ic = int ./cont
			else
				int
			end
			resolvedSelection[output_name] = ResolvedStatistik(
				() -> [
					interpolate_spectrum(x_spec, y_spec, ic, x[ix], y[iy]) for (ix, iy) in zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr)
				], 
				"intensity",
				spec[name]["wavelength"],
				"λ [Å]"
			)
		end
	end

	add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectra4MOST3", 
		"spectrum 6200-6790Å (4MOST)"
	)
	#=add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectraAPOGEESDSSV", 
		"spectrum 15140-17000Å (APOGEE & SDSSV)"
	)=#
	add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectraGaiaESOHR10", 
		"spectrum 5339-5619Å (GaiaESO HR10)"
	)
	add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectraGaiaRVS", 
		"spectrum 8470-8740Å (GaiaRVS)"
	)
	add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectraGaiaRVSCATriplet", 
		"spectrum 8480-8560Å (GaiaRVS Ca-triplet)"
	)
	add_spectrum_key!(
		resolvedSelection, 
		monitoring[itimeSurface], 
		"resolvedSpectraHalpha", 
		"spectrum Hα"
	)
	
	
	#==== Optical surface ====#
	resolvedSelection["optical surface - temperature vs. time"] = ResolvedStatistik(
		time_point(topticalsurfaces["Tplane"]),
		"temperature [K]",  
		time ./60,
		"time [min]"
	)
	resolvedSelection["optical surface - density vs. time"] = ResolvedStatistik(
		time_point(topticalsurfaces["lnDplane"]),
		"log10 density [g /cm3]",
		time ./60,
		"time [min]"
	)
	if haskey(topticalsurfaces, "fluxplane")
		resolvedSelection["optical surface - flux vs. time"] = ResolvedStatistik(
			time_point(topticalsurfaces["fluxplane"]),
			"F [erg /s /cm2]",
			time ./60,
		"time [min]"
		)
	end
	resolvedSelection["optical surface - vertical velocity vs. time"] = ResolvedStatistik(
		time_point(topticalsurfaces["uzplane"], x->x/1e5),
		"vz [km /s]", 
		time ./60,
		"time [min]"
	)
	#=if haskey(topticalsurfaces, "dtplane")
		resolvedSelection["optical surface - timestep"] = ResolvedStatistik(
			[
				topticalsurfaces["dtplane"][itimeSurface][ix, iy] for (ix, iy) in zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr)
			],
			L"\rm timestep\ [s]",
			time
		)
	end=#
	if haskey(topticalsurfaces, "qrplane")
		resolvedSelection["optical surface - heating vs. time"] = ResolvedStatistik(
			time_point(topticalsurfaces["qrplane"]),
			"Qr [erg /s /cm3]",
			time ./60,
		"time [min]"
		)
	end



	#==== upper surface ====#
	resolvedSelection["upper surface - temperature vs. time"] = ResolvedStatistik(
		time_point(tuppersurfaces["Tplane"]),
		"temperature [K]",  
		time ./60,
		"time [min]"
	)
	resolvedSelection["upper surface - density vs. time"] = ResolvedStatistik(
		time_point(tuppersurfaces["lnDplane"]),
		"log10 density [g /cm3]",
		time ./60,
		"time [min]"
	)
	if haskey(tuppersurfaces, "fluxplane")
		resolvedSelection["upper surface - flux vs. time"] = ResolvedStatistik(
			time_point(tuppersurfaces["fluxplane"]),
			"F [erg /s /cm2]",
			time ./60,
		"time [min]"
		)
	end
	resolvedSelection["upper surface - vertical velocity vs. time"] = ResolvedStatistik(
		time_point(tuppersurfaces["uzplane"], x->x/1e5),
		"vz [km /s]", 
		time ./60,
		"time [min]"
	)
	#=if haskey(tuppersurfaces, "dtplane")
		resolvedSelection["optical surface - timestep"] = ResolvedStatistik(
			[
				tuppersurfaces["dtplane"][itimeSurface][ix, iy] for (ix, iy) in zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr)
			],
			L"\rm timestep\ [s]",
			time
		)
	end=#
	if haskey(tuppersurfaces, "qrplane")
		resolvedSelection["upper surface - heating vs. time"] = ResolvedStatistik(
			time_point(tuppersurfaces["qrplane"]),
			"Qr [erg /s /cm3]",
			time ./60,
		"time [min]"
		)
	end

	

	

	if length(ix_pickOptVSurf_arr) > 0
		# check if there is any spectrum available
		if haskey(monitoring[itimeSurface], "resolvedSpectraHalpha")
			md"""
				Pick spatially resolved quantity to plot
			$(@bind resolved_choice Select(sort(keys(resolvedSelection)|>collect), default="spectrum Hα"))
			"""
		else
			md"""
				Pick spatially resolved quantity to plot
			$(@bind resolved_choice Select(sort(keys(resolvedSelection)|>collect), default="optical surface - temperature vs. time"))
			"""
		end
	end
end

# ╔═╡ 1691b6fd-8d7d-46b6-a692-5b0f5a986e91
begin
	"""
	Convolve a spectrum from its resolution to the given resolution.
	"""
	function convolve(wavelength::Vector, spectrum::Vector; R=nothing, Δλ=nothing)
	    @assert !isnothing(Δλ) | !isnothing(R)
	
	    Δλ = if isnothing(Δλ)
	        mean(wavelength) / R
	    else
	        Δλ
	    end
	    convolve(wavelength, spectrum, Δλ)
	end
	
	"""
	Convolve a spectrum from its resolution to the given resolution.
	"""
	function convolve(wavelength::Vector, spectrum::Vector, Δλ_goal)
	    Δλ          = wavelength[2] - wavelength[1]
	    kernel      = ImageFiltering.reflect(ImageFiltering.Kernel.gaussian( (Δλ_goal/Δλ *0.4246609,) ))
	    spec_conv   = ImageFiltering.imfilter(spectrum, kernel)
	    err_conv    = ImageFiltering.imfilter(error, kernel)
	    
	    (spec_conv, err_conv)
	end
end;

# ╔═╡ 63ec2623-eb88-43f1-a037-b86c013181f2
begin
	struct Atom
		name
		lines
	end
	struct Line
		λ
		info
	end
	struct LineList
		atoms
	end

	linelist = nothing
	pickOptVSurf
	redo
	removeSpectra2
	
	
	function readLinelist(path)
		lines = []
		open(path, "r") do f
			append!(lines, readlines(f))
		end

		i = 1
		read_data = false
		read_header = false
		atoms = Dict()
		order = []
		currentName = Ref("name")
		for line in lines
			ls = split(line, " ", keepempty=false)
			if (occursin("'", ls[1]) && !read_header)
				read_data=false
				read_header=true
				continue
			end
			if (occursin("'", ls[1]) && read_header)
				read_data=true
				read_header=false
				name = split(
				split(line, "'", keepempty=false)|>first, " ", keepempty=false
				)
				name = join([name[1], name[2]], " ")
				if !haskey(atoms, name)
					atoms[name] = Atom(name, [])
					append!(order, [name])
				end
				currentName[] = name
				continue
			end
			if read_data
				λ = parse(Float64, ls[1])
				l = Line(λ, ls)
				append!(atoms[currentName[]].lines, [l])
			end
		end
		
		LineList([atoms[name] for name in order])
	end
	if length(ix_pickOptVSurf_arr) > 0
		if occursin("spectrum", resolved_choice)
			linelist_paths = MUST.glob("*", "$(datafolder)/linelist")
			if length(linelist_paths) > 0
				linelists = readLinelist.(linelist_paths)
				linelist = LineList(vcat(getfield.(linelists, :atoms)...))
		
				atoms = getfield.(linelist.atoms, :name)
		
				#@show atoms
				
				# make clickable option available
				md"""
					Show selected lines from a linelist
				$(@bind lines_to_plot MultiCheckBox(atoms))
				"""
			end
		end
	end
end

# ╔═╡ 2dd31efe-7ca9-46d9-9f01-188dcf688b4f


# ╔═╡ c61e23a0-9b4b-492a-99f2-093f293d1d59
if length(ix_pickOptVSurf_arr) > 0
	pickOptVSurf
	redo
	removeSpectra2

	xlab = [@sprintf("%.3f", xi) for xi in x]
	ylab = [@sprintf("%.3f", yi) for yi in y]
	
	md"""
		Modify curve appearance
	$(@bind plot_names set_selected_names(zip(ix_pickOptVSurf_arr, iy_pickOptVSurf_arr), xlab, ylab))
	"""
end

# ╔═╡ 9db76078-3d67-47bd-bbb6-b16035a56912
yloc = if length(ix_pickOptVSurf_arr) > 0
	resolvedSelection[resolved_choice].yAxis()
else
	nothing
end;

# ╔═╡ e849b063-7500-44a8-8bfe-a757a661b60e
let
	pickOptVSurf
	redo
	removeSpectra2
	
	if length(iy_pickOptVSurf_arr) > 0
		#@show ix_pickOptVSurf, iy_pickOptVSurf
		#@show surfacesSliceSelection[surface_choice].zAxis[ix_pickOptVSurf, iy_pickOptVSurf]

		p = defaultPlot(scatter)
		xi = resolvedSelection[resolved_choice].xAxis
		
		# add linelist if wanted
		annotation = []
		if occursin("spectrum", resolved_choice)
			if length(lines_to_plot) > 0
				atoms = getfield.(linelist.atoms, :name)
				for atomname in lines_to_plot
					iatom = findfirst(atoms .== atomname)
					for line in linelist.atoms[iatom].lines
						if (line.λ > first(xi)) && (line.λ < last(xi))
							add_vline!(
								p, line.λ, 
								line_color="white",
								line_width=1,
								opacity=0.3
							)
							append!(
								annotation,
								[attr(
									showarrow=false,
									x=line.λ,
									y=1.075,
									yref="paper",
									text=atomname,
									xanchor="center",
									opacity=0.3
								)]
							)
						end
					end
				end
			end
		end

		i = 1
		for i in eachindex(iy_pickOptVSurf_arr)
			# spectrum here
			yi = yloc[i]
			add_trace!(
				p,
				scatter(
					x=xi, y=yi, name=plot_names[i], mode="lines"
				)
			)

			i +=1
		end

		relayout!(
			p, 
			xaxis_title=resolvedSelection[resolved_choice].xlabel, 
			yaxis_title=resolvedSelection[resolved_choice].label, 
			title=resolved_choice,
			legend=attr(yanchor="top", y=-0.1, xanchor="left", x=0.7),
			annotations=annotation
		)
		
		p
	else
		nothing
	end
end

# ╔═╡ e754eb4b-9304-440b-a306-6efa4f04b40b
if length(iy_pickOptVSurf_arr) > 0 
	let
		xlab = resolvedSelection[resolved_choice].xlabel
		ylab = resolvedSelection[resolved_choice].label
	
		data = JSON.json(
			Dict(
				xlab=>resolvedSelection[resolved_choice].xAxis,
				Dict(
					"$(ylab)-$(pl)"=>yloc[i] for (i, pl) in enumerate(plot_names)
				)...
			)
		)

		@info "You can download the data in the JSON format by clicking the button above."
		DownloadButton(data, replace(
			resolved_choice, 
			' '=>'_',
			'('=>"",
			')'=>""
		)*".json")
	end
end

# ╔═╡ Cell order:
# ╟─c7dc3b15-6555-4824-872a-d487fe5145ea
# ╟─2f1edd2a-b56e-11ee-29e7-c353938e7088
# ╟─6754b2c3-d205-4a12-88b3-53fe62c5637f
# ╟─41f0864e-26ee-46a6-b4ab-c401a4712941
# ╟─c596d1b3-32c5-4651-a8c1-3100fcd6cd59
# ╟─8af0c339-237c-42ca-bda5-0b0e44f11c30
# ╟─3c05fe1e-1c19-4c30-b34e-a874f63b91bc
# ╟─66e8e497-c3d4-48e7-9c2a-9de4d52b3dd4
# ╟─431bee5c-bfc1-4f53-8170-3ebee1277a2d
# ╟─ae5ab0a3-89db-4f1b-8e36-027e84d65c45
# ╟─8252972b-e496-4f13-b5d6-39e441fb9b84
# ╟─b76f0efd-6b2f-47d3-84e7-ca9fb866292e
# ╟─4eef1012-c8cd-4579-9a7e-ac2b4d4c7af6
# ╟─33c1da97-0760-495e-abd5-65531d5e1170
# ╟─8d6d674b-153b-4357-9f2d-c4e3cb05d059
# ╟─e938c3cf-e473-44ee-978c-4c0135856cab
# ╟─e61334ad-c74c-4334-87a3-46c15f622ce6
# ╟─b9a721cf-46ef-4e3c-a37c-8b35653e31cb
# ╟─a34793ae-db12-4dac-b2f8-348a88092815
# ╟─27885149-1f0f-4c75-ae98-96abd9509ddc
# ╟─333cbf0a-5a20-4612-8aa0-299e1ff40359
# ╟─1ef75139-4ee0-44fc-a731-a341ccd48e74
# ╟─aafbb147-e6ed-448f-bf4c-92f6c88d5199
# ╟─8a764271-7957-4a47-a979-0340a2ae7408
# ╟─ab72c112-5f17-4e71-a49e-5b2ffe2c6c20
# ╟─ada8d5a9-9a30-43a6-af66-050bcd000952
# ╟─2182de0b-bb9a-4fb2-b9df-45ef4ecd084e
# ╟─de173b50-fc57-4ac1-bd82-079921a29da5
# ╟─0639ce7d-955d-448f-84a0-353dfd4e93a3
# ╟─f1e4cbf4-0adb-4156-b696-5ea3801c6ae2
# ╟─585a336e-9c90-4b1b-93ac-8d81ff35d49d
# ╟─c82a3030-2bba-49ee-9eb6-cf07e901438d
# ╟─0584c2c0-de6e-43a9-a7cc-3dfc0742de21
# ╟─c1ca0584-92cd-497b-a8bd-8ca50ac84631
# ╟─3b070fdb-12c0-4ee8-b2dd-6d69fc972b81
# ╟─4eb0810c-b724-4284-8e97-5eb0b3d357ed
# ╟─0eb83755-4f86-4771-bcab-c7bbba096fb2
# ╟─1691b6fd-8d7d-46b6-a692-5b0f5a986e91
# ╟─63ec2623-eb88-43f1-a037-b86c013181f2
# ╟─2dd31efe-7ca9-46d9-9f01-188dcf688b4f
# ╟─c61e23a0-9b4b-492a-99f2-093f293d1d59
# ╟─9db76078-3d67-47bd-bbb6-b16035a56912
# ╟─e849b063-7500-44a8-8bfe-a757a661b60e
# ╟─e754eb4b-9304-440b-a306-6efa4f04b40b