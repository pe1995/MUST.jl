### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 5831f154-6c01-11ee-2397-ed194f5bc82f
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using CSV
	using MUST  # For the models
	using TSO   # For the EoS
	using DelimitedFiles
	using DataFrames
	using Glob
end

# ╔═╡ 7ba977c0-bc0d-4d71-9994-7f76ecc740a8
md"# Grid: From Stagger to Dispatch
This notebook intents to capture important quantities from the original Stagger grid, such that its parameters can be used to derive new input for other models in dispatch. For this purpose, we construct a automatic table here that stores all of this information, as well as average models, that can be used for interpolation later without involving the actual models at a later stage."

# ╔═╡ d1f0b74f-887c-44e6-935f-fc1ba94f1a68
md"## Location of Stagger grid
First we need to locate the Stagger models"

# ╔═╡ 95bce920-6553-427f-a741-46d1614e656c
begin
	folder = "/mnt/beegfs/gemini/groups/bergemann/users/shared-storage/bergemann-data/Stagger_zaz/grid/"
	listOfFolders = [f[1:end-1] for f in glob("*/", folder) if isdir(f)]
	listOfFolderNames = [f[first(findlast("/", f))+1:end] for f in listOfFolders]
	FolderNamesMask   = [f[1]=='t' for f in listOfFolderNames]
	listOfFolders     = listOfFolders[FolderNamesMask]
	listOfFolderNames = listOfFolderNames[FolderNamesMask];
end

# ╔═╡ a36f74bb-e049-4a78-a98e-f3d6b5bf01e8
md"## Store information
We construct a dataframe, to which we can store the information needed, as well as what model it is we used."

# ╔═╡ 941dfc62-8211-4c72-b766-67b540777f6a
begin
	meshFiles = ["" for _ in listOfFolders]
	snapshots = ["" for _ in listOfFolders]
	isPresent  = trues(length(meshFiles))
	
	for (i, f) in enumerate(listOfFolders)
	    mesh_file = glob("*.msh", f)
	    
	    mesh_file = if length(mesh_file) == 0
	        isPresent[i] = false
	        nothing
	    elseif length(mesh_file) > 1
	        @warn "More than one mesh file present in $(f)"
	        fist(mesh_file)
	    else
	        first(mesh_file)
	    end
	
	    meshFiles[i] = isnothing(mesh_file) ?  "" : basename(mesh_file)
	
	
	
	    snap_file = glob("*.aux", f)
	    
	    snap_file = if length(snap_file) == 0
	        isPresent[i] = false
	        nothing
	    elseif length(snap_file) > 1
	        rand(snap_file)
	    else
	        fist(snap_file)
	    end
	
	    snapshots[i] = isnothing(snap_file) ? 
		"" : 
		basename(snap_file)[1:findlast('.', snap_file |> basename)-1]
	end
	
	
	info = DataFrame(Dict(
		"name"=>listOfFolderNames[isPresent], 
		"folder"=>listOfFolders[isPresent], 
		"mesh"=>meshFiles[isPresent], 
		"snapshot"=>snapshots[isPresent])
	)
end

# ╔═╡ ae0eab4b-22c0-4147-b949-fa4fd84d687b
md"## Extract properties of selected models"

# ╔═╡ 4d1dcc2b-e22d-4f8b-837c-62a181f1d7a0
"""
	model_properties(name, path, snap_name)

Load the given Stagger snap with name `snap_name`, saved at path, as a `MUST.Box` and compute its resolution, as well as other parameters. Furthermore compute average models and store them under the same as the snap. Return a summary. 

Example:
```julia
for i in 1:nrow(info)
    @info "processing \$(info[i, "name"])"
    model_properties(info[i, "name"], info[i, "folder"], info[i, "snapshot"])
end
```
"""
function model_properties(name, path, snap_name)
    ## we get the temperature and logg from the name only at this stage
    teff =  parse(Float64, name[2:findfirst('g', name)-1]) .*100
    logg =  parse(Float64, name[findfirst('g', name)+1:findfirst('m', name)-1]) ./10
    feh  = -parse(Float64, name[findfirst('m', name)+1:end]) ./10

    b = MUST.Box(MUST.StaggerSnap(snap_name, path))
    T = MUST.plane_statistic(MUST.mean, b, :T)
    d = MUST.plane_statistic(MUST.mean, b, :d)
    z = MUST.axis(b, :z)

	# save the box in the multi format (for Nick)
	MUST.multiBox(b, "$(name)_stagger")

    resolution(model) = minimum(abs.(diff(MUST.axis(model, :z))))
    resolution_horizontal(model) = minimum(abs.(diff(MUST.axis(model, :x))))
	Δt(R, u, c=1.0) = c * R / u
    scaling(δt, goal_scale=1e-3) = exp10.(round(log10(δt / goal_scale), sigdigits=1))

    vmax = max(abs(maximum(b[:ux])), 
				abs(maximum(b[:uy])), 
				abs(maximum(b[:uz])))

    s = scaling(Δt(resolution(b), vmax, 0.1))

    #initial_model = SimpleInitialModel(teff, logg, z, T, log.(d))
    
    open("$(snap_name)_av.dat", "w") do f
        writedlm(f, [z T log.(d)])
    end

    ## save return the additional column values
    return (teff, logg, feh, 
            minimum(b[:x]), maximum(b[:x]), 
            minimum(b[:y]), maximum(b[:y]), 
            minimum(b[:z]), maximum(b[:z]), 
            min(minimum(b[:ux]), minimum(b[:uy]), minimum(b[:uz])),
            max(maximum(b[:ux]), maximum(b[:uy]), maximum(b[:uz])),
            s,
            resolution_horizontal(b),
            "$(snap_name)_av.dat")
end

# ╔═╡ 1086a2da-df30-4f07-863d-100590ce42f7
compute = true

# ╔═╡ f748b6fd-52b4-4367-a81f-5015246ea653
begin
	teff, logg, feh, mi_x, ma_x, mi_y, ma_y, mi_z, ma_z, vmin, vmax, tscale, hres, av_path = [], [], [], [], [], [], [], [], [], [], [], [], [], []
		
	if compute
		for i in 1:nrow(info)
		    @info "processing $(info[i, "name"])"
		    res = model_properties(
				info[i, "name"], 
				info[i, "folder"], 
				info[i, "snapshot"]
			)
		
		    append!(teff,    [res[1]])
		    append!(logg,    [res[2]])
		    append!(feh,     [res[3]])
		    append!(mi_x,    [res[4]])
		    append!(ma_x,    [res[5]])
		    append!(mi_y,    [res[6]])
		    append!(ma_y,    [res[7]])
		    append!(mi_z,    [res[8]])
		    append!(ma_z,    [res[9]])
		    append!(vmin,    [res[10]])
		    append!(vmax,    [res[11]])
		    append!(tscale,  [res[12]])
		    append!(hres,    [res[13]])
		    append!(av_path, [res[14]])
		end
	end
end

# ╔═╡ 00936c88-71ca-421a-80cc-c5c9b36a3ed4
md"Save in the `info` dataframe"

# ╔═╡ b99bf324-8f64-4beb-901e-d52c8ed6eea6
begin
	if compute
		info[!, "teff"]    = teff
		info[!, "logg"]    = logg
		info[!, "feh"]     = feh
		info[!, "mi_x"]    = mi_x
		info[!, "ma_x"]    = ma_x
		info[!, "mi_y"]    = mi_y
		info[!, "ma_y"]    = ma_y
		info[!, "mi_z"]    = mi_z
		info[!, "ma_z"]    = ma_z
		info[!, "vmin"]    = vmin
		info[!, "vmax"]    = vmax
		info[!, "tscale"]  = tscale
		info[!, "hres"]    = hres
		info[!, "av_path"] = av_path

		CSV.write("stagger_grid.mgrid", info)
		info
	end
end

# ╔═╡ 9cdefadb-0714-45e2-a538-b9a301c44035
grid = MUST.StaggerGrid("stagger_grid.mgrid")

# ╔═╡ 32523905-6be4-4d05-82a8-90120ef35d14
exampleRow = 6

# ╔═╡ 0198207a-cf6e-4765-8ac1-d35c4fb9a073
exampleVmax = grid.info[exampleRow, "vmax"]

# ╔═╡ f863519e-4621-472c-8f2d-3d84c6067b5f
exampleVmin = grid.info[exampleRow, "vmin"]

# ╔═╡ 9cef7095-eb84-414b-827c-8443d95f17cc
absoluteVmax = max(abs(exampleVmax), abs(exampleVmin))

# ╔═╡ f093783c-fb4c-4a26-9217-bb516e665d67
Δt(R, u, c=1.0) = c * R / u

# ╔═╡ 31247cf4-f6e7-4f54-b289-252843e25805
MUST.roundto(Δt(2.6e8, absoluteVmax, 0.75), 0.25)

# ╔═╡ Cell order:
# ╟─7ba977c0-bc0d-4d71-9994-7f76ecc740a8
# ╠═5831f154-6c01-11ee-2397-ed194f5bc82f
# ╟─d1f0b74f-887c-44e6-935f-fc1ba94f1a68
# ╠═95bce920-6553-427f-a741-46d1614e656c
# ╟─a36f74bb-e049-4a78-a98e-f3d6b5bf01e8
# ╟─941dfc62-8211-4c72-b766-67b540777f6a
# ╟─ae0eab4b-22c0-4147-b949-fa4fd84d687b
# ╟─4d1dcc2b-e22d-4f8b-837c-62a181f1d7a0
# ╠═1086a2da-df30-4f07-863d-100590ce42f7
# ╠═f748b6fd-52b4-4367-a81f-5015246ea653
# ╟─00936c88-71ca-421a-80cc-c5c9b36a3ed4
# ╠═b99bf324-8f64-4beb-901e-d52c8ed6eea6
# ╟─9cdefadb-0714-45e2-a538-b9a301c44035
# ╠═32523905-6be4-4d05-82a8-90120ef35d14
# ╠═0198207a-cf6e-4765-8ac1-d35c4fb9a073
# ╠═f863519e-4621-472c-8f2d-3d84c6067b5f
# ╠═9cef7095-eb84-414b-827c-8443d95f17cc
# ╠═f093783c-fb4c-4a26-9217-bb516e665d67
# ╠═31247cf4-f6e7-4f54-b289-252843e25805
