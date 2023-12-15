### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 8b1da31e-9b23-11ee-372a-7b8499963aa0
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using CSV
	using MUST  # For the models
	using TSO   # For the EoS
	using DelimitedFiles
	using DataFrames
	using Glob
	using Plots
end

# ╔═╡ 24d25bd6-3f93-49c6-8632-fd2222b545f2
md"# Example EoS
To compare the rosseland opacity from within the snapshot to the one from our opacities, we can load an example EoS. We can use this one for solar metallicity only."

# ╔═╡ 80d15e89-bc8b-44c8-9044-f4850c04bc57
eos = reload(
	SqEoS, "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6/combined_ross_eos_magg22.hdf5"
)

# ╔═╡ e42a941e-d24a-4962-8ed7-2b05add167f5
md"# Find Snapshots"

# ╔═╡ 9e105452-0791-4bab-9187-3e847ff7b926
begin
	folder = "/ptmp/peitner/Stagger_grid/"
	listOfFolders = [f[1:end-1] for f in glob("*/", folder) if isdir(f)]
	listOfFolderNames = [f[first(findlast("/", f))+1:end] for f in listOfFolders]
	FolderNamesMask   = [f[1]=='t' for f in listOfFolderNames]
	listOfFolders     = listOfFolders[FolderNamesMask]
	listOfFolderNames = listOfFolderNames[FolderNamesMask]
end

# ╔═╡ 32ea5bca-653e-4c7f-94a1-25de7a5b135e
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

# ╔═╡ afb08d1f-5ff4-467f-8d12-3927139b51bf
md"# Load example"

# ╔═╡ 4848a1c4-b827-4142-9d69-a092a53aea31
exampleId = 117

# ╔═╡ a9261b05-da84-4d21-b2b6-7a99fa88fc58
pathExample = info[exampleId, "folder"]

# ╔═╡ de19d6ad-a6cd-41d9-82c1-b3cef2955313
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
    logg =  try
		parse(Float64, name[findfirst('g', name)+1:findfirst('m', name)-1]) ./10
	catch
		parse(Float64, name[findfirst('g', name)+1:findfirst('p', name)-1]) ./10
	end
    feh  = try
		-parse(Float64, name[findfirst('m', name)+1:end]) ./10
	catch
		parse(Float64, name[findfirst('p', name)+1:end]) ./10
	end

    b = MUST.Box(MUST.StaggerSnap(snap_name, path))
    T = MUST.plane_statistic(MUST.mean, b, :T)
    d = MUST.plane_statistic(MUST.mean, b, :d)
    z = MUST.axis(b, :z)

	τ = MUST.optical_depth(b, opacity=:ross, density=:d)
	MUST.add!(b, τ, :τ_ross)
	bτ = MUST.height_scale_fast(b, :τ_ross)

	Tτ = MUST.plane_statistic(MUST.mean, bτ, :T)
    dτ = MUST.plane_statistic(MUST.mean, bτ, :d)
    zτ = MUST.plane_statistic(MUST.mean, bτ, :z)
    τ = MUST.axis(bτ, :τ_ross, 3)

	T, d, z, Tτ, dτ, zτ, τ 
end

# ╔═╡ 2dce8118-286b-4fb4-b293-ed69b3cac1a4
T, d, z, Tτ, dτ, zτ, τ = model_properties(
	info[exampleId, "name"], info[exampleId, "folder"], info[exampleId, "snapshot"]
)

# ╔═╡ 65fbee5f-2ba3-4dd2-b011-0e84657fae41
begin
	plot(zτ, Tτ)
	plot!(z, T)
end

# ╔═╡ 70b83ca5-b928-42dd-81ac-a1485a7d293c
begin
	plot(log10.(τ), Tτ)
end

# ╔═╡ Cell order:
# ╠═8b1da31e-9b23-11ee-372a-7b8499963aa0
# ╟─24d25bd6-3f93-49c6-8632-fd2222b545f2
# ╟─80d15e89-bc8b-44c8-9044-f4850c04bc57
# ╟─e42a941e-d24a-4962-8ed7-2b05add167f5
# ╟─9e105452-0791-4bab-9187-3e847ff7b926
# ╟─32ea5bca-653e-4c7f-94a1-25de7a5b135e
# ╟─afb08d1f-5ff4-467f-8d12-3927139b51bf
# ╠═4848a1c4-b827-4142-9d69-a092a53aea31
# ╠═a9261b05-da84-4d21-b2b6-7a99fa88fc58
# ╠═de19d6ad-a6cd-41d9-82c1-b3cef2955313
# ╠═2dce8118-286b-4fb4-b293-ed69b3cac1a4
# ╠═65fbee5f-2ba3-4dd2-b011-0e84657fae41
# ╠═70b83ca5-b928-42dd-81ac-a1485a7d293c
