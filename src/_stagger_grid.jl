# ============================================================================
# Type definitions
# ============================================================================

struct Atmos1DGrid{DF<:DataFrame} <:AbstractMUSTGrid
    name ::String
	path ::String
    info ::DF
end

# ============================================================================
# Constructors
# ============================================================================

Atmos1DGrid(path) = begin
    info = CSV.read(path, DataFrame)
    name = basename(path)[1:findlast('.', basename(path))-1]
	path = abspath(path)

    Atmos1DGrid(name, path, info)
end

StaggerGrid(args...; kwargs...) = Atmos1DGrid(args...; kwargs...)

# ============================================================================
# Checking paths
# ============================================================================

check_av_model_path(g::Atmos1DGrid) = begin
	ga = deepcopy(g)
	absolute_path!(ga)

	# check if av_models dir exists
	d = dirname(ga["av_path", 1])
	gridname = basename(dirname(d))
	if (!isdir(d)) && (isfile(d*".tar.gz"))
		_extract_tar_in_dir(dirname(d), name="av_models.tar.gz")
		@info "av_models unpacked."
	elseif !isdir(d) && (gridname=="MARCS_2.0")
		destination = d*".tar.gz"
		Downloads.download("https://keeper.mpdl.mpg.de/f/de46886c5929462ea1bd/?dl=1", destination)

		_extract_tar_in_dir(dirname(d), name="av_models.tar.gz")

		@info "av_models unpacked."
	elseif !isdir(d) && (gridname=="MARCS_2.1")
		destination = d*".tar.gz"
		Downloads.download("https://keeper.mpdl.mpg.de/f/46be98bda0974df1bb89/?dl=1", destination)

		_extract_tar_in_dir(dirname(d), name="av_models.tar.gz")

		@info "av_models unpacked."
	elseif !isdir(d)
		@warn "$(d) does not exist, and there is no $(d*".tar.gz") to extract."
	end
end

_extract_tar_in_dir(d; name="av_models.tar.gz") = begin
	c = ["-xzvf", name]
	currd = pwd()
	cd(d)
	run(`tar $(c)`)
	cd(currd)
end

# ============================================================================
# Getting info from the grid
# ============================================================================

get_mask(g::Atmos1DGrid; paras...) = begin
	grid = g.info
	mask = trues(nrow(grid))
	for (pname, val) in paras
		mask .= mask .& (grid[!, String(pname)] .== val)
	end

	mask
end

subselect(g::Atmos1DGrid; paras) = begin
	mask = get_mask(g; paras...)
	info_new = deepcopy(g.info[mask, :])
	Atmos1DGrid(g.name, g.path, info_new)
end

# ============================================================================
# Saving
# ============================================================================

save(grid::Atmos1DGrid, path) = CSV.write(path, grid.info)
save(grid::Atmos1DGrid) = CSV.write(grid.path, grid.info)

# ============================================================================
# Convenience
# ============================================================================

Base.getindex(g::Atmos1DGrid, k::String, i=!) = g.info[i, k]
Base.:+(g1::Atmos1DGrid, g2::Atmos1DGrid) = Atmos1DGrid(g1.name, g1.path, vcat(g1.info, g2.info))
Base.dirname(g::Atmos1DGrid) = dirname(g.path)

# ============================================================================
# Paths
# ============================================================================

absolute_path!(grid; from=@__FILE__) = begin
    if "av_path" in names(grid.info)
		grid.info[!, "av_path"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "av_path"]]
	end
	if "avo_path" in names(grid.info)
    	grid.info[!, "avo_path"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "avo_path"]]
    end
	if "eos_root" in names(grid.info)
		grid.info[!, "eos_root"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "eos_root"]]
	end
	if "matching_eos" in names(grid.info)
		grid.info[!, "matching_eos"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_eos"]]
	end
	if "matching_opa" in names(grid.info)
		grid.info[!, "matching_opa"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_opa"]]
	end
	if "matching_sopa" in names(grid.info)
		grid.info[!, "matching_sopa"] .= [isabspath(p) ? p : abspath.(joinpath.(dirname(from), p)) for p in grid.info[!, "matching_sopa"]]
	end
end

relative_path(relative_to, abs_path) = relpath(abs_path, dirname(relative_to))

relative_path!(grid; to=@__FILE__) = begin
    if "av_path" in names(grid.info)
		grid.info[!, "av_path"] .= relative_path.(to, grid.info[!, "av_path"])
	end
	if "avo_path" in names(grid.info)
    	grid.info[!, "avo_path"] .= relative_path.(to, grid.info[!, "avo_path"])
    end
	if "eos_root" in names(grid.info)
		grid.info[!, "eos_root"] .= relative_path.(to, grid.info[!, "eos_root"])
	end
	if "matching_eos" in names(grid.info)
		grid.info[!, "matching_eos"] .= relative_path.(to, grid.info[!, "matching_eos"])
	end
	if "matching_opa" in names(grid.info)
		grid.info[!, "matching_opa"] .= relative_path.(to, grid.info[!, "matching_opa"])
	end
	if "matching_sopa" in names(grid.info)
		grid.info[!, "matching_sopa"] .= relative_path.(to, grid.info[!, "matching_sopa"])
	end
end

# ============================================================================
# Replace Path
# ============================================================================

function replace_relative_path!(grid; kwargs...)
	# make paths absolute
	absolute_path!(grid)

	# replace paths with new ones
	replace_absolute_path!(grid; kwargs...)

	# make paths relative again
	relative_path!(grid)

	grid
end

function replace_absolute_path!(grid; kwargs...)
	# paths are already absolute
	for (k, v) in kwargs
		grid.info[k] = v
	end

	grid
end

# ============================================================================
# Interpolation
# ============================================================================

function interpolate_quantity(grid::Atmos1DGrid, what; teff, logg, feh, method="linear")
	feh_gr  = grid.info[!, "feh"]
	femask = feh_gr .≈ feh_gr[argmin(abs.(feh .- feh_gr))]

	logg_gr = grid.info[femask, "logg"]
	teff_gr = grid.info[femask, "teff"]
	what_gr = grid.info[femask, what]

	a = pyconvert(
		Any,
		first(
			scipy_interpolate.griddata(
				(logg_gr, teff_gr), what_gr, ([logg], [teff]),
				method=method
			)
		),
	)

	if isnan(a)
		pyconvert(
			Any,
			first(
				scipy_interpolate.griddata(
					(logg_gr, teff_gr), what_gr, ([logg], [teff]),
					method="nearest"
				)
			),
		)
	else
		a
	end
end

# ============================================================================
# Running
# ============================================================================

allowed_namelists(grid::Atmos1DGrid) = grid.info[!,"namelist_name"]

stage_namelists(grid::Atmos1DGrid, folder="run_grid"; clean_logs=true, clean_namelists=false) = begin
	nml_names = basename.(grid["namelist_name"])
	nml_path = grid["namelist_name"]
	folder_run = @in_dispatch(folder)

	if !isdir(folder_run)
		mkdir(folder_run)
	end

	# clean up dispatch from logs and namelists
	glob("*", folder_run) .|> rm

	for (i, n) in enumerate(nml_names)
		@info "copy $(nml_path[i]) to $(joinpath(folder_run, n))."
		cp(nml_path[i], joinpath(folder_run, n))
	end

	if clean_namelists
		glob("grid_*.nml", @in_dispatch("")) .|> rm
		for (i, n) in enumerate(nml_names)
			cp(joinpath(folder_run, n), joinpath(@in_dispatch(""), n))
		end
	end
	if clean_logs
		glob("grid_*.log", @in_dispatch("")) .|> rm
		glob("grid_*.err", @in_dispatch("")) .|> rm
	end
end

# ============================================================================
# Information of grid models from Name
# ============================================================================
