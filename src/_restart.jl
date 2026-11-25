# ============================================================================
# DISPATCH restarts from boxes
# ============================================================================

const restartGroup = "restart"

"""
    save_restart(box_path, rundir, iout)

Store the DISPATCH namelists of snapshot `iout` (snapshot.nml, patch namelists,
params.nml, rebinned_table.txt) in the group `restart` of the box file, so that
`write_restart` can rebuild a snapshot after the native output is deleted.
"""
function save_restart(box_path, rundir, iout)
	datadir = joinpath(rundir, @sprintf("%05d", iout))
	files = Dict(f=>joinpath(datadir, f) for f in readdir(datadir) if endswith(f, ".nml"))
	for f in ["params.nml", "rebinned_table.txt"]
		isfile(joinpath(rundir, f)) && (files[f] = joinpath(rundir, f))
	end

	HDF5.h5open(box_path, "r+") do fid
		haskey(fid, restartGroup) && HDF5.delete_object(fid, restartGroup)
		g = HDF5.create_group(fid, restartGroup)
		HDF5.attributes(g)["iout"] = iout
		for (name, path) in files
			g[name] = read(path, String)
		end
	end

	box_path
end

has_restart(box_path) = HDF5.h5open(fid->haskey(fid, restartGroup), box_path, "r")

"""
    write_restart(box_path, restart_dir)

Write a DISPATCH snapshot without guard zones from a box holding a `restart`
group. DISPATCH restarts from it with
`restart_params data='<dirname(restart_dir)>/', run='<basename(restart_dir)>', snapshot=<iout>`.
"""
function write_restart(box_path, restart_dir)
	files, iout = HDF5.h5open(box_path, "r") do fid
		g = fid[restartGroup]
		Dict(k=>read(g[k]) for k in keys(g)), read(HDF5.attributes(g)["iout"])
	end

	datadir = joinpath(restart_dir, @sprintf("%05d", iout))
	mkpath(datadir)
	for (name, text) in files
		dir = name in ("params.nml", "rebinned_table.txt") ? restart_dir : datadir
		write(joinpath(dir, name), text)
	end

	snapshot_nml, _, _, idxs = _snapshot_variables([joinpath(datadir, "snapshot.nml")])
	rank_files = sort([joinpath(datadir, f) for f in keys(files) if occursin("_patches.nml", f)])
	if length(rank_files) == 0
		rank_files = sort([joinpath(datadir, f) for f in keys(files) if occursin("snapshot_", f)])
	end
	rank_nmls = FreeNamelist.(rank_files, :patch_nml)
	patchMeta, patchDataFiles, _ = _patchmeta(datadir, rank_nmls, snapshot_nml)
	write(
		joinpath(datadir, "snapshot.nml"),
		replace(files["snapshot.nml"], r"(GUARD_ZONES\s*=\s*)T"i => s"\1F")
	)

	b = Box(box_path)
	global_x, global_y, global_z, patch_range = _build_cube_direct(patchMeta)
	if length.((global_x, global_y, global_z)) != size(b.data[:d]) || (b.z[1, 1, end] < b.z[1, 1, 1])
		error("box $(box_path) does not match the patch layout of its restart group")
	end

	units = StandardUnits(restart_dir)
	conversion = Dict(
		:d=>(r)->b.data[:d][r...] ./ units.d,
		:e=>(r)->b.data[:e][r...] ./ units.e,
		:px=>(r)->b.data[:d][r...] .* b.data[:ux][r...] ./ (units.d * units.u),
		:py=>(r)->b.data[:d][r...] .* b.data[:uy][r...] ./ (units.d * units.u),
		:pz=>(r)->b.data[:d][r...] .* b.data[:uz][r...] ./ (units.d * units.u),
		:bx=>(r)->b.data[:bx][r...] ./ units.B,
		:by=>(r)->b.data[:by][r...] ./ units.B,
		:bz=>(r)->b.data[:bz][r...] ./ units.B
	)

	mv = nmlValue(snapshot_nml, "mv")
	n = patchMeta[1].n
	variables = Dict(
		idxs[k]=>Symbol(lowercase(String(k)))
		for k in keys(idxs) if (idxs[k] < mv) && haskey(conversion, Symbol(lowercase(String(k))))
	)
	for iv in 0:mv-1
		haskey(variables, iv) ||
			@warn "variable $(iv) of snapshot $(iout) is not in the box; written as zeros"
	end

	for file in unique(patchDataFiles)
		open(file, "w") do io
			for i in findall(==(file), patchDataFiles)
				r = [patch_range[j, 1, i]:patch_range[j, 2, i] for j in 1:3]
				seek(io, (patchMeta[i].rec - 1) * mv * prod(n) * 4)
				for iv in 0:mv-1
					v = get(variables, iv, :none)
					data = haskey(conversion, v) ? conversion[v](r) : zeros(n...)
					write(io, Float32.(data))
				end
			end
		end
	end

	restart_dir
end
