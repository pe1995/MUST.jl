include("_patch_meta.jl")
include("_box_quants.jl")


"""
    _build_cube_direct(patches)

Build 3D cube from a list of patches. Prepare axis ranges
"""
function _build_cube_direct(patches)
	# (x, y, z), (li, ui), (patches...)
	patch_range = zeros(Int, 3, 2, length(patches))

	global_x = []
	global_y = []
	global_z = []
	for i in eachindex(patches)
		_add_if_not_exists!(global_x, patches[i].xi)
		_add_if_not_exists!(global_y, patches[i].yi)
		_add_if_not_exists!(global_z, patches[i].zi)
	end

	global_x .= sort(global_x)
	global_y .= sort(global_y)
	global_z .= sort(global_z)

	p = zeros(3)
	for i in eachindex(patches)
		p[1] = first(patches[i].xi)
		p[2] = first(patches[i].yi)
		p[3] = first(patches[i].zi)
		patch_range[:, 1, i] .= _find_in_meshgrid(p, global_x, global_y, global_z)

		p[1] = last(patches[i].xi)
		p[2] = last(patches[i].yi)
		p[3] = last(patches[i].zi)
		patch_range[:, 2, i] .= _find_in_meshgrid(p, global_x, global_y, global_z)
	end
	
	global_x, global_y, global_z, patch_range
end


"""
    _var_from_patch_direct(var, fname, shp, off, li, ui, idxd)

Read data from the binary file as a mmap. Directly call it into RAM.
This is always okay because it is done patchwise.
"""
function _var_from_patch_direct(var, fname, shp, off, li, ui, idxd)
    varidx = idxd[String(var)] + 1
    Mmap.mmap(
		fname, 
		Array{Float32, length(shp)}, 
		shp, 
		off[varidx]
	)[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
end





"""
    prepare_memory(dtype, shape, vars; folder, number, mmap=true)  
    
Prepare big arrays for the entire cube. If you don't have enough RAM use 
`mmap=true`. In this case there will be HDF5 backed file and the data
will be written piece-wise. Note that this is much slower!
"""
function prepare_memory(dtype, shape, vars; folder, number, mmap=false)
    if mmap
        _prepare_memory_mmap(dtype, shape, vars; folder, number)
    else
        _prepare_memory_array(dtype, shape, vars; folder, number)
    end
end

function _prepare_memory_mmap(dtype, shape, vars; folder, number)
    path = joinpath(folder, "box_sn$(number).hdf5")
    fid = HDF5.h5open(path, "w")

    for v in vars
        dset = HDF5.create_dataset(
            fid, 
            v, 
            datatype(dtype), 
            dataspace(shape...), alloc_time = HDF5.H5D_ALLOC_TIME_EARLY
        )
    end
    HDF5.create_dataset(
        fid, 
        "x", 
        datatype(dtype), 
        dataspace(shape...), alloc_time = HDF5.H5D_ALLOC_TIME_EARLY
    )
    HDF5.create_dataset(
        fid, 
        "y", 
        datatype(dtype), 
        dataspace(shape...), alloc_time = HDF5.H5D_ALLOC_TIME_EARLY
    )
    HDF5.create_dataset(
        fid, 
        "z", 
        datatype(dtype), 
        dataspace(shape...), alloc_time = HDF5.H5D_ALLOC_TIME_EARLY
    )

    p = AtmosphericParameters(
        Base.convert(dtype, -99.0), 
        Base.convert(dtype, -99.0), 
        Base.convert(dtype, -99.0), 
        Dict{Symbol, dtype}()
    )
    save(p, fid)
    
    fid
end

function _prepare_memory_array(dtype, shape, vars; folder, number)
    fid = Dict()
    
    for v in vars
        fid[String(v)] = zeros(dtype, shape...)
    end
    fid["x"] = zeros(dtype, shape...)
    fid["y"] = zeros(dtype, shape...)
    fid["z"] = zeros(dtype, shape...)
    
    fid
end






"""
    _check_files(iout, data, run)

Check basic folders/ files exist for the given run.
"""
function _check_files(iout, data, run, rundir)
    run = run=="" ? splitpath(rundir)[end] : run
	rundir = isnothing(rundir) ? joinpath(data, run) : rundir
	if !isdir(rundir)
		error("directory '$(rundir)' does not exist")
	end

	datadir = joinpath(rundir, @sprintf("%05d", iout))
	if !isdir(datadir)
		error("directory '$(datadir)' does not exist")
	end

	file = joinpath(rundir, "params.nml")
    params_list = FreeNamelist(file)
	files = [f for f in glob("*", datadir) if occursin("snapshot.nml", f)]
	
	if length(files)==0
		error("directory '$(datadir)' has no snapshot.nml file")
	end

	run, rundir, datadir, params_list, files
end





"""
    _snapshot_variables(files)

Extract available variables and their index from the namelists.
"""
function _snapshot_variables(files)
	snapshot_nmls = []
	idxs = Dict()
	for f in files
		nml_list = FreeNamelist(f)
		append!(snapshot_nmls, [nml_list])
		
		if :idx_nml in nmlFields(nml_list)
			idx_dict = nmlField(nml_list, :idx_nml)

			for key in keys(idx_dict)	
				if idx_dict[key] >= 0
					idxs[key] = idx_dict[key] 
				end
			end
		end
	end
	variables = keys(idxs) |> collect
	idxVars = [idxs[k] for k in variables]
	masks = sortperm(idxVars)
	variables = variables[masks]
	idxVars = idxVars[masks]
	variablesSym = Symbol.(variables)
	snapshot_nml = first(snapshot_nmls)

	snapshot_nml, variablesSym, idxVars, idxs
end





"""
    _squaregaseos(run)

Load the MUST EoS from the given run. If applicabale, please
consider using the TSO EoS instead, which is faster and more convenient.
"""
function _squaregaseos(run; inputNamelist=@in_dispatch(run), eos_path=nothing)
	eos_path = if isnothing(eos_path)
		if isfile(inputNamelist*".nml")
			inputNamelist = StellarNamelist(inputNamelist*".nml")    
			@in_dispatch(replace(
				nmlField(inputNamelist, :eos_params)["table_loc"], 
				"'"=>""
			))
		else
			nothing
		end
	else
		eos_path
	end
	if (!isnothing(eos_path)) && isdir(eos_path)
		SquareGasEOS(eos_path), [:T, :kr, :Pg, :Ne]
	else
		@warn "No EoS found!"
		nothing, []
	end
end






"""
    _patchmeta(datadir, rank_nmls, snapshot_nml)

Collect meta data from all ranks.
"""
function _patchmeta(datadir, rank_nmls, snapshot_nml)
	patchMeta = []
	patchDataFiles = []
	patchAuxFiles = []
	auxnames = []
	save_aux = true
	for rank in eachindex(rank_nmls)
		rank_nml = rank_nmls[rank]
		filename = joinpath(datadir, "snapshot_"*@sprintf("%05d", rank-1)*".dat")
		
		for p in eachindex(rank_nml)
			pnml = rank_nml[p]
			id = nmlValue(pnml, "id")
			
			auxname = joinpath(datadir, @sprintf("%05d", id)*".aux")
			if save_aux
				aux = readaux(auxname)
				append!(auxnames, [a.name for a in aux])
				save_aux = false
			end

			append!(patchMeta, [PatchMeta(pnml, snapshot_nml)])
			append!(patchDataFiles, [filename])
			append!(patchAuxFiles, [auxname])
		end
	end

	patchMeta, patchDataFiles, patchAuxFiles, auxnames
end

"""
    _patchdata!(args...)

Assemble all data from all patches and save them in the big storage arrays.
"""
function _patchdata!(temp_storage, r, patchMeta, patch_range, variablesSym, patchDataFiles, idxs, auxnames, patchAuxFiles, variablesAndAux, variablesDerived, quantities, eos_quantities, eos_sq, units, converters, l_conv, fid, lookup_generator=nothing)
	# Lookup function can be generated from the input generator
	# this is a useful hook to use EoS outside of MUST (e.g. TSO)
	# without needing to include TSO inside MUST
	lookup_function, lkp = if isnothing(lookup_generator)
		lookup_function = Dict(
			q=>(d, ee) -> lookup(
				eos_sq, Symbol(q), d, ee;
				to_log=false
			) for q in eos_quantities
		)

		lkp = (f, d, ee) -> f(d, ee)

		lookup_function, lkp
	else
		lookup_function =  Dict(
			q=>lookup_generator(eos_sq, Symbol(q)) |> first for q in eos_quantities
		)
		
		lkp = lookup_generator(eos_sq, Symbol(first(eos_quantities))) |> last

		lookup_function, lkp
	end
	
	for (i, patch) in enumerate(patchMeta)
		r .= patch_range[:, :, i]
		li = patch.li
		ui = patch.ui
		for (j, var) in enumerate(variablesSym)
			temp_storage[var] .= _var_from_patch_direct(
				var,
				patchDataFiles[i],
				Tuple(patch.ncell),
				patch.offset,
				li,
				ui,
				idxs
			) 
		end

		xx, yy, zz = meshgrid(patch.xi, patch.yi, patch.zi)
		
		if length(auxnames) > 0
			auxname = patchAuxFiles[i]
			aux = readaux(auxname)

			for a in aux
				asym = Symbol(a.name)
				if asym in variablesAndAux
					temp_storage[asym] .= a.data[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] 
				end
			end
		end

		# Check if there are derived variables that we need to add
		for d in variablesDerived
			if !(String(d.name) in auxnames)
				temp_storage[d.name] .= d.recipe(; 
					temp_storage...
				)
			end
		end

		# Check if there are convertsion that we need to apply
		for d in quantities
			temp_storage[d.name] .*= converters[d.name]
		end

		xx .*= l_conv(units)
		yy .*= l_conv(units)
		zz .*= l_conv(units)

		for q in eos_quantities
			temp_storage[q] .= lkp(
				lookup_function[q],
				log.(temp_storage[:d]), log.(temp_storage[:ee])
			)
		end

		for d in quantities
			var = String(d.name)
			_save_to_structure!(fid[var], r, temp_storage[d.name])
		end

		for q in eos_quantities
			_save_to_structure!(fid[String(q)], r, temp_storage[q])
		end

		_save_to_structure!(fid["x"], r, xx)
		_save_to_structure!(fid["y"], r, yy)
		_save_to_structure!(fid["z"], r, zz)
	end
end







@inline function _save_to_structure!(fid, r, arr)
    fid[r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] = arr
end
@inline function _save_to_structure!(fid::AbstractArray, r, arr)
    fid[r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] .= arr
end







function _save_box(number, fid, time, logg, folder)
    HDF5.delete_object(fid, "time")
    HDF5.delete_object(fid, "logg")
    fid["time"] = time
    fid["logg"] = logg

	# see if we can get teff
	teff = _get_teff_from_fid(fid)

	HDF5.delete_object(fid, "teff")
    fid["teff"] = teff
    
    close(fid)
end

function _save_box(number, fid::Dict, time, logg, folder)
    dtype = eltype(fid["x"])

	# see if we can get teff
	teff = _get_teff_from_fid(fid)

    p = AtmosphericParameters(
        Base.convert(dtype, time), 
        Base.convert(dtype, teff), 
        Base.convert(dtype, logg), 
        Dict{Symbol, dtype}()
    )

    data = Dict(Symbol(p)=>v for (p, v) in fid if !(p in ["x", "y", "z"]))
    b = Box(
        fid["x"],
        fid["y"],
        fid["z"],
        data,
        p
    )

    save(b, number, folder=folder)
end

function _to_box(number, fid::Dict, time, logg, folder)
    dtype = eltype(fid["x"])

	# see if we can get teff
	teff = _get_teff_from_fid(fid)

    p = AtmosphericParameters(
        Base.convert(dtype, time), 
        Base.convert(dtype, -99.0), 
        Base.convert(dtype, logg), 
        Dict{Symbol, dtype}()
    )

    data = Dict(Symbol(p)=>v for (p, v) in fid if !(p in ["x", "y", "z"]))
    Box(
        fid["x"],
        fid["y"],
        fid["z"],
        data,
        p
    )
end

_get_teff_from_fid(fid, flux_name="flux") = begin
	if flux_name in keys(fid)
		(mean(fid[flux_name][:, :, end]) /σ_S) ^0.25
	else
		-99.0
	end
end





"""
    Box(iout::Int; run="", data=@in_dispatch("data"), rundir=nothing, quantities=defaultQuantities, eos_reader=_squaregaseos, lookup_generator=nothing, mmap=false)

Read and convert data from the given run, saved in the given data folder. This function entirely relies 
on julia functions and does not use the dispatch python interface. 
Quantities by default are d, e, ux, uy, uz, ee, dt, flux, and whatever is available in the EoS. To add more quantites see the `MUST.Quantity` type.
__Note__: TSO hooks are included to use the TSO EoS interface directly by passing the suitable `eos_reader` and `lookup_generator`.
This is implemented for you and ready to be used in the ingredient `convert2must.jl`. Without using this ingredient, the dafault setup will be used,
which is the `MUST` EoS reader.
"""
function Box(iout::Int; run="", data=@in_dispatch("data"), rundir=nothing, quantities=defaultQuantities, eos_path=nothing, eos_reader=(x)->_squaregaseos(x, eos_path=eos_path), lookup_generator=nothing, mmap=false, save_snapshot=false)	
	# check the inputs
    run, rundir, datadir, params_list, files = _check_files(iout, data, run, rundir)

	mmap = if !save_snapshot 
		false
	else
		mmap
	end

	# load the units and prepare the converters
	units = StandardUnits(rundir)
	l_conv = standardConversion(:l)
	converters = Dict(
		q.name=>q.conversion(units)
		for q in quantities
	)

	# add parameter groups from data/run/NNNNN/*snapshot.nml
	snapshot_nml, variablesSym, idxVars, idxs = _snapshot_variables(files)

	# rank namelists
    rank_files = [f for f in glob("*", datadir) if occursin("_patches.nml", f)]
	rank_nmls = FreeNamelist.(rank_files, :patch_nml)

	# general quantities
    logg = Base.convert(Float32, log10.(abs.(nmlValue(params_list, "g_cgs"))))
	time = Base.convert(Float32, nmlValue(snapshot_nml, "time")) * units.t

    # Use the Squaregas EOS 
    eos_sq, eos_quantities = eos_reader(run)
	
	# First we read the meta data so that we can build the data cube
	patchMeta, patchDataFiles, patchAuxFiles, auxnames = _patchmeta(
		datadir, rank_nmls, snapshot_nml
	)

	# decide if there are quantities to skip because there is no aux data for them
	# and they can not be derived otherwise
	variableMaskOk = [ (!q.derived)|(q.name in Symbol.(auxnames))|(!isnothing(q.recipe)) for q in quantities]
	#if count(.!variableMaskOk) > 0
	#	@warn "Exluding $(varnames(quantities)[.!variableMaskOk]) because they are derived, not available as aux parameters and no recipe to derive them was given."
	#end
	quantities = quantities[variableMaskOk]

	# build the 3D cube from patch meta data
	global_x, global_y, global_z, patch_range = _build_cube_direct(patchMeta)
	shape = length.([global_x, global_y, global_z])
		
	# prepare the HDF5 file for the final box
	variablesAndAux = [varnames(quantities)..., eos_quantities...]
	variablesDerived = derived(quantities)
	fid = prepare_memory(
		Float32, shape, String.(variablesAndAux); folder=rundir, number=iout,
		mmap=mmap
	)
	
	# After building the cube we read the data
	# the data can directly be added to the hdf5 box
	r = zeros(Int, 3, 2)
	li = patchMeta[1].li
	ui = patchMeta[1].ui
	n = ui .- (li .- 1)
	temp_storage = Dict{Symbol, Any}(
		k=>zeros(Float32, n...) 
		for k in unique([variablesSym..., variablesAndAux...])
	)
	_patchdata!(
		temp_storage, 
		r,
		patchMeta, 
		patch_range, 
		variablesSym, 
		patchDataFiles, 
		idxs, 
		auxnames, 
		patchAuxFiles, 
		variablesAndAux, 
		variablesDerived, 
		quantities, 
		eos_quantities, 
		eos_sq,
		units,
		converters,
		l_conv,
		fid,
		lookup_generator
	)

	if save_snapshot
		_save_box(iout, fid, time, logg, rundir)
		# return the mmaped box
		Box(
			"box_sn$(iout)",
			folder=rundir
		)
	else
		_to_box(iout, fid, time, logg, rundir)
	end
end
