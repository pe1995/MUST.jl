### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 847ff48c-f0df-11ee-2af0-ff2f7b38a5b4
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PythonPlot
	using TSO
end

# ╔═╡ e610e24b-5ae7-4650-bd5d-9a7929ba263a
plt = matplotlib.pyplot

# ╔═╡ 18e400aa-1f9d-4890-8ac4-75a2782f89f6
FortranFiles = MUST.FortranFiles

# ╔═╡ 810d479a-6a8b-4858-a207-6941e52b12e7
read = MUST.read

# ╔═╡ 00b86f70-8b57-49b3-bcec-bff646ace479
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ 06c56269-64c7-48da-9b1e-f1ff67e6696e
begin
	function parse_multiple(s::AbstractString)
		sComponents = []
		
		# split at commas
		sComma = strip.(split(s, ",", keepempty=false))
	
		# split for multiplication signs
		sMulti = [strip.(split(si, "*", keepempty=false)) for si in sComma]
	
		# flatten the list
		for component in sMulti
			c = if length(component) > 1
				repeat([last(component)], parse(Int, first(component)))
			else
				component
			end
	
			append!(sComponents, c)
		end
	
		# try to convert to number
		sComponentsParsed = []
		for component in sComponents
			c_new = if occursin(".", component)
				try 
					parse(Float64, component)
				catch
					component
				end
			else
				try 
					parse(Int, component)
				catch
					component
				end
			end
	
			append!(sComponentsParsed, [c_new])
		end
			
		sComponentsParsed
	end
	parse_multiple(s::AbstractArray) = s
end

# ╔═╡ a393d4e7-f72f-40fd-91f6-29f8009ea272
begin
	function read_stepwise(f, dtype)
		rec = FortranFiles.Record(f)
		res = []
		while true
			try
				a = read(rec, dtype)
				append!(res, [a])
			catch
				break
			end
		end
	
		close(rec)
		res
	end
	
	function read_string(f)
		str_arr = read_stepwise(f, FortranFiles.FString{1})
		str_arr = String[FortranFiles.trimstring(a) for a in str_arr]
		join(str_arr)
	end
	
	function read_numbers(f, dtype=Int32)
		a = dtype[read_stepwise(f, dtype)...]
		if length(a) == 1
			first(a)
		else
			a
		end
	end
end

# ╔═╡ d5aaa674-7af2-426a-9cc1-5c42f2c8248b
begin
	struct AuxPatch
		version
		id
		name
		rnk
		shp
		tp
		data
	end
	
	function readaux(fname)
		auxs = []
		if isfile(fname)
			f = FortranFiles.FortranFile(fname)
			version, id = read(f, (Int32, 2))
			while true
				try			
					name = read_string(f)
					rnk = read_numbers(f)
					shp = read_numbers(f)
					tp = read_string(f)
			
					data = if first(tp) == 'r'
						a = zeros(Float32, shp...)
						read(f, a)
						a
					elseif first(tp) == 'd'
						a = zeros(Float64, shp...)
						read(f, a)
						a
					else
						a = zeros(Int32, shp...)
						read(f, a)
						a
					end
	
					append!(auxs, [AuxPatch(version, id, name, rnk, shp, tp, data)])
				catch
					break
				end
			end
			close(f)
		end
		auxs
	end
end

# ╔═╡ 872a2102-da21-4357-954a-fa25152395e8
function _add_axes(llc, ng, gn, ds)
	first = llc .- ng .* ds
	n = gn
   
	x=first[1] .+ ds[1] .* range(0, n[1]-1)
	y=first[2] .+ ds[2] .* range(0, n[2]-1)
	z=first[3] .+ ds[3] .* range(0, n[3]-1)

	x, y, z
end

# ╔═╡ 229350b2-7594-4c56-b0ca-85c5113a0a54
begin
	struct PatchMeta
		id
		mv
		rec
		nv
		ncell
		li
		ui
		gn
		n
		ng
		offset
		s
		pos
		llc
		urc
		ds
		x
		y 
		z
		xi
		yi 
		zi
	end
	
	function PatchMeta(patch, snap, noaux=false)
		guards = MUST.nmlValue(snap, "guard_zones")
		n = parse_multiple(MUST.nmlValue(snap, "n"))
		
		li, ui = if !guards
	        [1, 1, 1], n
		else
			(
				parse_multiple(MUST.nmlValue(snap, "li")),
				parse_multiple(MUST.nmlValue(snap, "ui"))
			)
		end
	
		mv = MUST.nmlValue(snap, "mv")
		nv = MUST.nmlValue(snap, "nv")
		rec = MUST.nmlValue(patch, "record")
		ncell = parse_multiple(MUST.nmlValue(patch, "ncell"))
		gn = parse_multiple(MUST.nmlValue(snap, "gn"))
		ng = parse_multiple(MUST.nmlValue(snap, "ng"))
	    offset = [iv + (rec-1) * mv for iv in 0:mv-1] .* (4*prod(gn))
		id = MUST.nmlValue(patch, "id")
	
		# for axes
		s = parse_multiple(MUST.nmlValue(patch, "size"))
		pos = parse_multiple(MUST.nmlValue(patch, "position"))
		llc= pos .- s./2.0
		urc= pos .+ s./2.0
	
		ds = parse_multiple(MUST.nmlValue(patch, "ds"))
	
		x, y, z = _add_axes(llc, ng, gn, ds)
		xi, yi, zi = x[li[1]:ui[1]], y[li[2]:ui[2]], z[li[3]:ui[3]]
	
		PatchMeta(
			id,
			mv,
			rec,
			nv,
			ncell,
			li,
			ui,
			gn,
			n,
			ng,
			offset,
			s,
			pos,
			llc,
			urc,
			ds,
			x,
			y,
			z,
			xi|>collect,
			yi|>collect,
			zi|>collect
		)
	end
end

# ╔═╡ 0b0d9a9e-2251-4f6e-a761-2b247ea20471
function _build_cube(patches)
	# (x, y, z), (li, ui), (patches...)
	patch_range = zeros(Int, 3, 2, length(patches))

	global_x = []
	global_y = []
	global_z = []
	for i in eachindex(patches)
		MUST._add_if_not_exists!(global_x, patches[i].xi)
		MUST._add_if_not_exists!(global_y, patches[i].yi)
		MUST._add_if_not_exists!(global_z, patches[i].zi)
	end

	global_x .= sort(global_x)
	global_y .= sort(global_y)
	global_z .= sort(global_z)

	p = zeros(3)
	for i in eachindex(patches)
		p[1] = first(patches[i].xi)
		p[2] = first(patches[i].yi)
		p[3] = first(patches[i].zi)
		patch_range[:, 1, i] .= MUST._find_in_meshgrid(p, global_x, global_y, global_z)

		p[1] = last(patches[i].xi)
		p[2] = last(patches[i].yi)
		p[3] = last(patches[i].zi)
		patch_range[:, 2, i] .= MUST._find_in_meshgrid(p, global_x, global_y, global_z)
	end
	
	global_x, global_y, global_z, patch_range
end

# ╔═╡ f0c9dc10-8448-4e2c-9aca-317ee396ea92
function _var_from_patch_direct(var, fname, shp, off, li, ui, idxd)
    varidx = idxd[String(var)] + 1
    MUST.Mmap.mmap(
		fname, 
		Array{Float32, length(shp)}, 
		shp, 
		off[varidx]
	)[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
end

# ╔═╡ 1f02c11c-e8f5-45a9-b043-7e883fdbafd0
begin
	function prepare_memory(dtype, shape, vars; folder, number, mmap=true)
		if mmap
			prepare_memory_mmap(dtype, shape, vars; folder, number)
		else
			prepare_memory_array(dtype, shape, vars; folder, number)
		end
	end
	
	function prepare_memory_mmap(dtype, shape, vars; folder, number)
		path = joinpath(folder, "box_$(number).hdf5")
		fid = MUST.HDF5.h5open(path, "w")
	
		for v in vars
			dset = MUST.HDF5.create_dataset(
				fid, 
				v, 
				MUST.datatype(dtype), 
				MUST.dataspace(shape...), alloc_time = MUST.HDF5.H5D_ALLOC_TIME_EARLY
			)
		end
		MUST.HDF5.create_dataset(
			fid, 
			"x", 
			MUST.datatype(dtype), 
			MUST.dataspace(shape...), alloc_time = MUST.HDF5.H5D_ALLOC_TIME_EARLY
		)
		MUST.HDF5.create_dataset(
			fid, 
			"y", 
			MUST.datatype(dtype), 
			MUST.dataspace(shape...), alloc_time = MUST.HDF5.H5D_ALLOC_TIME_EARLY
		)
		MUST.HDF5.create_dataset(
			fid, 
			"z", 
			MUST.datatype(dtype), 
			MUST.dataspace(shape...), alloc_time = MUST.HDF5.H5D_ALLOC_TIME_EARLY
		)
	
		p = MUST.AtmosphericParameters(
			Base.convert(dtype, -99.0), 
			Base.convert(dtype, -99.0), 
			Base.convert(dtype, -99.0), 
			Dict{Symbol, dtype}()
		)
		MUST.save(p, fid)
		
		fid
	end
	
	function prepare_memory_array(dtype, shape, vars; folder, number)
		fid = Dict()
		
		for v in vars
			fid[String(v)] = zeros(dtype, shape...)
		end
		fid["x"] = zeros(dtype, shape...)
		fid["y"] = zeros(dtype, shape...)
		fid["z"] = zeros(dtype, shape...)
		
		fid
	end
end

# ╔═╡ 9c90d6a9-3249-49ed-ae1d-398157922bdd
function StandardUnits(name::String; folder=MUST.@in_dispatch("data"))
	rundir = joinpath(folder, name)
	file = joinpath(rundir, "params.nml")
    ps = MUST.FreeNamelist(file)

	scalingp = MUST.nmlField(ps, Symbol("scaling_params"))
    pnames = keys(scalingp) |> collect
	
    l_name = "l_cgs" 
    d_name = "d_cgs" 
    v_name = "v_cgs" 
    t_name = "t_cgs" 

    l = scalingp[l_name]
    d = scalingp[d_name]

	if t_name in pnames
        t = scalingp[t_name]
        v = l / t
	elseif v_name in pnames
        v = scalingp[v_name]
        t = l / v
	else
		error("Either t or v scaling needs to be provided")
    end

    MUST.StandardUnits(l=l,d=d,t=t,u=v)
end

# ╔═╡ 99464f53-968b-466a-a357-b72ec132f135
begin
	struct Quantity
		name
		recipe
		conversion
		derived
	end

	standardConversion(name) = begin
		(units) -> getfield(units, name)
	end

	Quantity(name) = Quantity(name, identity, standardConversion(name), false)
	Quantity(name, recipe) = Quantity(name, recipe, standardConversion(name), true)
	Quantity(name, recipe, conversion) = Quantity(name, recipe, conversion, true)

	derived(q) = [qi for qi in q if qi.derived]
	nonderived(q) = [qi for qi in q if !qi.derived]
	varnames(q) = [qi.name for qi in q]
end

# ╔═╡ 5816cfd3-4c3b-4320-b4ea-1185daae848e
begin
	dQ = Quantity(:d)
	eQ = Quantity(:e)
	uxQ = Quantity(:ux, (; d, px, kwargs...)->(px ./ d), standardConversion(:u))
	uyQ = Quantity(:uy, (; d, py, kwargs...)->(py ./ d), standardConversion(:u))
	uzQ = Quantity(:uz, (; d, pz, kwargs...)->(pz ./ d), standardConversion(:u))
	eeQ = Quantity(:ee, (; d, e, kwargs...)->(e ./ d), standardConversion(:ee))
	dtQ = Quantity(:dt_rt, identity, standardConversion(:t))
	fluxQ = Quantity(:flux, identity, standardConversion(:flux))
	
	defaultQuantities = [
		dQ,
		eQ,
		uxQ,
		uyQ,
		uzQ,
		eeQ,
		dtQ,
		fluxQ
	]
end

# ╔═╡ e0a79867-5207-4247-843f-62837596b2aa
function _check_files(iout, data, run)
	rundir = joinpath(data, run)
	if !isdir(rundir)
		error("directory '$(rundir)' does not exist")
	end

	datadir = joinpath(rundir, MUST.@sprintf("%05d", iout))
	if !isdir(datadir)
		error("directory '$(datadir)' does not exist")
	end

	file = joinpath(rundir, "params.nml")
    params_list = MUST.FreeNamelist(file)
	files = [f for f in MUST.glob("*", datadir) if occursin("snapshot.nml", f)]
	
	if length(files)==0
		error("directory '$(datadir)' has no snapshot.nml file")
	end

	rundir, datadir, params_list, files
end

# ╔═╡ 1adf1f53-f90b-4c9d-a670-998d432e19dd
function _snapshot_variables(files)
	snapshot_nmls = []
	idxs = Dict()
	for f in files
		nml_list = MUST.FreeNamelist(f)
		append!(snapshot_nmls, [nml_list])
		
		if :idx_nml in MUST.nmlFields(nml_list)
			idx_dict = MUST.nmlField(nml_list, :idx_nml)

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

# ╔═╡ 40eb7042-f7af-4520-b1bc-51998ba27e20
function _squaregaseos(run)
	inputNamelist = @in_dispatch run
	eos_path = if isfile(inputNamelist*".nml")
    	inputNamelist = MUST.StellarNamelist(inputNamelist*".nml")    
		replace(
			MUST.nmlField(inputNamelist, :eos_params)["table_loc"], 
			"'"=>""
		)
	else
		nothing
	end
	if (!isnothing(eos_path)) && isdir(@in_dispatch(eos_path))
		MUST.SquareGasEOS(@in_dispatch(eos_path)), [:T, :kr, :Pg, :Ne]
	else
		@warn "No EoS found!"
		nothing, []
	end
end

# ╔═╡ 9a9aab09-58db-49b5-9821-acf4a0260698
function _sqeos(run)
	inputNamelist = @in_dispatch run
	eos_path = if isfile(inputNamelist*".nml")
    	inputNamelist = MUST.StellarNamelist(inputNamelist*".nml")    
		p = replace(
			MUST.nmlField(inputNamelist, :eos_params)["table_loc"], 
			"'"=>""
		)
		joinpath(p, "eos.hdf5")
	else
		nothing
	end
	if (!isnothing(eos_path)) && isfile(@in_dispatch(eos_path))
		TSO.reload(SqEoS, @in_dispatch(eos_path)), [:T, :kr, :Pg, :Ne]
	else
		@warn "No EoS found!"
		nothing, []
	end
end

# ╔═╡ 2591487d-12a1-4324-85ce-3bcbbdaa734f
function _lookup_generator(eos, q)
	q = if q == :T
		:lnT
	elseif q == :Pg
		:lnPg
	elseif q == :Ne
		:lnNe
	elseif q == :kr
		:lnRoss
	else
		q
	end
		
	TSO.lookup_function(TSO.@axed(eos), q), (args...) -> exp.(TSO.lookup(args...))
end

# ╔═╡ be2e26c4-0a15-4108-a44a-2624973fe6ef
function _patchmeta(datadir, rank_nmls, snapshot_nml)
	patchMeta = []
	patchDataFiles = []
	patchAuxFiles = []
	auxnames = []
	save_aux = true
	for rank in eachindex(rank_nmls)
		rank_nml = rank_nmls[rank]
		filename = joinpath(datadir, "snapshot_"*MUST.@sprintf("%05d", rank-1)*".dat")
		
		for p in eachindex(rank_nml)
			pnml = rank_nml[p]
			id = MUST.nmlValue(pnml, "id")
			
			auxname = joinpath(datadir, MUST.@sprintf("%05d", id)*".aux")
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

# ╔═╡ a79da710-d267-4e88-8e38-1d07bea6ca53
to = TSO.TimerOutput()

# ╔═╡ c41689fe-2c35-456c-8ef3-0716b43dca25
begin
	@inline function _save_to_structure!(fid, r, arr)
		fid[r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] = arr
	end
	@inline function _save_to_structure!(fid::AbstractArray, r, arr)
		fid[r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] .= arr
	end
end

# ╔═╡ 78b4b739-6e6b-4621-a125-8538c2000c61
function _patchdata!(temp_storage, r, patchMeta, patch_range, variablesSym, patchDataFiles, idxs, auxnames, patchAuxFiles, variablesAndAux, variablesDerived, quantities, eos_quantities, eos_sq, units, converters, l_conv, fid, lookup_generator=nothing)
	# Lookup function can be generated from the input generator
	# this is a useful hook to use EoS outside of MUST (e.g. TSO)
	# without needing to include TSO inside MUST
	lookup_function, lkp = if isnothing(lookup_generator)
		lookup_function = Dict(
			q=>(d, ee) -> MUST.lookup(
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
		TSO.@timeit to "_var_from_patch_direct" for (j, var) in enumerate(variablesSym)
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

		TSO.@timeit to "meshgrid" xx, yy, zz = MUST.meshgrid(patch.xi, patch.yi, patch.zi)
		
		TSO.@timeit to "aux" if length(auxnames) > 0
			auxname = patchAuxFiles[i]
			aux = readaux(auxname)

			for a in aux
				if Symbol(a.name) in variablesAndAux
					temp_storage[Symbol(a.name)] .= a.data[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] 
				end
			end
		end

		# Check if there are derived variables that we need to add
		TSO.@timeit to "recipe" for d in variablesDerived
			if !(String(d.name) in auxnames)
				temp_storage[d.name] .= d.recipe(; 
					temp_storage...
				)
			end
		end

		# Check if there are convertsion that we need to apply
		TSO.@timeit to "converters" for d in quantities
			temp_storage[d.name] .*= converters[d.name]
		end

		xx .*= l_conv(units)
		yy .*= l_conv(units)
		zz .*= l_conv(units)

		TSO.@timeit to "lookup" for q in eos_quantities
			temp_storage[q] .= lkp(
				lookup_function[q],
				log.(temp_storage[:d]), log.(temp_storage[:ee])
			)
		end

		TSO.@timeit to "fid1" for d in quantities
			var = String(d.name)
			_save_to_structure!(fid[var], r, temp_storage[d.name])
		end

		TSO.@timeit to "fid2" for q in eos_quantities
			_save_to_structure!(fid[String(q)], r, temp_storage[q])
		end

		_save_to_structure!(fid["x"], r, xx)
		_save_to_structure!(fid["y"], r, yy)
		_save_to_structure!(fid["z"], r, zz)
	end
end

# ╔═╡ 01745045-8419-4d58-990d-3572ad1f47c2
begin
	function _save_box(number, fid, time, logg, folder)
		MUST.HDF5.delete_object(fid, "time")
		MUST.HDF5.delete_object(fid, "logg")
		fid["time"] = time
		fid["logg"] = logg
		
		close(fid)
	end
	function _save_box(number, fid::Dict, time, logg, folder)
		dtype = eltype(fid["x"])
		p = MUST.AtmosphericParameters(
			Base.convert(dtype, time), 
			Base.convert(dtype, -99.0), 
			Base.convert(dtype, logg), 
			Dict{Symbol, dtype}()
		)
	
		data = Dict(Symbol(p)=>v for (p, v) in fid if !(p in ["x", "y", "z"]))
		b = MUST.Box(
			fid["x"],
			fid["y"],
			fid["z"],
			data,
			p
		)
	
		MUST.save(b, number, folder=folder)
	end
end

# ╔═╡ 486a3d26-cafb-4dbd-b735-c0a3c89c02d9
function snapshot(iout; run="", data=MUST.@in_dispatch("data"), quantities=defaultQuantities, eos_reader=_squaregaseos, lookup_generator=nothing)	
	# check the inputs
    TSO.@timeit to "_check_files" rundir, datadir, params_list, files = _check_files(iout, data, run)

	# load the units and prepare the converters
	units = StandardUnits(run; folder=data)
	l_conv = standardConversion(:l)
	converters = Dict(
		q.name=>q.conversion(units)
		for q in quantities
	)

	# add parameter groups from data/run/NNNNN/*snapshot.nml
	TSO.@timeit to "_snapshot_variables" snapshot_nml, variablesSym, idxVars, idxs = _snapshot_variables(files)

	# rank namelists
    rank_files = [f for f in MUST.glob("*", datadir) if occursin("_patches.nml", f)]
	rank_nmls = MUST.FreeNamelist.(rank_files, :patch_nml)

	# general quantities
    logg = Base.convert(Float32, log10.(abs.(MUST.nmlValue(params_list, "g_cgs"))))
	time = Base.convert(Float32, MUST.nmlValue(snapshot_nml, "time")) * units.t
	
    # Use the Squaregas EOS 
    eos_sq, eos_quantities = eos_reader(run)
	
	# First we read the meta data so that we can build the data cube
	TSO.@timeit to "_patchmeta" patchMeta, patchDataFiles, patchAuxFiles, auxnames = _patchmeta(
		datadir, rank_nmls, snapshot_nml
	)

	# build the 3D cube from patch meta data
	TSO.@timeit to "_build_cube" global_x, global_y, global_z, patch_range = _build_cube(patchMeta)
	shape = length.([global_x, global_y, global_z])
		
	# prepare the HDF5 file for the final box
	variablesAndAux = [varnames(quantities)..., eos_quantities...]
	variablesDerived = derived(quantities)
	TSO.@timeit to "prepare_memory" fid = prepare_memory(
		Float32, shape, String.(variablesAndAux); folder=rundir, number=iout,
		mmap=false
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
	TSO.@timeit to "_patchdata" _patchdata!(
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

	_save_box(iout, fid, time, logg, rundir)
	
	# return the mmaped box
	MUST.Box(
		"box_sn$(iout)",
		folder=rundir
	)
end

# ╔═╡ 35102657-6a9a-4fbb-82b9-24f5289946c4
StandardUnits("grid_t5777g44m00_deep")

# ╔═╡ 71163085-0adb-49a3-8ce5-028d18c8f5d7
begin
	TSO.reset_timer!(to)
	b = snapshot(
		10, 
		run="grid_t5777g44m00_deep"
	)
	show(to)
end

# ╔═╡ ef238fdf-31eb-42c5-b91b-134196d5cf8f
b.parameter.time

# ╔═╡ 41c7c4ce-ff05-43ca-adc0-f1a4771c9844
b.parameter.logg

# ╔═╡ b3f76f05-57c7-4890-92a0-2f268f72ec5f
let
	plt.close()

	plt.plot(profile(MUST.mean, b, :z, :dt_rt)...)

	gcf()
end

# ╔═╡ 51ab063f-a103-4b77-80f6-b7b5601df785


# ╔═╡ 04701ec7-6831-4617-b81c-1518264f7aaa
md"# Test within MUST"

# ╔═╡ 19fc41ec-abcb-4c1d-9a6b-e16221ad7244
c2m = MUST.ingredients("convert2must.jl")

# ╔═╡ 14672d91-6b5d-4a0c-a2c2-e52150f766cb
eos = TSO.reload(SqEoS, @in_dispatch("input_data/grd/DIS_MARCS_E_t5777g44m00_v0.5.1/eos.hdf5"))

# ╔═╡ ff186ff3-c7de-4f93-a9b2-0331bcf2652b
opa = TSO.reload(SqOpacity, @in_dispatch("input_data/grd/DIS_MARCS_E_t5777g44m00_v0.5.1/binned_opacities.hdf5"))

# ╔═╡ b3f3425b-3ccb-4de5-b32f-a698df810ce3
b2, bt2 = MUST.Box(
	"grid_t5777g44m00_deep",
	10, 
	data=@in_dispatch("data"),
	legacy=false
)

# ╔═╡ 1da12e54-263d-41af-a760-3da1af97f9f3
let
	plt.close()

	plt.plot(profile(MUST.mean, bt2, :z, :dt_rt)...)

	gcf()
end

# ╔═╡ Cell order:
# ╠═847ff48c-f0df-11ee-2af0-ff2f7b38a5b4
# ╠═e610e24b-5ae7-4650-bd5d-9a7929ba263a
# ╠═18e400aa-1f9d-4890-8ac4-75a2782f89f6
# ╠═810d479a-6a8b-4858-a207-6941e52b12e7
# ╠═00b86f70-8b57-49b3-bcec-bff646ace479
# ╠═06c56269-64c7-48da-9b1e-f1ff67e6696e
# ╠═a393d4e7-f72f-40fd-91f6-29f8009ea272
# ╠═d5aaa674-7af2-426a-9cc1-5c42f2c8248b
# ╠═872a2102-da21-4357-954a-fa25152395e8
# ╠═229350b2-7594-4c56-b0ca-85c5113a0a54
# ╠═0b0d9a9e-2251-4f6e-a761-2b247ea20471
# ╠═f0c9dc10-8448-4e2c-9aca-317ee396ea92
# ╠═1f02c11c-e8f5-45a9-b043-7e883fdbafd0
# ╠═9c90d6a9-3249-49ed-ae1d-398157922bdd
# ╠═99464f53-968b-466a-a357-b72ec132f135
# ╠═5816cfd3-4c3b-4320-b4ea-1185daae848e
# ╠═e0a79867-5207-4247-843f-62837596b2aa
# ╠═1adf1f53-f90b-4c9d-a670-998d432e19dd
# ╠═40eb7042-f7af-4520-b1bc-51998ba27e20
# ╠═9a9aab09-58db-49b5-9821-acf4a0260698
# ╠═2591487d-12a1-4324-85ce-3bcbbdaa734f
# ╠═be2e26c4-0a15-4108-a44a-2624973fe6ef
# ╠═a79da710-d267-4e88-8e38-1d07bea6ca53
# ╠═78b4b739-6e6b-4621-a125-8538c2000c61
# ╠═c41689fe-2c35-456c-8ef3-0716b43dca25
# ╠═01745045-8419-4d58-990d-3572ad1f47c2
# ╠═486a3d26-cafb-4dbd-b735-c0a3c89c02d9
# ╠═35102657-6a9a-4fbb-82b9-24f5289946c4
# ╠═71163085-0adb-49a3-8ce5-028d18c8f5d7
# ╠═ef238fdf-31eb-42c5-b91b-134196d5cf8f
# ╠═41c7c4ce-ff05-43ca-adc0-f1a4771c9844
# ╠═b3f76f05-57c7-4890-92a0-2f268f72ec5f
# ╟─51ab063f-a103-4b77-80f6-b7b5601df785
# ╟─04701ec7-6831-4617-b81c-1518264f7aaa
# ╠═19fc41ec-abcb-4c1d-9a6b-e16221ad7244
# ╠═14672d91-6b5d-4a0c-a2c2-e52150f766cb
# ╠═ff186ff3-c7de-4f93-a9b2-0331bcf2652b
# ╠═b3f3425b-3ccb-4de5-b32f-a698df810ce3
# ╠═1da12e54-263d-41af-a760-3da1af97f9f3
