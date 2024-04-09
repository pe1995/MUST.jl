mutable struct AtmosphericParameters{T<:AbstractFloat} 
    time             ::T
    teff             ::T
    logg             ::T
    composition      ::Dict{Symbol,T}
end

mutable struct Space <:AbstractSpace
    data             ::Dict{Symbol,Vector{<:Union{Float32, Float64, Int16, Int32, Int64}}}
    patch_dimensions ::Array{Int,2}
    parameter        ::AtmosphericParameters
end

mutable struct Box <:AbstractSpace
    x                ::Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}
    y                ::Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}
    z                ::Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}
    data             ::Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}
    parameter        ::AtmosphericParameters
end



#==== Functionality ====#

AtmosphericParameters() = AtmosphericParameters(-99.0, -99.0, -99.0, Dict{Symbol,Float64}())

"""
Atmospheric parameters of a dispatch snapshot.
Creates a link between the patch time, teff as well as physical parameters.
"""
function set!(atmp::AtmosphericParameters, teff_function, ini_namelist)
    teff    = teff_function(atmp.time)
    logg    = log10(ini_namelist.stellar_params["g_cgs"])

    epn = MUST.@in_dispatch(String(strip(ini_namelist.eos_params["table_loc"])[2:end-1]))
    eos_params_path = joinpath(epn, "tabparam.in")
    
    composition = try
        composition_from_eos(eos_params_path)
    catch
        Dict{Symbol,typeof(logg)}()
    end

    atmp.teff = teff
    atmp.logg = logg
    atmp.composition = composition

    nothing
end

read_teff(path) = begin
    path = joinpath(path)
    if !ispath(path) @warn "$(path) does not exist."
        return nothing
    end
    readdlm(path)
end 

teff_interpolated(path) = begin
    teff = read_teff(path)
    return isnothing(teff) ? 
        (args...)->0.0 :
        LinearInterpolation(teff[:,1], teff[:,2], extrapolation_bc=Flat())
end

function composition_from_eos(path)
    content = read_eos_params(path)
    k       = [keys(content)...]
    cel     = findfirst(occursin.("cel",k))
    cel     = k[cel]
    el      = split(content[cel][2:end-1], " ")
    el      = [String(s) for s in el if length(s)>0]
    
    cel     = findfirst(occursin.("abund",k))
    cel     = k[cel]
    vl      = String.(split(content[cel][2:end-1], " "))
    vl      = [s for s in vl if length(s)>0]
    vl      = parse.(Float64, vl)
    Dict(Symbol(c)=>v for (c,v) in zip(el,vl))
end






function _var_from_patch(var, fname, shp, off, li, ui, idxd; density=:d)
    if var==:ux
        varidx = pyconvert(Any, idxd[String(:px)]) + 1
        didx = pyconvert(Any, idxd[String(density)]) + 1

        d = mmap(fname, Array{Float32, length(shp)}, shp, off[didx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
        mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] ./ d
    elseif var==:uy
        varidx = pyconvert(Any, idxd[String(:py)]) + 1
        didx = pyconvert(Any, idxd[String(density)]) + 1

        d = mmap(fname, Array{Float32, length(shp)}, shp, off[didx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
        mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] ./ d
    elseif var==:uz
        varidx = pyconvert(Any, idxd[String(:pz)]) + 1
        didx = pyconvert(Any, idxd[String(density)]) + 1

        d = mmap(fname, Array{Float32, length(shp)}, shp, off[didx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
        mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] ./ d
    elseif var==:ee
        varidx = pyconvert(Any, idxd[String(:e)]) + 1
        didx = pyconvert(Any, idxd[String(density)]) + 1

        d = mmap(fname, Array{Float32, length(shp)}, shp, off[didx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
        mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]] ./ d
    else
        varidx = pyconvert(Any, idxd[String(var)]) + 1
        mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
    end
end

"""
Create a Space object from a Dispatch simulation snapshot.

Parameters
----------
snapshot: PyObject
    A dispatch snapshot.
quantites: Symbols
    Variables to extract from the snapshot

Returns
-------
MUST.Space type object

"""
function Space(snapshot::Py, quantities::Symbol...; density=:d, use_numpy=false) 
    qs = Dict{Symbol,Vector{Float32}}(q=>Float32[] for q in quantities)
    qs[:x] = Float32[]; qs[:y] = Float32[]; qs[:z] = Float32[]
    qs[:i_patch] = Int[]
    patch_dimensions = zeros(Int,(length(snapshot.patches),3))
    
    for (i,patch) in enumerate(snapshot.patches)
        # grid of this patch relative to the global grid
        #patch_size = (length(patch.xi),length(patch.yi),length(patch.zi))
        x = pyconvert.(Float32, patch.xi)
        y = pyconvert.(Float32, patch.yi) 
        z = pyconvert.(Float32, patch.zi)
        patch_size = (length(x), length(y), length(z))
        patch_dimensions[i,1] = length(patch.xi)
        patch_dimensions[i,2] = length(patch.yi)
        patch_dimensions[i,3] = length(patch.zi)

        # New: Dont rely on numpy for mmaps
        fname = pyconvert(Any, patch.filename)
		shp = Tuple(pyconvert.(Any, patch.ncell))
		off = pyconvert(Vector{Int}, patch.offset)
		li = pyconvert(Vector{Int}, patch.li)
		ui = pyconvert(Vector{Int}, patch.ui)
        idxd = patch.idx.__dict__
        
        # extract the quantites for this patch
        for (iq, q) in enumerate(quantities)            
            q_matrix = use_numpy ? 
                pyconvert(Array{Float32, 3}, numpy.array(patch.var(String(q)))) :
                _var_from_patch(q, fname, shp, off, li, ui, idxd, density=density) 
            
            q_array = zeros(Float32, prod(patch_size))
            coords = zeros(Float32, (prod(patch_size),3))
            j = 1
            @inbounds for iz in 1:patch_size[3]
                @inbounds for iy in 1:patch_size[2]
                    @inbounds for ix in 1:patch_size[1]
                        q_array[j]  = q_matrix[ix, iy, iz]
                        coords[j,1] = x[ix]
                        coords[j,2] = y[iy]
                        coords[j,3] = z[iz]
                        j += 1
                    end
                end
            end
            append!(qs[q], q_array)
            if iq==1
                append!(qs[:x],coords[:,1])
                append!(qs[:y],coords[:,2])
                append!(qs[:z],coords[:,3])
                append!(qs[:i_patch],[i for _ in 1:length(q_array)])
            end
        end
    end

    time = pyconvert(Any, snapshot.nml_list["snapshot_nml"]["time"])
    Space(qs, patch_dimensions, AtmosphericParameters(time, 
                                                        Base.convert(typeof(time), -99.0), 
                                                        Base.convert(typeof(time), -99.0), 
                                                        Dict{Symbol, typeof(time)}())
        )
end

"""
Create a Space object from a Legacy Stagger Snapshot.

Parameters
----------
snapshot: PyObject
    A dispatch snapshot.
quantites: Symbols
    Variables to extract from the snapshot

Returns
-------
MUST.Space type object

"""
function Space(snapshot::MUST.StaggerLegacySnap, quantities::Symbol...)
    snap = snapshot.snap
    qs = Dict{Symbol,Vector{<:Union{Float32, Float64, Int16, Int32, Int64}}}(q=>Float32[] for q in quantities)
    qs[:x] = Float32[]; qs[:y] = Float32[]; qs[:z] = Float32[]
    qs[:i_patch] = Int[]
    patch_dimensions = zeros(Int, (1,3))
    patch_dimensions[1,1] = snap.nx*snap.ny*snap.nz
    patch_dimensions[1,2] = snap.nx*snap.ny*snap.nz
    patch_dimensions[1,3] = snap.nx*snap.ny*snap.nz

    x,y,z = snap.read_mem("xx"), snap.read_mem("yy"), snap.read_mem("zz")

    qm = zeros(Float32, (snap.nx, snap.ny, snap.nz))
    for (i,q) in enumerate(quantities)
        qm = snap.read_mem(String(q))
        qs[q] = zeros(Float32, snap.nx*snap.ny*snap.nz)
        if i==1
            coords = zeros(Float32, (snap.nx*snap.ny*snap.nz,3))
        end
        c=1
        @inbounds for k in 1:snap.nz
        @inbounds for j in 1:snap.ny
        @inbounds for n in 1:snap.nx
            qs[q][c] = qm[n,j,k]
            if i==1
                coords[c,1] = x[n]
                coords[c,2] = y[j]
                coords[c,3] = z[k]
            end
            c+=1
        end
        end
        end
        if i==1
            append!(qs[:x],coords[:,1])
            append!(qs[:y],coords[:,2])
            append!(qs[:z],coords[:,3])
            append!(qs[:i_patch], Int[1 for _ in 1:size(coords,1)])
        end
    end
    Space(qs, patch_dimensions, AtmosphericParameters(-99.0, -99.0, -99.0, Dict{Symbol, Float64}()))
end

"""
Interpolate a Space object on an regular grid and save as Arrays in the Box type.
Use this if you want 3D data cubes.

Parameters
----------
x: Float Vector
    A dispatch snapshot.
y: Float Vector
    A dispatch snapshot.
z: Float Vector
    A dispatch snapshot.
s: Space
    The space object which is the Basis of the Box

Returns
-------
MUST.Box type object

"""
function Box(s::Space, x::Vector{T}, y::Vector{T}, z::Vector{T}) where {T<:AbstractFloat}
    # filter the space to match the 3d cube
    region_mask = (;row...) ->  (minimum(x) .<= row[:x] .<= maximum(x)) .&
                                (minimum(y) .<= row[:y] .<= maximum(y)) .&
                                (minimum(z) .<= row[:z] .<= maximum(z))
    s_masked = filter(region_mask, s)

    # 3D Mesh of the new coordinate Box
    #x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z, indexing="ij")
    #x_grid = pyconvert(Array{eltype(x)}, x_grid)
    #y_grid = pyconvert(Array{eltype(x)}, y_grid)
    #z_grid = pyconvert(Array{eltype(x)}, z_grid)
    x_grid, y_grid, z_grid = meshgrid(x, y, z)


    # result dict
    results::Dict{Symbol,Array{T,3}} = Dict(q=>similar(x_grid) for q in keys(s.data) if !(q in [:x,:y,:z,:i_patch]))

    mask_data = []
    grid_data = []
    if length(x) > 1
        append!(mask_data, [s_masked.data[:x]])
        append!(grid_data, [x_grid])
    end
    if length(y) > 1
        append!(mask_data, [s_masked.data[:y]])
        append!(grid_data, [y_grid])
    end
    if length(z) > 1
        append!(mask_data, [s_masked.data[:z]])
        append!(grid_data, [z_grid])
    end


    # Interpolate the Space onto the given grid
    for quantity in Symbol.(keys(s.data))
        quantity in [:x,:y,:z,:i_patch] ? continue : nothing
        results[quantity][1:length(x),1:length(y),1:length(z)] = pyconvert(Array{eltype(x_grid)}, scipy_interpolate.griddata(Tuple(mask_data), 
                                                                                            s_masked.data[quantity], 
                                                                                            Tuple(grid_data)))
    end

    Box(x_grid[2:end-1,2:end-1,2:end-1],
        y_grid[2:end-1,2:end-1,2:end-1],
        z_grid[2:end-1,2:end-1,2:end-1],
        Dict(q=>results[q][2:end-1,2:end-1,2:end-1] for q in keys(results)), s.parameter)
end

uniform_grid(s::Space, N::Int,axis=:x) = begin
    data = s.data[axis]
    maxd = maximum(data)
    mind = minimum(data)
    extd = maxd-mind

    mind += extd*0.001
    maxd -= extd*0.001
    extd = maxd-mind
    step = extd/N
    Vector(mind:max(step,eps(1.0)):maxd)
end

Box(s::Space, N::Int, args...; kwargs...) = begin Box(s,
                            uniform_grid(s, N, :x), 
                            uniform_grid(s, N, :y), 
                            uniform_grid(s, N, :z), args...; kwargs...)
end

"""
Interpolate a Space object on an regular grid and save as new Space object.

Parameters
----------
s: Space
    The space object which is the basis of the Box
x: Float Vector
    A dispatch snapshot.
y: Float Vector
    A dispatch snapshot.
z: Float Vector
    A dispatch snapshot.

Returns
-------
MUST.Space type object
"""
function spacebox(s::MUST.Space,x::Vector{T}, y::Vector{T}, z::Vector{T}, args...; kwargs...) where {T<:AbstractFloat}
    box = Box(s, x, y, z, args...; kwargs...)
    Space(box)
end

spacebox(s::Space, N::Int, args...; kwargs...) = begin spacebox(s,
                                                    uniform_grid(s, N, :x), 
                                                    uniform_grid(s, N, :y), 
                                                    uniform_grid(s, N, :z), args...; kwargs...)
end

"""
Convert Space to Box, assuming that the space has already uniform grid (No interpolation)
"""
function Box(s::Space)   
    # create x,y,z from unique s.data -> meshgrid x,y,z -> pick data from dataframe where coordinates match
    x = sort(unique(s.data[:x]))
    y = sort(unique(s.data[:y]))
    z = sort(unique(s.data[:z]))

    #x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z, indexing="ij")
    #x_grid = pyconvert(Array{eltype(x)}, x_grid)
    #y_grid = pyconvert(Array{eltype(x)}, y_grid)
    #z_grid = pyconvert(Array{eltype(x)}, z_grid)
    x_grid, y_grid, z_grid = meshgrid(x, y, z)

    results::Dict{Symbol, Array{<:Union{Float32, Float64, Int16, Int32, Int64}, 3}} = Dict(q=>zeros(typeof(s.data[q][1]),size(x_grid)) 
                                                                                        for q in keys(s.data) if !(q in [:x,:y,:z]))
    for q in keys(results)
        q == :i_patch ? continue : nothing
        results[q] .= Base.convert(typeof(results[q][1,1,1]), NaN)
    end

    p = zeros(3) 
    skip_q = [:x, :y, :z]
    @inbounds for irow in eachindex(s.data[:x])
        p[1] = s.data[:x][irow]
        p[2] = s.data[:y][irow] 
        p[3] = s.data[:z][irow] 
        ix,iy,iz = _find_in_meshgrid(p, x, y, z)
        for quantity in keys(s.data)
            quantity in skip_q ? continue : nothing
            
            # go though the data points and save the quantity
            results[quantity][ix,iy,iz] = s.data[quantity][irow]
        end
    end

    Box(x_grid, y_grid, z_grid, results, s.parameter)
end 

function Space(box::MUST.Box)
    q_dict = Dict(q => reshape(box.data[q], :) for q in keys(box.data))
    q_dict[:x] = reshape(box.x, :); q_dict[:y] = reshape(box.x, :); q_dict[:z] = reshape(box.x, :)
    q_dict[:i_patch] = ones(length(q_dict[:x]))
    
    p_dim = ones(1,3)
    p_dim[1,1] = length(x); p_dim[1,2] = length(y); p_dim[1,3] = length(z)
    Space(q_dict, p_dim, box.parameter)
end

@inline _find_in_meshgrid(p, x, y, z) = begin
    ix = findfirst(i->isapprox(p[1], i), x)
    iy = findfirst(i->isapprox(p[2], i), y)
    iz = findfirst(i->isapprox(p[3], i), z)
    (ix,iy,iz)
end   



#= Boxes from snapshots wihout Spaces =#

_add_if_not_exists!(arr, data) = begin
	@inbounds for entry in data
		if !(entry in arr)
			append!(arr, entry)
		end
	end
end

function _build_cube(patch_data)
	# (x, y, z), (li, ui), (patches...)
	patch_range = zeros(Int, 3, 2, length(patch_data))

	global_x = []
	global_y = []
	global_z = []
	@inbounds for i in eachindex(patch_data)
		_add_if_not_exists!(global_x, patch_data[i][:x])
		_add_if_not_exists!(global_y, patch_data[i][:y])
		_add_if_not_exists!(global_z, patch_data[i][:z])
	end

	global_x .= sort(global_x)
	global_y .= sort(global_y)
	global_z .= sort(global_z)

	p = zeros(3)
	@inbounds for i in eachindex(patch_data)
		p[1] = first(patch_data[i][:x])
		p[2] = first(patch_data[i][:y])
		p[3] = first(patch_data[i][:z])
		patch_range[:, 1, i] .= _find_in_meshgrid(p, global_x, global_y, global_z)

		p[1] = last(patch_data[i][:x])
		p[2] = last(patch_data[i][:y])
		p[3] = last(patch_data[i][:z])
		patch_range[:, 2, i] .= _find_in_meshgrid(p, global_x, global_y, global_z)
	end
	
	global_x, global_y, global_z, patch_range
end

"""
	Box(snap::Py)

Convert a snapshot to a `Box`. Assumes that the patch arangement is cubic, i.e. there is no local mesh refinement. If this is not true, convert to `Space` first, and the interpolate to `Box`.
"""
function Box(snap::Py, quantities::Symbol...; density=:d, use_mmap=false)
	# first we loop through and get the sizes of all patches
	patch_data = []
	@inbounds for (i, patch) in enumerate(snap.patches)
		x = pyconvert(Array{Float32}, patch.xi)
        y = pyconvert(Array{Float32}, patch.yi) 
        z = pyconvert(Array{Float32}, patch.zi)

        fname = pyconvert(Any, patch.filename)
		shp = Tuple(pyconvert.(Any, patch.ncell))
		off = pyconvert(Vector{Int}, patch.offset)
		li = pyconvert(Vector{Int}, patch.li)
		ui = pyconvert(Vector{Int}, patch.ui)
        idxd = patch.idx.__dict__
        aux_vars = patch.aux.vars
       
		meta = Dict(
			:x=>x,
			:y=>y,
			:z=>z,
			:fname=>fname,
			:shp=>shp,
            :off=>off,
            :li=>li,
            :ui=>ui,
            :idxd=>idxd,
            :aux_vars=>aux_vars
		)
		append!(patch_data, [meta])
	end

	x, y, z, patch_range = _build_cube(patch_data)
    q_and_aux = [quantities..., Symbol.(pyconvert(Vector{String}, first(patch_data)[:aux_vars].keys()))...]

	# now we create the data arrays
	data = if !use_mmap
		Dict{Symbol,Array{Float32,3}}(
			q=>Array{Float32, 3}(undef, length(x), length(y), length(z)) 
			for q in q_and_aux
		)
	else
		pname = tempname(pwd())
		io = open(pname, "w+")
		d = mmap(io, Array{Float32, 4}, (length(x), length(y), length(z), length(quantities)))
		Dict{Symbol,Array{Float32,3}}(
			q=>@view d[:, :, :, i] 
			for (i, q) in enumerate(quantq_and_auxities)
		)
	end    

	# and we loop through the patches, read the mmaps, and fill them in the data arrays
	r = zeros(Int, 3, 2)
	@inbounds for (i, pd) in enumerate(patch_data)
		r .= patch_range[:, :, i]
		@inbounds for (j, q) in enumerate(quantities)
			data[q][r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] .= _var_from_patch(
				q, 
				pd[:fname], 
				pd[:shp], 
				pd[:off], 
				pd[:li], 
				pd[:ui], 
				pd[:idxd], 
				density=density
			)
		end

        # add aux variables
        av = pyconvert(Vector{String}, pd[:aux_vars].keys())
        li = pd[:li]
        ui = pd[:ui]
        for (j, aux) in enumerate(av)
            data[Symbol(aux)][r[1,1]:r[1,2], r[2,1]:r[2,2], r[3,1]:r[3,2]] .= pyconvert(
                Array{Float32}, pd[:aux_vars][aux]["v"]
            )[li[1]:ui[1],li[2]:ui[2],li[3]:ui[3]]
        end
	end

	time = pyconvert(Any, snap.nml_list["snapshot_nml"]["time"])
	xx, yy, zz = meshgrid(x, y, z)
	Box(
		xx, yy, zz, data,
		AtmosphericParameters(
			time, 
			Base.convert(typeof(time), -99.0), 
			Base.convert(typeof(time), -99.0), 
			Dict{Symbol, typeof(time)}()
		)
	)
end







"""
Filter the Space object.

Parameters
----------
f: Function that accepts any Key of the space dictionary 
    (NOTE: Takes the DataFrame Row for now)
s: Space
    The space object to filter

Returns
-------
MUST.Space type object

"""
function Base.filter(f, s::MUST.Space)
    mask = f(;s.data...)
	new_frame = Dict(q=>s.data[q][mask] for q in keys(s.data))
	Space(new_frame, s.patch_dimensions, s.parameter)
end

function Base.filter!(f, s::MUST.Space)
    mask = f(;s.data...)
    for q in keys(s.data)
        s.data[q] = s.data[q][mask]
	end
end

Base.getindex(s::T,key) where {T<:MUST.Space} = s.data[key]
Base.getindex(s::T,key) where {T<:MUST.Box}   = key in keys(s.data) ? s.data[key] : getfield(s,key)

Base.setindex!(s::T, val, key) where {T<:MUST.Space} = begin
    s.data[key] = val
end

Base.setindex!(s::T, val, key) where {T<:MUST.Box} = begin
    if (key in fieldnames(typeof(s)))
        setfield!(s, key, val)
    else
        s.data[key] = val
    end
end

Base.length(s::MUST.Space) = length(s.data[:x])
Base.length(s::MUST.Box)   = length(s.x)
Base.keys(s::T) where {T<:Union{MUST.Box,MUST.Space}} = keys(s.data)

"""
Convert to given units using a AtmosUnits type.

Parameters
----------
s: Space
    The space object
u: AstroUnits
    Object containting the conversion factor
params: Symbols
    Map between Space entry and unit that should be applied
    
Returns
-------
nothing

"""
function convert!(s::AbstractSpace, u::AtmosUnits; params...)
    for (s_para, u_para) in params
        if s_para == :time
            s.parameter.time = s.parameter.time * getfield(u, u_para)
        else
            if !(s_para in keys(s.data)) & !(s_para in [:x, :y, :z])
                continue
            else
                s[s_para] .= s[s_para] .* getfield(u, u_para)
            end
        end
    end
    nothing
end

function convert(s::AbstractSpace, u::AtmosUnits; params...)
    s_new = deepcopy(s)
    convert!(s_new, u; params...)
    s_new
end

"""
Add the given field from the EOS.

Parameters
----------
s: Space
    The space object
eos: Dispatch EOS
    EOS to convert the quantites
    
Returns
-------
nothing

"""
function add_from_EOS!(s::SB, eos::AbstractEOS, quantity::Symbol;
                        EOS_paras=(:d,:ee), 
                        convert_to=Base.identity, switch_name_to=nothing) where {SB<:Union{Space, Box}}
    qframe = MUST.add_from_EOS(s, eos, quantity; 
                                EOS_paras =EOS_paras, 
                                convert_to=convert_to)
    isnothing(switch_name_to) ? s.data[quantity] = qframe : s.data[switch_name_to] = qframe
    nothing
end

function add_from_EOS(s::MUST.Space, eos::AbstractEOS, quantity::Symbol;
                        EOS_paras=(:d,:ee), 
                        convert_to=Base.identity)
    convert_to.(lookup(eos, String(quantity), s[EOS_paras[1]], s[EOS_paras[2]]))
end

function add_from_EOS(b::MUST.Box, eos::AbstractEOS, quantity::Symbol;
                EOS_paras=(:d,:ee), 
                convert_to=Base.identity)
    r = similar(b[EOS_paras[1]]) 
    r .= convert_to.(
        lookup(
            eos, 
            String(quantity), 
            b[EOS_paras[1]],
            b[EOS_paras[2]]
            )
        )

    #=@inbounds for j in axes(b.z, 2)
        @inbounds for i in axes(b.z ,1)
            r[i, j, :] = convert_to.(lookup(eos, String(quantity), 
                            view(b[EOS_paras[1]], i, j, :), 
                            view(b[EOS_paras[2]], i, j, :)))
        end
    end=#

    r
end


"""
Add the given field from the given array.

Parameters
----------
s: Space
    The space object
field_name: Symbol
    field that should be added
arr:
    array that should be added. Has to match in size
    
Returns
-------
nothing

"""
function add!(s::S, arr::Vector{T}, name::Symbol) where {T<:AbstractFloat, S<:MUST.Space}
    @assert length(arr) == length(s)
    s.data[name] = Base.convert(Vector{typeof(s.data[:x][1])},arr)
    nothing
end 

function add!(s::S, arr::Array{T,3}, name::Symbol) where {T<:AbstractFloat, S<:MUST.Box}
    @assert length(arr) == length(s)
    s.data[name] = Base.convert(Array{typeof(s.x[1,1,1]),3},arr)
    nothing
end 

function add!(s::MUST.Space, b::MUST.Box, name::Symbol)
    p = zeros(3)
    x,y,z = b.x[:,1,1],b.y[1,:,1],b.z[1,1,:]

    res = rand(typeof(b.data[name][1,1,1]), length(s))
    for i in 1:length(s)
        p[1] = s.data[:x][i]
        p[2] = s.data[:y][i]
        p[3] = s.data[:z][i]
        ix,iy,iz = _find_in_meshgrid(p, x, y, z)
        res[i] = b.data[name][ix,iy,iz]
    end
    add!(s, res, name)
    nothing
end 

"""
    add!(s; kwargs...)

Add given arrays to given space.
"""
add!(s; kwargs...) = begin
    for (k, v) in kwargs
        add!(s, v, k)
    end

    nothing
end

"""
Compute the horizonal average of a Box object.

Parameters
----------
s: Box
    The space object
q: Symbol
    Quantity to compute average
    
Returns
-------
nothing

"""
plane_average(b::MUST.Box, q::Symbol) = plane_statistic(Statistics.mean, b, q)

"""
Compute the horizonal statistic of a Box object.

Parameters
----------
s: Box
    The space object
q: Symbol
    Quantity to compute average
    
Returns
-------
nothing

"""
function plane_statistic(stat::F, b::MUST.Box, q::Symbol) where {F<:Function}
    av = zeros(typeof(b.z[1,1,1]), size(b.z,3))
    for i in 1:size(b.z,3)
        av[i] = stat(b[q][:,:,i])
    end
    av
end

function _check_path(folder, filename)
    if !isnothing(filename)
        file    = isnothing(filename) ? "space_"*randstring(6)*".hdf5" : filename*".hdf5"
        path    = isnothing(folder) ? @in_dispatch(file) : joinpath(folder, file)
    else
        is_free = false
        while !is_free
            file    = "space_"*randstring(6)*".hdf5" 
            path    = isnothing(folder) ? @in_dispatch(file) : joinpath(folder, file)
            is_free = !isfile(path)
            !is_free ? println("$(path) is not free.") : nothing
        end
    end
    path
end

"""
Save the content of a Space object as HDF5 file.

Parameters
----------
s: Space object
    The space object to be saved
    
Returns
-------
nothing

"""
function save(s::MUST.Space; folder=nothing, name=nothing)
    path = _check_path(folder, name)
    fid  = HDF5.h5open(path, "w")
    for q in keys(s.data)
        fid[String(q)] = s.data[q]
    end
    fid["patch_dimensions"] = s.patch_dimensions

    # Save the parameters
    save(s.parameter, fid)

    close(fid)
    path
end

function save(s::MUST.Box; folder=nothing, name=nothing)
    path = _check_path(folder, name)
    fid  = HDF5.h5open(path, "w")
    for q in keys(s.data)
        fid[String(q)] = s.data[q]
    end
    fid["x"] = s.x
    fid["y"] = s.y
    fid["z"] = s.z

    # Save the parameters
    save(s.parameter, fid)

    close(fid)
    path
end

"""
    save(s::Box, number::Int)

Identify on which scale the box is and save it with the correct name.
"""
save(s::Box, number::Int; folder=nothing) = begin
    name = is_geo_scale(s) ? "box_sn$(number)" : "box_tau_sn$(number)"

    save(s, name=name, folder=folder)
end

is_geo_scale(s::Box) = begin
    all_equal = trues(size(s.z, 3))
    zref = s.z[1, 1, :]
    for iy in axes(s.y, 2)
        for ix in axes(s.x, 1)
            all_equal .= all_equal .& (zref .≈ s.z[ix, iy, :])
        end
    end

    all(all_equal)
end

function save(p::AtmosphericParameters, fid)
    eles = [keys(p.composition)...]
    vals = eltype(values(p.composition))[p.composition[e] for e in eles]
    fid["time"] = p.time
    fid["teff"] = p.teff
    fid["logg"] = p.logg
    fid["composition_e"] = String[String(e) for e in eles]
    fid["composition_v"] = vals
end

"""
Load the content of a Space object from HDF5 file. The data will be loaded as a memmap.

Parameters
----------
filename: String
    Name of the saved memmap
folder: 
    Folder where it has been saved.
    
Returns
-------
Space object

"""
function Space(name::String; folder::F=nothing) where {F<:Union{String,Nothing}}
    if !isfile(name)
        path = isnothing(folder) ? @in_dispatch(name*".hdf5") : joinpath(folder, name*".hdf5")
        @assert isfile(path) 
    else
        path = name
    end

    fid = HDF5.h5open(path, "r")
    res = Dict{Symbol,Vector{<:Union{Float32, Float64, Int16, Int32, Int64}}}()

    aux_fieldnames = _get_para_fieldnames(AtmosphericParameters)
    append!(aux_fieldnames, ["patch_dimensions"])

    for q in keys(fid)
        q in aux_fieldnames ? continue : nothing
        res[Symbol(q)] = HDF5.readmmap(fid[q])
    end
    patch_dimensions = HDF5.readmmap(fid["patch_dimensions"])
    params           = _read_params(AtmosphericParameters, fid)
    
    close(fid)

    Space(res, patch_dimensions, params)
end

function Box(name::String; folder::F=nothing) where {F<:Union{String,Nothing}}
    if !isfile(name)
        path = isnothing(folder) ? @in_dispatch(name*".hdf5") : joinpath(folder, name*".hdf5")
        if !isfile(path) 
            error("$(path) does not exist.")
        end

        path
    else
        path = name
    end

    fid = HDF5.h5open(path, "r")
    res = Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}()

    aux_fieldnames = _get_para_fieldnames(AtmosphericParameters)
    append!(aux_fieldnames, ["x","y","z"])

    for q in keys(fid)
        q in aux_fieldnames ? continue : nothing
        res[Symbol(q)] = HDF5.readmmap(fid[q])
    end
    x = HDF5.readmmap(fid["x"])
    y = HDF5.readmmap(fid["y"])
    z = HDF5.readmmap(fid["z"])

    params = _read_params(AtmosphericParameters, fid)

    close(fid)

    Box(x, y, z, res, params)
end

_get_para_fieldnames(t::Type{AtmosphericParameters}) = ["time", "teff", "logg", "composition_e", "composition_v", "parameter_type"]
_read_params(T::Type{AtmosphericParameters}, res) = begin 
    T(HDF5.read(res["time"]), HDF5.read(res["teff"]), HDF5.read(res["logg"]), 
                                                Dict(Symbol(e) => v for (e,v) in zip(HDF5.read(res["composition_e"]), HDF5.read(res["composition_v"]))) )

end





"""
Reduce a Box object to one single plane by applying the function f to each column in x,y

Parameters
----------
b: Box
    Box that should be reduced
f: Function
    taking 1-D array of a column and return a single value, which is the desired height index. 
    If index is false, then not the index of the desired height coordinate is returned, but the z value itself.
    The oter quantites are then interpolated to this z value.
    The function has to accept all columns as kwargs!
index: Bool
    
Returns
-------
Box object

Example
-------
f(; tau, columns...) = argmin(tau)
b_new = reduce_by_column(b, f)

"""
function reduce_by_column(f::F, b::MUST.Box; index=false) where {F<:Function}
    shape    = (size(b.x,1),size(b.x,2),1)
    data_new = Dict(q=>rand(typeof(b.data[q][1,1,1]), shape) for q in keys(b.data))
    x_new = rand(typeof(b.x[1,1,1]), shape)
    y_new = rand(typeof(b.x[1,1,1]), shape)
    z_new = rand(typeof(b.x[1,1,1]), shape)

    for j in 1:shape[2]
        for i in 1:shape[1]
            zv = @view b.z[i,j,:]
            xv = @view b.x[i,j,:] 
            yv = @view b.y[i,j,:] 
            qv = Dict(q => @view(b.data[q][i,j,:]) for q in keys(b.data))

            mm = (.!(isnan.(zv))) .& (.!(isnan.(xv))) .& (.!(isnan.(yv))) 
            for q in keys(qv)
                mm = mm .& .!(isnan.(qv[q]))
            end

            qv = Dict(q => @view(b.data[q][i,j,mm]) for q in keys(b.data))
            zv = @view b.z[i,j,mm]
            xv = @view b.x[i,j,mm] 
            yv = @view b.y[i,j,mm]

            if index
                idx = length(zv) == 0 ? 1 : f(; qv..., x=xv, y=yv, z=zv) 
                mm  = length(zv) == 0 ? [true for _ in 1:length(b.z[i,j,:])] : mm

                for q in keys(b.data)
                    data_new[q][i,j,1] = b.data[q][i,j,mm][idx]
                end
                x_new[i,j,1] = b.x[i,j,mm][idx]
                y_new[i,j,1] = b.y[i,j,mm][idx]
                z_new[i,j,1] = b.z[i,j,mm][idx]
            else
                if length(zv) > 1 
                    z_int = f(; qv..., x=xv, y=yv, z=zv) 

                    szv = sortperm(zv)
                    for q in keys(b.data)
                        data_new[q][i,j,1] = LinearInterpolation(zv[szv], qv[q][szv], extrapolation_bc=Flat())(z_int)
                    end
                    x_new[i,j,1] = LinearInterpolation(zv[szv], xv[szv], extrapolation_bc=Flat())(z_int)
                    y_new[i,j,1] = LinearInterpolation(zv[szv], yv[szv], extrapolation_bc=Flat())(z_int)
                    z_new[i,j,1] = z_int
                else
                    for q in keys(b.data)
                        data_new[q][i,j,1] = b.data[q][i,j,1]
                    end
                    x_new[i,j,1] = b.x[i,j,1]
                    y_new[i,j,1] = b.y[i,j,1]
                    z_new[i,j,1] = b.z[i,j,1]
                end
            end
        end
    end

    Box(x_new, y_new, z_new, data_new, b.parameter)
end

"""
Linear Interpolate the columns of the Box as a function of height. For every column a Function is returned.
"""
function interpolate_by_column(b::MUST.Box)
    shape    = (size(b.x,1),size(b.x,2),1)
    data_new = Dict{Symbol, Matrix{Interpolations.Extrapolation}}(q => Matrix{Any}(undef,shape[1:2]...) for q in keys(b.data))
    
    data_new[:x] = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)
    data_new[:y] = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)

    @inbounds for j in 1:shape[2]
        @inbounds for i in 1:shape[1]
            zv = @view b.z[i,j,:]
            xv = @view b.x[i,j,:] 
            yv = @view b.y[i,j,:] 
            qv = Dict(q => @view(b.data[q][i,j,:]) for q in keys(b.data))

            mm = (.!(isnan.(zv))) .& (.!(isnan.(xv))) .& (.!(isnan.(yv))) 
            for q in keys(qv)
                mm = mm .& .!(isnan.(qv[q]))
            end

            qv = Dict(q => @view(b.data[q][i,j,mm]) for q in keys(b.data))
            zv = @view b.z[i,j,mm]
            xv = @view b.x[i,j,mm] 
            yv = @view b.y[i,j,mm]

            szv = sortperm(zv)
            for q in keys(b.data)
                data_new[q][i,j] = LinearInterpolation(Interpolations.deduplicate_knots!(zv[szv]), qv[q][szv], extrapolation_bc=Flat())
            end
            data_new[:x][i,j] = LinearInterpolation(Interpolations.deduplicate_knots!(zv[szv]), xv[szv], extrapolation_bc=Flat())
            data_new[:y][i,j] = LinearInterpolation(Interpolations.deduplicate_knots!(zv[szv]), yv[szv], extrapolation_bc=Flat())
        end
    end

    data_new
end

"""
Linear Interpolate the columns of the Box as a function of height. For every column a Function is returned.
"""
function interpolate_by_column(b::MUST.Box, val; logspace=true)
    shape    = (size(b.x,1),size(b.x,2),1)
    data_new = Dict{Symbol, Matrix{Interpolations.Extrapolation}}(q => Matrix{Any}(undef,shape[1:2]...) for q in keys(b.data))

    data_new[:x] = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)
    data_new[:y] = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)
    data_new[:z] = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)

    @inbounds for j in 1:shape[2]
        @inbounds for i in 1:shape[1]
            zv = @view b.z[i,j,:]
            xv = @view b.x[i,j,:] 
            yv = @view b.y[i,j,:] 
            qv = Dict(q => @view(b.data[q][i,j,:]) for q in keys(b.data))

            mm = (.!(isnan.(zv))) .& (.!(isnan.(xv))) .& (.!(isnan.(yv))) 
            for q in keys(qv)
                mm = mm .& .!(isnan.(qv[q]))
            end

            qv = Dict(q => @view(b.data[q][i,j,mm]) for q in keys(b.data))
            zv = @view b.z[i,j,mm]
            xv = @view b.x[i,j,mm] 
            yv = @view b.y[i,j,mm]

            szv = sortperm(qv[val])
            for q in keys(b.data)
                if q == val
                    continue
                end
                data_new[q][i,j] = logspace ? 
                    LinearInterpolation(
                        Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), 
                        qv[q][szv], 
                        extrapolation_bc=Flat()
                    ) : 
                    LinearInterpolation(
                        Interpolations.deduplicate_knots!(qv[val][szv]), 
                        qv[q][szv], 
                        extrapolation_bc=Flat()
                    ) 
            end

            data_new[:z][i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), zv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv], zv[szv]), extrapolation_bc=Flat()) 

            data_new[:x][i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), xv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv], xv[szv]), extrapolation_bc=Flat()) 
            data_new[:y][i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), yv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv], yv[szv]), extrapolation_bc=Flat()) 
        end
    end

    data_new
end

"""
Linear Interpolate the columns of the Box as a function of height. For every column a Function is returned.
"""
function interpolate_by_column(b::MUST.Box, val, qin; logspace=true)
    shape    = (size(b.x,1),size(b.x,2),1)
    data_new = Matrix{Any}(undef,shape[1:2]...)
    data_z   = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)
    data_x   = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)
    data_y   = Matrix{Interpolations.Extrapolation}(undef,shape[1:2]...)


    names = [val, qin]

    @inbounds for j in 1:shape[2]
        @inbounds for i in 1:shape[1]
            zv = @view b.z[i,j,:]
            xv = @view b.x[i,j,:] 
            yv = @view b.y[i,j,:] 
            qv = Dict(q => @view(b.data[q][i,j,:]) for q in names)

            mm = (.!(isnan.(zv))) .& (.!(isnan.(xv))) .& (.!(isnan.(yv))) 
            for q in keys(qv)
                mm = mm .& .!(isnan.(qv[q]))
            end

            zv = @view b.z[i,j,mm]
            xv = @view b.x[i,j,mm] 
            yv = @view b.y[i,j,mm]
            qv  = Dict(q => @view(b.data[q][i,j,mm]) for q in names)
            szv = sortperm(qv[val])
            for q in names
                if q == val
                    continue
                end
                data_new[i,j] = logspace ? 
                    LinearInterpolation(
                        Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), 
                        qv[q][szv], 
                        extrapolation_bc=Flat()
                    ) : 
                    LinearInterpolation(
                        Interpolations.deduplicate_knots!(qv[val][szv]), 
                        qv[q][szv], 
                        extrapolation_bc=Flat()
                    ) 
            end

            data_z[i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), zv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv]), zv[szv], extrapolation_bc=Flat())
                            
            data_x[i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), xv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv]), xv[szv], extrapolation_bc=Flat()) 
            data_y[i,j] = logspace ? 
                            LinearInterpolation(Interpolations.deduplicate_knots!(log.(10, qv[val][szv])), yv[szv], extrapolation_bc=Flat()) : 
                            LinearInterpolation(Interpolations.deduplicate_knots!(qv[val][szv]), yv[szv], extrapolation_bc=Flat()) 
        end
    end

    Dict(qin=>data_new, :x=>data_x, :y=>data_y, :z=>data_z)
end

"""
Apply the function f to every column of the Box and save the result in res.
"""
function apply_by_column!(f::Function, res::T, b::MUST.Box; check_nan=false) where {T<:AbstractArray}
    
    #zv = zeros(size(b.z,3))
    #xv = zeros(size(b.z,3))
    #yv = zeros(size(b.z,3))
    #qv = Dict(q => zeros(size(b.data[q],3)) for q in keys(b.data))

    qv = Dict{Symbol, AbstractVector{eltype(b.x)}}()
    zv = AbstractVector{eltype(b.x)}
    xv = AbstractVector{eltype(b.x)}
    yv = AbstractVector{eltype(b.x)}

    if check_nan
        @inbounds for j in 1:size(b.x,2)
            @inbounds for i in 1:size(b.x,1)
                zv = @view b.z[i,j,:]
                xv = @view b.x[i,j,:] 
                yv = @view b.y[i,j,:] 
                qv = Dict(q => @view(b.data[q][i,j,:]) for q in keys(b.data))

                mm = (.!(isnan.(zv))) .& (.!(isnan.(xv))) .& (.!(isnan.(yv))) 
                for q in keys(qv)
                    mm = mm .& .!(isnan.(qv[q]))
                end

                qv = Dict(q => @view(b.data[q][i,j,mm]) for q in keys(b.data))
                zv = @view b.z[i,j,mm]
                xv = @view b.x[i,j,mm] 
                yv = @view b.y[i,j,mm]

                res[i,j] = f(; qv..., x=xv, y=yv, z=zv) 
            end
        end
    else
        @inbounds for j in 1:size(b.x,2)
            @inbounds for i in 1:size(b.x,1)
                zv = @view b.z[i,j,:]
                xv = @view b.x[i,j,:] 
                yv = @view b.y[i,j,:] 
                for q in keys(b.data)
                    qv[q] = @view b.data[q][i,j,:]
                end
                res[i,j] = f(; qv..., x=xv, y=yv, z=zv) 
            end
        end
    end
end

height_where(; kwargs...) = begin
    @assert length(keys(kwargs)) == 1
    col = keys(kwargs)[1]
    val = values(kwargs)[1]

    f_mod(; kw...) = begin
        res = kw[col]
        mask  = (.!isnan.(res)) .&  (.!isnan.(kw[:z]))
        smask = sortperm(res[mask])
        umask = uniqueidx(res[mask][smask])
        LinearInterpolation(res[mask][smask][umask], kw[:z][mask][smask][umask], extrapolation_bc=Flat())(val)
    end
    return f_mod
end

"""
Reduce a Box object to one single plane by applying the function f to each plane x,y in z

Parameters
----------
b: Box
    Box that should be reduced
f: Function
    taking 3-D array of a plane and return a single value, which is the desired height index. 
    If index is false, then not the index of the desired height coordinate is returned, but the z value itself.
    The oter quantites are then interpolated to this z value.
    The function has to accept all columns as kwargs!
index: Bool
    
Returns
-------
Box object

Example
-------
f(; tau, columns...) = argmin(tau)
b_new = reduce_by_plane(f, b)

"""
function reduce_by_plane(f::F, b::MUST.Box; index=false) where {F<:Function}
    shape    = (size(b.x,1),size(b.x,2),1)
    data_new = Dict(q=>rand(typeof(b.data[q][1,1,1]), shape) for q in keys(b.data))
    x_new = rand(typeof(b.x[1,1,1]), shape)
    y_new = rand(typeof(b.x[1,1,1]), shape)
    z_new = rand(typeof(b.x[1,1,1]), shape)

    if index
        idx = size(b.z,3) == 0 ? 1 : f(; b.data..., x=b.x, y=b.y, z=b.z) 
        
        for q in keys(b.data)
            data_new[q][:,:,1] .= b.data[q][:,:,idx]
        end
        x_new[:,:,1] .= b.x[:,:,idx]
        y_new[:,:,1] .= b.y[:,:,idx]
        z_new[:,:,1] .= b.z[:,:,idx]
    else
        if size(b.z,3) > 1 
            z_int = f(; b.data..., x=b.x, y=b.y, z=b.z) 

            for j in 1:shape[2]
            for i in 1:shape[1]
                szv = (.!isnan.(b.z[i,j,:])) .& (.!isnan.(b.x[i,j,:])) .& (.!isnan.(b.y[i,j,:]))
                for q in keys(b.data)
                    szv_loc = szv .& (.!isnan.(b.data[q][i,j,:]))
                    zv  = @view b.z[i,j,szv_loc]
                    qv  = @view b.data[q][i,j,szv_loc]
                    zsm = sortperm(zv)
                    data_new[q][i,j,1] = LinearInterpolation(zv[zsm], qv[zsm], extrapolation_bc=Flat())(z_int)
                end

                x_new[i,j,1] = b.x[i,j,1]
                y_new[i,j,1] = b.y[i,j,1]
                z_new[i,j,1] = z_int
            end
            end

        else
            for q in keys(b.data)
                data_new[q][i,j,1] = b.data[q][i,j,1]
            end
            x_new[i,j,1] = b.x[i,j,1]
            y_new[i,j,1] = b.y[i,j,1]
            z_new[i,j,1] = b.z[i,j,1]
        end
    end

    Box(x_new, y_new, z_new, data_new, b.parameter)
end

"""
Proxy for reduce_by_plane which returns the plane closest to value.
The function f has to be a map 2D-Plane -> AbstractFloat, which is compated to value. If
remove_nan is true, NaN will be removed beforehand. Be aware that in that case the array
passed to f will be 1D instead.
"""
function reduce_by_value(f::F, b::MUST.Box; remove_nan=true, kwargs...) where {F<:Function}
    @assert length(kwargs) == 1
    f_mod(; kw...) = begin
        res = zeros(size(kw[:z],3))
        for plane in 1:size(kw[:z],3)
            res[plane] = !remove_nan ? f( kw[keys(kwargs)[1]][:,:,plane] ) : 
                                       f( filter(!isnan,kw[keys(kwargs)[1]][:,:,plane]) )
        end
        #argmin(abs.(res .- kwargs[1]))
        mask = .!isnan.(res) .&  .!isnan.(kw[:z][1,1,:])
        smask = sortperm(res[mask])
        LinearInterpolation(res[mask][smask], kw[:z][1,1,mask][smask], extrapolation_bc=Flat())(kwargs[1])
    end

    reduce_by_plane(f_mod, b; index=false)
end


"""
    interpolate_to(box, v; logspace=true, kwargs...)

Interpolate the `Box` of values v to the plane given in kwargs. Specify `logspace=true`
if the given quantity is in log.
"""
function interpolate_to(box, v::Symbol; logspace=true, kwargs...)
    z     = first(keys(kwargs))
    z_val = first(values(kwargs))

    ips = interpolate_by_column(box, z, v, logspace=logspace)
    ip  = ips[v]

    col_new = zeros(eltype(box.x), size(box.x)[1:2]..., 1)
    col_z   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)
    col_x   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)
    col_y   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)

    @inbounds for j in axes(box.x, 2)
        @inbounds for i in axes(box.x, 1)
            col_new[i, j, 1] = ip[i, j](z_val)
            col_z[  i, j, 1] = ips[:z][i, j](z_val)
            col_x[  i, j, 1] = ips[:x][i, j](z_val)
            col_y[  i, j, 1] = ips[:y][i, j](z_val)
        end
    end

    Box(col_x, col_y, col_z, Dict{Symbol,Array{eltype(box.x),3}}(v=>col_new), deepcopy(box.parameter))
end

"""
    interpolate_to(box; logspace=true, kwargs...)

Interpolate the `Box` to the plane given in kwargs. Specify `logspace=true`
if the given quantity is in log.
"""
function interpolate_to(box; logspace=true, kwargs...)
    z     = first(keys(kwargs))
    z_val = first(values(kwargs))

    ips = interpolate_by_column(box, z, logspace=logspace)

    cols_new = Dict(p=>zeros(eltype(box.x), size(box.x)[1:2]..., 1) for p in keys(b.data))
    col_z   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)
    col_x   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)
    col_y   = zeros(eltype(box.x), size(box.x)[1:2]..., 1)

    @inbounds for j in axes(box.x, 2)
        @inbounds for i in axes(box.x, 1)
            for p in keys(b.data)
                cols_new[p][i, j, 1] = ips[p][i, j](z_val)
            end
            col_z[  i, j, 1] = ips[:z][i, j](z_val)
            col_x[  i, j, 1] = ips[:x][i, j](z_val)
            col_y[  i, j, 1] = ips[:y][i, j](z_val)
        end
    end

    Box(col_x, col_y, col_z, Dict{Symbol,Array{eltype(box.x),3}}(v=>cols_new[v] for v in keys(b.data)), deepcopy(box.parameter))
end






"""
Switch the height scale of a box. Creates a new box by reducing the old box by column for every new height value.
"""
function height_scale_slow(b::MUST.Box, new_scale, limits=nothing, logspace=true)
    N_points = size(b.z, 3)
    TT = eltype(b.z)

    if isnothing(limits)
        # create a new height scale 
        min_plane = MUST.plane_statistic(minimum, b, new_scale)
        max_plane = MUST.plane_statistic(maximum, b, new_scale)
        low_lim   = maximum(min_plane)
        high_lim  = minimum(max_plane)
    else
        low_lim  = limits[1]
        high_lim = limits[2]
    end

    if !logspace
        h_scale = range( TT(low_lim), TT(high_lim); length=N_points)
    else
        h_scale = Base.convert.(TT, 10.0 .^ range( log(10, TT(low_lim)), log(10, TT(high_lim)); length=N_points))
    end

    new_box = deepcopy(b)
    interp  = MUST.interpolate_by_column(b)

    # Save loop allocations
    Nx = size(b.x,1)
    Ny = size(b.x,2)
    #xv = zeros(eltype(b.x), N_points); yv=zeros(eltype(b.x), N_points); zv=zeros(eltype(b.x), N_points)
    #qv = Dict(q => zeros(eltype(b.x), size(b.data[q],3)) for q in keys(b.data))
    z_int = zeros(eltype(b.x), Nx, Ny)

    # Step by step construct a box object from individual planes
    @inbounds for i in 1:N_points
        @inline scale_column(; z, kwargs...) = linear_interpolate(kwargs[new_scale], z, h_scale[i])

        apply_by_column!(scale_column, z_int, b; check_nan=false)
        
        @inbounds for j in 1:Ny
            @inbounds for k in 1:Nx
                for q in keys(new_box.data)
                    new_box[q][k,j,i] = interp[q][k,j](z_int[k,j])
                end
                new_box.x[k,j,i] = interp[:x][k,j](z_int[k,j])
                new_box.y[k,j,i] = interp[:y][k,j](z_int[k,j])
                new_box.z[k,j,i] = z_int[k,j]
            end
        end
    end 

    new_box
end


"""
Switch the height scale of a box. Creates a new box by reducing the old box by column for every new height value.
"""
function height_scale(b::MUST.Box, new_scale, limits=nothing, logspace=true)
    N_points = size(b.z, 3)
    TT = eltype(b.z)

    if isnothing(limits)
        # create a new height scale 
        min_plane = MUST.plane_statistic(minimum, b, new_scale)
        max_plane = MUST.plane_statistic(maximum, b, new_scale)
        low_lim   = maximum(min_plane)
        high_lim  = minimum(max_plane)
    else
        low_lim  = limits[1]
        high_lim = limits[2]
    end

    if !logspace
        h_scale = range( TT(low_lim), TT(high_lim); length=N_points)
    else
        h_scale = Base.convert.(TT, 10.0 .^ range( log(10, TT(low_lim)), log(10, TT(high_lim)); length=N_points))
    end

    new_box = deepcopy(b)
    interp  = MUST.interpolate_by_column(b, new_scale, logspace=logspace)

    # Save loop allocations
    Nx = size(b.x,1)
    Ny = size(b.x,2)

    # Step by step construct a box object from individual planes
    @inbounds for i in 1:N_points
        lh = log(10, h_scale[i])
        @inbounds for j in 1:Ny
            @inbounds for k in 1:Nx
                for q in keys(new_box.data)
                    if q == new_scale
                        continue
                    end
                    
                    new_box[q][k,j,i] = interp[q][k,j](lh)
                end
                new_box.z[k,j,i] = interp[:z][k,j](lh)
                new_box[new_scale][k,j,i] = h_scale[i]
            end
        end
    end 

    new_box
end

"""
Switch the height scale of a box. Creates a new box by reducing the old box by column for every new height value.
"""
function height_scale_fast(b::MUST.Box, new_scale; logspace=true, dont_log=[:T, :ux, :uy, :uz, :ee, :e])
    N_points = size(b.z, 3)
    TT = eltype(b.z)

    # create a new height scale 
    min_plane = MUST.plane_statistic(minimum, b, new_scale)
    max_plane = MUST.plane_statistic(maximum, b, new_scale)
    low_lim   = maximum(min_plane)
    high_lim  = minimum(max_plane)

    if !logspace
        h_scale = range( TT(low_lim), TT(high_lim); length=N_points)
    else
        h_scale = Base.convert.(TT, 10.0 .^ range( log(10, TT(low_lim)), log(10, TT(high_lim)); length=N_points))
    end

    new_box = deepcopy(b)
    sorting, interps = interpolators(log.(10, b[new_scale]), log.(10, h_scale))

    # Save loop allocations
    Nx = size(b.x, 1)
    Ny = size(b.x, 2)

    # Step by step construct a box object from individual planes
    @inbounds for j in 1:Ny
        @inbounds for k in 1:Nx
            for q in keys(new_box.data)
                if q == new_scale
                    continue
                end
                @inbounds new_box.data[q][k,j,:] .= if all(b.data[q][k, j, :] .> 0.0) & !(q in dont_log)
                    exp.(
                        gevaluate!(
                            interps[k, j], 
                            @views log.(b.data[q][k, j, sorting[k, j, :]])
                        )
                    )
                else
                    gevaluate!(
                        interps[k, j], 
                        @views b.data[q][k, j, sorting[k, j, :]]
                    )
                end
            end
            @inbounds new_box.z[k, j, :] .= gevaluate!(
                interps[k, j], 
                @views b.z[k, j, sorting[k, j, :]]
            )

            @inbounds new_box.data[new_scale][k, j, :] .= h_scale
        end
    end
    new_box
end

function interpolators(old_cube, new_scale)
    sorting = zeros(Int, size(old_cube)...)
    interps = Matrix{GridInterpolation}(undef, size(old_cube)[1:2]...)
    g_new = Grid(new_scale)
    g_old = similar(old_cube, size(old_cube, 3))
    for j in axes(old_cube, 2)
        for i in axes(old_cube, 1)
            @inbounds g_old .= old_cube[i, j, :]
            @inbounds sorting[i, j, :] .= sortperm(g_old)
            @inbounds interps[i, j] = ginterpolate(Grid(g_old[sorting[i, j, :]]), g_new)
        end
    end
        
    sorting, interps
end


@inline function linear_interpolate(x::A, y::B, x0::T) where {T<:AbstractFloat, A<:AbstractArray{T,1}, B<:AbstractArray{T,1}}
    i_min::Int64 = argmin(abs.(x .- x0))
    i_min = if i_min == length(x)
        i_min -1
    elseif i_min == 1
        i_min +1
    else
        i_min
    end

    i_min2::Int64 = i_min
    if x[i_min] > x0
        # If i_min is larger than the desired value, pick the smaller neighbor
        i_min2 = argmin([x[i_min-1], x[i_min+1]])
    elseif x[i_min] < x0
        # If i_min is smaller than the desired value, pick the larger neighbor
        i_min2 = argmax([x[i_min-1], x[i_min+1]])
    else 
        return y[i_min]
    end

    m::T = (y[i_min+(-1)^i_min2] - y[i_min]) / (x[i_min+(-1)^i_min2] - x[i_min])
    return m * x0 + (y[i_min] - m* x[i_min])
end








#=== Specific functions ===#

optical_depth(ρ::Vector{T}, κ::Vector{T2}, z::Vector{T3}) where {T, T2, T3} = begin
    τ_ross = zeros(T2, length(ρ)) 
    ρκ = ρ .* κ
    for j in eachindex(z)
        if j==1 
            τ_ross[1] = 0 + (z[2] - z[1]) * 0.5 * (ρκ[j])
        else
            τ_ross[j] = τ_ross[j-1] + (z[j] - z[j-1]) * 0.5 * (ρκ[j] + ρκ[j-1])
        end
    end

    τ_ross
end

"""
Compute the optical depth from opacity+density or rk
"""
optical_depth(b::MUST.Box; kwargs... ) = begin
    if :rk in keys(kwargs)
        return optical_depth(b, b.data[kwargs[:rk]])
    else
        return optical_depth(b, b.data[kwargs[:opacity]], b.data[kwargs[:density]])
    end
end

optical_depth(b::MUST.Box, k, rho) = optical_depth(b, k .* rho)

function optical_depth(b::MUST.Box, rk)
    #=τ(z,ix,iy) = begin
        mask   = b.z[ix,iy,:] .>= z
        count(mask) < 2 ? 
            0.0 : 
            integrate(b.z[ix,iy,mask],rk[ix,iy,mask])
    end=#

    Nx = size(b.x,1)
    Ny = size(b.x,2)
    Nz = size(b.x,3)
    optical_depth = zeros(Nx,Ny,Nz)

    @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            z = b.z[i,j,:]
            @inbounds for k in Nz:-1:1
                optical_depth[i,j,k] = if k==Nz
                    0 + abs(z[Nz-1] - z[Nz]) * 0.5 * (rk[i, j, k])
                else
                    optical_depth[i,j,k+1] + abs(z[k] - z[k+1]) * 0.5 * (rk[i, j, k] + rk[i, j, k+1])
                end
            end
        end
    end

    #=for k in Nz-1:-1:1
        for j in 1:Ny
            for i in 1:Nx
                z = b.z[i,j,k]
                optical_depth[i,j,k] = τ(z, i, j)
            end
        end
    end=#

    # Always compute the log, makes the extrapolation easier
    #optical_depth = log10.(optical_depth)

    # Last point extrapolate
    #m = (optical_depth[:,:,end-2] .- optical_depth[:,:,end-1]) ./ (b.z[:,:,end-2] .- b.z[:,:,end-1])
    #a = optical_depth[:,:,end-1] .- m .* b.z[:,:,end-1]
    #optical_depth[:,:,end] = m .* b.z[:,:,end] .+ a

    optical_depth
end

function Boxes(folder::String; snaps=nothing)
    snaps = isnothing(snaps) ? Base.Colon() : snaps
    content_of_folder = glob("*/", folder)
    snapshots = sort(MUST.list_of_snapshots(content_of_folder))[snaps]
    
    boxes     = []
    boxesT    = []

    for (i_s,snap) in enumerate(snapshots)
        try
            append!(boxes, [MUST.Box("box_sn$(snapshots[i_s])", folder=folder)])
        catch
            nothing
        end

        try
            append!(boxesT, [MUST.Box("box_tau_sn$(snapshots[i_s])", folder=folder)])
        catch
            nothing
        end

    end

    return if length(boxes) == 1
        boxes[1], boxesT[1]
    else
        boxes, boxesT
    end
end


"""
    converted_snapshots(folder)

Return a list of already converted snapshots in the convert.jl
naming convention.
"""
function converted_snapshots(folder)
	files_converted = glob("*.hdf5", folder)
	snaps = Dict()
	for file in files_converted
		snname = basename(file)

		
		if occursin("tau", snname) 
			continue 
		end 

        if occursin("tav", snname) 
			continue 
		end 
		
		snid   = parse(Int, snname[last(findfirst("sn", snname))+1:end-5])
		is_τ = isfile(joinpath(folder,"box_tau_sn$(snid).hdf5"))

		snname  = String(split(snname, ".hdf5") |> first)
		sntname = "box_tau_sn$(snid)"

		if is_τ 
			snaps[snid] = (snname, sntname) 
		else
			snaps[snid] = (snname, nothing) 
		end
	end

	snaps["folder"] = folder
	snaps
end

"""
Sort converted snapshots.
"""
list_snapshots(snapshots) = sort([k for k in keys(snapshots) if k != "folder"])

"""
    pick_snapshot(snapshots, i)

Pick the ith snapshot from the list of available snapshots. 
If ``i == :recent``, pick the most recent available snapshot. 
If ``i == :time_average``, pick the time average, if available.
"""
function pick_snapshot(snapshots, i; skip_last_if_missing=true, verbose=0)
	i = if i == :recent 
        if ((isnothing(last(snapshots[last(list_snapshots(snapshots))]))) & 
                        skip_last_if_missing)
            list_snapshots(snapshots)[end-1]
        else 
            last(list_snapshots(snapshots))
        end
    elseif i <0
        list_snapshots(snapshots)[end+i]
	else 
		i
	end

	snap = if i == :time_average
        ["box_tav", "box_tau_tav"]
        isfile(joinpath(snapshots["folder"], "box_tau_tav.hdf5")) ? 
            ["box_tav", "box_tau_tav"] :
            ["box_tav", nothing]
    else   
        snapshots[i]
    end

	if isnothing(last(snap))
		verbose>0 && @info "snapshot $(i) loaded."
		MUST.Box(first(snap), folder=snapshots["folder"]), nothing
	else
		verbose>0 && @info "snapshot $(i) + τ-shot loaded."
		MUST.Box(first(snap), folder=snapshots["folder"]), 
		MUST.Box(last(snap), folder=snapshots["folder"])
	end
end

"""
    pick_snapshot(folder, i)

Pick ith snapshots from all available snapshots in ``folder``.
"""
pick_snapshot(folder::String, i; kwargs...) = pick_snapshot(converted_snapshots(folder), i; kwargs...)











## functions between different boxes

"""
Compute the statistic between every box in boxes array. Need to be on the same grid for this!
""" 
function time_statistic(f::Function, boxes)
    ### Loop through the box and average point by point
    ### storage arrays
    var = zeros(eltype(first(boxes).z), length(boxes))
    mat_var = zeros(eltype(first(boxes).z), size(first(boxes).z))
    storage = deepcopy(first(boxes))
    x,y,z = deepcopy(first(boxes).x), deepcopy(first(boxes).y), deepcopy(first(boxes).z)

    for variable in keys(storage.data) ### Loop through the variables
        variable == :i_patch && continue
        @assert all([ all(size(b.data[variable]) .== size(first(boxes).data[variable])) for b in boxes])
    end

    for variable in keys(storage.data) ### Loop through the variables
        variable == :i_patch && continue
        @inbounds for k in axes(z, 3)
            @inbounds for j in axes(y, 2)
                @inbounds for i in axes(x, 1) ### Loop through the box
                    for b in eachindex(boxes)
                        var[b] = boxes[b].data[variable][i,j,k] 
                    end
                    #### Save the statistic of those variables
                    mat_var[i,j,k] = f(var)
                end
            end
        end
        storage.data[variable] .= mat_var
    end

    storage
end








## Helper functions

axis(b::Box, which, dimension::Int=0) = begin
    if which in [:x, :y, :z]
        dim = findfirst(which .== [:x, :y, :z])
    else
        dim = dimension
    end

    return if dim==0
        sb = size(b)
        view(b[which], [sb[i] > 1 ? Base.Colon() : 1 for i in 1:length(sb)]...)
    else
        view(b[which], [i==dim ? Base.Colon() : 1 for i in 1:3]...)
    end
end

Base.axes(b::Box, args...; kwargs...) = axes(b.z, args...; kwargs...)
Base.size(b::Box, args...; kwargs...) = size(b.z, args...; kwargs...)

closest(a, v) = argmin(abs.(a .- v))
    

is_log(x) = begin
	sx = String(x)
	
	xnew, f = if occursin("log10", sx)
		Symbol(sx[findfirst("log10", sx)[end]+1:end]), log10
	elseif occursin("log", sx)
		Symbol(sx[findfirst("log", sx)[end]+1:end]), log
	else
		x, identity
	end

	xnew, f
end

"""
    profile(f, model, [x=:z, y=:T])

Compute the 1D profile from the 3D Box object.
The `f` can be any funtion that is applied plane-wise to the data.
`log10` or `log` can be specified in the variable name to automatically
return log quantities for convenience when plotting.
"""
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if ndims(axis(model, xs)) > 1
		logx.(axis(model, xs, 3)), logy.(plane_statistic(f, model, ys)) 
	else
		logx.(axis(model, xs)), logy.(plane_statistic(f, model, ys))
	end
end



"""
	time_average_profiles(model_folder)

Time average the ```plane_statistics``` of all snapshots located in that folder.
"""
function time_average_profile(f, model_folder, args...; hscale=:τ, kwargs...)
	x, y = [], []
    which = (hscale==:z) ? first : last 
	for i in list_snapshots(converted_snapshots(model_folder))
		snap = which(pick_snapshot(model_folder, i))
		xi, yi = profile(f, snap, args...; kwargs...)

		append!(x, [xi])
		append!(y, [yi])
	end

	# interpolate to common x axis
	x_common = range(
		maximum(minimum.(x)), 
		minimum(maximum.(x)), 
		length=length(x |> first)
	)
	y_common = []
	for i in eachindex(x)
		m = sortperm(x[i])
		append!(y_common, [linear_interpolation(x[i][m], y[i][m]).(x_common)])
	end
	
	# compute the mean and the std
	y_mean = zeros(length(x_common))
	y_std  = zeros(length(x_common))
	yi     = zeros(length(x))

	for i in eachindex(x_common)
		for j in eachindex(yi)
			yi[j] = y_common[j][i]
		end
		y_mean[i] = MUST.mean(yi)
		y_std[i]  = MUST.std(yi)
	end

	x_common, y_mean, y_std
end



mesh(m::Box) = meshgrid(
	axis(m, :x) ./1e8, 		
	axis(m, :y) ./1e8, 								
	axis(m, :z) ./1e8
)






#= Flipping the atmos in the right direction, and also make the z axis negative below =#

"""
    flip!(box::Box)

Flip the box object and possibly reverse the sign of the height scale (only geometrical).
It used the density to determine where the bottom is.
"""
function flip!(box::Box; density=:d, uz=:uz)
    # Determine if it is oriented from bottom to top
    d = box[density][1, 1, :]
    is_upside_down = first(d) < last(d)
    
    if is_upside_down
        for (k, v) in box.data
            reverse!(box.data[k]; dims=3)
        end

        reverse!(box.x, dims=3)
        reverse!(box.y, dims=3)
        reverse!(box.z, dims=3)
    end

    # Now it is from bottom to top, so the first value in z should be the smallest value
    if box.z[1, 1, 1] > box.z[1, 1, end]
        box.z .*= -1
        box.data[uz] .*= -1
    end

    box
end








#= NEW: Convert from Dispatch entirely in Julia =#
include("_box_direct.jl")