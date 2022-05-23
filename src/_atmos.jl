abstract type AbstractSpace end

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
"""
Atmospheric parameters of a dispatch snapshot.
Creates a link between the patch time, teff as well as physical parameters.
"""
function set!(atmp::AtmosphericParameters, teff_function, ini_namelist)
    teff    = teff_function(atmp.time)
    logg    = log10(ini_namelist.stellar_params["g_cgs"])

    epn = MUST.@in_dispatch(String(strip(ini_namelist.eos_params["table_loc"])[2:end-1]))
    eos_params_path = joinpath(epn, "tabparam.in")
    
    composition = composition_from_eos(eos_params_path)

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
    data = readdlm(path)
end 

teff_interpolated(path) = begin
    teff = read_teff(path)
    LinearInterpolation(teff[:,1], teff[:,2], extrapolation_bc=Flat())
end

function composition_from_eos(path)
    content = read_eos_params(path)
    el      = split(content["cel"][2:end-1], " ")
    el      = [String(s) for s in el if length(s)>0]
    vl      = String.(split(content["abund"][2:end-1], " "))
    vl      = [s for s in vl if length(s)>0]
    vl      = parse.(Float64, vl)
    Dict(Symbol(c)=>v for (c,v) in zip(el,vl))
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
function Space(snapshot::PyCall.PyObject, quantities::Symbol...) 
    qs = Dict{Symbol,Vector{Float32}}(q=>Float32[] for q in quantities)
    qs[:x] = Float32[]; qs[:y] = Float32[]; qs[:z] = Float32[]
    qs[:i_patch] = Int[]
    patch_dimensions = zeros(Int,(length(snapshot.patches),3))
    
    for (i,patch) in enumerate(snapshot.patches)
        
        # grid of this patch relative to the global grid
        patch_size = (length(patch.xi),length(patch.yi),length(patch.zi))
        x = patch.xi; y = patch.yi; z = patch.zi
        patch_dimensions[i,1] = length(patch.xi)
        patch_dimensions[i,2] = length(patch.yi)
        patch_dimensions[i,3] = length(patch.zi)
        
        # extract the quantites for this patch
        for (iq,q) in enumerate(quantities)
            q_matrix = patch.var(String(q))
            q_array  = zeros(Float32, prod(patch_size))
            coords   = zeros(Float32, (prod(patch_size),3))
            j = 1
            for iz in 1:patch_size[3]
                for iy in 1:patch_size[3]
                    for ix in 1:patch_size[3]
                        q_array[j]  = q_matrix[ix,iy,iz]
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

    time = snapshot.nml_list["snapshot_nml"]["time"]
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
    region_mask = (row) ->  (minimum(x) <= row.x <= maximum(x)) &
                            (minimum(y) <= row.y <= maximum(y)) &
                            (minimum(z) <= row.z <= maximum(z))
    s_masked = filter(region_mask,s)

    # 3D Mesh of the new coordinate Box
    x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z, indexing="ij")

    # result dict
    results::Dict{Symbol,Array{T,3}} = Dict(q=>similar(x_grid) for q in keys(s.data) if !(q in [:x,:y,:z,:i_patch]))

    mask_data = []
    grid_data = []
    if length(x) > 1
        append!(mask_data, [s_masked.data.x])
        append!(grid_data, [x_grid])
    end
    if length(y) > 1
        append!(mask_data, [s_masked.data.y])
        append!(grid_data, [y_grid])
    end
    if length(z) > 1
        append!(mask_data, [s_masked.data.z])
        append!(grid_data, [z_grid])
    end


    # Interpolate the Space onto the given grid
    for quantity in Symbol.(names(s.data))
        quantity in [:x,:y,:z,:i_patch] ? continue : nothing
        results[quantity][1:length(x),1:length(y),1:length(z)] = scipy_interpolate.griddata(Tuple(mask_data), 
                                                                                            s_masked.data[quantity], 
                                                                                            Tuple(grid_data))
    end

    Box(x_grid[2:end-1,2:end-1,2:end-1],
        y_grid[2:end-1,2:end-1,2:end-1],
        z_grid[2:end-1,2:end-1,2:end-1],
        Dict(q=>results[q][2:end-1,2:end-1,2:end-1] for q in keys(results)), s.parameter)
end

@inline uniform_grid(s::Space, N::Int,axis=:x) = begin
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

    x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z, indexing="ij")
    x_grid = Base.convert(Array{typeof(x[1])}, x_grid)
    y_grid = Base.convert(Array{typeof(x[1])}, y_grid)
    z_grid = Base.convert(Array{typeof(x[1])}, z_grid)

    results::Dict{Symbol, Array{<:Union{Float32, Float64, Int16, Int32, Int64}, 3}} = Dict(q=>zeros(typeof(s.data[q][1]),size(x_grid)) 
                                                                                        for q in keys(s.data) if !(q in [:x,:y,:z]))
    for q in keys(results)
        q == :i_patch ? continue : nothing
        results[q] .= Base.convert(typeof(results[q][1,1,1]), NaN)
    end

    p        = zeros(3) 
    skip_q   = [:x,:y,:z]
    for irow in 1:length(s.data[:x])
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
    ix = findfirst(p[1] .≈ x)
    iy = findfirst(p[2] .≈ y)
    iz = findfirst(p[3] .≈ z)
    (ix,iy,iz)
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
function convert!(s::Space, u::AtmosUnits; params...)
    for (s_para, u_para) in params
        s.data[s_para] .= s.data[s_para] .* getfield(u,u_para)
    end
    nothing
end

function convert(s::Space, u::AtmosUnits; params...)
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
function add_from_EOS!(s::MUST.Space, eos, quantity::Symbol;
                        EOS_paras=(:d,:ee), EOS_log=(true,true),
                        convert_to=Base.identity, switch_name_to=nothing)
    qframe = MUST.add_from_EOS(s, eos, quantity; 
                                EOS_paras =EOS_paras, 
                                EOS_log   =EOS_log, 
                                convert_to=convert_to)
    isnothing(switch_name_to) ? s.data[quantity] = qframe : s.data[switch_name_to] = qframe
    nothing
end
function add_from_EOS(s::MUST.Space, eos, quantity::Symbol;
                        EOS_paras=(:d,:ee), EOS_log=(true,true),
                        convert_to=Base.identity)
    f1 = EOS_log[1] ? log : Base.identity
    f2 = EOS_log[2] ? log : Base.identity
    convert_to.(eos.lookup(String(quantity), f1.(s[EOS_paras[1]]), f2.(s[EOS_paras[2]])))
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
        av[i] = stat(b.data[q][:,:,i])
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

function save(p::AtmosphericParameters, fid)
    eles = [keys(p.composition)...]
    vals = [p.composition[e] for e in eles]
    fid["time"] = p.time
    fid["teff"] = p.teff
    fid["logg"] = p.logg
    fid["composition_e"] = String.(eles)
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
        @assert isfile(path) 
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
_read_params(T::Type{AtmosphericParameters}, res) = T(HDF5.read(res["time"]), HDF5.read(res["teff"]), HDF5.read(res["logg"]), 
                                                Dict(Symbol(e) => v for (e,v) in zip(HDF5.read(res["composition_e"]), HDF5.read(res["composition_v"]))) )

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

            mm = (.!(isnan.(zv))) .& (.!(isnan.(zv))) .& (.!(isnan.(zv))) 
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
b_new = reduce_by_plane(b, f)

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
