abstract type AbstractSpace end

mutable struct Space <:AbstractSpace
    data             ::Dict{Symbol,Vector{Any}}
    patch_dimensions ::Array{Int,2}
end

mutable struct Box{T<:AbstractFloat} <:AbstractSpace
    x                ::Array{T,3}
    y                ::Array{T,3}
    z                ::Array{T,3}
    data             ::Dict{Symbol,Array{T,3}}
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

    Space(qs, patch_dimensions)
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
    qs = Dict{Symbol,Vector{Float32}}(q=>Float32[] for q in quantities)
    qs[:x] = Float32[]; qs[:y] = Float32[]; qs[:z] = Float32[]
    qs[:i_patch] = Int[]
    patch_dimensions = zeros(Int, (1,3))
    patch_dimensions[1,1] = snap.nx*snap.ny*snap.nz
    patch_dimensions[1,2] = snap.nx*snap.ny*snap.nz
    patch_dimensions[1,3] = snap.nx*snap.ny*snap.nz

    qm = zeros(Float32, (snap.nx, snap.ny, snap.nz))
    for (i,q) in enumerate(quantities)
        qm = snap.read_mem(String(q))
        qs[q] = zeros(Float32, snap.nx*snap.ny*snap.nz)
        if i==1
            coords = zeros(Float32, (snap.nx*snap.ny*snap.nz,3))
        end
        c=1
        for k in 1:snap.nz
        for j in 1:snap.ny
        for n in 1:snap.nx
            qs[q][c] = qm[n,j,k]
            if i==1
                coords[c,1] = snap.xx[n]
                coords[c,2] = snap.yy[j]
                coords[c,3] = snap.zz[k]
            end
            c+=1
        end
        end
        end
        if i==1
            append!(qs[:x],coords[:,1])
            append!(qs[:y],coords[:,2])
            append!(qs[:z],coords[:,3])
            append!(qs[:i_patch],[1 for _ in size(coords,1)])
        end
    end
    Space(qs, patch_dimensions)
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
        Dict(q=>results[q][2:end-1,2:end-1,2:end-1] for q in keys(results)))
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
    box    = Box(s, x, y, z, args...; kwargs...)
    
    q_dict = Dict(q => reshape(box.data[q], :) for q in keys(box.data))
    q_dict[:x] = reshape(box.x, :); q_dict[:y] = reshape(box.x, :); q_dict[:z] = reshape(box.x, :)
    q_dict[:i_patch] = ones(length(q_dict[:x]))
    
    p_dim = ones(1,3)
    p_dim[1,1] = length(x); p_dim[1,2] = length(y); p_dim[1,3] = length(z)
    
    Space(q_dict,p_dim)
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
    x = unique(s.data[:x])
    y = unique(s.data[:y])
    z = unique(s.data[:z])

    x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z, indexing="ij")
    results::Dict{Symbol,Array{typeof(x[1]),3}} = Dict(q=>similar(x_grid) for q in keys(s.data) if !(q in [:x,:y,:z,:i_patch]))

    for quantity in keys(s.data)
        quantity in [:x,:y,:z,:i_patch] ? continue : nothing
        # go though the data points and save the quantity
        p        = zeros(3) 
        for irow in 1:length(s.data[quantity])
            p[1] = s.data[:x][irow]
            p[2] = s.data[:y][irow] 
            p[3] = s.data[:z][irow] 
            ix,iy,iz = _find_in_meshgrid(p, x, y, z)
            results[quantity][ix,iy,iz] = s.data[quantity][irow]
        end
    end

    Box(x_grid, y_grid, z_grid, results)
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
    mask = f(s.data)
	new_frame = Dict(q=>s.data[q][mask] for q in keys(s.data))
	Space(new_frame, s.patch_dimensions)
end

function Base.filter!(f, s::MUST.Space)
    mask = f(s.data)
    for q in keys(s.data)
        s.data[q] = s.data[q][mask]
	end
end

Base.getindex(s::T,key) where {T<:Union{MUST.Box,MUST.Space}} = s.data[key]

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
                        convert_to=identity)
    qframe = MUST.add_from_EOS(s, eos, quantity; 
                                EOS_paras =EOS_paras, 
                                EOS_log   =EOS_log, 
                                convert_to=convert_to)
    s.data[quantity] = qframe
    nothing
end
function add_from_EOS(s::MUST.Space, eos, quantity::Symbol;
    EOS_paras=(:d,:ee), EOS_log=(true,true),
    convert_to=identity)
    f1 = EOS_log[1] ? log : indentity
    f2 = EOS_log[2] ? log : indentity
    convert_to.(eos.lookup.(String(quantity), f1.(s[EOS_paras[1]]), f2.(s[EOS_paras[2]])))
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
function plane_average(b::MUST.Box, q::Symbol)
    plane_statistic(Statistics.mean, b, q)
end

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