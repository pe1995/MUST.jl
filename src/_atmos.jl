abstract type AbstractSpace end

mutable struct Space <:AbstractSpace
    data             ::DataFrames.DataFrame
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
        for q in quantities
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
            append!(qs[:x],coords[:,1])
            append!(qs[:y],coords[:,2])
            append!(qs[:z],coords[:,3])
            append!(qs[:i_patch],[i for _ in 1:length(q_array)])
        end
    end

    Space(DataFrame(qs), patch_dimensions)
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
    region_mask = (row) ->  (minimum(x) < row.x < maximum(x)) &
                            (minimum(y) < row.y < maximum(y)) &
                            (minimum(z) < row.x < maximum(z))
    s_masked = filter(region_mask,s)

    # 3D Mesh of the new coordinate Box
    x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z)

    # result dict
    results::Dict{Symbol,Array{T,3}} = Dict(q=>similar(x_grid) for q in Symbol.(names(s.data)) if !(q in [:x,:y,:z,:i_patch]))

    # Interpolate the Space onto the given grid
    for quantity in Symbol.(names(s.data))
        quantity in [:x,:y,:z,:i_patch] ? continue : nothing
        results[quantity] .= scipy_interpolate.griddata((s_masked.data.x, s_masked.data.y, s_masked.data.z), 
                                                            s_masked.data[:,quantity], 
                                                            (x_grid, y_grid, z_grid))
    end

    Box(x_grid,y_grid,z_grid,results)
end

@inline uniform_grid(s::Space, N::Int,axis=:x) = begin
    data = getfield(s,axis)
    maxd = maximum(data)
    mind = minimum(data)
    extd = maxd-mind

    mind += extd*0.01
    maxd -= extd*0.01
    extd = maxd-mind
    step = extd/N
    Vector(mind:step:maxd)
end
Box(s::Space, N::Int, args...; kwargs...) = Box(s,
                            uniform_grid(s, N, :x), 
                            uniform_grid(s, N, :y), 
                            uniform_grid(s, N, :z), args...; kwargs...)

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
    
    Space(DataFrame(q_dict),p_dim)
end
spacebox(s::Space, N::Int, args...; kwargs...) = spacebox(s,
                                                    uniform_grid(s, N, :x), 
                                                    uniform_grid(s, N, :y), 
                                                    uniform_grid(s, N, :z), args...; kwargs...)

"""
Convert Space to Box, assuming that the space has already uniform grid (No interpolation)
"""
function spacebox(s::Space)    
    # create x,y,z from unique s.data -> meshgrid x,y,z -> pick data from dataframe where coordinates match
    x = unique(s.data.x)
    y = unique(s.data.y)
    z = unique(s.data.z)

    x_grid, y_grid, z_grid = numpy.meshgrid(x, y, z)
    results::Dict{Symbol,Array{T,3}} = Dict(q=>similar(x_grid) for q in Symbol.(names(s.data)) if !(q in [:x,:y,:z,:i_patch]))

    for quantity in Symbol.(names(s.data))
        quantity in [:x,:y,:z,:i_patch] ? continue : nothing
        # go though the data points and save the quantity
        for row in s.data
            p        = (row.x, row.y, row.z)
            ix,iy,iz = _find_in_meshgrid(p, x, y, z)
            results[quantity][ix,iy,iz] = getfield(row,quantity)
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
	new_frame = filter(f, s.data)
	Space(new_frame, s.patch_dimensions)
end

"""
Filter the Space object.

Parameters
----------
s: Space
    The space object to filter
    
Returns
-------
nothing

"""
function Base.filter!(s::MUST.Space, args...; kwargs...)
	new_frame = filter(s.data,args...; kwargs...)
	s.data    = new_frame
    nothing
end

Base.getindex(s::T,key) where {T<:MUST.Space} = s.data[:,key]
Base.getindex(s::T,key) where {T<:MUST.Box}   = s.data[key]