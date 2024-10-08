#= MUST Grids =#

"""
    RegularBoxAxis(values)

Regular Box Axis constructed from 1D array of axis nodes.
It memorizes if the axis needs to be permuted in order to interpolate its data.
"""
struct RegularBoxAxis{T} <: AbstractBoxAxis
    nodes    ::Vector{T}
    sorted   ::Bool
    sortmask ::Vector{Int64}
end

"""
    RegularBoxGrid(axis...)

Regular Box grid constructed from regular Box axes.
"""
struct RegularBoxGrid{B<:AbstractBoxAxis} <: AbstractBoxGrid
    axes::Vector{B}
end




#= Constructors =#

RegularBoxGrid(axes::AbstractBoxAxis...) = RegularBoxGrid([axes...])
RegularBoxGrid(axes::AbstractArray{T,1}...) where {T<:AbstractFloat} = begin
    RegularBoxGrid(RegularBoxAxis.(axes))
end
#=Axis(values) = if is_uniform(values)
    #@show flipped(values)
    RegularBoxAxis(flipped(values))
else
    #@warn "Non-uniform spacing detected."
    #@show flipped(values)
    RegularBoxAxis(flipped(values))
end
=#
Axis(values) = RegularBoxAxis(sort(values), issorted(values), sortperm(values))

"""
    Grid(ax::AbstractVector...) 

Construct a N dimensional interpolation grid for N given arrays.
"""
Grid(ax::AbstractVector...) = RegularBoxGrid(Axis.(ax)...)






#= Utilities =#

permutation(a::RegularBoxAxis) = a.sortmask

"""
    is_uniform(values)

Check if the given axis is regular.
"""
is_uniform(values) = begin
    diff1 = values[2] - values[1]
    all(diff(values) .≈ diff1)
end

flipped(values) = begin
    if issorted(values) 
        values
    else
        #@warn "One Axis is not sorted!"
        if issorted(reverse(values))
            reverse(values)
        else
            error("One axis is not sorted, reversing it did not fix the issue.")
        end
    end
end

nodes(a::AbstractBoxAxis) = a.nodes
nodes(g::RegularBoxGrid) = [nodes(a) for a in axis(g)]
axis(g::RegularBoxGrid) = g.axes
 