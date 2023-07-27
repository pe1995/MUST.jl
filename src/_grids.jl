#= MUST Grids =#

"""
    RegularBoxAxis(values)

Regular Box Axis constructed from 1D array of axis nodes.
"""
struct RegularBoxAxis{T} <: AbstractBoxAxis
    nodes::Vector{T}
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
RegularBoxAxis(values)  = if is_uniform(values)
    RegularBoxAxis(flipped(values))
else
    @warn "Non-uniform spacing detected."
    RegularBoxAxis(flipped(values))
end

Grid(ax::AbstractVector...) = RegularBoxGrid(RegularBoxAxis.(ax)...)





#= Utilities =#

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
        @warn "One Axis is not sorted!"
        if issorted(reverse(values))
            reverse(values)
        else
            error("One axis is not sorted, reversing it did not fix the issue.")
        end
    end
end

nodes(axis::AbstractBoxAxis) = axis.nodes
nodes(axis::RegularBoxGrid)  = [nodes(a) for a in axis]
axis(g::RegularBoxGrid)      = g.axes
 