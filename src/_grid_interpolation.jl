#= Methods for Grid interpolation =#
#= 
    provide capabilities of interpolating 
    each axis of a grid from the current axis to a new axis. 
    For each axis, it creates a AxisInterpolation object, that stores 
    buffers and interpolation weights, so that the integration 
    of muliple variables can be done in no time after the first call.

    For grids, the axis interpolation is called one step at a time,
    and the returned values are further interpolated in the other directions.
=#

struct AxisInterpolation{T<:AbstractFloat, I<:Integer, 
                        A1<:AbstractBoxAxis, A2<:AbstractBoxAxis}
    old_axis::A1
    new_axis::A2
    weights::Array{T, 2}
    indices::Array{I, 1}
end

struct GridInterpolation{A<:AxisInterpolation}
    axes_interpolation::Vector{A}
    buffers
end





#= Constructors =#

"""
    AxisInterpolation(old_axis, new_axis)

Compute weights and indices corresponding to the 
interpolation from the old_axis to the new_axis and
store the values in an AxisInterpolation object.
"""
AxisInterpolation(old_axis::AbstractBoxAxis, new_axis::AbstractBoxAxis) = begin
    w, i = interpolation_weights(nodes(old_axis), nodes(new_axis))

    AxisInterpolation(old_axis, new_axis, w, i)
end


"""
    GridInterpolation(old_axis, new_axis)

Construct a grid interpolator from a list of axis interpolation objects.
"""
GridInterpolation(from::Vector{<:AbstractBoxAxis}, to::Vector{<:AbstractBoxAxis}) = begin
    interpolators = AxisInterpolation.(from, to)
    buffers = []

    s_old = (length(nodes(f)) for f in from) |> collect
    s_new = (size(interpolators[j].weights, 1) for j in eachindex(interpolators)) |> collect
    for i in eachindex(interpolators)
        #i == length(interpolators) ? break : nothing

        s = (j <= i ? s_new[j] : s_old[j] for j in eachindex(interpolators))
        values_int = similar(interpolators[i].weights, s...)

        append!(buffers, [values_int])
    end

    GridInterpolation(interpolators, buffers)
end

GridInterpolation(from::AbstractBoxGrid, to::AbstractBoxGrid) = begin
    GridInterpolation(axis(from), axis(to))
end






#= Interpolation Utilities =#

"""
    interpolation_weights(in, out)

Compute the linear interpolation weights corresponding
to the 2 arrays.
"""
function interpolation_weights(ingrid, outgrid) 
	weights = similar(outgrid, length(outgrid), 2)
	indices = zeros(Int, length(outgrid))

	for i in eachindex(outgrid)
		if outgrid[i] < first(ingrid)
			indices[i] = 1
			weights[i, 1] = 1.0
			weights[i, 2] = 0.0
		elseif outgrid[i] >= last(ingrid)
			indices[i] = length(ingrid)-1
			weights[i, 1] = 0.0
			weights[i, 2] = 1.0
		else
			indices[i] = findfirst(x->x>outgrid[i], ingrid) -1	

			x0 = ingrid[indices[i]]
			x1 = ingrid[indices[i]+1]
			
			weights[i, 1] = (x1-outgrid[i]) / (x1-x0)
			weights[i, 2] = (outgrid[i]-x0) / (x1-x0)
		end
	end

	weights, indices
end

"""
    interpolate_axis(values, weights, indices)

Apply weights to values at indices in order to interpolate
from an old to a new grid.
"""
interpolate_axis(values, weights, indices) = begin
	values_new = similar(values, size(weights, 1))
    interpolate_axis!(values_new, values, weights, indices)

    values_new
end

"""
    interpolate_axis!(values_new, values, weights, indices)

Apply weights to values at indices in order to interpolate
from an old to a new grid.
"""
function interpolate_axis!(values_new, values, weights, indices)	
	for i in eachindex(values_new)
		values_new[i] = weights[i, 1] * values[indices[i]] + 
					weights[i, 2] * values[indices[i]+1]
	end

	values_new
end

"""
    interpolate_grid!(grid, values)

Interpolate the values between two grids consisting of any number
of axis as given by the grid object. The interpolation
is done piece-wise linear in the order given when the grid was constructed.
Intermediate results are saved and in the buffers, which will be re-used
upon a second visit of this function. A copy of the result is returned, unless
requested otherwise.
"""
function interpolate_grid!(ip::GridInterpolation, values; return_copy=true)
    #idx = zeros(Int, length(ip.axes_interpolation))
    buffers = [deepcopy(values), ip.buffers...]

    for i in eachindex(ip.axes_interpolation)
        # For each interpolation axis (x, y, z, etc.)
        dim_selection = (j==i ? 1 : Base.Colon() for j in eachindex(ip.axes_interpolation)) |> collect
        c = CartesianIndices(view(ip.buffers[i], dim_selection...))
        for ci in c
            # Now we loop through all other dimensions, and pick the axis
            # we are currently interpolating
            idx = _fill_index(Base.Colon(), ci, i)
            interpolate_axis!(  @view(buffers[i+1][idx...]),
                                @view(buffers[i][idx...]), 
                                ip.axes_interpolation[i].weights,
                                ip.axes_interpolation[i].indices)
        end
    end


    return_copy ? deepcopy(buffers[end]) : buffers[end]
end    

function _fill_index(to_fill, cartesian, at)
    new_index = [Base.Colon(), Tuple(cartesian)...]
    i = 1
    for j in eachindex(new_index)
        new_index[j] = if j==at
            to_fill
        else
            i += 1
            cartesian[i-1]
        end
    end
    new_index
end





#= Interface =#

"""
    ginterpolate(grid, target_grid)

Create an interpolation object and compute interpolation weights, 
that interpolate from the grid to the target grid.
"""
ginterpolate(args...; kwargs...) = GridInterpolation(args...; kwargs...)

"""
    gevaluate!(ip::GridInterpolation, values::Array)

Evaluate the interpolation weights for all values in `values`.
Effectively interpolates these values to the new grid that was
given to the interpolator `ip`.
"""
gevaluate!(args...; kwargs...) = interpolate_grid!(args...; kwargs...)
gevaluate(grid_in, grid_out, values, args...; kwargs...)  = begin
    ip = ginterpolate(grid_in, grid_out)
    gevaluate!(ip, values, args...; kwargs...)
end
