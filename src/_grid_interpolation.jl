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
    method::Symbol
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
AxisInterpolation(old_axis::AbstractBoxAxis, new_axis::AbstractBoxAxis; method=:linear) = begin
    w, i = interpolation_weights(nodes(old_axis), nodes(new_axis))

    AxisInterpolation(old_axis, new_axis, w, i, method)
end


"""
    GridInterpolation(old_axis, new_axis)

Construct a grid interpolator from a list of axis interpolation objects.
"""
GridInterpolation(from::Vector{<:AbstractBoxAxis}, to::Vector{<:AbstractBoxAxis}; method=:linear) = begin
    interpolators = AxisInterpolation.(from, to; method=method)

    initial_size = [length(nodes(interpolators[i].old_axis)) for i in eachindex(interpolators)]
    buffers = [similar(interpolators[1].weights, initial_size...)]

    s_old = (length(nodes(f)) for f in from) |> collect
    s_new = (size(interpolators[j].weights, 1) for j in eachindex(interpolators)) |> collect
    for i in eachindex(interpolators)
        s = (j <= i ? s_new[j] : s_old[j] for j in eachindex(interpolators))
        values_int = similar(interpolators[i].weights, s...)

        append!(buffers, [values_int])
    end

    GridInterpolation(interpolators, buffers)
end

GridInterpolation(from::AbstractBoxGrid, to::AbstractBoxGrid; method=:linear) = begin
    GridInterpolation(axis(from), axis(to); method=method)
end






#= General convenience =#

method(ip::AxisInterpolation) = ip.method
method(ip::GridInterpolation) = [i.method for i in ip.axes_interpolation]

fromaxis(ip::AxisInterpolation) = ip.old_axis
toaxis(ip::AxisInterpolation) = ip.new_axis







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
		x0, x1 = if outgrid[i] < first(ingrid)
            # Extrapolate below the first point
            indices[i] = 1
            x0, x1 = ingrid[1], ingrid[2]

			#weights[i, 1] = 1.0
			#weights[i, 2] = 0.0
		elseif outgrid[i] >= last(ingrid)
			# Extrapolate above the last point
            indices[i] = length(ingrid) - 1
            x0, x1 = ingrid[end-1], ingrid[end]

			#weights[i, 1] = 0.0
			#weights[i, 2] = 1.0
		else
			indices[i] = findfirst(x->x>outgrid[i], ingrid) -1	

			x0, x1 = ingrid[indices[i]], ingrid[indices[i]+1]
			
			#weights[i, 1] = (x1-outgrid[i]) / (x1-x0)
			#weights[i, 2] = (outgrid[i]-x0) / (x1-x0)
		end

        weights[i, 1] = (x1-outgrid[i]) / (x1-x0)
	    weights[i, 2] = (outgrid[i]-x0) / (x1-x0)
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
		@inbounds values_new[i] = weights[i, 1] * values[indices[i]] + 
					weights[i, 2] * values[indices[i]+1]
	end

	values_new
end

"""
    pchip_mono8!(newy, yy, newx, xx)

Monotonic piecewise cubic hermite interpolation
Coded by Richard Hoppe, July 2023
"""
function pchip_mono8!(newy, yy, newx, xx)
    # real(8), intent(in) :: xx(:)            # X values of the data points
    # real,    intent(in) :: yy(:)            # Y values of the data points
    # real(8), intent(in) :: newx(:)          # X values to interpolate at
    # real,    intent(out) :: newy(:)         # Interpolated Y values
    # real(8), allocatable :: dydx(:)         # derivative
    # real(8), allocatable :: dx(:)
    # real(8), allocatable :: dy(:)
    # integer :: n                 ! Number of data points
    # integer :: n_interp          ! Number of interpolation points
    # integer :: j                 ! Loop index for interpolation points
    # integer :: i                 ! Loop index for data points
    # real(8) :: d1, d2, d, alpha, t1, t2, t3, p, a
    n = length(xx)
    n_interp = length(newx)

    dydx = fill!(similar(xx, n), 0.0)
    dx   = similar(xx, n-1)
    dy   = similar(xx, n-1)

    dx = xx[2:n] .- xx[1:n-1]
    dy = yy[2:n] .- yy[1:n-1]

    @inbounds for i in 2:n-1
        if (dy[i-1] == 0) | (dy[i] == 0)
            continue
        end

        d1 = dy[i-1]/dx[i-1]
        d2 = dy[i]/dx[i]

        if (sign(d1) == sign(d2))
            # harmonic mean derivative
            alpha = (1. + dx[i-1]/(dx[i] + dx[i-1])) / 3.0
            dydx[i] = 1. / ((1. - alpha)/d2 + alpha/d1)
        end
    end

    # boundary slopes from quadratic polynomial using 2 points and 1 dydx
    a = (dy[1] - dx[1]*dydx[2]) / (xx[2]^2 - xx[1]^2 - 2. * xx[2]*dx[1])
    dydx[1] = dydx[2] - 2. * a *dx[1]

    a = (dy[n-1] - dx[n-1]*dydx[n-1]) / (xx[n]^2 - xx[n-1]^2 - 2. * xx[n-1]*dx[n-1])
    dydx[n] = dydx[n-1] + 2. * a *dx[n-1]

    i = 1
    @inbounds for j in 1:n_interp
        # Find the interval containing newx(j)
        @inbounds for k in i:n-1
            if (newx[j] <= xx[k + 1])
                i = k
                break #! Found the interval
            end
        end
        i = min(i, n-1)
        
        # Calculate the interpolation parameter
        t1 = (newx[j] - xx[i]) / dx[i]

        t2 = t1*t1
        t3 = t2*t1
        p = 3.0*t2 - 2.0*t3

        newy[j] = (1.0 - p)*yy[i] + p*yy[i+1] + dx[i]*((t3 - 2.0*t2 + t1)*dydx[i] + (t3-t2)*dydx[i+1]) 
    end
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
    buffers = ip.buffers
    buffers[1] .= values

    for i in eachindex(ip.axes_interpolation)
        # For each interpolation axis (x, y, z, etc.)
        dim_selection = (j==i ? 1 : Base.Colon() for j in eachindex(ip.axes_interpolation)) |> collect
        c = CartesianIndices(view(ip.buffers[i], dim_selection...))
        for ci in c
            # Now we loop through all other dimensions, and pick the axis
            # we are currently interpolating
            #idx = _fill_index(Base.Colon(), ci, i)
            idx_from = _fill_index(permutation(fromaxis(ip.axes_interpolation[i])), ci, i)
            idx_to = _fill_index(permutation(toaxis(ip.axes_interpolation[i])), ci, i)
            if method(ip.axes_interpolation[i]) == :linear
                interpolate_axis!(  
                    @view(buffers[i+1][idx_to...]),
                    @view(buffers[i][idx_from...]),
                    ip.axes_interpolation[i].weights,
                    ip.axes_interpolation[i].indices
                )
            elseif method(ip.axes_interpolation[i]) == :pchip
                pchip_mono8!(  
                    @view(buffers[i+1][idx_to...]),
                    @view(buffers[i][idx_from...]),
                    nodes(ip.axes_interpolation[i].new_axis),
                    nodes(ip.axes_interpolation[i].old_axis)
                )
            else
                error("Interpolation method not recognized. Pick :linear or :pchip")
            end
        end
    end


    return_copy ? deepcopy(buffers[end]) : buffers[end]
end    

function _fill_index(to_fill, cartesian, at)
    new_index = Any[Base.Colon(), Tuple(cartesian)...]
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