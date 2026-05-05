# ============================================================================= 
# Methods for Grid interpolation 
# ============================================================================= 
# provide capabilities of interpolating 
# each axis of a grid from the current axis to a new axis. 
# For each axis, it creates a AxisInterpolation object, that stores 
# buffers and interpolation weights, so that the integration 
# of muliple variables can be done in no time after the first call.
# For grids, the axis interpolation is called one step at a time,
# and the returned values are further interpolated in the other directions.
# ============================================================================= 

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

# ============================================================================= 
# Constructors 
# ============================================================================= 

"""
    AxisInterpolation(old_axis, new_axis)

Compute weights and indices corresponding to the 
interpolation from the old_axis to the new_axis and
store the values in an AxisInterpolation object.
"""
AxisInterpolation(old_axis::AbstractBoxAxis, new_axis::AbstractBoxAxis; method=:linear, extrapolation=:linear) = begin
    w, i = interpolation_weights(nodes(old_axis), nodes(new_axis); extrapolation=extrapolation)

    AxisInterpolation(old_axis, new_axis, w, i, method)
end

"""
    GridInterpolation(old_axis, new_axis)

Construct a grid interpolator from a list of axis interpolation objects.
"""
GridInterpolation(from::Vector{<:AbstractBoxAxis}, to::Vector{<:AbstractBoxAxis}; method=:linear, extrapolation=:linear) = begin
    interpolators = AxisInterpolation.(from, to; method=method, extrapolation=extrapolation)

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

GridInterpolation(from::AbstractBoxGrid, to::AbstractBoxGrid; method=:linear, extrapolation=:linear) = begin
    GridInterpolation(axis(from), axis(to); method=method, extrapolation=extrapolation)
end

# ============================================================================= 
# General convenience 
# ============================================================================= 

method(ip::AxisInterpolation) = ip.method
method(ip::GridInterpolation) = [i.method for i in ip.axes_interpolation]

fromaxis(ip::AxisInterpolation) = ip.old_axis
toaxis(ip::AxisInterpolation) = ip.new_axis

# ============================================================================= 
# Interpolation Utilities 
# ============================================================================= 

"""
    interpolation_weights(in, out)

Compute the linear interpolation weights corresponding
to the 2 arrays.
"""
function interpolation_weights(ingrid, outgrid; extrapolation=:linear) 
	weights = similar(outgrid, length(outgrid), 2)
	indices = zeros(Int, length(outgrid))
    outtype = eltype(outgrid)
    nantype = Base.convert(outtype, NaN)

	for i in eachindex(outgrid)
		x0, x1 = if outgrid[i] < first(ingrid)
            indices[i] = 1
            x0, x1 = if extrapolation==:linear
                ingrid[1], ingrid[2]
            elseif extrapolation==:NaN
                nantype, nantype
            elseif extrapolation==:constant
                outgrid[i], outgrid[i]+1
            else
                error("Unknown extrapolation method: $extrapolation")
            end
		elseif outgrid[i] > last(ingrid)
            indices[i] = length(ingrid) - 1
            x0, x1 = if extrapolation==:linear
                ingrid[end-1], ingrid[end]
            elseif extrapolation==:NaN
                nantype, nantype
            elseif extrapolation==:constant
                outgrid[i]-1, outgrid[i]
            else
                error("Unknown extrapolation method: $extrapolation")
            end
		else
			indices[i] = min(searchsortedlast(ingrid, outgrid[i]), length(ingrid) - 1)
			x0, x1 = ingrid[indices[i]], ingrid[indices[i]+1]
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
	@inbounds @simd for i in eachindex(values_new)
		values_new[i] = weights[i, 1] * values[indices[i]] + 
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
    @inline _dx(idx) = xx[idx+1] - xx[idx]
    @inline _dy(idx) = yy[idx+1] - yy[idx]

    @inbounds for i in 2:n-1
        if (_dy(i-1) == 0) || (_dy(i) == 0)
            continue
        end

        d1 = _dy(i-1)/_dx(i-1)
        d2 = _dy(i)/_dx(i)

        if (sign(d1) == sign(d2))
            # harmonic mean derivative
            alpha = (1. + _dx(i-1)/(_dx(i) + _dx(i-1))) / 3.0
            dydx[i] = 1. / ((1. - alpha)/d2 + alpha/d1)
        end
    end

    # boundary slopes from quadratic polynomial using 2 points and 1 dydx
    a = (_dy(1) - _dx(1)*dydx[2]) / (xx[2]^2 - xx[1]^2 - 2. * xx[2]*_dx(1))
    dydx[1] = dydx[2] - 2. * a *_dx(1)

    a = (_dy(n-1) - _dx(n-1)*dydx[n-1]) / (xx[n]^2 - xx[n-1]^2 - 2. * xx[n-1]*_dx(n-1))
    dydx[n] = dydx[n-1] + 2. * a *_dx(n-1)

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
        t1 = (newx[j] - xx[i]) / _dx(i)

        t2 = t1*t1
        t3 = t2*t1
        p = 3.0*t2 - 2.0*t3

        newy[j] = (1.0 - p)*yy[i] + p*yy[i+1] + _dx(i)*((t3 - 2.0*t2 + t1)*dydx[i] + (t3-t2)*dydx[i+1]) 
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
function _interpolate_grid_dim!(buf_out, buf_in, axis_ip, i::Int, n_dims::Int)
    dim_selection = if n_dims == 1
        (i == 1 ? 1 : Colon(),)
    elseif n_dims == 2
        (i == 1 ? 1 : Colon(), i == 2 ? 1 : Colon())
    elseif n_dims == 3
        (i == 1 ? 1 : Colon(), i == 2 ? 1 : Colon(), i == 3 ? 1 : Colon())
    else
        Tuple(j == i ? 1 : Colon() for j in 1:n_dims)
    end
    c = CartesianIndices(view(buf_in, dim_selection...))
    
    from_perm = permutation(fromaxis(axis_ip))
    to_perm = permutation(toaxis(axis_ip))
    m = method(axis_ip)

    if m == :linear
        w = axis_ip.weights
        ind = axis_ip.indices
        for ci in c
            idx_from = _fill_index(from_perm, ci, i)
            idx_to = _fill_index(to_perm, ci, i)
            interpolate_axis!(  
                view(buf_out, idx_to...),
                view(buf_in, idx_from...),
                w,
                ind
            )
        end
    elseif m == :pchip
        nodes_new = nodes(axis_ip.new_axis)
        nodes_old = nodes(axis_ip.old_axis)
        for ci in c
            idx_from = _fill_index(from_perm, ci, i)
            idx_to = _fill_index(to_perm, ci, i)
            pchip_mono8!(  
                view(buf_out, idx_to...),
                view(buf_in, idx_from...),
                nodes_new,
                nodes_old
            )
        end
    else
        error("Interpolation method not recognized. Pick :linear or :pchip")
    end
end

function interpolate_grid!(ip::GridInterpolation, values; return_copy=true)
    buffers = ip.buffers
    buffers[1] .= values
    n_dims = length(ip.axes_interpolation)

    for i in 1:n_dims
        _interpolate_grid_dim!(buffers[i+1], buffers[i], ip.axes_interpolation[i], i, n_dims)
    end

    return_copy ? deepcopy(buffers[end]) : buffers[end]
end    

@inline function _fill_index(to_fill, cartesian::CartesianIndex{N}, at) where N
    ntuple(Val(N+1)) do j
        if j == at
            to_fill
        elseif j < at
            cartesian[j]
        else
            cartesian[j-1]
        end
    end
end

# ============================================================================= 
# Interface 
# ============================================================================= 

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