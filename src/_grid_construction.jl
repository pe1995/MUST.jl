#= Constructing grids from existing Boxes =#

"""
    Grid(box::Box)

Create Axis and Grid objects from a Box object.
The correct axis (1D) are picked if the box is regular,
i.e. the MUST.axis(:x), MUST.axis(:y) and MUST.axis(:z)
function return 1D arrays of the cube axis.
"""
Grid(b::Box; order=[1, 2, 3]) = begin
    @assert check_uniform(b)

    x, y, z = axis(b, :x)[:], axis(b, :y)[:], axis(b, :z)[:]

    # Convert the axis to a higher precision
    #x = Base.convert.(Float64, x)
    #y = Base.convert.(Float64, y)
    #z = Base.convert.(Float64, z)

    axes = [Axis(x), Axis(y), Axis(z)][order]

    RegularBoxGrid(axes)
end

Grid(b::Box, zscale::Symbol; order=[1, 2, 3], scale=log) = begin
    #@assert check_uniform(b)

    x, y, z = axis(b, :x)[:], axis(b, :y)[:], scale.(axis(b, zscale, 3)[:])

    # Convert the axis to a higher precision
    #x = Base.convert.(Float64, x)
    #y = Base.convert.(Float64, y)
    #z = Base.convert.(Float64, z)

    axes = [Axis(x), Axis(y), Axis(z)][order]

    RegularBoxGrid(axes)
end

#Grid(x, y, z) = RegularBoxGrid(
#    RegularBoxAxis(x), 
#    RegularBoxAxis(y), 
#    RegularBoxAxis(z)
#)


#= Utilities =#

check_uniform(b::Box) = begin
    x, y, z = axis(b, :x), axis(b, :y), axis(b, :z)
    x = Base.convert.(Float64, x)
    y = Base.convert.(Float64, y)
    z = Base.convert.(Float64, z)
  
    for k in axes(b.z, 3)
        for j in axes(b.y, 2)
            if !all(x .≈ @view(b.x[:, j, k]))
                return false
            end
        end
    end

    for k in axes(b.z, 3)
        for i in axes(b.y, 2)
            if !all(y .≈ @view(b.y[i, :, k]))
                return false
            end
        end
    end

    for j in axes(b.y, 2)
        for i in axes(b.x, 1)
            if !all(z .≈ @view(b.z[i, j, :]))
                return false
            end
        end
    end

    true
end

"""
    scale_axis(axis, factor)

Scale the number of points on an axis by the given factor or to a specific number.
"""
function scale_axis(axis; factor=nothing, N=nothing)
    if isnothing(factor) & isnothing(N)
        error("N or factor required.")
    end

    len = if isnothing(N)
        Int(ceil(length(axis).*factor))
    else
        N
    end

	axis_new = range(first(axis), last(axis), length=len)
	Base.convert.(eltype(axis), axis_new)
end





#= Scaling a box in resolution =#

function gresample(b::Box; 
        nx=size(b, 1), ny=size(b, 2), nz=size(b, 3), 
        method=:linear, order=[1, 2, 3], 
        dont_log=[:T, :ux, :uy, :uz, :ee, :e])
	if (nx==size(b, 1)) & (ny==size(b, 2)) & (nz==size(b, 3))
		@warn "Size of new box = size of old box."
	end

    b = deepcopy(b)
    natural_order = sortperm(order)

	# The coordinate grid of the input box
	grid = Grid(b, order=order)

	# build the new axis
	x_new = isnothing(nx) ? deepcopy(nodes(grid.axes[1])) : scale_axis(axis(b, :x), N=nx)
	y_new = isnothing(ny) ? deepcopy(nodes(grid.axes[2])) : scale_axis(axis(b, :y), N=ny)
	z_new = isnothing(nz) ? deepcopy(nodes(grid.axes[3])) : scale_axis(axis(b, :z), N=nz)
	target_grid = Grid([x_new, y_new, z_new][order]...)

    #@show nodes(axis(target_grid)[3])
    #@show nodes(axis(grid)[3])


	# interpolate 
	ip = ginterpolate(grid, target_grid, method=method)

	# compute the interpolate quantities
	data_new = Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}()
    d_tmp = permutedims(first(values(b.data)), order)

    TN = eltype(b.x)

	for f in keys(b.data)
		d_tmp .= if all(b.data[f] .> 0.0) & (eltype(b.data[f]) <: AbstractFloat) & !(f in dont_log)
            #@info "logging $(f)"
			permutedims(log.(b.data[f]), order)
		else
			permutedims(b.data[f], order)
		end
		
        # interpolate
		data_new[f] = Base.convert.(TN, permutedims(gevaluate!(ip, d_tmp), natural_order))

        # apply exp again if needed
		data_new[f] .= if all(b.data[f] .> 0.0) & 
							(eltype(b.data[f]) <: AbstractFloat) & !(f in dont_log)
			exp.(data_new[f])
		else
			data_new[f]
		end
	end

	xx, yy, zz = meshgrid(nodes(target_grid)...)
    xx = Base.convert.(TN, xx)
    yy = Base.convert.(TN, yy)
    zz = Base.convert.(TN, zz)

	Box(xx, yy, zz, data_new, deepcopy(b.parameter))
end


#=
function gresample(b::Box, zscale::Symbol; 
                    nx=size(b, 1), ny=size(b, 2), nz=size(b, 3), 
                    method=:linear, order=[1, 2, 3], scale=:log)
	if (nx==size(b, 1)) & (ny==size(b, 2)) & (nz==size(b, 3))
		@warn "Size of new box = size of old box."
	end

    b = deepcopy(b)

    fs = scale == :log ? log : identity
    fs_back = scale == :log ? exp : identity

    natural_order = sortperm(order)

	# The coordinate grid of the input box
	grid = Grid(b, zscale, order=order, scale=fs)

	# build the new axis
	x_new = scale_axis(axis(b, :x), N=nx)
	y_new = scale_axis(axis(b, :y), N=ny)
	z_new = scale_axis(fs.(axis(b, zscale, 3)), N=nz)
	target_grid = Grid([x_new, y_new, z_new][order]...)

    @show grid.axes[3].nodes
    @show target_grid.axes[3].nodes


	# interpolate 
	ip = ginterpolate(grid, target_grid, method=method)

	# compute the interpolate quantities
	data_new = Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}()
    d_tmp = permutedims(first(values(b.data)), order)

    TN = eltype(b.x)
	for f in keys(b.data)
        if f == zscale
            continue
        end

		d_tmp .= if all(b.data[f] .> 0.0) & (eltype(b.data[f]) <: AbstractFloat)
			permutedims(log.(b.data[f]), order)
		else
			permutedims(b.data[f], order)
		end
		
        # interpolate
		data_new[f] = Base.convert.(TN, permutedims(gevaluate!(ip, d_tmp), natural_order))

        # apply exp again if needed
		data_new[f] .= if all(b.data[f] .> 0.0) & 
							(eltype(b.data[f]) <: AbstractFloat)
			exp.(data_new[f])
		else
			data_new[f]
		end
	end

	xx, yy, zz = meshgrid(x_new, y_new, z_new)
    xx = Base.convert.(TN, xx)
    yy = Base.convert.(TN, yy)
    zz = Base.convert.(TN, zz)

    data_new[zscale] = fs_back.(zz)
    d_tmp .= permutedims(b.z, order)
    z_real = Base.convert.(TN, permutedims(gevaluate!(ip, d_tmp), natural_order))

	Box(xx, yy, z_real, data_new, deepcopy(b.parameter))
end
=#