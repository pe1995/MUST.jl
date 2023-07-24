#= Constructing grids from existing Boxes =#

"""
    Grid(box::Box)

Create Axis and Grid objects from a Box object.
The correct axis (1D) are picked if the box is regular,
i.e. the MUST.axis(:x), MUST.axis(:y) and MUST.axis(:z)
function return 1D arrays of the cube axis.
"""
Grid(b::Box) = begin
    @assert check_uniform(b)

    x, y, z = axis(b, :x)[:], axis(b, :y)[:], axis(b, :z)[:]

    # Convert the axis to a higher precision
    x = Base.convert.(Float64, x)
    y = Base.convert.(Float64, y)
    z = Base.convert.(Float64, z)


    axes = [RegularBoxAxis(x), RegularBoxAxis(y), RegularBoxAxis(z)]
    RegularBoxGrid(axes)
end

Grid(x, y, z) = RegularBoxGrid(
    RegularBoxAxis(x), 
    RegularBoxAxis(y), 
    RegularBoxAxis(z)
)




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

function gresample(b::Box; nx=size(b, 1), ny=size(b, 2), nz=size(b, 3))
	if (nx==size(b, 1)) & (ny==size(b, 2)) & (nz==size(b, 3))
		@warn "Size of new box = size of old box."
		deepcopy(b)
	end

	# The coordinate grid of the input box
	grid = Grid(b)

	# build the new axis
	x_new = scale_axis(axis(b, :x), N=nx)
	y_new = scale_axis(axis(b, :y), N=ny)
	z_new = scale_axis(axis(b, :z), N=nz)
	target_grid = Grid(x_new, y_new, z_new)

	# interpolate 
	ip = ginterpolate(grid, target_grid)

	# compute the interpolate quantities
	data_new = Dict{Symbol,Array{<:Union{Float32, Float64, Int16, Int32, Int64},3}}()

    TN = eltype(b.x)
	for f in keys(b.data)
		d = if all(b.data[f] .> 0.0) & (eltype(b.data[f]) <: AbstractFloat)
			log.(b.data[f])
		else
			b.data[f]
		end
		
        # interpolate
		data_new[f] = Base.convert.(TN, gevaluate!(ip, d))

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

	Box(xx, yy, zz, data_new, deepcopy(b.parameter))
end