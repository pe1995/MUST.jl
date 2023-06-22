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

    for j in axes(b.y, 3)
        for i in axes(b.x, 2)
            if !all(z .≈ @view(b.z[i, j, :]))
                return false
            end
        end
    end

    true
end