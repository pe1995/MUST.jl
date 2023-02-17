"""
Convert Box objects to legacy bifrost/stagger format.
"""
function legacyBox(b::Box, name::String)
    # For legacy code, we need .snap, .aux and .mesh files
    # We will only save quantities that are needed for e.g. Tabgen, the rest will be filled
    # with 0 values.
    snap_file = "$(name).snap"
    aux_file  = "$(name).aux"
    mesh_file = "$(name).msh"

    _write_mesh(b, mesh_file)
    _write_snap(b, snap_file)
    _write_aux(b, aux_file)

    nothing
end

"""
Create a Space object from a Legacy Stagger Snapshot.

Parameters
----------
snapshot: StaggerSnap
    A Stagger snapshot.
quantites: Symbols
    Variables to extract from the snapshot

Returns
-------
MUST.Space type object

"""
function Space(snap::MUST.StaggerSnap, quantities::Symbol...)
    qs = Dict{Symbol,Vector{<:Union{Float32, Float64, Int16, Int32, Int64}}}(q=>Float32[] for q in quantities)
    qs[:x] = Float32[]; qs[:y] = Float32[]; qs[:z] = Float32[]
    qs[:i_patch] = Int[]
    patch_dimensions = zeros(Int, (1,3))
    patch_dimensions[1,1] = length(snap[:x])
    patch_dimensions[1,2] = length(snap[:y])
    patch_dimensions[1,3] = length(snap[:z])

    x, y, z = snap[:x], snap[:y], snap[:z]
    nx, ny, nz = length(x), length(y), length(z)

    qm = zeros(Float32, (nx, ny, nz))
    for (i,q) in enumerate(quantities)
        qm .= snap[q]
        qs[q] = zeros(Float32, nx*ny*nz)
        if i==1
            coords = zeros(Float32, (nx*ny*nz,3))
        end
        c=1
        @inbounds for k in 1:nz
        @inbounds for j in 1:ny
        @inbounds for n in 1:nx
            qs[q][c] = qm[n,j,k]
            if i==1
                coords[c,1] = x[n]
                coords[c,2] = y[j]
                coords[c,3] = z[k]
            end
            c+=1
        end
        end
        end
        if i==1
            append!(qs[:x],coords[:,1])
            append!(qs[:y],coords[:,2])
            append!(qs[:z],coords[:,3])
            append!(qs[:i_patch], Int[1 for _ in 1:size(coords,1)])
        end
    end
    Space(qs, patch_dimensions, AtmosphericParameters(-99.0, -99.0, -99.0, Dict{Symbol, Float64}()))
end


"""
Convert legacy snapshot to Box object.
"""
function Box(s::StaggerSnap; units=StaggerCGS(), eos=nothing, gridded=true)
    stagger_box = if !gridded
        @warn "Data is assumed to be unstructured! Interpolation will be very slow for the whole cube."
        ## First create a space object, because the snapshot is not spaced equally
        stagger_space = Space(s, :r, :px, :py, :pz, :temp, :pp, :ross, :tau, :ee, :e)

        ## Now we interpolate this space to a box object with equal size 
        new_x = uniform_grid(stagger_space, length(s[:x]), :x)
        new_y = uniform_grid(stagger_space, length(s[:y]), :y)
        new_z = uniform_grid(stagger_space, length(s[:z]), :z)
        Box(stagger_space, new_x, new_y, new_z)
    else
        uniform!(s)
        new_x = s[:x]
        new_y = s[:y]
        new_z = s[:z]
        s
    end

    @assert all(new_x[2:end] .- new_x[1:end-1] .≈ new_x[2] - new_x[1])
    @assert all(new_y[2:end] .- new_y[1:end-1] .≈ new_y[2] - new_y[1])
    @assert all(new_z[2:end] .- new_z[1:end-1] .≈ new_z[2] - new_z[1])

    ## The snapshot is already ordered in the normal way
    d  = stagger_box[:r]
    ux = stagger_box[:px] ./ d
    uy = stagger_box[:py] ./ d
    uz = -1 .* stagger_box[:pz] ./ d
    
    data = Dict{Symbol, Array{Float32, 3}}( :d=>d, :pp=>stagger_box[:pp],   :T=>stagger_box[:temp], 
                                            :ross=>stagger_box[:ross] ./ d, :τ=>stagger_box[:tau],
                                            :px=>stagger_box[:px],          :py=>stagger_box[:py], :pz=>-1 .*stagger_box[:pz],
                                            :ux=>ux,                        :uy=>uy,               :uz=>uz,
                                            :ee=>stagger_box[:ee],          :e=>stagger_box[:e])

    xx,yy,zz = meshgrid(reverse(new_x), reverse(new_y), reverse(-1 .* new_z))

    for key in keys(data)
        data[key] = reverse(data[key])
    end

    b = Box(xx, yy, zz, data, AtmosphericParameters())

    # A conversion to cgs has to be done
    # Apply the conversion
    convert!(b, units;  d=:d, pp=:p, ross=:k,
                        px=:pm, py=:pm, pz=:pm,
                        ux=:u, uy=:u, uz=:u, 
                        ee=:ee, e=:e,
                        x=:l, y=:l, z=:l)

    if !isnothing(eos)
        @info "Recomputing Energy from EoS."
        # We compute the energy cube, which is needed to get the gas pressure from the EOS
        e_cube_stagger = MUST.bisect_energy(eos, size(b[:d]); d=b[:d], T=b[:T])
        b.data[:ee] = e_cube_stagger
        b.data[:e]  = e_cube_stagger .* b[:d]
        
        pg = MUST.lookup(eos, :Pg, b[:d], e_cube_stagger)
        MUST.add!(b, pg, :Pg);


        MUST.add_from_EOS!(b, eos, :kr)
        τ  = MUST.optical_depth(b, opacity=:kr, density=:d)
        MUST.add!(b, τ,  :τ_ross)
    end

    τ2 = MUST.optical_depth(b, opacity=:ross, density=:d)
    MUST.add!(b, τ2, :τ_ross_stag)

    b
end

"""
    uniform!(snap)

Interpolate a Stagger snapshot column wise to new equidistant grid. It assumes
that the data is gridded and already unifrom in x and y.
It will interpolate the cube column by column in z and keep the dimensions.
"""
function uniform!(s::StaggerSnap)
    T     = eltype(valtype(s.data))
    new_z = range(T(minimum(s[:z])), T(maximum(s[:z])), length=length(s[:z])) |> collect
    old_z = s[:z]
    v     = similar(old_z)

    @show size(v) size(old_z)

    for j in eachindex(s[:y])
        for i in eachindex(s[:x])
            for q in keys(s.data)
                v .= view(s.data[q], i, j, :)
                s.data[q][i, j, :] .= linear_interpolation(old_z, v).(new_z)
            end
        end
    end

    s.mesh.y .= new_z 
end

function _write_aux(b, name)
    n = prod(size(b.x))
    f = FortranFile(name, "w", access="direct", recl=n*4)

    units = StaggerCGS()

    tmp = reverse(Base.convert.(Float32, b[:pg] ./ units.p))
    FortranFiles.write(f, tmp, rec=1)   

    tmp = reverse(Base.convert.(Float32, b[:T]))
    FortranFiles.write(f, tmp, rec=2);

    close(f)

    nothing
end

function _write_mesh(b, name)
    mx, my, mz = Base.convert.(Int32, size(b.x)) # Size of the Box
    dummy_x = zeros(Float32, mx)              # Fill not needed arrays with 0
    dummy_y = zeros(Float32, my)
    dummy_z = zeros(Float32, mz)

    units = StaggerCGS()

    open(name, "w") do io
        println(io,@sprintf "%d" mx)
        println(io,join([@sprintf "%e" Base.convert(Float32, x/units.l) for x in axis(b, :x)], " "))
        println(io,join([@sprintf "%e" x for x in dummy_x], " "))
        println(io,join([@sprintf "%e" x for x in dummy_x], " "))
        println(io,join([@sprintf "%e" x for x in dummy_x], " "))
        println(io,@sprintf "%d" my)
        println(io,join([@sprintf "%e" Base.convert(Float32, x/units.l) for x in axis(b, :y)], " "))
        println(io,join([@sprintf "%e" x for x in dummy_y], " "))
        println(io,join([@sprintf "%e" x for x in dummy_y], " "))
        println(io,join([@sprintf "%e" x for x in dummy_y], " "))
        println(io,@sprintf "%d"  mz)
        println(io,join([@sprintf "%e" Base.convert(Float32, x/units.l) for x in axis(b, :z)], " "))
        println(io,join([@sprintf "%e" x for x in dummy_z], " "))
        println(io,join([@sprintf "%e" x for x in dummy_z], " "))
        println(io,join([@sprintf "%e" x for x in dummy_z], " "))
    end
end

function _write_snap(b, name)
    n = prod(size(b.x))
    f = FortranFile(name, "w", access="direct", recl=n*4);
    
    units = StaggerCGS()

    tmp = reverse(Base.convert.(Float32, b[:d] ./ units.d))  # r
    FortranFiles.write(f, Base.convert.(Float32, tmp), rec=1)  

    tmp = reverse(Base.convert.(Float32, b[:ux] ./ units.u .* b[:d] ./ units.d)) # px
    FortranFiles.write(f, Base.convert.(Float32, tmp), rec=2)
    
    tmp = reverse(Base.convert.(Float32, b[:uy] ./ units.u .* b[:d] ./ units.d)) # py
    FortranFiles.write(f, Base.convert.(Float32, tmp), rec=3)
    
    tmp = reverse(Base.convert.(Float32, b[:uz] ./ units.u .* b[:d] ./ units.d)) # pz
    FortranFiles.write(f, Base.convert.(Float32, tmp), rec=4)
    
    tmp = reverse(Base.convert.(Float32, b[:e] ./ units.e)) # e
    FortranFiles.write(f, Base.convert.(Float32, tmp), rec=5)
    
    close(f)

    nothing
end