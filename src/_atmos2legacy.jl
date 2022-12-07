"""Convert Box objects to legacy bifrost/stagger format."""
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

"""Convert legacy snapshot to Box object."""
function Box(s::StaggerSnap; units=StaggerCGS(), eos=nothing)
    # The snapshot is already ordered in the normal way
    d  = s[:r]
    ux = s[:px] ./ d
    uy = s[:py] ./ d
    uz = -1 .* s[:pz] ./ d
    
    data = Dict{Symbol, Array{Float32, 3}}( :d=>d, :pp=>s[:pp], :T=>s[:temp], 
                                            :ross=>s[:ross] ./ d, :τ=>s[:tau],
                                            :px=>s[:px],     :py=>s[:py], :pz=>-1 .*s[:pz],
                                            :ux=>ux,         :uy=>uy,     :uz=>uz,
                                            :ee=>s[:ee],     :e=>s[:e])

    xx,yy,zz = meshgrid(reverse(s[:x]), reverse(s[:y]), reverse(-1 .* s[:z]))

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