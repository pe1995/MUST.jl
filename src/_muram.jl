
"""
    MURaMBox(path)

Read a Muram 3D data cube from path in the NetCDF format and convert it to a Box.
"""
function MURaMBox(model_path)
    U  = ncread(model_path, "U")
    W  = ncread(model_path, "W")
    T  = ncread(model_path, "T")
    P  = ncread(model_path, "P")
    By = ncread(model_path, "By")
    V  = ncread(model_path, "V")
    Bx = ncread(model_path, "Bx")
    R  = ncread(model_path, "R")
    E  = ncread(model_path, "E")
    Bz = ncread(model_path, "Bz")

    s = size(U);

    dx = ncgetatt(model_path, "Global", "dx")
    dy = ncgetatt(model_path, "Global", "dy")
    dz = ncgetatt(model_path, "Global", "dz")

    time  = ncgetatt(model_path, "Global", "time")
    paras = AtmosphericParameters() 
    paras.time = time;

    x = range(-(s[1]-1) * dx/2,   step=dx, length=s[1])
    y = range(-(s[2]-1) * dy/2,   step=dy, length=s[2])
    z = -range(-(s[3]-1) * dz*1/2, step=dz, length=s[3])


    xx, yy, zz = meshgrid(x, y, z)


    data = Dict{Symbol, Array{eltype(U), 3}}(   :ux  => U ,
                                                :uz  => W ,
                                                :T   => T ,
                                                :P   => P ,
                                                :By  => By,
                                                :uy  => V ,
                                                :Bx  => Bx,
                                                :d   => R ,
                                                :E   => E ,
                                                :Bz  => Bz )

    Box(xx, yy, zz, data, paras)
end


"""
    MURaMBox(path, eos)

Read a Muram 3D data cube from path in the NetCDF format and convert it to a Box.
Compute the rosseland optical depth scale and interpolate the cube to it. Return both.
Can either be computed with the MUST eos interface or with the more convenient TSO.

(exp.(lookup(eos, :lnRoss, log.(b[:d]), log.(b[:T]))) for TSO.jl)
"""
function MURaMBox(model_path, opacity)
    ## The default box without rosseland
    b = MURaMBox(model_path)

    ## lookup the rosseland opacity
    b.data[:kr] = opacity

    ## optical depth form this
    τ = MUST.optical_depth(b, opacity=:kr, density=:d)
    b.data[:τ_ross] = τ

    ## Interpolate cube
    b_t = MUST.height_scale(b, :τ_ross)

    return b, b_t
end