
function MURaMBox(model_path)
    ds = NCDataset(model_path, "r")
    U  = ds["U"][:,:,:]
    W  = ds["W"][:,:,:]
    T  = ds["T"][:,:,:]
    P  = ds["P"][:,:,:]
    By = ds["By"][:,:,:]
    V  = ds["V"][:,:,:]
    Bx = ds["Bx"][:,:,:]
    R  = ds["R"][:,:,:]
    E  = haskey(ds, "E") ? ds["E"][:,:,:] : nothing
    Bz = ds["Bz"][:,:,:]
    tau = haskey(ds, "tau") ? ds["tau"][:,:,:] : nothing

    s = size(U);

    dx = ds.attrib["dx"]
    dy = ds.attrib["dy"]
    dz = ds.attrib["dz"]

    time  = ds.attrib["time"]
    paras = AtmosphericParameters()
    paras.time = time;

    x = range(-(s[1]-1) * dx/2,   step=dx, length=s[1])
    y = range(-(s[2]-1) * dy/2,   step=dy, length=s[2])
    z = -range(-(s[3]-1) * dz*1/2, step=dz, length=s[3])

    xx, yy, zz = meshgrid(x, y, z)

    data = Dict{Symbol, Array{eltype(U), 3}}(
        :ux  => U ,
        :uz  => W ,
        :T   => T ,
        :P   => P ,
        :By  => By,
        :uy  => V ,
        :Bx  => Bx,
        :d   => R ,
        :Bz  => Bz
    )

    if !isnothing(E)
        data[:E] = E
    end

    if !isnothing(tau)
        data[:tau] = tau
    end

    Box(xx, yy, zz, data, paras)
end

function MURaMBox_optical_depth(b, opacity)
    ## lookup the rosseland opacity
    b.data[:kr] = opacity

    ## optical depth form this
    τ = optical_depth(b, opacity=:kr, density=:d)
    b.data[:τ_ross] = τ

    ## Interpolate cube
    height_scale(b, :τ_ross)
end
