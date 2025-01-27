"""
    CO5BOLDBox(model_path)

Read CO5BOLD model from path. Store as `MUST.Box` type.
"""
function CO5BOLDBox(model_path)
    f = open(model_path, "r")

    time = read(f, Float32)
    nx, ny, nz = read(f, Float32), read(f, Float32), read(f, Float32)

    x = zeros(Float32, nx)
    y = zeros(Float32, ny)
    z = zeros(Float32, nz) 

    read!(f, x)
    read!(f, y)
    read!(f, z)

    ltemp = zeros(Float32, nx, ny, nz)
    lrho = zeros(Float32, nx, ny, nz)
    lpe = zeros(Float32, nx, ny, nz) 
    px = zeros(Float32, nx, ny, nz)
    py = zeros(Float32, nx, ny, nz)
    pz = zeros(Float32, nx, ny, nz) 

    read!(f, ltemp)
    read!(f, lrho)
    read!(f, lpe)
    read!(f, px)
    read!(f, py)
    read!(f, pz)

    close(f)

    xx, yy, zz = meshgrid(x, y, z)
    paras = AtmosphericParameters() 
    paras.time = time

    temp = exp.(ltemp)
    rho = exp.(lrho)

    MUST.Box(
        xx, yy, zz,
        Dict{Symbol, Any}(
            :ux  => px ./rho,
            :uy  => py ./rho,
            :uz  => pz ./rho,
            :T   => temp,
            :Pe  => pe,
            :d   => rho,
        ),
        paras
    )
end