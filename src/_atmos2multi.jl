function multiBox(b::Box, eos::AbstractEOS, name::String; downsample_xy=1, downsamlpe_z=1)
    mesh_path = name*".mesh"
    atmos_path = name*".atmos"

    _write_mesh_multi(b, mesh_path; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
    _write_atmos_multi(b, atmos_path, eos; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
end

function _write_mesh_multi(b, path; downsample_xy=1, downsamlpe_z=1)
    x = axis(b, :x)[1:downsample_xy:end]
    y = axis(b, :y)[1:downsample_xy:end]
    z = axis(b, :z)[1:downsamlpe_z:end]
    n = (length(x), length(y), length(z))
    open(path, "w") do f
        write(f, join(n, " ")*"\n")
        write(f, join(x, " ")*"\n")
        write(f, join(y, " ")*"\n")
        write(f, join(z, " ")*"\n")
    end
end

function _write_atmos_multi(b, path, eos; downsample_xy=1, downsamlpe_z=1)
    s = size(b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end])
    res::Array{Float32, 4} = zeros(Float32, s..., 6)

    res[:, :, :, 1] .= Base.convert.(Float32, lookup(eos, :Ne, b[:d][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end], b[:ee][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]))
    res[:, :, :, 2] .= b[:T][ 1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]
    res[:, :, :, 3] .= b[:ux][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 4] .= b[:uy][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 5] .= b[:uz][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 6] .= b[:d][ 1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]
    
    reverse!(res, dims=3)

    f = open(path, "w")
    write(f, res)
    close(f)
end


## Reading

function multiBox(name::String)
    mesh_path = name*".mesh"
    atmos_path = name*".atmos"

    n, x, y, z = _read_mesh_multi(mesh_path)
    ne, T, ux, uy, uz, d = _read_atmos_multi(atmos_path, n)

    data = Dict{Symbol, typeof(ne)}(:Ne=>ne, :T=>T, :ux=>ux, :uy=>uy, :uz=>uz, :d=>d)

    xx, yy, zz = meshgrid(x, y, z)

    Box(xx, yy, zz, data, AtmosphericParameters())
end

function _read_mesh_multi(path)
    open(path, "r") do f
        n = parse.(Int, split(strip(readline(f)), " "))
        x = parse.(Float32, split(strip(readline(f)), " "))
        y = parse.(Float32, split(strip(readline(f)), " "))
        z = parse.(Float32, split(strip(readline(f)), " "))

        n, x, y, reverse(z)
    end
end

function _read_atmos_multi(path, n)    
    r  = zeros(Float32, n..., 6)
    ne = zeros(Float32, n...)
    T  = zeros(Float32, n...)
    px = zeros(Float32, n...)
    py = zeros(Float32, n...)
    pz = zeros(Float32, n...)
    d  = zeros(Float32, n...)

    reverse(read!(path, r), dims=3)

    ne .= r[:, :, :, 1]
    T  .= r[:, :, :, 2]
    px .= r[:, :, :, 3]
    py .= r[:, :, :, 4]
    pz .= r[:, :, :, 5]
    d  .= r[:, :, :, 6]

    ne, T, px, py, pz, d
end


