function multiBox(b::Box, eos::AbstractEOS, name::String; downsample_xy=1, downsamlpe_z=1)
    mesh_path = name*".mesh"
    atmos_path = name*".atmos"

    _write_mesh_multi(b, mesh_path; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
    _write_atmos_multi(b, atmos_path, eos; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
end

function multiBox(b::Box, name::String; downsample_xy=1, downsamlpe_z=1)
    mesh_path = name*".mesh"
    atmos_path = name*".atmos"

    _write_mesh_multi(b, mesh_path; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
    _write_atmos_multi(b, atmos_path, nothing; downsample_xy=downsample_xy, downsamlpe_z=downsamlpe_z)
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

function _write_atmos_multi(b, path, eos=nothing; downsample_xy=1, downsamlpe_z=1)
    s = size(b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end])
    res::Array{Float32, 4} = zeros(Float32, s..., 6)

    ee = similar(res, size(b[:d])[1:3]...)
    ne = similar(res, size(b[:d])[1:3]...)

    if !haskey(b.data, :ne)
        @info "recomputing Ne"
        for k in axes(res, 3)
            for j in axes(res,2 )
                for i in axes(res, 1)
                    ee[i, j, k] = MUST.bisect(eos, ee=limits(eos)[3:4], d=b[:d][i, j, k], T=b[:T][i, j, k])
                    ne[i, j, k] = Base.convert(Float32, lookup(eos, :Ne, b[:d][i, j, k], ee[i, j, k]))
                    if ne[i, j, k] == 0.0
                        @info "Ne is 0, i,j,k,d,T,ee: $(i),$(j),$(k),$(b[:d][i, j, k]),$(b[:T][i, j, k]),$(ee[i, j, k])"
                    end
                end
            end
        end
    else
        ne .= b[:ne]
    end

    res[:, :, :, 1] .= ne[ 1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]
    res[:, :, :, 2] .= b[:T][ 1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]
    res[:, :, :, 3] .= b[:ux][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 4] .= b[:uy][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 5] .= b[:uz][1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end] #.* b[:d][1:downsample_xy:end,1:downsample_xy:end,1:downsamlpe_z:end]
    res[:, :, :, 6] .= b[:d][ 1:downsample_xy:end, 1:downsample_xy:end, 1:downsamlpe_z:end]
    
    #@show minimum(res[:,:,:,1])

    #reverse!(res, dims=3)

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

    data = Dict{Symbol, Any}(:Ne=>ne, :T=>T, :ux=>ux, :uy=>uy, :uz=>uz, :d=>d)

    xx, yy, zz = meshgrid(x, y, z)

    Box(xx, yy, zz, data, AtmosphericParameters())
end

function _read_mesh_multi(path)
    open(path, "r") do f
        n = parse.(Int, split(strip(readline(f)), " "))
        x = parse.(Float32, split(strip(readline(f)), " "))
        y = parse.(Float32, split(strip(readline(f)), " "))
        z = parse.(Float32, split(strip(readline(f)), " "))

        n, x, y, z
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

    read!(path, r)

    ne .= r[:, :, :, 1]
    T  .= r[:, :, :, 2]
    px .= r[:, :, :, 3]
    py .= r[:, :, :, 4]
    pz .= r[:, :, :, 5]
    d  .= r[:, :, :, 6]

    ne, T, px, py, pz, d
end


