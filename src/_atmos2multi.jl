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


#= Reading =#

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









#= converting Box to average 3D model for M3D =#

"""
    save_text_m3d(z, T, ρ, f_new; header=nothing, vmic=zeros(length(z)))

Save the model in a format readable by M3D.
"""
save_text_m3d(f_new, z, ρ, T; header=nothing, vmic=zeros(length(z)), pe=zeros(length(z))) = begin
    open(f_new, "w") do f
        h = isnothing(header) ?  "Average 3D" : header
		write(f, h*"\n")
		write(f, "$(length(z))\n")
		
		for i in eachindex(z)
			line = @sprintf "%15.6E  %10.1f  %14.4E  %14.4E  %5.2f\n" z[i] T[i] pe[i] ρ[i] vmic[i]
			write(f, line)
		end
	end
end

"""
    save_text_m1d(τ, T, ne, f_new; logg, header=nothing, vmic=zeros(length(z)))

Save the model in a format readable by M3D.
"""
save_text_m1d(f_new, τ, T, ne; logg, header=nothing, vmic=zeros(length(τ)), v=zeros(length(τ))) = begin
    open(f_new, "w") do f
        h = isnothing(header) ? "Average 3D" : header
		write(f, "  "*h*"\n")
		write(f, "  TAU\n")
		write(f, "* surface gravity log(g)\n")
		write(f, "  $(logg)\n")
		write(f, "* Number of depth points\n")
		write(f, "  $(length(τ))\n")
		write(f, "* log10 tau (500nm), temperature [K], electron density [g cm-3], v (always 0), microturbulence [km s-1]\n")
		
		for i in eachindex(τ)
			line = @sprintf "  %15.6E  %15.6E  %15.6E  %.1f  %15.6E\n" τ[i] T[i] ne[i] v[i] vmic[i]
			write(f, line)
		end
	end
end

"""
    save_text_m1d(τ, T, ne, f_new; logg, header=nothing, vmic=zeros(length(z)))

Save the model in a format readable by M3D.
"""
save_text_m1d_dscale(f_new, τ; header=nothing, kwargs...) = begin
    open(f_new, "w") do f
        h = isnothing(header) ? "Average 3D" : header
		write(f, "  "*h*"\n")
		write(f, "  TAU\n")
		write(f, "  $(length(τ))     $(first(τ))\n")		
		for i in eachindex(τ)
			line = @sprintf "  %15.6E\n" τ[i]
			write(f, line)
		end
	end
end


"""
    save_average_m3d(b::Box, scale, f_new; header=nothing)

Save the average geometical model in a format readable by M3D.
"""
save_average_m3d(b::Box, scale, f_new; kwargs...) = begin
    z, ρ, T, ne, pe, vmic = average!(b, scale=scale)

    # convert vmic to km/s
    vmic = vmic ./ 1e5
    save_text_m3d(f_new, z, ρ, T; vmic=vmic, pe=pe, kwargs...) 
end

save_tau500_average_m3d(b::Box, f_new; kwargs...) = save_average_m3d(b, :log10τ500, f_new; kwargs...)
save_tauross_average_m3d(b::Box, f_new; kwargs...) = save_average_m3d(b, :log10τ_ross, f_new; kwargs...)
save_geo_average_m3d(b::Box, f_new; kwargs...) = save_average_m3d(b, :z, f_new; kwargs...)


"""
    save_geo_average_m3d(b::Box, f_new; header=nothing)

Save the average geometical model in a format readable by M3D.
"""
save_tau_average_m1d(b::Box, f_new; kwargs...) = begin
    logτ, ρ, T, ne, pe, vmic = average!(b, scale=:log10τ500)

    filename = basename(f_new)
    dir = dirname(f_new)
    f_new_atmos = joinpath(dir, "atmos.$(filename)")
    f_new_dscale = joinpath(dir, "dscale.$(filename)")

    # convert vmic to km/s
    vmic = vmic ./ 1e5
    save_text_m1d(f_new_atmos, logτ, T, ne; vmic=vmic, kwargs...) 
    save_text_m1d_dscale(f_new_dscale, logτ; kwargs...) 
end




