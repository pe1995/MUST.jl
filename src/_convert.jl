function snapshotBox_py(
    number;
    folder,
    optical_depth_scale=true,
    save_snapshot=true,
    add_selection=false,
    to_multi=false,
    is_box=true,
    use_mmap=false,
    convert_units=Dict(
        :dt_rt=>:t,
        :qr=>:qr,
        :pg=>:p,
        :flux=>:flux
    )
    )
    folder = if !isdir(folder)
        @in_dispatch folder
    else
        folder
    end

    nml_name = @in_dispatch splitpath(folder)[end]

    nml = StellarNamelist(nml_name*".nml")

    logg = log10.(abs.(nml.stellar_params["g_cgs"]))

    eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
    eos_sq = SquareGasEOS(@in_dispatch(eos_path))

    save_info = ispath(joinpath(folder, "teff.dat"))
    teff_func = if save_info
        teff_interpolated(joinpath(folder, "teff.dat"))
    else
        nothing
    end

    content_of_folder = glob("*/", folder)
    snapshots = list_of_snapshots(content_of_folder)
    snapshot_order = sortperm(snapshots)

    snapshots = snapshots[snapshot_order]
    snapshots_dir = content_of_folder[snapshot_order]

    @assert number in snapshots

    try
        rm.(glob("*.sr", snapshots_dir[findfirst(number .== snapshots)]))

        snap = dispatch.snapshot(number, data=folder)

        units = StandardUnits(snap)

        s = if !is_box
            Space(snap, :d, :ee, :ux, :uy, :uz, :e)
        else
            Box(snap, :d, :ee, :ux, :uy, :uz, :e, use_mmap=use_mmap)
        end

        convert!(s, units; d=:d, ee=:ee, e=:e,
                            ux=:u, uy=:u, uz=:u,
                            x=:l, y=:l, z=:l, time=:t, convert_units...)

        add_from_EOS!(s, eos_sq, :T)
        add_from_EOS!(s, eos_sq, :kr)
        add_from_EOS!(s, eos_sq, :Pg)
        add_from_EOS!(s, eos_sq, :Ne)

        if save_info
            @info "setting teff for sn $(number)"
            set!(s.parameter, teff_func, nml)
        end

        b_s = if !is_box
            b = Box(
                s,
                MUST.uniform_grid(s, Int(ceil(length(s.data[:x].^(1/3)))), :x),
                MUST.uniform_grid(s, Int(ceil(length(s.data[:x].^(1/3)))), :y),
                MUST.uniform_grid(s, Int(ceil(length(s.data[:x].^(1/3)))), :z)
            )
            s = nothing
            b
        else
            s
        end

        τ = optical_depth(b_s, opacity=:kr, density=:d)
        add!(b_s, τ, :τ_ross)

        b_s.parameter.logg = logg

        merge!(b_s.meta, Dict(:snapshot=>number, :kind=>"geo", :run=>last(splitpath(abspath(folder)))))
        box_path = save_snapshot ? save(b_s; name="box_sn$(number)", folder=folder) : nothing

        b_τ = if optical_depth_scale
            b_τ = height_scale_fast(b_s, :τ_ross)
            merge!(b_τ.meta, b_s.meta, Dict(:kind=>"tau"))
            save_snapshot && save!(box_path; tau=b_τ)

            b_τ
        else
            nothing
        end

        if to_multi
            b_s.data[:ne] = b_s.data[:Ne]
            multiBox(b_s, joinpath(folder, "m3dis_$(number)"))
        end

        rm.(glob("*.sr", snapshots_dir[findfirst(number .== snapshots)]))

        b_s, b_τ
    catch
        @warn "snapshot $(number) could not be loaded."
        nothing, nothing
    end
end

function snapshotBox_jl(
    number;
    folder,
    optical_depth_scale=true,
    to_multi=false,
    use_mmap=false,
    save_snapshot=true,
    save_folder=nothing,
    quantites=defaultQuantities,
    eos_path=nothing,
    eos_reader=(x)->_squaregaseos(x, eos_path=eos_path, inputNamelist=joinpath(folder, "params")),
    lookup_generator=nothing,
    kwargs...)
    folder = if !isdir(folder)
        @in_dispatch folder
    else
        folder
    end
    out_folder = isnothing(save_folder) ? folder : save_folder
    b_s, b_τ = begin
        b_s = @optionalTiming boxingTime MUST.Box(
            number,
            rundir=folder,
            quantities=quantites,
            mmap=use_mmap,
            eos_reader=eos_reader,
            lookup_generator=lookup_generator,
            save_snapshot=false
        )

        τ = @optionalTiming opticalDepthTime optical_depth(b_s, opacity=:kr, density=:d)
        add!(b_s, τ, :τ_ross)
        b_s.meta[:snapshot] = number
        b_s.meta[:kind] = "geo"
        b_s.meta[:run] = last(splitpath(abspath(folder)))
        if save_snapshot
            box_path = save(b_s; name="box_sn$(number)", folder=out_folder)
            try
                save_restart(box_path, folder, number)
            catch e
                @warn "restart data of snapshot $(number) could not be saved." exception=e
            end
        end

        b_τ = if optical_depth_scale
            b_τ = @optionalTiming heightScaleTime height_scale(b_s, :τ_ross, extrapolation=:constant)
            merge!(b_τ.meta, b_s.meta, Dict(:kind=>"tau"))
            save_snapshot && save!(box_path; tau=b_τ)

            b_τ
        else
            nothing
        end

        b_s, b_τ
    end

    if !isnothing(b_s)
        @optionalTiming multiTime if to_multi
            try
                snapshotsresample = MUST.gresample(
                    b_s,
                    nz=size(b_s, 3) * 2 - 1,
                    nx=size(b_s, 1),
                    ny=size(b_s, 2),
                )
                snapshotsresample.data[:ne] = snapshotsresample.data[:Ne]
                multiBox(snapshotsresample, joinpath(out_folder, "m3dis_$(number)"))
            catch
                @warn "snapshot $(number) could not be converted to M3D format."
                nothing
            end
        end
    end

    b_s, b_τ
end

"""
    snapshotBox(number; folder, optical_depth_scale=true, save_snapshot=true, add_selection=false, to_multi=true)

Convert the snapshot with number `number` from dispatch to a Space object.
In the current version, this relies on the dispatch python package.
This may change in the future.

# Examples
```julia
MUST.@import_dispatch "/path/to/dispatch2"
box, boxτ = snapshotBox(10, folder=@in_dispatch("data/mydispatchrun"))
```
"""
snapshotBox(args...; legacy=false, kwargs...) = legacy ? snapshotBox_py(args...; kwargs...) : snapshotBox_jl(args...; kwargs...)

const boxingTime = Ref(false)
const opticalDepthTime = Ref(false)
const heightScaleTime = Ref(false)
const multiTime = Ref(false)
