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
function snapshotBox(
    number; 
    folder, 
    optical_depth_scale=true, 
    save_snapshot=true, 
    add_selection=false, 
    to_multi=false,
    is_box=true, 
    use_mmap=false
    )
    folder = if !isdir(folder)
        #@warn "Given folder does not exist. Trying to add dispatch location"
        @in_dispatch folder
    else
        folder
    end

    # Name of the namelist of the current folder
    nml_name = @in_dispatch splitpath(folder)[end]

    # Init namelist
    nml = StellarNamelist(nml_name*".nml")

    # read logg from the namelist
    logg = log10.(abs.(nml.stellar_params["g_cgs"]))

    # Use the Squaregas EOS 
    eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
    eos_sq = SquareGasEOS(@in_dispatch(eos_path))

    save_info = ispath(joinpath(folder, "teff.dat")) 
    teff_func = if save_info
        # Read and interpolate teff data
        teff_interpolated(joinpath(folder, "teff.dat"))
    else
        nothing
    end

    content_of_folder = glob("*/", folder)
    snapshots = sort(list_of_snapshots(content_of_folder))

    @assert number in snapshots

    try
        # The dispatch snapshot object (Python)
        snap = dispatch.snapshot(number, data=folder)

        # Units for conversion to CGS
        units = StandardUnits(snap)

        # Convert its content to pure Julia
        s, add_selection = if add_selection
            try
                !is_box ? 
                    (Space(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg), true) :
                    (Box(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg, use_mmap=use_mmap), true)
            catch
                !is_box ?
                    (Space(snap, :d, :ee, :ux, :uy, :uz, :e), false) :
                    (Box(snap, :d, :ee, :ux, :uy, :uz, :e, use_mmap=use_mmap), false)
            end
        else
            !is_box ?
                (Space(snap, :d, :ee, :ux, :uy, :uz, :e), false) :
                (Box(snap, :d, :ee, :ux, :uy, :uz, :e, use_mmap=use_mmap), false)
        end

        # Apply the conversion
        if add_selection
            convert!(s, units; d=:d, ee=:ee, e=:e, 
                                qr=:qr, pg=:p,
                                ux=:u, uy=:u, uz=:u,
                                x=:l, y=:l, z=:l, time=:t)
        else
            convert!(s, units; d=:d, ee=:ee, e=:e, 
                                ux=:u, uy=:u, uz=:u,
                                x=:l, y=:l, z=:l, time=:t)
        end

        # Add additional columns already in CGS after converting
        add_from_EOS!(s, eos_sq, :T)
        add_from_EOS!(s, eos_sq, :kr)
        add_from_EOS!(s, eos_sq, :Ne)

        # Check if Teff should be added
        if save_info
            @info "setting teff for sn $(number)"
            set!(s.parameter, teff_func, nml)
        end

        # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
        b_s = if !is_box 
            b = Box(s) 
            s = nothing
            b
        else
            s
        end

        # Add the optical depth
        τ = optical_depth(b_s, opacity=:kr, density=:d)
        add!(b_s, τ, :τ_ross)

        # add logg
        b_s.parameter.logg = logg

        # First save
        save_snapshot && save(b_s; name="box_sn$(number)", folder=folder)

        # Convert the height scale from cm to optical depth
        b_τ = if optical_depth_scale
            b_τ = height_scale_fast(b_s, :τ_ross)
            save_snapshot && save(b_τ; name="box_tau_sn$(number)", folder=folder)

            b_τ 
        else
            nothing
        end

        # convert to Multi format and save in the same folder, keep the number of the snapshot!
        if to_multi
            b_s.data[:ne] = b_s.data[:Ne]
            multiBox(b_s, joinpath(folder, "m3dis_$(number)"))
        end

        b_s, b_τ
    catch
        @warn "snapshot $(number) could not be loaded."
        nothing, nothing
    end
end