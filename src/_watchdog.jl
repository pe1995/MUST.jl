"""
    WatchDog

Monitor progress of completed or running DISPATCH simulations.
Save statistics of snapshots available during computation. Administer
time dependent output during runtime.
"""
mutable struct WatchDog
    name
    folder
    snapshots
    snapshotsCompleted
    refinement
    currentSnapshot
    functions
end

WatchDog(name; folder=@in_dispatch("data/"), functions...) = WatchDog(name, joinpath(folder, name), [], [], 0, 0, functions)


"""
    monitor(w::WatchDog; timeout=2*60*60, check_every=5, delay=0, snapshotbuffer=1, save_box=false, reverse=false, keeplast=-1, batch=10)

Monitor the progress of the linked simulation. Compute statistics and save 
depth profiles. Check for new snapshots after `check_every` seconds.
Cancel the monitoring if `timeout` seconds have passed without
finding a new snapshot.
"""
function monitor(w::WatchDog; timeout=2*60*60, check_every=5, delay=0, snapshotbuffer=1, save_box=false, reverse=false, keeplast=-1, batch=10, spectra=-1, onlyrefinement=false)
    time_start = time()
    time_current = time()
    time_passed_since(t_ref) = time() - t_ref

    # keep looping until timeout
    while true
        # list snapshots
        updatesnaps!(w)

        # check how many snapshots we have comleted
        # delete the snapshots if they are older than keeplast
        if keeplast > 0
            n_snapsCompleted = length(w.snapshotsCompleted)
            if n_snapsCompleted > keeplast
                snaps2remove = w.snapshotsCompleted[1:end-keeplast]
                for (i, snap2remove) in enumerate(snaps2remove)
                   deleteSnapshot(w, snap2remove)
                end
            end
        end

        # check if there is a new snapshot
        enum = enumerate(w.snapshots)
        enum = reverse ? Iterators.reverse(enum) : enum
        ibatch = 0
        for (i, snapf) in enum
            # for safety reasons, we dont convert the last N snaps
            # this avoids them being read while not written properly
            if i >= length(w.snapshots)-snapshotbuffer
                continue
            end

            if !(snapf in w.snapshotsCompleted)
                # give the snapshot some time to be written (optional)
                sleep(delay)

                # Check how detailed the output should be
                w.currentSnapshot = snapf
                w.refinement = 0
                if spectra > 0
                    # Increase refinemet if spectra are requested
                    w.refinement = (((snapf % spectra) == 0) && snapf>0) ? 1 : 0
                end

                # save the snapshot to disk if we have a refinement level
                save_refinement = if w.refinement > 0
                    true
                else
                    save_box
                end

                # skip if only a refinement is requested
                if onlyrefinement
                    if w.refinement == 0
                        @info "Skipping snapshot $(snapf) because no refinement is requested."
                        continue
                    end
                end

                success = try
                    @info "Analysing snapshot $(snapf)..."
                    analyse(w, snapf, save_box=save_refinement)
                    true
                catch
                    @info "...snapshot $(snapf) failed."
                    false
                end

                if success
                    @info "...snapshot $(snapf) done."
                    ibatch += 1
                    time_current = time()
                end
            end

            if ibatch>=batch
                # if more than `batch` snapshots have been converted
                # we exit the conversion loop and check for new snapshots
                # this is usefull, because this provides more recent 
                # snapshots quicker
                break
            end
        end

        # check how much time has passed
        # if a snapshot was found the timer is reset
        if time_passed_since(time_current) > timeout
            @info "Watchdog monitoring end. Total time: $(time_passed_since(time_start)) s."
            break
        else
            @info "time passed (since last conversion): $(time_passed_since(time_current))"
            @info "time passed (total): $(time_passed_since(time_start))"
        end 

        sleep(check_every)
    end
end







#= WatchDog analysis of snapshot =#

"""
    analyse(w::WatchDog, snapshot; save_box=false)

Monitor variables of given snapshot. Save relevant quantities in one single 
HDF5 file, so that they can later be accessed fast and easy.
You can pass any number of functions to the WatchDog, that will receive the watchdog and the boxes.
They can compute whatever, but they need to return a Dict with whatever they computed
and the corresponding data.
"""
function analyse(w::WatchDog, snapshot; save_box=false)
    b, bτ = snapshotBox(
        snapshot; 
        folder=w.folder, 
        optical_depth_scale=true, 
        save_snapshot=save_box, 
        to_multi=false,
        is_box=true, 
        use_mmap=false,
        legacy=false
    )

    monitoring = Dict()
    for (fname, f) in w.functions
        monitoring[fname] = f(w, b, bτ)
    end

    b = nothing
    bτ = nothing
    Base.GC.gc()

    monitoring["general"] = Dict(
        "name" => w.name,
        "folder" => w.folder,
        "snapshot" => snapshot
    )

    # save the data
    save(w, monitoring)
end







#= Default monitoring functions =#

atmosphericParameters(w::WatchDog, b, bτ) = begin
    Dict(
        "teff" => b.parameter.teff,
        "logg" => b.parameter.logg,
        "time" => b.parameter.time
    )
end



_geostatistic(f, b) = begin
    variables = keys(b.data) |> collect
    
    results = Dict()
    for v in variables
        x, y = profile(f, b, :z, v) 
        results[:z] = x
        results[v] = y
    end

    results
end

# quantiles
geometrical15thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.15), b)
geometrical30thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.30), b)
geometrical45thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.45), b)
geometrical60thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.60), b)
geometrical75thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.75), b)
geometrical90thQuantile(w::WatchDog, b, bτ) = _geostatistic(x->quantile(reshape(x, :), 0.90), b)

geometricalAverages(w::WatchDog, b, bτ) = _geostatistic(mean, b)
geometricalMinimum(w::WatchDog, b, bτ) = _geostatistic(minimum, b)
geometricalMaximum(w::WatchDog, b, bτ) = _geostatistic(maximum, b)
geometricalRMS(w::WatchDog, b, bτ) = _geostatistic(x->sqrt(mean(x .^2)), b)
geometricalStd(w::WatchDog, b, bτ) = _geostatistic(std, b)

geoMassFlux(w::WatchDog, b, bτ) = begin
    add!(b, ρuz = b[:d] .* b[:uz])
    massFlux = profile(mean, b, :z, :ρuz)
    _, d = profile(mean, b, :z, :d)


    Dict(
        "massFlux" => massFlux[2] ./ d,
        "z" => massFlux[1]
    )
end

upperBoundarySurface(w::WatchDog, b, bτ) = begin
    uzplane = b[:uz][:, :, end]
    uxplane = b[:ux][:, :, end]
    uyplane = b[:uy][:, :, end]
    Tplane = b[:T][:, :, end]
    Dplane = log.(b[:d][:, :, end])
    dtplane = haskey(b.data, :dt_rt) ? b[:dt_rt][:, :, end] : nothing
    fluxplane = haskey(b.data, :flux) ? b[:flux][:, :, end] : nothing
    qrplane = haskey(b.data, :qr) ? b[:qr][:, :, end] : nothing

    d = Dict(
        "uzplane" => uzplane,
        "uxplane" => uxplane,
        "uyplane" => uyplane,
        "Tplane" => Tplane,
        "lnDplane" => Dplane,
        "x" => b.x[:, :, end],
        "y" => b.y[:, :, end]
    )

    if !isnothing(dtplane)
        d["dtplane"] = dtplane[:, :]
    end
    if !isnothing(fluxplane)
        d["fluxplane"] = fluxplane[:, :]
    end
    if !isnothing(qrplane)
        d["qrplane"] = qrplane[:, :]
    end

    d
end

lowerBoundarySurface(w::WatchDog, b, bτ) = begin
    uzplane = b[:uz][:, :, 1]
    uxplane = b[:ux][:, :, 1]
    uyplane = b[:uy][:, :, 1]
    Tplane = b[:T][:, :, 1]
    Dplane = log.(b[:d][:, :, 1])
    dtplane = haskey(b.data, :dt_rt) ? b[:dt_rt][:, :, 1] : nothing
    fluxplane = haskey(b.data, :flux) ? b[:flux][:, :, 1] : nothing
    qrplane = haskey(b.data, :qr) ? b[:qr][:, :, 1] : nothing

    d = Dict(
        "uzplane" => uzplane,
        "uxplane" => uxplane,
        "uyplane" => uyplane,
        "Tplane" => Tplane,
        "lnDplane" => Dplane,
        "x" => b.x[:, :, 1],
        "y" => b.y[:, :, 1]
    )

    if !isnothing(dtplane)
        d["dtplane"] = dtplane[:, :]
    end
    if !isnothing(fluxplane)
        d["fluxplane"] = fluxplane[:, :]
    end
    if !isnothing(qrplane)
        d["qrplane"] = qrplane[:, :]
    end

    d
end

minimumTempSurface(w::WatchDog, b, bτ) = begin
    # plane of minimum temperature
    iplane = argmin(b[:T])[3]

    uzplane = b[:uz][:, :, iplane]
    uxplane = b[:ux][:, :, iplane]
    uyplane = b[:uy][:, :, iplane]

    Tplane = b[:T][:, :, iplane]
    Dplane = log.(b[:d][:, :, iplane])
    dtplane = haskey(b.data, :dt_rt) ? b[:dt_rt][:, :, iplane] : nothing
    fluxplane = haskey(b.data, :flux) ? b[:flux][:, :, iplane] : nothing

    d = Dict(
        "uzplane" => uzplane,
        "uxplane" => uxplane,
        "uyplane" => uyplane,
        "Tplane" => Tplane,
        "lnDplane" => Dplane,
        "x" => b.x[:, :, iplane],
        "y" => b.y[:, :, iplane]
    )
    if !isnothing(dtplane)
        d["dtplane"] = dtplane
    end
    if !isnothing(fluxplane)
        d["fluxplane"] = fluxplane
    end

    d
end

centerVerticalCut(w::WatchDog, b, bτ) = begin
    ixd = floor(Int, size(b, 1) /2)

    uzplane = b[:uz][ixd, :, :]
    uxplane = b[:ux][ixd, :, :]
    uyplane = b[:uy][ixd, :, :]

    Tplane = b[:T][ixd, :, :]
    Dplane = log.(b[:d][ixd, :, :])
    dtplane = haskey(b.data, :dt_rt) ? b[:dt_rt][ixd, :, :] : nothing
    fluxplane = haskey(b.data, :flux) ? b[:flux][ixd, :, :] : nothing
    
    d = Dict(
        "uzplane" => uzplane,
        "uxplane" => uxplane,
        "uyplane" => uyplane,
        "Tplane" => Tplane,
        "lnDplane" => Dplane,
        "y" => b.y[ixd, :, :],
        "z" => b.z[ixd, :, :]
    )
    if !isnothing(dtplane)
        d["dtplane"] = dtplane
    end
    if !isnothing(fluxplane)
        d["fluxplane"] = fluxplane
    end

    d
end





_optstatistic(f, bτ) = begin
    variables = keys(bτ.data) |> collect
    
    results = Dict()
    for v in variables
        x, y = profile(f, bτ, :log10τ_ross, v) 
        results[:log10τ_ross] = x
        results[v] = y
    end

    results
end

opticalAverages(w::WatchDog, b, bτ) = _optstatistic(mean, bτ)
opticalMaximum(w::WatchDog, b, bτ) = _optstatistic(maximum, bτ)
opticalMinimum(w::WatchDog, b, bτ) = _optstatistic(minimum, bτ)
opticalRMS(w::WatchDog, b, bτ) = _optstatistic(x->sqrt(mean(x .^2)), bτ) 
opticalStd(w::WatchDog, b, bτ) = _optstatistic(std, bτ)

optMassFlux(w::WatchDog, b, bτ) = begin
    add!(bτ, ρuz = bτ[:d] .* bτ[:uz])
    massFlux = profile(mean, bτ, :log10τ_ross, :ρuz)
    _, d = profile(mean, bτ, :log10τ_ross, :d)

    Dict(
        "massFlux" => massFlux[2] ./d,
        "log10τ_ross" => massFlux[1]
    )
end

opticalSurfaces(w::WatchDog, b, bτ) = begin
    add!(b, lnd=log.(b[:d]))
    uzplane = interpolate_to(b, :uz; logspace=true, τ_ross=0.0)[:uz]
    uxplane = interpolate_to(b, :ux; logspace=true, τ_ross=0.0)[:ux]
    uyplane = interpolate_to(b, :uy; logspace=true, τ_ross=0.0)[:uy]
    Tplane = interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
    Dplane = interpolate_to(b, :lnd; logspace=true, τ_ross=0.0)[:lnd]
    dtplane = haskey(b.data, :dt_rt) ? interpolate_to(b, :dt_rt; logspace=true, τ_ross=0.0)[:dt_rt] : nothing
    fluxplane = haskey(b.data, :flux) ? interpolate_to(b, :flux; logspace=true, τ_ross=0.0)[:flux] : nothing
    qrplane = haskey(b.data, :qr) ? interpolate_to(b, :qr; logspace=true, τ_ross=0.0)[:qr] : nothing

    d = Dict(
        "uzplane" => uzplane[:, :, 1],
        "uxplane" => uxplane[:, :, 1],
        "uyplane" => uyplane[:, :, 1],
        "Tplane" => Tplane[:, :, 1],
        "lnDplane" => Dplane[:, :, 1],
        "x" => b.x[:, :, 1],
        "y" => b.y[:, :, 1]
    )

    if !isnothing(dtplane)
        d["dtplane"] = dtplane[:, :, 1]
    end
    if !isnothing(fluxplane)
        d["fluxplane"] = fluxplane[:, :, 1]
    end
    if !isnothing(qrplane)
        d["qrplane"] = qrplane[:, :, 1]
    end

    d
end





#= Spectrum synthesis =#

resolvedSpectra(lambda_min, lambda_max, w, b, bτ; R=SPECTRUM_RESOLUTION[], Δλ=nothing) = begin
    dλ = if isnothing(Δλ)
        (lambda_max+lambda_min)/2.0/R
    else
        Δλ
    end
    d = Dict{Any, Any}()
    
    # only compute spectra if the WatchDog says its time
    if w.refinement > 0
        try
            λs = lambda_min
            λe = lambda_max
            Δλ = dλ
            name = w.name
            run = w.folder
            isnap = w.currentSnapshot
            _, _, feh = parametersFromName(name)
            multi_name = joinpath(run, "m3dis_$(isnap)")

            # linelists relative to M3D 
            linelists = String[
                "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata",
                "input_multi3d/LINE-LISTS/nlte_ges_linelist_jmg04sep2023_I_II"
            ]

            # make sure, that M3D is setup correctly
            if !isdir(@in_m3dis("input_multi3d/MUST"))
                @warn "There is not yet a MUST directory in M3D. This is needed to keep the paths short!"
                symlink(@in_dispatch("data"), @in_m3dis("input_multi3d/MUST"))
            end

            linelists = if !isdir(@in_m3dis("input_multi3d/LINE-LISTS"))
                @warn "The default linelists can not be found! Provide the following: $(linelists)"
                []
            else
                linelists
            end
            
            # downsample to make it faster
            bDown = gresample(b; nx=SPECTRUM_DOWNSAMPLING[], ny=SPECTRUM_DOWNSAMPLING[], nz=size(b, 3) * 2 - 1)
            bDown.data[:ne] = pop!(bDown.data, :Ne)
            multiBox(bDown, joinpath(run, "m3dis_$(isnap)"))

            # setup for M3D
            spectrum_namelist = Dict(
                :model_folder=>joinpath("input_multi3d","MUST", name),
                :linelist=>nothing,
                :absmet=>nothing,
                :linelist_params=>(:line_lists=>linelists,),
                :atom_params=>(:atom_file=>"",),
                :spectrum_params=>(:daa=>Δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>false),
                :atmos_params=>(
                    :dims=>30, 
                    :atmos_format=>"must",
                    :use_density=>true, 
                    :use_ne=>false,
                    :FeH=>feh
                ),
                :m3d_params=>(
                    :save_resolved=>true,
                    :long_scheme=>"disk_center"
                ),
            )
            m3dis_kwargs = Dict(:threads=>SPECTRUM_THREADS[])

            @info "Spectrum synthesis between $(lambda_min)Å and $(lambda_max)Å (Δλ=$(dλ)) requested for snapshot $(isnap)."
            # Running M3D
            result = spectrum(
                "m3dis_$(isnap)"; 
                name=name, 
                NLTE=false, 
                slurm=false, 
                namelist_kwargs=spectrum_namelist, 
                m3dis_kwargs=m3dis_kwargs
            )

            # get the output arrays
            i = pyconvert(Array, result.i3)
            c = pyconvert(Array, result.c3)
            l = pyconvert(Array, result.lam)

            d["intensity"] = i 
            d["continuum"] = c
            d["wavelength"] = l
            d["x"] = bDown.x[:, :, end]
            d["y"] = bDown.y[:, :, end]

            # delete the M3D output folder
            mn = join(["m3dis_$(isnap)", name], "_")
            output_folder = @in_m3dis("data/$(mn)")
            if isdir(output_folder)
                rm(output_folder, recursive=true)
            end
        catch
            @warn "Spectrum synthesis between $(lambda_min)Å and $(lambda_max)Å (Δλ=$(dλ)) failes for snapshot $(isnap)."
            #error("Spectrum synthesis failed.")
        end
    end

    d
end


resolvedSpectra4MOST3(w::WatchDog, b, bτ) = resolvedSpectra(6200.0, 6790.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])
resolvedSpectraAPOGEESDSSV(w::WatchDog, b, bτ) = resolvedSpectra(15140.0, 17000.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])
resolvedSpectraGaiaESOHR10(w::WatchDog, b, bτ) = resolvedSpectra(5339.0, 5619.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])
resolvedSpectraGaiaRVS(w::WatchDog, b, bτ) = resolvedSpectra(8470.0, 8740.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])
resolvedSpectraGaiaRVSCATriplet(w::WatchDog, b, bτ) = resolvedSpectra(8480.0, 8560.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])
resolvedSpectraHalpha(w::WatchDog, b, bτ) = resolvedSpectra(6562.8-20.0, 6562.8+20.0, w, b, bτ; Δλ=SPECTRUM_RESOLUTION[])







"""
    defaultMonitoring(name; additional_functions...) 

Create a WatchDog with all default monitoring functions. Additional functions
can be passed if wanted.
"""
defaultWatchDog(name; folder=@in_dispatch("data/"), additional_functions...) = WatchDog(
    name;
    folder=folder,
    atmosphericParameters = atmosphericParameters,
    geometricalAverages = geometricalAverages,
    geometricalRMS = geometricalRMS,
    geometricalStd = geometricalStd, 
    opticalAverages = opticalAverages,
    opticalRMS = opticalRMS,
    opticalStd = opticalStd,
    optMassFlux = optMassFlux,
    geoMassFlux = geoMassFlux,
    minimumTempSurface = minimumTempSurface,
    opticalSurfaces = opticalSurfaces,
    opticalMaximum = opticalMaximum,
    opticalMinimum = opticalMinimum,
    geometricalMinimum = geometricalMinimum,
    geometricalMaximum = geometricalMaximum,
    upperBoundarySurface = upperBoundarySurface,
    lowerBoundarySurface = lowerBoundarySurface,
    geometrical15thQuantile = geometrical15thQuantile,
    geometrical30thQuantile = geometrical30thQuantile,
    geometrical45thQuantile = geometrical45thQuantile,
    geometrical60thQuantile = geometrical60thQuantile,
    geometrical75thQuantile = geometrical75thQuantile,
    geometrical90thQuantile = geometrical90thQuantile,
    centerVerticalCut = centerVerticalCut,
    #resolvedSpectraAPOGEESDSSV = resolvedSpectraAPOGEESDSSV,
    resolvedSpectraGaiaESOHR10 = resolvedSpectraGaiaESOHR10,
    resolvedSpectraGaiaRVS = resolvedSpectraGaiaRVS,
    resolvedSpectraGaiaRVSCATriplet = resolvedSpectraGaiaRVSCATriplet,
    resolvedSpectraHalpha = resolvedSpectraHalpha,
    additional_functions...
)







#= Utility functions =#

updatesnaps!(w::WatchDog) = begin
    w.snapshots = sort(list_of_snapshots(w.folder))
    w.snapshotsCompleted = snapshotnumber.(availableSnaps(w))
end

snapshotnumber(snap) = parse(Int, split(split(snap, "_", keepempty=false) |> last, ".hdf5", keepempty=false) |> first)

monitoringPath(w) = joinpath(w.folder, "monitoring")
monitoringPath(w, snap) = joinpath(w.folder, "monitoring", "snap_$(snap).hdf5")
availableSnaps(w) = begin
    l = glob("snap*", monitoringPath(w))
    ln = snapshotnumber.(l)
    m = sortperm(ln)

    l[m]
end






"""
    deleteSnapshot(w, snap)

Delete the snapshot with number snap from the monitored run.
The function removes the snapshot data to save disk space, but
does not remove it from the monitoring. This allows to check the progress
without saving every snapshot.
"""
deleteSnapshot(w, snap) = begin
    # check if monitoring really exists
    if !(snap in w.snapshotsCompleted)
        @warn "Deleting requested for snapshot $(snap), which is not monitored yet."
        return
    end

    _, _, datadir, _, _ = _check_files(snap, nothing, w.name, w.folder, check_existing=false)
    if isdir(datadir)
        rm(datadir, recursive=true)
        @info "Snapshot $(snap) removed from $(datadir)."
    end
end

"""
    deleteMonitoring(w, snap)

Delete the snapshot with number snap from the monitoring.
This is usefull if you want to re-do the analysis of the given snapshot.
"""
deleteMonitoring(w, snap) = begin
    # check if monitoring really exists
    if !(snap in w.snapshotsCompleted)
        @warn "Deleting requested for snapshot $(snap), which is not monitored yet."
        return
    end

   # remove the monitoring
   rm(monitoringPath(w, snap))
end







add_to_hdf5!(fid, fname, val)       = fid["$(fname)"] = val
add_to_hdf5!(fid, fname, val::Bool) = fid["$(fname)"] = Int(val)

save(w::WatchDog, monitoring) = begin
    # check if folder exists
    if !isdir(monitoringPath(w))
        mkdir(monitoringPath(w))
    end

    isnap = monitoring["general"]["snapshot"]

    # within the same folder we create the monitoring
    fid = HDF5.h5open(monitoringPath(w, isnap), "w")

    # for each function we called, we add a group to the monitoring
    for (fname, data) in monitoring
        create_group(fid, "$(fname)")

        # all variables in the group we add to the monitoring
        for (quantity, var) in data
            add_to_hdf5!(fid["$(fname)"], quantity, var)
        end
    end

    close(fid)
end





reload(s::Type{S}, name, snap; folder=@in_dispatch("data/"), mmap=false, groups=nothing) where {S<:WatchDog} = begin
    reload(s(name; folder=folder), snap, mmap=mmap, groups=groups)
end

reload!(s::Type{S}, name; folder=@in_dispatch("data/"), mmap=false, asDict=false, groups=nothing) where {S<:WatchDog} = begin
    reload!(s(name; folder=folder), mmap=mmap, asDict=asDict, groups=groups)
end

reload(w::S, snap; mmap=false, groups=nothing) where {S<:WatchDog} = begin
    fid = HDF5.h5open(monitoringPath(w, snap), "r")
    fvals = Dict()

    for gname in keys(fid)
        # check if this group should be loaded
        if !isnothing(groups)
            if !(gname in groups)
                continue
            end
        end

        # read from file
        fvals[gname] = Dict()
        group = fid[gname]
        for dname in keys(group)
            fvals[gname][dname] = if (mmap && 
                                    HDF5.ismmappable(group[dname]) && 
                                    (HDF5.get_jl_type(group[dname])<:AbstractArray))
                HDF5.readmmap(group[dname])
            else
                HDF5.read(group[dname])
            end
        end
    end

    close(fid)

    fvals
end

"""
    reload!(w::S; mmap=false, asDict=false, groups=nothing, lastN=:all)

Load the monitoring of the given WatchDog. Optionally specify a list of groups you
want to load, or how many snapshots you want to read from the end. This is 
usefull when you only need a specific field of the last couple of snapshots.
"""
reload!(w::S; mmap=false, asDict=false, groups=nothing, lastN=:all) where {S<:WatchDog} = begin
    listOfSnaps = availableSnaps(w)
    w.snapshotsCompleted = []
    firstSnap = lastN == :all ? 1 : max(length(listOfSnaps) - (lastN-1), 1)
    if !asDict
        l  = []
        for snap in listOfSnaps[firstSnap:end]
            try
                si = reload(w, snapshotnumber(snap); mmap=mmap, groups=groups)
                append!(l, [si])
                append!(w.snapshotsCompleted, [snapshotnumber(snap)])
            catch
                nothing
            end
        end

        l
    else
        l = Dict()
        for snap in listOfSnaps[firstSnap:end]
            try
                l[snap] = reload(w, snapshotnumber(snap); mmap=mmap, groups=groups)
                append!(w.snapshotsCompleted, [snapshotnumber(snap)])
            catch
                nothing
            end
        end

        l
    end
end









#= Convenience functions =#

"""
    timeevolution(m, group, field) 

Return an array containing the time evolution from all available snapshots in `m` 
for the given topic `group` and variable `field` in that topic.
"""
timeevolution(m, group, field) = [
    moni[group][field] 
    for moni in m
]

"""
    timeevolution(m, group) 

Return an Dictionary containing arrays of the time evolution from all available 
snapshots in `m` for the given topic `group` with all variables in that topic as dictionary entried.
"""
timeevolution(m, group) = begin
    mg = m[1][group]
    Dict(
        k => [moni[group][k] for moni in m]
        for k in keys(mg)
    )
end








#= GLobal parameters =#

const SPECTRUM_THREADS = Ref(2)
const SPECTRUM_RESOLUTION = Ref(0.05)
const SPECTRUM_DOWNSAMPLING = Ref(40)

