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
    functions
end

WatchDog(name; folder=@in_dispatch("data/"), functions...) = WatchDog(name, joinpath(folder, name), [], [], functions)


"""
    monitor(w::WatchDog)

Monitor the progress of the linked simulation. Compute statistics and save 
depth profiles. Check for new snapshots after `check_every` seconds.
Cancel the monitoring if `timeout` seconds have passed without
finding a new snapshot.
"""
function monitor(w::WatchDog; timeout=2*60*60, check_every=5, delay=0, snapshotbuffer=1, save_box=false)
    time_start = time()
    time_current = time()
    time_passed_since(t_ref) = time() - t_ref

    # keep looping until timeout
    while true
        # list snapshots
        updatesnaps!(w)

        # check if there is a new snapshot
        for (i, snapf) in enumerate(w.snapshots)
            # for safety reasons, we dont convert the last N snaps
            # this avoids them being read while not written properly
            if i >= length(w.snapshots)-snapshotbuffer
                continue
            end

            if !(snapf in w.snapshotsCompleted)
                # give the snapshot some time to be written (optional)
                sleep(delay)

                success = try
                    @info "Analysing snapshot $(snapf)..."
                    analyse(w, snapf, save_box=save_box)
                    true
                catch
                    @info "...snapshot $(snapf) failed."
                    false
                end

                if success
                    @info "...snapshot $(snapf) done."
                    time_current = time()
                end
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
        w.folder, 
        optical_depth_scale=true, 
        save_snapshot=save_box, 
        to_multi=false,
        is_box=true, 
        use_mmap=false
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
    Tplane = b[:T][:, :, end]
    Dplane = log.(b[:d][:, :, end])

    Dict(
        "uzplane" => uzplane[:, :],
        "Tplane" => Tplane[:, :],
        "lnDplane" => Dplane[:, :],
        "x" => b.x[:, :, end],
        "y" => b.y[:, :, end]
    )
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
    Tplane = interpolate_to(b, :T; logspace=true, τ_ross=0.0)[:T]
    Dplane = interpolate_to(b, :lnd; logspace=true, τ_ross=0.0)[:lnd]

    Dict(
        "uzplane" => uzplane[:, :, 1],
        "Tplane" => Tplane[:, :, 1],
        "lnDplane" => Dplane[:, :, 1],
        "x" => b.x[:, :, 1],
        "y" => b.y[:, :, 1]
    )
end






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
    opticalSurfaces = opticalSurfaces,
    opticalMaximum = opticalMaximum,
    opticalMinimum = opticalMinimum,
    geometricalMinimum = geometricalMinimum,
    geometricalMaximum = geometricalMaximum,
    upperBoundarySurface=upperBoundarySurface,
    geometrical15thQuantile=geometrical15thQuantile,
    geometrical30thQuantile=geometrical30thQuantile,
    geometrical45thQuantile=geometrical45thQuantile,
    geometrical60thQuantile=geometrical60thQuantile,
    geometrical75thQuantile=geometrical75thQuantile,
    geometrical90thQuantile=geometrical90thQuantile,
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




reload(s::Type{S}, name, snap; folder=@in_dispatch("data/"), mmap=false) where {S<:WatchDog} = begin
    reload(s(name; folder=folder), snap, mmap=mmap)
end

reload(s::Type{S}, name; folder=@in_dispatch("data/"), mmap=false, asDict=false) where {S<:WatchDog} = begin
    reload(s(name; folder=folder), mmap=mmap, asDict=asDict)
end

reload(w::S, snap; mmap=false) where {S<:WatchDog} = begin
    fid = HDF5.h5open(monitoringPath(w, snap), "r")
    fvals = Dict()

    for gname in keys(fid)
        fvals[gname] = Dict()
        group = fid[gname]
        for dname in keys(group)
            fvals[gname][dname] = mmap ? HDF5.readmmap(group[dname]) : HDF5.read(group[dname])
        end
    end

    close(fid)

    fvals
end

reload(w::S; mmap=false, asDict=false) where {S<:WatchDog} = begin
    listOfSnaps = availableSnaps(w)
    
    if !asDict
        l  = []
        for snap in listOfSnaps
            try
                si = reload(w, snapshotnumber(snap); mmap=mmap)
                append!(l, [si])
            catch
                nothing
            end
        end

        l
    else
        l = Dict()
        for snap in listOfSnaps
            try
                l[snap] = reload(w, snapshotnumber(snap); mmap=mmap)
            catch
                nothing
            end
        end

        l
    end
end
