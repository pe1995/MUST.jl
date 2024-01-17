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

WatchDog(name; folder=@in_dispatch("data/$(name)"), functions...) = WatchDog(name, folder, [], [], functions)


"""
    monitor(w::WatchDog)

Monitor the progress of the linked simulation. Compute statistics and save 
depth profiles. Check for new snapshots after `check_every` seconds.
Cancel the monitoring if `timeout` seconds have passed without
finding a new snapshot.
"""
function monitor(w::WatchDog; timeout=2*60*60, check_every=30)
    time_start = time()
    time_current = time()
    time_passed_since(t_ref) = time() - t_ref

    # keep looping until timeout
    while true
        # list snapshots
        updatesnaps!(w)

        # check if there is a new snapshot
        for (i, snapf) in enumerate(w.snapshots)
            if !(snapf in w.snapshotsCompleted)
                success = try
                    @info "Analysing snapshot $(snapf)..."
                    analyse(w, snapf)
                    true
                catch
                    @info "... snapshot $(snapf) failed."
                    false
                end

                if success
                    @info "... snapshot $(snapf) done."
                    time_current = time()
                    @show time_passed_since(time_current)
                    @show time_passed_since(time_start)
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
        save_snapshot=false, 
        add_selection=false, 
        to_multi=false,
        is_box=true, 
        use_mmap=false
    )

    monitoring = Dict()
    for (fname, f) in w.functions
        monitoring[fname] = f(w, b, bτ)
    end

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

geometricalAverages(w::WatchDog, b, bτ) = _geostatistic(mean, b)
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
        "y" => b.x[:, :, 1]
    )
end






"""
    defaultMonitoring(name; additional_functions...) 

Create a WatchDog with all default monitoring functions. Additional functions
can be passed if wanted.
"""
defaultWatchDog(name; folder=@in_dispatch("data/$(name)"), additional_functions...) = WatchDog(
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

reload(s::Type{S}, name::String, snap; folder=@in_dispatch("data/$(name)"), mmap=false) where {S<:WatchDog} = begin
    w = s(name)

    fid = HDF5.h5open(monitoringPath(w, snap), "w")
    fvals = Dict()

    for gname in keys(fid)
        fvals[gname] = Dict()
        for dname in keys(fid[gname])
            fvals[gname][dname] = mmap ? HDF5.readmmap(fid[gname][dname]) : HDF5.read(fid[gname][dname])
        end
    end

    close(fid)

    fvals
end

reload(s::Type{S}, name::String; folder=@in_dispatch("data/$(name)"), mmap=false, asDict=false) where {S<:WatchDog} = begin
    w = s(name)
    listOfSnaps = availableSnaps(w)
    
    if !asDict
        [reload(s, name, snapshotnumber(snap); folder=folder, mmap=mmap) for snap in listOfSnaps]
    else
        Dict(
            snap => reload(s, name, snap; folder=folder, mmap=mmap) 
            for snap in listOfSnaps
        )
    end
end

