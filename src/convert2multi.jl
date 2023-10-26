using MUST
using TSO

function snaps2multi(folder, snaps...; 
						label, eos=nothing,
						n_horizontal=nothing, n_vertical=nothing, outfolder="", method=:linear, name="m3dis")
    i = 0
    for snap in snaps
        snapshot, snapshot_τ = try
            pick_snapshot(folder, snap)
        catch
            continue
        end

        if !isnothing(eos)
            aos = @axed eos
            ne = lookup(
                aos, :lnNe, 
                log.(snapshot[:d]), 
                log.(snapshot[:ee])
            ) 
            snapshot.data[:ne] = exp.(ne)
        elseif !haskey(snapshot.data, :ne)
            @warn "No electron density (:ne) found! Filling it with 1.0..."
            snapshot.data[:ne] = fill!(similar(snapshot.data[:d]), 1.0)
        end

        snapshotsresample = gresample(
            snapshot, 
            nz=n_vertical, 
            nx=n_horizontal, 
            ny=n_horizontal,
            method=method
        )

        if length(outfolder) > 1
            (!isdir(outfolder)) && mkdir(outfolder)
        end

        i += 1
        output_name = joinpath(outfolder, join([name, "$(label)_$(i)"], "_"))
    
        MUST.multiBox(
            snapshotsresample, 
            output_name
        )
    end
end

function snaps2multi(snaps::MUST.Box...; 
                    label, eos=nothing, n_horizontal=nothing, n_vertical=nothing, outfolder="", name="m3dis")
    labels = if typeof(label) <: AbstractVector
        label
    else
        ["$(label)_$(j)" for (j, _) in enumerate(snaps)]
    end

    for (i, snap) in enumerate(snaps)
        snapshot = snap
        snapshotsresample = gresample(
            snapshot, 
            nz=n_vertical, 
            nx=n_horizontal, 
            ny=n_horizontal
        )

        if length(outfolder) > 1
            (!isdir(outfolder)) && mkdir(outfolder)
        end
        output_name = joinpath(outfolder, join([name, "$(labels[i])"], "_"))
       
        if !isnothing(eos)
            aos = @axed eos
            ne = lookup(
                aos, :lnNe, 
                log.(snapshotsresample[:d]), 
                log.(snapshotsresample[:ee])
            ) 
            snapshotsresample.data[:ne] = exp.(ne)
        elseif !haskey(snapshotsresample.data, :ne)
            @warn "No electron density (:ne) found! Filling it with 1.0..."
            snapshotsresample.data[:ne] = fill!(similar(snapshotsresample.data[:d]), 1.0)
        end
            
        MUST.multiBox(
            snapshotsresample, 
            output_name
        )
    end
end
