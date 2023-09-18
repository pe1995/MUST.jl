using MUST
using TSO

function snaps2multi(folder, snaps...; 
						eos, label, 
						n_horizontal=nothing, n_vertical=nothing, outfolder="", method=:linear, name="m3dis")
    i = 0
    for snap in snaps
        snapshot, snapshot_τ = try
            pick_snapshot(folder, snap)
        catch
            continue
        end

        aos = @axed eos
        ne = lookup(
            aos, :lnNe, 
            log.(snapshot[:d]), 
            log.(snapshot[:ee])
        ) 
        snapshot.data[:ne] = exp.(ne)

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

        #@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
    end
end

function snaps2multi(snaps::MUST.Box...; 
                    eos, label, n_horizontal=nothing, n_vertical=nothing, outfolder="", name="m3dis")
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
        aos = @axed eos

        if is_internal_energy(aos)
            ne = lookup(
                aos, :lnNe, 
                log.(snapshotsresample[:d]), log.(snapshotsresample[:ee])
            ) 
        else
            ne = lookup(
                aos, :lnNe, 
                log.(snapshotsresample[:d]), log.(snapshotsresample[:T])
            ) 
        end
        
        snapshotsresample.data[:ne] = exp.(ne)
    
        MUST.multiBox(
            snapshotsresample, 
            output_name
        )

        #@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
    end
end
