using MUST
using TSO

function snaps2multi(folder, snaps...; 
						eos, label, 
						n_horizontal=10, n_vertical=280, outfolder="", method=:linear)
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
        output_name = joinpath(outfolder, "m3dis_sun_$(label)_$(i)")
    
        MUST.multiBox(
            snapshotsresample, 
            output_name
        )

        #@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
    end
end

function snaps2multi(snaps::MUST.Box...; 
                    eos, label, n_horizontal=10, n_vertical=280, outfolder="")
    i = 0
    for snap in snaps
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

        i += 1
        output_name = joinpath(outfolder, "m3dis_sun_$(label)_$(i)")
        aos = @axed eos

        ne = lookup(
            aos, :lnNe, 
            log.(snapshotsresample[:d]), log.(snapshotsresample[:ee])
        ) 
        
        snapshotsresample.data[:ne] = exp.(ne)
    
        MUST.multiBox(
            snapshotsresample, 
            output_name
        )

        #@info "New size: $(size(@view(ne[1:downsample_xy:end, 1:downsample_xy:end, :])))"
    end
end
