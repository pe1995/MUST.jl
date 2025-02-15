using Pkg; Pkg.activate("."); #Pkg.update(); 
using MUST
using Glob
using Distributed, SlurmClusterManager#, ParallelDataTransfer
PythonCall = MUST.PythonCall

if "SLURM_NTASKS" in keys(ENV)
    addprocs(SlurmManager())
    for i in workers()
        host,pid = fetch(@spawnat i (gethostname(), getpid()))
        @info "Worker $(i) is running on $(host) with ID $(pid)" 
    end
else
    @warn "No Slurm environment detected. Using default addprocs."
    addprocs(1)
end

@everywhere begin
    using Pkg; Pkg.activate("."); #Pkg.update(); 
    using MUST
    using Glob
    using Distributed
    using Interpolations

    PythonCall = MUST.PythonCall
end




@everywhere begin
    #MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
    #MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" EOS select

    MUST.@import_dispatch "/home/eitner/shared/model_grid/dispatch2"
    MUST.@import_dispatch "/home/eitner/shared/model_grid/dispatch2" EOS select

    add_selection = false
end

folder = MUST.@in_dispatch ARGS[1]
if length(ARGS) > 1
    snapshots = sort(parse.(Int, ARGS[2:end]))
else
    content_of_folder = glob("*/", folder)
    snapshots  = sort(MUST.list_of_snapshots(content_of_folder))
    snapshots = snapshots[end-1:end-1]
end



# Name of the namelist of the current folder
nml_name = MUST.@in_dispatch splitpath(folder)[end]

# Init namelist
const nml = MUST.StellarNamelist(nml_name*".nml")

# Use the new Squaregas EOS 
eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
const eos_sq = MUST.SquareGasEOS(MUST.@in_dispatch(eos_path))

# Effective temperature present?
save_info = ispath(joinpath(folder, "teff.dat")) 
if save_info
    # Read and interpolate teff data
    teff_func = MUST.teff_interpolated(joinpath(folder, "teff.dat"))
    
    # Send the data to the worksers
    MUST.sendsync(workers(), teff=teff_func)
end

# Send the data to the workers
MUST.sendsync(workers(), folder=folder, do_teff=save_info, ini_nml=nml, eos=eos_sq)

@everywhere function convert_snapshots(snapshots)
    for i_s in eachindex(snapshots)
        @info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
        try
            # The dispatch snapshot object (Python)
            snap = dispatch.snapshot(snapshots[i_s], data=folder)

            # Units for conversion to CGS
            units = MUST.StaggerCGS(snap)

            # Convert its content to pure Julia
            s = if add_selection
                MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg)
            else
                MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)
            end

            # Apply the conversion
            if add_selection
                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    qr=:qr, pg=:p,
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l, time=:t)
            else
                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l, time=:t)
            end

            # Add additional columns already in CGS after converting
            MUST.add_from_EOS!(s, eos, :T)
            MUST.add_from_EOS!(s, eos, :kr)
            MUST.add_from_EOS!(s, eos, :Ne)

            # Check if Teff should be added
            if do_teff
                @info "setting teff for sn $(snapshots[i_s])"
                MUST.set!(s.parameter, teff, ini_nml)
            end

            # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
            b_s = MUST.Box(s)
            s = nothing

            # convert to Multi format and save in the same folder, keep the number of the snapshot!
            b_s.data[:ne] = b_s.data[:Ne]
            MUST.multiBox(b_s, joinpath(folder, "m3dis_$(snapshots[i_s])"))
        catch
            @warn "snapshot $(snapshots[i_s]) could not be loaded."
            continue
        end
    end
    nothing
end

futures = Distributed.Future[]
wrk     = sort(workers())
for (i,split) in enumerate(MUST.split_similar(snapshots, nworkers()))

    # Send the split to the worker
    MUST.sendsync(wrk[i], split_l=split)

    # Execute the function convert_snapshots on the passed split
    append!(futures, [@spawnat wrk[i] convert_snapshots(split_l)])
end

# Wait for the execution to finish
MUST.wait(futures)




