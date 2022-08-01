using Pkg; Pkg.activate("."); #Pkg.update(); 
using MUST
using Glob
using PyCall
using Distributed, SlurmClusterManager#, ParallelDataTransfer

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
    using PyCall
    using Distributed
    using Interpolations
end

@everywhere begin
    MUST.@import_dispatch "../../../dispatch2_clean/dispatch2"
    MUST.@import_dispatch "../../../dispatch2_clean/dispatch2" EOS select

    # Different EOS
    const eos_legacy = MUST.@legacyPythonEOS
end

folder = MUST.@in_dispatch ARGS[1]
if length(ARGS) > 1
    snapshots = sort(parse.(Int, ARGS[2:end]))
else
    content_of_folder = glob("*/", folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))
end

# Name of the namelist of the current folder
nml_name = MUST.@in_dispatch splitpath(folder)[end]

# Init namelist
const nml = MUST.StellarNamelist(nml_name*".nml")

# Effective temperature present?
save_info = ispath(joinpath(folder, "teff.dat")) 
if save_info
    # Read and interpolate teff data
    teff_func = MUST.teff_interpolated(joinpath(folder, "teff.dat"))
    
    # Send the data to the worksers
    MUST.sendsync(workers(), teff=teff_func)
end

# Send the data to the workers
MUST.sendsync(workers(), folder=folder, do_teff=save_info, ini_nml=nml)#, eos=eos_sq)

@everywhere function convert_snapshots(snapshots)
    for i_s in eachindex(snapshots)
        @info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
        try
            # The dispatch snapshot object (Python)
            snap = dispatch.snapshot(snapshots[i_s], data=folder)

            # Convert its content to pure Julia
            s = MUST.Space(snap, :d, :ee, :uz, :e, :tt)

            # Units for conversion to CGS
            units = MUST.StaggerCGS(snap)

            #MUST.add_from_EOS!(s, eos_legacy, :T)
            # Check if Teff should be added
            if do_teff
                @info "setting teff for sn $(snapshots[i_s])"
                MUST.set!(s.parameter, teff, ini_nml)
            end

            # Apply the conversion
            MUST.convert!(s, units; d=:d,   ee=:ee, uz=:u,  e=:e, 
                                    x=:l,   y=:l,   z=:l)

            # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
            b_s = MUST.Box(s)

            # First save
            MUST.save(b_s; name="box_sn$(snapshots[i_s])",   folder=folder)

        catch
            @warn "snapshot $(snapshots[i_s]) could not be loaded."
            continue
        end
    end
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




