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
    addprocs(4)
end

@everywhere begin
    using Pkg; Pkg.activate("."); #Pkg.update(); 
    using MUST
    using Glob
    using PyCall
    using Distributed
end

@everywhere begin
    MUST.@import_dispatch "../../../dispatch2"
    MUST.@import_dispatch "../../../dispatch2" EOS select

    # Different EOS
    #const eos        = dispatch.EOS.bifrost(top=MUST.dispatch_location, dir="data/eos/stagger/");
    #const eos_legacy = dispatch.EOS.stagger(top=MUST.dispatch_location, file="table.dat");

    # Use the new Squaregas EOS 
    eos = MUST.SquareGasEOS(MUST.@in_dispatch("input_data/square_gas_1"));
end

folder = MUST.@in_dispatch ARGS[1]
if length(ARGS) > 1
    snapshots = sort(parse.(Int, ARGS[2:end]))
else
    content_of_folder = glob("*/", folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))
end

# Send the data to the worksers
MUST.sendsync(workers(), folder=folder)

@everywhere function convert_snapshots(snapshots)
    for i_s in 1:length(snapshots)
        @info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
        try
            # The dispatch snapshot object (Python)
            snap = dispatch.snapshot(snapshots[i_s], data=folder)

            # Convert its content to pure Julia
            s = MUST.Space(snap, :d, :ee, :uz, :e)

            # Units for conversion to CGS
            units = MUST.StaggerCGS(snap)

            # Apply the conversion
            MUST.convert!(s, units; d=:d,   ee=:ee, uz=:u,  e=:e, 
                                    x=:l,   y=:l,   z=:l)

            # Add additional columns already in CGS after converting
            MUST.add_from_EOS!(s, eos, :T)
            MUST.add_from_EOS!(s, eos, :kr)

            # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
            b_s = MUST.Box(s)

            # Add the optical depth
            τ = MUST.optical_depth(b_s, opacity=:kr, density=:d)
            MUST.add!(b_s, τ, :τ_ross)

            # First save
            MUST.save(b_s; name="box_sn$(snapshots[i_s])",   folder=folder)

            # NEW: Convert the height scale from cm to optical depth
            b_s = MUST.height_scale(b_s, :τ_ross);

            # Write to HDF5 file. Can easily be read as a Memory map later with the build in functions
            #MUST.save(s;   name="space_sn$(snapshots[i_s])", folder=folder)
            MUST.save(b_s; name="box_tau_sn$(snapshots[i_s])",   folder=folder)
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

    @show split wrk[i]

    # Execute the function convert_snapshots on the passed split
    append!(futures, [@spawnat wrk[i] convert_snapshots(split_l)])
end

# Wait for the execution to finish
MUST.wait(futures)




