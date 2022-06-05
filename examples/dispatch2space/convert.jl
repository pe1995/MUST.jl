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
    using Interpolations
end

@everywhere begin
    MUST.@import_dispatch "../../../dispatch2"
    MUST.@import_dispatch "../../../dispatch2" EOS select

    # Different EOS
    #const eos        = dispatch.EOS.bifrost(top=MUST.dispatch_location, dir="data/eos/stagger/");
    #const eos_legacy = dispatch.EOS.stagger(top=MUST.dispatch_location, file="table.dat");
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

# Use the new Squaregas EOS 
eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
@show eos_path
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

            # Check if Teff should be added
            if do_teff
                @info "setting teff for sn $(snapshots[i_s])"
                MUST.set!(s.parameter, teff, ini_nml)
            end

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

    # Execute the function convert_snapshots on the passed split
    append!(futures, [@spawnat wrk[i] convert_snapshots(split_l)])
end

# Wait for the execution to finish
MUST.wait(futures)




