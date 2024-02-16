using Pkg; Pkg.activate("."); 
using MUST
using Distributed, SlurmClusterManager
using Glob
using ProgressMeter
PythonCall = MUST.PythonCall

# add the workers either with slurm or interactive
if "SLURM_NTASKS" in keys(ENV)
    addprocs(SlurmManager())
    for i in workers()
        host,pid = fetch(@spawnat i (gethostname(), getpid()))
        @info "Worker $(i) is running on $(host) with ID $(pid)" 
    end
else
    @warn "No Slurm environment detected. Using default addprocs."
    addprocs(3)
end

@everywhere begin
    using Pkg; Pkg.activate(".")
    using MUST
    using Distributed
    using Glob
    using ProgressMeter
    PythonCall = MUST.PythonCall
end

# make dispatch available on all workers
@everywhere begin
    MUST.@import_dispatch "../../../dispatch2"
end

# select the last N snapshots
folder = MUST.@in_dispatch ARGS[1]
content_of_folder = glob("*/", folder)
snapshots = sort(MUST.list_of_snapshots(content_of_folder))

snapshots = if length(ARGS) > 1
    istart = length(snapshots) - 1 - parse(Int, ARGS[2])
     snapshots[istart:end-2]
else
    snapshots[1:end-2]
end

# make this available to all workers
MUST.sendsync(workers(), dir=folder)

# skip already converted snapshots
csnaps = MUST.list_snapshots(MUST.converted_snapshots(folder))
mask = [!(s in csnaps) for s in snapshots]
snapshots = snapshots[mask]

println("")
println("====================================================================")
if count(.!mask) > 0
    @info "$(count(.!mask)) snapshots skipped because they are already converted."
end
@info "conversion started..."

# convert the snapshots
@showprogress pmap(snapshots) do snap
    MUST.snapshotBox(
        snap; 
        folder=dir, 
        optical_depth_scale=true, 
        save_snapshot=true, 
        add_selection=false, 
        to_multi=false,
        is_box=true, 
        use_mmap=false
    )

    nothing
end
@info "...conversion finished."
println("====================================================================")
