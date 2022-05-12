"""Submit dispatch job without waiting."""
function srun_dispatch(nml_name, threads::Int, memMB::Int, timeout::String; wait=false)
    ddir    = MUST.@in_dispatch("")
    command = `srun -N 1 -n 1 -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./dispatch.x $(nml_name)`
    #run(command, wait=wait)
    run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"), 
                          stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
end

"""
Run a MUSTGrid from start to end.
    TODO: This can dispatch on RestartMUSTGrid to start the tasks 
    async. To do this one should add a flag, that decides if 
    one should wait for the code to exit. If not, then in the RestartMUSTGrid
    each namelist can check Base.process_exited() before running it. For this the 
    RestartMUSTGrid shoudl get the running process list.
"""
function run!(grid::AbstractMUSTGrid; threads::Int=12, memMB::Int=1333, timeout="00:40:00",submission_pause=1)
    allowed_nmls = allowed_namelists(grid)

    results = []
    for name in allowed_nmls
        append!(results, [srun_dispatch(name, threads, memMB, timeout)])
        
        @info "$(name) submitted."
        sleep(submission_pause)
    end

    # wait for completion
    for (i,r)  in enumerate(results)
        s = success(r)
        @info "$(allowed_nmls[i]) finished with success status $(s)."
    end

    # check success
    grid.info[!,"$(grid.name)_success"] = check_success(grid)
end

allowed_namelists(grid::AbstractMUSTGrid) = grid.info[!,"$(grid.name)_name"]
allowed_namelists(grid::RestartMUSTGrid)  = grid.info[grid.info[!,Symbol("$(grid.from_name)_success")],"$(grid.name)_name"]

"""
Generic check_success.
    Check if the simulation has ended based on the expected number of snapshots
"""
check_success(grid::AbstractMUSTGrid) = begin
    nml_names   = grid.info[!,"$(grid.name)_name"]
    nml_folders = [MUST.@in_dispatch(joinpath("data",split(n, ".nml")[1]))  for n in nml_names] 
    nmls        = [MUST.StellarNamelist(MUST.@in_dispatch(n)) for n in nml_names]

    suc = [check_success(grid, n, f) for (n,f) in zip(nmls,nml_folders)]
end

function check_success(grid::AbstractMUSTGrid, nml::AbstractNamelist, output_folder::String)
    # Check if the simulation reached close to the end
    end_point = floor(nml.io_params["end_time"] / nml.io_params["out_time"])

    content_of_folder = glob("*/", output_folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))

    length(snapshots) -2 >= end_point
end

function check_success(grid::RestartMUSTGrid, nml::AbstractNamelist, output_folder::String)
    nml_name  = split(output_folder, "/")[end]*".nml"
    idx       = findfirst(grid.info[!,"$(grid.name)_name"] .== nml_name)

    if !grid.info[idx,"$(grid.from_name)_success"] 
        return false
    end

    i_restart = 0.0
    try
        i_restart = first(get_restart_snap_nml(nml).snapshot_nml["time"])
    catch
        @warn "Restart Snapshot not found for $(nml_name) (after checking for $(grid.from_name) success)."
        return false
    end

    # Check if the simulation reached close to the end
    end_point = floor((nml.io_params["end_time"]-i_restart) / nml.io_params["out_time"])

    content_of_folder = glob("*/", output_folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))

    length(snapshots) -2 >= end_point
end

slurm_setup() = begin
    @assert "SLURM_CPUS_PER_TASK" in keys(ENV)
    @assert "SLURM_NTASKS_PER_NODE" in keys(ENV)
    @assert "SLURM_MEM_PER_NODE" in keys(ENV)

    threads, tasks, mem = parse.(Int, [ENV["SLURM_CPUS_PER_TASK"], ENV["SLURM_NTASKS_PER_NODE"], ENV["SLURM_MEM_PER_NODE"]])
    mem = floor(Int, mem / tasks)
    (threads, tasks, mem)
end
