"""Submit dispatch job without waiting."""
function srun_dispatch(nml_name, threads, memMB; wait=false)
    ddir    = MUST.@in_dispatch("")
    command = `srun -N 1 -n 1 -c $(threads) --mem-per-cpu=$(memMB)mb --exclusive -D $(ddir) ./dispatch.x $(nml_name)`
    run(command, wait=wait)
end

function run!(grid::AbstractMUSTGrid; threads=12, memMB=1333, submission_pause=5)
    allowed_nmls = allowed_namelists(grid)

    results = []
    for name in allowed_nmls
        append!(results, [srun_dispatch(name, threads, memMB)])
        
        @info "$(name) submitted."
        sleep(submission_pause)
    end

    # wait for completion
    for (i,r)  in enumerate(results)
        s = success(r)
        @info "$(nml_names[i]) finished with success status $(s)."
    end

    # check success
    grid.info[!,"$(grid.name)_success"] = check_success(grid)
end

allowed_namelists(grid::AbstractMUSTGrid) = grid.info[!,"$(grid.name)_name"]
allowed_namelists(grid::RestartMUSTGrid)  = grid.info[grid.info[Symbol("$(grid.from_name)_success")],"$(grid.name)_name"]

"""
Generic check_success.
    Check if the simulation has ended based on the expected number of snapshots
"""
check_success(grid::AbstractMUSTGrid) = begin
    nml_names   = grid.info[!,"$(grid.name)_name"]
    nml_folders = [MUST.@in_dispatch(joinpath("data",split(n, ".nml")[1]))  for n in nml_names] 
    nmls        = [MUST.StellarNamelist(MUST.@in_dispatch(n)) for n in nml_names]

    suc = [check_success(grid, n, f) for (n,f) in zip(nml_names,nml_folders)]
end

function check_success(grid::AbstractMUSTGrid, nml::AbstractNamelist, output_folder::String)
    # Check if the simulation reached close to the end
    end_point = floor(nml.io_params["end_time"] / nml.io_params["out_time"])

    content_of_folder = glob("*/", output_folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))

    @show length(snapshots) end_point
    length(snapshots) -2 >= end_point
end

function check_success(grid::RestartMUSTGrid, nml::AbstractNamelist, output_folder::String)
    i_restart = first(get_restart_snap_nml(nml).snapshot_nml["time"])

    # Check if the simulation reached close to the end
    end_point = floor((nml.io_params["end_time"]-i_restart) / nml.io_params["out_time"])

    content_of_folder = glob("*/", output_folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))

    length(snapshots) -2 >= end_point
end

slurm_setup() = begin
    threads, tasks, mem = parse.(Int, [ENV["SLURM_CPUS_PER_TASK"], ENV["SLURM_NTASKS_PER_NODE"], ENV["SLURM_MEM"]])
    mem = floor(mem / tasks / threads)
    (threads, tasks, mem)
end
