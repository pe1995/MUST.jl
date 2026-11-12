# ============================================================================
# Methods for submitting dispatch jobs via a slurm submission
# ============================================================================

function srun_dispatch(nml_name; threads::Int, memMB::Int, timeout::String, wait=false)
    ddir = @in_dispatch("")
    #run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"),

    exec_name = "$(nml_name).export.sh"
	exec_path = @in_dispatch("$(nml_name).export.sh")

	# Do an export script before
	open(exec_path, "w") do f
		write(f, "#!/bin/bash\n")
		write(f, "export KMP_NUM_THREADS=$(threads)\n")
		write(f, "export OMP_NUM_THREADS=$(threads)\n")
        write(f, "export KMP_STACKSIZE=256M\n")
        write(f, "export OMP_STACKSIZE=256M\n")
		write(f, "./dispatch.x $(nml_name)")
	end

	run(`chmod +x $(exec_path)`)
    command = `srun -N 1 -n 1 -J m3dis -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./$(exec_name)`

    run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"),
                          stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
end

function run_dispatch(nml_name; threads=70, wait=true, ddir=@in_dispatch(""))

    exec_name = "$(nml_name).export.sh"
    exec_path = @in_dispatch("$(nml_name).export.sh")

    # Do an export script before
    open(exec_path, "w") do f
    write(f, "#!/bin/bash\n")
    write(f, "export KMP_NUM_THREADS=$(threads)\n")
    write(f, "export OMP_NUM_THREADS=$(threads)\n")
    write(f, "export KMP_STACKSIZE=256M\n")
    write(f, "export OMP_STACKSIZE=256M\n")
    write(f, "./dispatch.x $(nml_name)")
    end

    run(`chmod +x $(exec_path)`)
    currdir = pwd()
    cd(ddir)
    command = `./$(exec_name)`
    r = run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"),
            stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
    cd(currdir)

    r
end

function srun_m3dis(nml_name;
					threads=70, memMB=100000, timeout="15:00:00", wait=false)
	ddir = @in_m3dis("")

	exec_name = "$(nml_name).export.sh"
	exec_path = @in_m3dis("$(nml_name).export.sh")

	# Do an export script before
	open(exec_path, "w") do f
		write(f, "#!/bin/bash\n")
		write(f, "export KMP_NUM_THREADS=$(threads)\n")
		write(f, "export OMP_NUM_THREADS=$(threads)\n")
        write(f, "export KMP_STACKSIZE=256M\n")
        write(f, "export OMP_STACKSIZE=256M\n")
		write(f, "./dispatch.x $(nml_name)")
	end

	run(`chmod +x $(exec_path)`)
    command = `srun -N 1 -n 1 -J m3dis -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./$(exec_name)`

    run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"),
                          stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
end

function run_m3dis(nml_name; threads=70, wait=true, ddir=@in_m3dis(""))

    exec_name = "$(nml_name).export.sh"
    exec_path = @in_m3dis("$(nml_name).export.sh")

    # Do an export script before
    open(exec_path, "w") do f
    write(f, "#!/bin/bash\n")
    write(f, "export KMP_NUM_THREADS=$(threads)\n")
    write(f, "export OMP_NUM_THREADS=$(threads)\n")
    write(f, "export KMP_STACKSIZE=256M\n")
    write(f, "export OMP_STACKSIZE=256M\n")
    write(f, "./dispatch.x $(nml_name)")
    end

    run(`chmod +x $(exec_path)`)
    currdir = pwd()
    cd(ddir)
    command = `./$(exec_name)`
    r = run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"),
            stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
    cd(currdir)

    r
end

# ============================================================================
# Grid related running
# ============================================================================

function run!(grid::AbstractMUSTGrid; threads::Int=12, memMB::Int=1333, timeout="00:40:00",submission_pause=1, slurm=true, use_status=true)
    allowed_nmls = allowed_namelists(grid)

    if !slurm
        @warn "Submitting jobs without memory management."
    end

    results = []
    for name in allowed_nmls
        r = if slurm
            srun_dispatch(name, threads=threads, memMB=memMB, timeout=timeout)
        else
            run_dispatch(name, threads=threads)
        end
        append!(results, [r])

        @info "$(name) submitted."
        sleep(submission_pause)
    end

    # wait for completion
    succ = []
    for (i,r)  in enumerate(results)
        s = success(r)
        append!(succ, [s])
        @info "$(allowed_nmls[i]) finished with success status $(s)."
    end

    # check success
    grid.info[!,"$(grid.name)_success"] = use_status ? succ : check_success(grid)
end

allowed_namelists(grid::AbstractMUSTGrid) = grid.info[!,"$(grid.name)_name"]
allowed_namelists(grid::RestartMUSTGrid)  = grid.info[grid.info[!,Symbol("$(grid.from_name)_success")],"$(grid.name)_name"]

check_success(grid::AbstractMUSTGrid) = begin
    nml_names   = allowed_namelists(grid)
    nml_folders = [@in_dispatch(joinpath("data", split(n, ".nml")[1]))  for n in nml_names]
    nmls        = [StellarNamelist(@in_dispatch(n)) for n in nml_names]

    [check_success(grid, n, f) for (n,f) in zip(nmls,nml_folders)]
end

function check_success(grid::AbstractMUSTGrid, nml::AbstractNamelist, output_folder::String)
    # Check if the simulation reached close to the end
    end_point = floor(nml.io_params["end_time"] / nml.io_params["out_time"])

    content_of_folder = glob("*/", output_folder)
    snapshots         = sort(list_of_snapshots(content_of_folder))

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
