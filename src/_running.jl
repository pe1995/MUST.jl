## Methods for submitting dispatch jobs via a slurm submission 

"""
    srun_dispatch(nml_name, threads::Int, memMB::Int, timeout::String; wait=false)

Submit dispatch job.
"""
function srun_dispatch(nml_name; threads::Int, memMB::Int, timeout::String, wait=false)
    ddir = @in_dispatch("")
    #command = `srun -N 1 -n 1 -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./dispatch.x $(nml_name)`
    #run(command, wait=wait)
    #run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"), 
    #                      stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)

    exec_name = "$(nml_name).export.sh"
	exec_path = @in_dispatch("$(nml_name).export.sh")
	
	# Do an export script before
	open(exec_path, "w") do f
		write(f, "#!/bin/bash\n")
		write(f, "source ~/.bashrc\n")
		write(f, "export KMP_NUM_THREADS=$(threads)\n")
		write(f, "export OMP_NUM_THREADS=$(threads)\n")
		write(f, "./dispatch.x $(nml_name)")
	end
	
	run(`chmod +x $(exec_path)`)
    command = `srun -N 1 -n 1 -J m3dis -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./$(exec_name)`
    
    run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"), 
                          stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
end

"""
    run_dispatch(nml_name; threads=70, wait=true)

Similar `to srun_dispatch()`, but running in the shell directly.
"""
function run_dispatch(nml_name; threads=70, wait=true, ddir=@in_dispatch(""))

    exec_name = "$(nml_name).export.sh"
    exec_path = @in_dispatch("$(nml_name).export.sh")

    # Do an export script before
    open(exec_path, "w") do f
    write(f, "#!/bin/bash\n")
    write(f, "source ~/.bashrc\n")
    write(f, "export KMP_NUM_THREADS=$(threads)\n")
    write(f, "export OMP_NUM_THREADS=$(threads)\n")
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








"""
	srun_m3dis(nml_name; threads=70, memMB=200000, timeout="05:00:00", wait=false)

Submit a job corresponding to the given namelist and slurm setup. Optinally wait for the job to finish. If this function is executed within a slurm environment, it will create a job step within that environment according to the given parameters. Note that in this case wait should be set to false, otherwise the current thread will be halted until the job step exits. If this function is executed within a e.g. notebook, you should set wait=true, so that it waits until it receivs a finishing signinal from the slurm process. 
The slurm job will be run inside the loaded Multi environment, so make sure you have a.) created the namelist, and b) loaded multi correctly using the `@import_m3dis` macro.
"""
function srun_m3dis(nml_name; 
					threads=70, memMB=100000, timeout="15:00:00", wait=false)
	ddir = @in_m3dis("")
	
	exec_name = "$(nml_name).export.sh"
	exec_path = @in_m3dis("$(nml_name).export.sh")
	
	# Do an export script before
	open(exec_path, "w") do f
		write(f, "#!/bin/bash\n")
		write(f, "source ~/.bashrc\n")
		write(f, "export KMP_NUM_THREADS=$(threads)\n")
		write(f, "export OMP_NUM_THREADS=$(threads)\n")
		write(f, "./dispatch.x $(nml_name)")
	end
	
	run(`chmod +x $(exec_path)`)
    command = `srun -N 1 -n 1 -J m3dis -c $(threads) --mem=$(memMB)mb --time=$(timeout) --exclusive -D $(ddir) ./$(exec_name)`
    
    run(pipeline(command, stdout=joinpath(ddir, "$(nml_name).log"), 
                          stderr=joinpath(ddir, "$(nml_name).err")), wait=wait)
end

"""
	run_m3dis(nml_name; threads=70, wait=true)

Similar `to srun_m3dis()`, but running in the shell directly.
"""
function run_m3dis(nml_name; threads=70, wait=true, ddir=@in_m3dis(""))

    exec_name = "$(nml_name).export.sh"
    exec_path = @in_m3dis("$(nml_name).export.sh")

    # Do an export script before
    open(exec_path, "w") do f
    write(f, "#!/bin/bash\n")
    write(f, "source ~/.bashrc\n")
    write(f, "export KMP_NUM_THREADS=$(threads)\n")
    write(f, "export OMP_NUM_THREADS=$(threads)\n")
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










#= Grid related running =#

"""
Run a MUSTGrid from start to end.
    TODO: This can dispatch on RestartMUSTGrid to start the tasks 
    async. To do this one should add a flag, that decides if 
    one should wait for the code to exit. If not, then in the RestartMUSTGrid
    each namelist can check Base.process_exited() before running it. For this the 
    RestartMUSTGrid shoudl get the running process list.
"""
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

"""
Generic check_success.
    Check if the simulation has ended based on the expected number of snapshots
"""
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

