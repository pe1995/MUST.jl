# Creating a grid of Namelists for execution
using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using MUST
MUST.@import_dispatch "../../../dispatch2/"


function create_namelists(ngrid::Int)
    central_seeds = Dict("tt_K" => 0.6204e4, "d_cgs" => 1.6397e-7, "ee_min" => 5.5, "ee0"=>5.3)
    relative_lims = Dict("tt_K" => 0.2, "d_cgs" => 0.2, "ee_min" => 0.05, "ee0" => 0.05)
    limits        = Dict(p => [central_seeds[p]-relative_lims[p]*central_seeds[p], central_seeds[p]+relative_lims[p]*central_seeds[p]] for p in keys(central_seeds))
    grid          = Dict(p => MUST.randrange(limits[p]..., Float32, ngrid) for p in keys(central_seeds))

    dummy_stellar = MUST.StellarNamelist(MUST.@in_dispatch("stellar_marcs.nml"))
    dummy_restart = MUST.StellarNamelist(MUST.@in_dispatch("restart_marcs.nml"))

    summary = Dict("phase1_name"=>[], "phase2_name"=>[], grid...)

    for i in 1:ngrid
        d = grid["d_cgs"][i]; t = grid["tt_K"][i]; em = grid["ee_min"][i]; e0 = grid["ee0"][i]
        
        # name of namelists
        id = "grid$(i)"
        name_phase1 = "$(id)_phase1.nml"
        name_phase2 = "$(id)_phase2.nml"

        # copy the dummy namelist
        phase1 = deepcopy(dummy_stellar)
        phase2 = deepcopy(dummy_restart)

        # Set the values
        phase1.stellar_params["tt_K"]   = t
        phase1.stellar_params["d_cgs"]  = d
        phase1.stellar_params["ee_min"] = em
        phase1.newton_params["ee0"]     = e0

        phase2.restart_params["run"]    = "'"* String(name_phase1[1:first(findlast(".nml", name_phase1))-1]) *"'"
        phase2.newton_params["ee0"]     = e0

        # write the namelists
        MUST.write(phase1, MUST.@in_dispatch(name_phase1))
        MUST.write(phase2, MUST.@in_dispatch(name_phase2))

        # save the summary
        append!(summary["phase1_name"], [name_phase1])
        append!(summary["phase2_name"], [name_phase2])
    end

    DataFrame(summary)
end

function run_namelists(summary; slurm_script="run_10threads.slurm", slurm_script_large="run_16threads.slurm",user="peitner", 
                                job_name="MUST", job_limit=300, interval=30)
    curr_dir = pwd()
    cd(MUST.@in_dispatch(""))

    slurm_ids = Int[]
    slurm_ids_phase2 = zeros(nrow(summary))
    completed_phase1 = []
    running_phase2   = 0
    success_phase1   = [false for _ in 1:nrow(summary)]
    success_phase2   = [false for _ in 1:nrow(summary)]

    for i in 1:nrow(summary)
        # Limit the total number of tasks
        limit_tasks(user, job_name, job_limit, interval)

        name    = summary[i,Symbol("phase1_name")]
        command = `sbatch $(slurm_script) $(name)`
        append!(slurm_ids, [ parse(Int, String( strip(split(read(command, String))[end]) ) ) ])
        @info "Submitted $(name) ( $(slurm_ids[i]) )"

        sleep(2)
    end

	
    CSV.write("summary.csv", summary)

    # Execute the second phase if one of the jobs was successful
    # Now that every task has been submitted this means that they are 
    # either present in the job_list or already done
    while running_phase2 < nrow(summary)
        jobs_running = running_tasks(user, job_name)
        @info "Running Namelists: \n"*
              "           total number of jobs : $(nrow(summary))\n"*
              "           jobs running         : $(length(jobs_running))\n"*
              "           completed Phase 1    : $(length(completed_phase1))\n"*
              "           started Phase 2      : $(running_phase2)\n"

        for i in 1:length(slurm_ids)
            (slurm_ids[i] in jobs_running)     ? continue : nothing                    # Skip if still running
            (slurm_ids[i] in completed_phase1) ? continue : nothing                    # Skip if already done
            
            success_phase1[i] = check_success(summary[i,Symbol("phase1_name")])        # it is finished, so check success
            append!(completed_phase1, slurm_ids[i])

            running_phase2 +=1

            if success_phase1[i] 
                # Limit the total number of tasks
            	limit_tasks(user, job_name, job_limit, interval)

                name    = summary[i,Symbol("phase2_name")]
                command = `sbatch $(slurm_script_large) $(name)`
                slurm_ids_phase2[i] = parse(Int, String( strip(split(read(command, String))[end]) ) )
                @info "Submitted $(name) ( $(slurm_ids_phase2[i]) )"
            else
                slurm_ids_phase2[i] = 0
                success_phase2[i]   = false
            end
        end

        sleep(interval)
    end

    # Wait until all jobs are submitted
    wait_for_completion(;user=user, job_name=job_name, interval=interval)

    # Check success of all
    success_phase2 = check_success.(summary[:,Symbol("phase2_name")])

    summary[!,:phase1_success] = success_phase1
    summary[!,:phase2_success] = success_phase2
    summary[!,:phase1_sid] = slurm_ids
    summary[!,:phase2_sid] = slurm_ids_phase2

    cd(curr_dir)
end

function check_success!(summary; phase)
    success = rand(Bool, nrow(summary))

    for row in 1:nrow(summary)
        name      = summary[row,Symbol("phase$(phase)_name")]
        nml       = MUST.StellarNamelist(MUST.@in_dispatch(name))
        
        name      = String(name[1:first(findlast(".", name))-1])
        output_at = MUST.@in_dispatch "data/$(name)"

        success[row] = !isdir(output_at) ? false : MUST.check_success(nml, output_at)
    end

    summary["phase$(phase)_success"] = success
end

function check_success(name)
    nml       = MUST.StellarNamelist(MUST.@in_dispatch(name))
    name      = String(name[1:first(findlast(".", name))-1])
    output_at = MUST.@in_dispatch "data/$(name)"
    !isdir(output_at) ? false : MUST.check_success(nml, output_at)
end

function limit_tasks(user="peitner", job_name="MUST", job_limit=300, interval=30)
    n_running = n_running_tasks(user, job_name)
    while n_running >= job_limit -1
        @warn "Too Many Jobs. Waiting for completion. ( $(n_running) / $(job_limit) )"
        sleep(interval)
    end
end

function wait_for_completion(; user="peitner", job_name="MUST", interval=30)
    n_total   = n_running_tasks(user, job_name)
    n_running = n_total
    elapsed   = 0
    while n_running > 0
        @info "All Jobs submitted. Waiting for completion: $(n_running) / $(n_total)     [ $(elapsed)s ]"
        
        sleep(interval)
        elapsed += interval
        
        n_running = n_running_tasks(user, job_name)
    end

    @info "All jobs finished."
end

function running_tasks(user="peitner", job_name="MUST")
    command = `squeue -u $(user) --name $(job_name)`
    jobs = split(Base.read(command, String), "\n")[2:end-1]
    jobs = length(jobs) == 0 ? Int[] : parse.(Int, String.( first.( split.( strip.(jobs), " "))))
end

n_running_tasks(args...; kwargs...) = length(running_tasks(args...; kwargs...))


#=== Execution ====#

# Number of namelists to create
#N_namelists = parse(Int,ARGS[1])

# Create the namelists
#summary = create_namelists(N_namelists)

# execute phase 1 and 2 using slurm
#run_namelists(summary) 

# Save the summary
#CSV.write("summary.csv", summary)

summary = DataFrame(CSV.File("summary.csv"))
summary[!,"phase2_success"] = check_success.(summary[:,"phase2_name"])

CSV.write("summary.csv", summary)