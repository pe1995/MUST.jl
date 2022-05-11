using Pkg; Pkg.activate(".")
using MUST
using CSV
using DataFrames
MUST.@import_dispatch "../../../dispatch2/"

function create_namelists(ngrid::Int)
    central_seeds = Dict("tt_K" => 11000, "d_cgs" => 1.6397e-7, "ee_min" => 5.5, "ee0"=>5.3)
    relative_lims = Dict("tt_K" => 0.2, "d_cgs" => 0.2, "ee_min" => 0.001, "ee0" => 0.001)
    limits        = Dict(p => [central_seeds[p]-relative_lims[p]*central_seeds[p], central_seeds[p]+relative_lims[p]*central_seeds[p]] for p in keys(central_seeds))
    grid          = Dict(p => MUST.randrange(limits[p]..., Float32, ngrid) for p in keys(central_seeds))

    dummy_stellar = MUST.StellarNamelist(MUST.@in_dispatch("stellar.nml"))
    dummy_restart = MUST.StellarNamelist(MUST.@in_dispatch("restart.nml"))

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

"""Submit dispatch job without waiting."""
function run_dispatch(nml_name, threads, memMB)
    ddir    = MUST.@in_dispatch("")
    command = `srun -N 1 -n 1 -c $(threads) --mem-per-cpu=$(memMB)mb --exclusive -D $(ddir) ./dispatch.x $(nml_name)`
    run(command, wait=false)
end

function run_namelists(nml_names; threads=12, memMB=1333)
    results = []
    for name in nml_names
        append!(results, [run_dispatch(name, threads, memMB)])
        
        @info "$(name) submitted."
        sleep(5)
    end

    # wait for completion
    for (i,r)  in enumerate(results)
        s = success(r)
        @info "$(nml_names[i]) finished with success status $(s)."
    end
end

# the job configurations
threads = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
tasks   = parse(Int, ENV["SLURM_NTASKS"]) - 1
mem     = parse(Int, ENV["SLURM_MEM_PER_CPU"]) 

summary = create_namelists(11)

run_namelists(summary[:,:phase1_name], threads=threads, memMB=mem)

CSV.write("summary.csv", summary)