using Pkg; Pkg.activate(".")
using MUST
using CSV
using DataFrames
MUST.@import_dispatch "../../../dispatch2/"

# How large should the grid be
ngrid = parse(Int, ARGS[1])

# Slurm setup
threads, tasks, mem = MUST.slurm_setup()
@info threads tasks mem

# One may choose any default namelist as a template for
# the others to copy. We make sure that the namelists are
# correctly setup for our needs. If nothing is given, stellar.nml and restart.nml are
# assumed to be the default setting
phase1_nml_template = MUST.StellarNamelist(MUST.@in_dispatch("stellar.nml"))
phase2_nml_template = MUST.StellarNamelist(MUST.@in_dispatch("restart.nml"))

MUST.set!(phase1_nml_template, io_params=("print_seconds" => 120,
                                          "print_every"   => 120,
                                          "end_time"      => 50.0,
                                          "out_time"      => 5.0))
MUST.set!(phase2_nml_template, io_params=("print_seconds" => 120,
                                          "print_every"   => 120,
                                          "end_time"      => 80.0,
                                          "out_time"      => 5.0),
                                restart_params=("snapshot" => 8,))

#========== PHASE 1 ==========#
# Create a random grid from central seeds
central_seeds = Dict{Symbol, Dict{String, Float32}}(   
                        :stellar_params => Dict{String, Float32}("tt_k" => 9000.0, "d_cgs" => log(1.6397e-7), "ee_min" => 5.5), 
                        :newton_params  => Dict{String, Float32}("ee0"=>5.3)
                    )

relative_lims = Dict{Symbol, Dict{String, Float32}}(   
                        :stellar_params => Dict{String, Float32}("tt_k" => 0.3, "d_cgs" => 0.1, "ee_min" => 0.0001), 
                        :newton_params  => Dict{String, Float32}("ee0"=>0.0001)
                    )

# The grid can be constructed using e.g. a RangeMUSTGrid
phase1_grid = MUST.RangeMUSTGrid(central_seeds, relative_lims)
MUST.randomgrid!(phase1_grid, ngrid)

# Anything that is log spaced can now be turned into linear again
MUST.modify!(phase1_grid, :stellar_params, "d_cgs") do x
    exp(x)
end

# Create the namelists of this grid (in the dispatch folder)
MUST.create_namelists!(phase1_grid, phase="phase1", default_nml=phase1_nml_template)

# Save the intermediate results
CSV.write("summary.csv", phase1_grid.info)

# Run the first phase and check for completion
MUST.run!(phase1_grid, threads=threads, memMB=mem, timeout="00:40:00")

# Save the intermediate results
CSV.write("summary.csv", phase1_grid.info)

@info "Phase 1 completed."

#========== PHASE 2 ==========#
# In the second phase take the first phase grid and add the restart information
# Also take over the grid parameters to the new grid. If you dont want the 
# new namelists to have the parameters of the old namelist you can leave out the last parameter.
phase2_grid = MUST.RestartMUSTGrid(phase1_grid, from_phase="phase1", grid=phase1_grid.grid)

# Also create the namelists for this grid
MUST.create_namelists!(phase2_grid, phase="phase2", default_nml=phase2_nml_template)

# Save the intermediate results
CSV.write("summary.csv", phase2_grid.info)

# Run the first phase and check for completion
MUST.run!(phase2_grid, threads=threads, memMB=mem, timeout="03:00:00")

@info "Phase 2 completed."

#========== COMPLETION ==========#
# One can add as many phases as needed.
# The last phase contains the info from all phases
CSV.write("summary.csv", phase2_grid.info)