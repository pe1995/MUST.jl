using Pkg; Pkg.activate(".")
using MUST
using CSV
using DataFrames
MUST.@import_dispatch "../../../dispatch2_clean/dispatch2/"

arguments   = MUST.MUSTGridArgs()
ngrid       = last(arguments[:index_range]) - (first(arguments[:index_range]) -1)   # How large should the grid be
start_index = first(arguments[:index_range])                                        # Where the namelist naming should start
                                                                                    #   this is important if more than 1 instance is running 
info_path   = "$(arguments[:info_path]).csv"                                        # Where the info should be stored    

# Slurm setup
threads, tasks, mem = MUST.slurm_setup()
@info threads tasks mem

# One may choose any default namelist as a template for
# the others to copy. We make sure that the namelists are
# correctly setup for our needs. If nothing is given, stellar.nml and restart.nml are
# assumed to be the default setting
phase1_nml_template = MUST.StellarNamelist(MUST.@in_dispatch("p10105.nml"))
MUST.set!(phase1_nml_template, io_params=("end_time" => 120.0,
                                          "out_time" => 5.0))

#========== PHASE 1 ==========#
# Create a random grid from central seeds
#=central_seeds = Dict{Symbol, Dict{String, Float32}}(   
                        :stellar_params => Dict{String, Float32}("tt_k" => 11000.0, "d_cgs" => log10(3.0445492234271924e-7)))

relative_lims = Dict{Symbol, Dict{String, Float32}}(   
                        :stellar_params => Dict{String, Float32}("tt_k" => 0.1, "d_cgs" => 0.005))
=#
central_seeds = Dict{Symbol, Dict{String, Float32}}(   
                            :stellar_params => Dict{String, Float32}("tt_k" => 12000.0, "d_cgs" => log10(4.5e-7)))
    
relative_lims = Dict{Symbol, Dict{String, Float32}}(   
                            :stellar_params => Dict{String, Float32}("tt_k" => 0.1, "d_cgs" => 0.07))

# The grid can be constructed using e.g. a RangeMUSTGrid
phase1_grid = MUST.RangeMUSTGrid(central_seeds, relative_lims, phase="phase1")
MUST.randomgrid!(phase1_grid, ngrid)

# Anything that is log spaced can now be turned into linear again
MUST.modify!(phase1_grid, :stellar_params, "d_cgs") do x
    10.0 .^(x)
end

# Create the namelists of this grid (in the dispatch folder)
MUST.create_namelists!(phase1_grid, default_nml=phase1_nml_template, start_index=start_index)

# Save the intermediate results
CSV.write(info_path, phase1_grid.info)

# Run the first phase and check for completion
MUST.run!(phase1_grid, threads=threads, memMB=mem, timeout="03:00:00")

# Save the intermediate results
CSV.write(info_path, phase1_grid.info)

@info "Phase 1 completed."