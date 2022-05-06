# Creating a grid of Namelists for execution
using MUST

central_seeds = Dict("tt_K" => 0.6204e4,       "d_cgs" => 1.6397e-7,        "ee_min" => 5.5)
limits        = Dict(p => [central_seeds[p]-0.1*central_seeds[p], central_seeds[p]+0.1*central_seeds[p]] for p in keys(central_seeds))
grid          = Dict(p => MUST.randrange(limits[p]) for p in keys(central_seeds))
