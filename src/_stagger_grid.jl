#= Type definitions =#

"""
Information of a part of the Stagger grid. Contains location of
snapshots and averages.
"""
struct StaggerGrid{DF<:DataFrame} <:AbstractMUSTGrid
    name :: String
    info :: DF
end



#= Constructors =#

"""
    StaggerGrid(path)

Load a already computed average Stagger grid. 
The save location and layout of the Stagger grid may
vary from case to case, so only the reading function is 
provided here. For an example case please go to the 
examples folder and look at the `initial_models/from_stagger.ipynb` notebook.

Methods for the interaction with the `TSO.jl` module are located in the same 
folder for now to keep the depencies separate. At a later stage, `TSO.jl` might
be added as a dependecy here.
"""
StaggerGrid(path) = begin
    info = CSV.read(path, DataFrame)
    name = basename(path)[1:findlast('.', basename(path))-1]

    StaggerGrid(name, info)
end



#= Saving =#

save(grid::StaggerGrid, path) = CSV.write(path, grid.info)