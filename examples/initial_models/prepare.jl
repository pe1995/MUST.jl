#=
Preparing initial models from Average 3D Stagger models.
This script contains the complete end-to-end procedure, which includes the following steps:

    (A) Read the prepared Stagger grid, which contains information about the size of the models and
        the initial profile. 
    (B) Use the initial profile to compute formation opacities and bin structure for each model.
        Each model will be added to the name of the mother opacity table, such that they all end
        up in different folders. If there is a prefered mother opacity table, then its path 
        should be given in the grids info field. If there is no such field the fallback table
        will be used.
    (C) The opacity table will be converted to Energy and upsampled if needed.
    (D) The input namelist for dispatch will be created from an example namelist. The resolution will be picked such that
        it agrees with the corresponding resolution of the Stagger model, which is given by the average model. Everything will be linked 
        to the correct folders and the namelist will be copied to dispatch. Everything should be ready now and documented in 
        the grid.info field.
=#

#= Packages =#

using Pkg; Pkg.activate(".")
using MUST
using TSO
using DataFrames
using DelimitedFiles
using Distributed
using SlurmClusterManager

MUST.@import_dispatch "../../../dispatch2"

if "SLURM_NTASKS" in keys(ENV)
    addprocs(SlurmManager())
    for i in workers()
        host,pid = fetch(@spawnat i (gethostname(), getpid()))
        @info "Worker $(i) is running on $(host) with ID $(pid)" 
    end
else
    @warn "No Slurm environment detected. Using default addprocs."
    addprocs(1)
end

@everywhere begin
    using Pkg; Pkg.activate(".")
    using MUST
    using TSO
    using DataFrames
    using DelimitedFiles
    using Distributed
end

@everywhere begin
    host = "gemini"

    if host == "raven"
        name_extension    = "DIS_MARCS"
        version           = "0.1"
        dispatch_location = "/u/peitner/DISPATCH/dispatch2/"
        initial_grid_path = "stagger_grid.mgrid"
        final_grid_path   = "dispatch_grid.mgrid"
        mother_table_path = "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6"
        extension         = "magg22"
        version           = "v0.1"
        Nbins             = 7
    elseif host == "gemini"
        name_extension    = "DIS_MARCS"
        version           = "0.1"
        dispatch_location = "/home/eitner/shared/model_grid/dispatch2"
        initial_grid_path = "stagger_grid.mgrid"
        final_grid_path   = "dispatch_grid.mgrid"
        mother_table_path = "/home/eitner/shared/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_v1.6"
        extension         = "magg22"
        version           = "v0.1"
        Nbins             = 7
    end

    MUST.@import_dispatch dispatch_location
end

# import the relevant functions
@everywhere prepare4dispatch = MUST.ingredients("prepare4dispatch.jl")


#=================== Step (A): The Initial Grid ==============================#
begin
    grid = MUST.StaggerGrid(initial_grid_path)

    ## Check for opacity table field
    if !("eos_root" in names(grid.info))
        @warn "The given grid does not contain information about the root EoS and opacity source. Adding the fallback table."
        grid.info[!, "eos_root"] = [mother_table_path for _ in 1:nrow(grid.info)]
    end
end

#=================== Step (B): Fomation opacities ============================#
begin
    for i in 1:nrow(grid.info)
        prepare4dispatch.formation_opacities(
            grid.info[i, "eos_root"], 
            grid.info[i, "av_path"], 
            grid.info[i, "name"];
            logg=grid.info[i, "logg"],
            extension=extension)
    end

    ## Do the binning in parallel across many workers, this is the most time consuming part
    args = [
        (
            grid.info[i, "eos_root"], 
            grid.info[i, "av_path"], 
            grid.info[i, "name"],
            grid.info[i, "logg"], 
            :kmeans,
            Nbins,
            extension
        ) for i in 1:nrow(grid.info)
    ]
    @everywhere bin_opacities(args) = prepare4dispatch.bin_opacities(args[1:3]...; logg=args[4], method=args[5], Nbins=args[6], extension=args[7])
    Distributed.pmap(bin_opacities, args)    
    
    ## Save the eos info
    grid.info[!, "name_extension"]   = [name_extension for _ in 1:nrow(grid.info)]
    grid.info[!, "version"]          = [version for _ in 1:nrow(grid.info)]
    grid.info[!, "binned_tables"]    = [TSO.join_full(name_extension, grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "binned_E_tables"]  = [TSO.join_full(name_extension, "E", grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "rad_bins"]         = [Nbins for _ in 1:nrow(grid.info)]


    ## compute the resolution and the rounded size of the box
    ## use the EoS that was just created for this
    prepare4dispatch.resolution!(grid, patch_size=20, cut_bottom=0.35)
end

#====================== Step (C): Conversion =================================#
begin
    prepare4dispatch.fromT_toE.(grid.info[!, "binned_tables"], grid.info[!, "binned_E_tables"], upsample=2048)

    ## Copy the average model in the same folder so that we can link it all to the right place
    for i in 1:nrow(grid.info)
        cp(grid.info[i, "av_path"], joinpath(grid.info[i, "binned_E_tables"], "inim.dat"), force=true)
    end
end

#=============== Step (D): Input namelists for dispatch ======================#
begin
    prepare4dispatch.create_namelist!(grid)

    for i in 1:nrow(grid.info)
        cp(grid.info[i, "namelist_name"], joinpath(grid.info[i, "binned_E_tables"], "ininml.dat"), force=true)
    end

    MUST.save(grid, final_grid_path)
end