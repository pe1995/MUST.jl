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
        dispatch_location = "/u/peitner/DISPATCH/dispatch2/"
        #initial_grid_path = "stagger_grid.mgrid"
        #final_grid_path   = "dispatch_grid.mgrid"
        initial_grid_path = "random_grid.mgrid"
        initial_cl_path   = "random_grid_avail.mgrid"
        initial_mod_path  = "random_grid_solar.mgrid"
        final_grid_path   = "random_models.mgrid"
        mother_table_path = "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6"
        extension         = "magg22"
        version           = "v0.1"
        Nbins             = 8
        clean             = true
        use_adiabat       = true
    elseif host == "gemini"
        name_extension    = "DIS_MARCS"
        dispatch_location = "/home/eitner/shared/model_grid/dispatch2"

        initial_grid_path = "stagger_grid.mgrid"
        initial_cl_path   = "stagger_grid_avail.mgrid"
        initial_mod_path  = "stagger_grid_solar.mgrid"
        final_grid_path   = "dispatch_grid.mgrid"
        #initial_grid_path = "random_setup.mgrid"
        #final_grid_path   = "random_grid.mgrid"
        #initial_grid_path = "node_setup.mgrid"
        #final_grid_path   = "node_grid.mgrid"

        #mother_table_path = "/home/eitner/shared/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_v1.6"
        #mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_asplund_m0_a0_v1.6"
        mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.6"
        #mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.5"
       
        extension         = "magg_m0_a0"
        version           = "v0.4"
        Nbins             = 10
        clean             = false
        use_adiabat       = false
    elseif host == "cloud"
        name_extension    = "DIS_MARCS"
        dispatch_location = "/home/ubuntu/DISPATCH/dispatch2"

        initial_grid_path = "random_grid.mgrid"
        initial_cl_path   = "random_grid_avail.mgrid"
        initial_mod_path  = "random_grid_solar.mgrid"
        final_grid_path   = "random_models.mgrid"
        #initial_grid_path = "random_setup.mgrid"
        #final_grid_path   = "random_grid.mgrid"
        #initial_grid_path = "node_setup.mgrid"
        #final_grid_path   = "node_grid.mgrid"

        mother_table_path = "/home/ubuntu/DISPATCH/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.7"
        extension         = "magg_m0_a0"
        version           = "v0.4"
        Nbins             = 8
        clean             = false
        use_adiabat       = true
    end

    MUST.@import_dispatch "../../../dispatch2"
    #MUST.@import_dispatch dispatch_location
end

# import the relevant functions
@everywhere prepare4dispatch = MUST.ingredients("prepare4dispatch.jl")


#=================== Step (A): The Initial Grid ==============================#
begin
    grid = MUST.StaggerGrid(initial_grid_path)
    deleteat!(grid.info, .!isfile.(grid["av_path"]))
    MUST.save(grid, initial_cl_path)

    deleteat!(grid.info, grid["feh"].!=0.0)
    MUST.save(grid, initial_mod_path)


    ## Check for opacity table field
    if !("eos_root" in names(grid.info))
        @warn "The given grid does not contain information about the root EoS and opacity source. Adding the fallback table."
        grid.info[!, "eos_root"] = [mother_table_path for _ in 1:nrow(grid.info)]
    end
end

#================ Step (B): Fomation opacities (parallel) ====================#
begin
    args = [
        (
            grid.info[i, "eos_root"], 
            grid.info[i, "av_path"], 
            grid.info[i, "name"],
            grid.info[i, "logg"],
            extension,
            :kmeans, 
            Nbins,
            false #i==1 ? true : false
        ) for i in 1:nrow(grid.info)
    ]

    ## To save disk space, the formation opacities can be deleted after binning
    ## For this it is required to run this sequentially
    @everywhere formation_and_bin(args) = begin
        eos_root = args[1]
        av_path = args[2]
        name = args[3]
        logg = args[4]
        extension = args[5]
        method = args[6]
        Nbins = args[7]
        do_ross = args[8]


        @info "Opacity Binning: $(name)"

        ## formation opacities
        prepare4dispatch.formation_opacities(
            eos_root, av_path, name; 
            logg=logg, extension=extension, do_ross=do_ross
        )

        qlim = round(
            prepare4dispatch.quadrantlimit(eos_root, name, extension=extension, λ_lim=5.0), 
            sigdigits=3
        )

        quadrants = [ 
            TSO.Quadrant((0.0, 4.0), (qlim, 4.5), 2, stripes=:κ),
            TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
            TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (-100, qlim), 4, stripes=:λ),
        ]
           
        ## bin opacities
        prepare4dispatch.bin_opacities(
            eos_root, av_path, name; 
            logg=logg, method=method, Nbins=Nbins, extension=extension,
            version=version,
            quadrants=quadrants
        )

        ## remove the formation opacities
        clean && prepare4dispatch.clean(eos_root, name, extension=extension)
    end

    ## End-to-end binning, with clean-up
    #Distributed.pmap(formation_and_bin, args)
    
    ## Save the eos info
    grid.info[!, "name_extension"]   = [name_extension for _ in 1:nrow(grid.info)]
    grid.info[!, "version"]          = [version for _ in 1:nrow(grid.info)]
    grid.info[!, "binned_tables"]    = [TSO.join_full(name_extension, grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "binned_E_tables"]  = [TSO.join_full(name_extension, "E", grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "rad_bins"]         = [Nbins for _ in 1:nrow(grid.info)]

    ## compute the resolution and the rounded size of the box
    ## use the EoS that was just created for this
    prepare4dispatch.resolution!(grid, patch_size=22, τ_up=-4.0, τ_surf=0.0, τ_down=6.0, scale_resolution=0.7)
end

#====================== Step (C): Conversion =================================#
begin
    prepare4dispatch.fromT_toE.(grid.info[!, "binned_tables"], grid.info[!, "binned_E_tables"], upsample=2048)

    ## Copy the average model in the same folder so that we can link it all to the right place
    for i in 1:nrow(grid.info)
        if !use_adiabat
            cp(grid.info[i, "av_path"], joinpath(grid.info[i, "binned_E_tables"], "inim.dat"), force=true)
        else
            prepare4dispatch.adiabat(
                joinpath(grid.info[i, "binned_E_tables"], "eos.hdf5"),
                grid.info[i, "av_path"], 
                grid.info[i, "logg"],
                saveat=joinpath(grid.info[i, "binned_E_tables"], "inim.dat"),
                ee_min=grid.info[i, "ee_min"],
                common_size=grid.info[i, "initial_model_size"]
            )
        end
    end
end

#=============== Step (D): Input namelists for dispatch ======================#
begin
    prepare4dispatch.create_namelist!(grid)

    for i in 1:nrow(grid.info)
        cp(grid.info[i, "namelist_name"], joinpath(grid.info[i, "binned_E_tables"], "ininml.dat"), force=true)
    end

    MUST.save(grid, final_grid_path)

    # Stage the grid for execution, possible remove other output
    MUST.stage_namelists(grid, clean_namelists=false, clean_logs=false)
end

#=============================================================================#
