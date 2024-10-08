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

# Input file with information needed for binning
input_file = ARGS[1]

if "SLURM_NTASKS" in keys(ENV)
    addprocs(SlurmManager())
    for i in workers()
        host,pid = fetch(@spawnat i (gethostname(), getpid()))
        @info "Worker $(i) is running on $(host) with ID $(pid)" 
    end
else
    @warn "No Slurm environment detected."
    #addprocs(2)
end

@everywhere begin
    using Pkg; Pkg.activate(".")
    using MUST
    using TSO
    using DataFrames
    using DelimitedFiles
    using Distributed
end

@everywhere include($(input_file))
@everywhere MUST.@import_dispatch dispatch_location
@everywhere prepare4dispatch = MUST.ingredients("prepare4dispatch.jl")


#===================== Step (-1): Add timings ================================#
begin
    if !("SLURM_NTASKS" in keys(ENV))
        @info "Activate timing..."
        TSO.activate_timing!.(TSO.timers)
        TSO.start_timing!()
    end
end


#=================== Step (A): The Initial Grid ==============================#
begin
    grid = MUST.StaggerGrid(initial_grid_path)
    deleteat!(grid.info, .!isfile.(grid["av_path"]))
    MUST.save(grid, initial_cl_path)

    # To save time, remove models now that are not interesting at the moment
    #deleteat!(grid.info, grid["feh"].!=0.0)
    #deleteat!(grid.info, (grid["teff"].==5000.0) .| ((grid["teff"].!=5777.0)))
    #deleteat!(grid.info, (grid["logg"].==4.0) .| ((grid["logg"].!=4.0)))
    MUST.save(grid, initial_mod_path)


    # Check for opacity table field
    if !("eos_root" in names(grid.info))
        @warn """The given grid does not contain information about the root EoS and opacity source. Adding the fallback table.
        This should not happen if you interpolated the models from a grid! Please make sure to provide the `eos_root` column in the grid!"""
        grid.info[!, "eos_root"] = [mother_table_path for _ in 1:nrow(grid.info)]
    end

    # within the EoS root folder, check if the desired EoS and opacities are available
    iwarn = false
    for i in 1:nrow(grid.info)
        if !isfile(grid.info[i, "eos_root"], eos_path)
            error("For grid entry $(i) there is no EoS at the given `eos_path` at the given `eos_root`!")
        end
        if !isfile(grid.info[i, "eos_root"], opa_path)
            error("For grid entry $(i) there is no Opacitiy table at the given `opa_path` at the given `eos_root`!")
        end
        if !iwarn
            if !isfile(grid.info[i, "eos_root"], sopa_path)
                @warn "For grid entry $(i) there is no scattering at the given `sopa_path` at the given `eos_root`!"
                global iwarn = true
                #global sopa_path = nothing
            end
        end
    end

    # save the paths for later
    if !("eos_original" in names(grid.info))
        grid.info[!, "eos_original"] = [eos_path for _ in 1:nrow(grid.info)]
    end
    if !("opa_original" in names(grid.info)) 
        grid.info[!, "opa_original"] = [opa_path for _ in 1:nrow(grid.info)]
    end
    if !("sopa_original" in names(grid.info))
        grid.info[!, "sopa_original"] = [sopa_path for _ in 1:nrow(grid.info)]
    end
    if (iwarn)
        grid.info[!, "sopa_original"] = ["" for _ in 1:nrow(grid.info)]
    end
end

#================ Step (B): Fomation opacities (parallel) ====================#
begin
    # compute the rosseland opacity if wanted
    do_ross = [(i==1) ? recompute_ross : false for i in 1:nrow(grid.info)]
    eos_final = [
        do_ross[i] ? 
        prepare4dispatch.compute_rosseland_opacities(
            grid.info[i, "eos_root"], 
            grid.info[i, "eos_original"], 
            grid.info[i, "opa_original"]
        ) : 
        grid.info[i, "eos_original"] 
        for i in 1:nrow(grid.info)
    ]
    grid.info[!, "eos_final"] = eos_final

    args = [
        (
            grid.info[i, "eos_root"], 
            grid.info[i, "av_path"], 
            grid.info[i, "name"],
            grid.info[i, "logg"],
            grid.info[i, "eos_final"],
            grid.info[i, "opa_original"],
            grid.info[i, "sopa_original"],
            :kmeans, 
            Nbins,
            (i==1) ? recompute_ross : false
        ) for i in 1:nrow(grid.info)
    ]

    # To save disk space, the formation opacities can be deleted after binning
    # For this it is required to run this sequentially
    @everywhere formation_and_bin(args) = begin
        eos_root = args[1]
        av_path = args[2]
        name = args[3]
        logg = args[4]
        eos_path = args[5]
        opa_path = args[6]
        sopa_path = args[7]
        method = args[8]
        Nbins = args[9]

        sopa_path = if (sopa_path=="")
            nothing
        else
            sopa_path
        end
        
        # formation opacities
        (!skip_formation) && prepare4dispatch.formation_opacities(
            name, eos_root, eos_path, opa_path, av_path; 
            logg=logg
        )

        if !skip_binning
            @info "Opacity Binning: $(name)"

            quadrants = make_quadrants(name, eos_root, opa_path)

            if Nbins != sum([q.nbins for q in quadrants])
                @warn "Given number of bins does not match the given binning! ($Nbins <-> $(sum([q.nbins for q in quadrants])))"
            end

            # bin opacities
            prepare4dispatch.bin_opacities(
                name, eos_root, eos_path, opa_path, av_path; 
                logg=logg, method=method, Nbins=Nbins, scattering_path=sopa_path,
                name_extension=name_extension,
                version=version,
                quadrants=quadrants
            )
        end

        # remove the formation opacities
        clean && prepare4dispatch.clean(eos_root, name)
    end

    # End-to-end binning, with clean-up
    Distributed.pmap(formation_and_bin, args)
    
    # Save the eos info
    grid.info[!, "name_extension"]   = [name_extension for _ in 1:nrow(grid.info)]
    grid.info[!, "version"]          = [version for _ in 1:nrow(grid.info)]
    grid.info[!, "binned_tables"]    = [TSO.join_full(name_extension, grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "binned_E_tables"]  = [TSO.join_full(name_extension, "E", grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "rad_bins"]         = [Nbins for _ in 1:nrow(grid.info)]

    # compute the resolution and the rounded size of the box
    # use the EoS that was just created for this
    prepare4dispatch.resolution!(
        grid, 
        patch_size=patch_size, 
        τ_up=τ_up, 
        τ_surf=τ_surf, 
        τ_down=τ_down, 
        scale_resolution=scale_resolution,
        τ_ee0=τ_ee0, 
        τ_eemin=τ_eemin, 
        τ_zee0=τ_zee0, 
        τ_rho0=τ_rho0,
        dxdz_max=dxdz_max
    )
end

#====================== Step (C): Conversion =================================#
begin
    prepare4dispatch.fromT_toE.(grid.info[!, "binned_tables"], grid.info[!, "binned_E_tables"], upsample=4096)

    # Copy the average model in the same folder so that we can link it all to the right place
    for i in 1:nrow(grid.info)
        if use_adiabat & use_avnewz
            error("You can not use the adiabat and the av model with new z scale as initial condition!")
        elseif use_adiabat
            prepare4dispatch.adiabat(
                joinpath(grid.info[i, "binned_E_tables"], "eos.hdf5"),
                grid.info[i, "av_path"], 
                grid.info[i, "logg"],
                saveat=joinpath(grid.info[i, "binned_E_tables"], "inim.dat"),
                ee_min=grid.info[i, "ee_min"],
                common_size=grid.info[i, "initial_model_size"]
            )
        elseif use_avnewz
            prepare4dispatch.new_zscale(
                joinpath(grid.info[i, "binned_tables"], "eos.hdf5"),
                grid.info[i, "av_path"], 
                grid.info[i, "logg"],
                saveat=joinpath(grid.info[i, "binned_E_tables"], "inim.dat"),
                common_size=grid.info[i, "initial_model_size"]
            )
        else
            cp(grid.info[i, "av_path"], joinpath(grid.info[i, "binned_E_tables"], "inim.dat"), force=true)
        end
    end

    if use_avnewz
        @info "Recomputing model dimensions based on new z scale."
        prepare4dispatch.resolution!(
            grid, 
            patch_size=patch_size, 
            τ_up=τ_up, 
            τ_surf=τ_surf, 
            τ_down=τ_down, 
            scale_resolution=scale_resolution,
            τ_ee0=τ_ee0, 
            τ_eemin=τ_eemin, 
            τ_zee0=τ_zee0, 
            τ_rho0=τ_rho0,
            dxdz_max=dxdz_max,
            use_inim=true
        )
    end
end

#=============== Step (D): Input namelists for dispatch ======================#
begin
    prepare4dispatch.create_namelist!(grid; namelist_kwargs...)

    for i in 1:nrow(grid.info)
        cp(grid.info[i, "namelist_name"], joinpath(grid.info[i, "binned_E_tables"], "ininml.dat"), force=true)
        cp(joinpath(grid.info[i, "binned_tables"], "bin_assignment.hdf5"), joinpath(grid.info[i, "binned_E_tables"], "bin_assignment.hdf5"), force=true)
    end

    MUST.save(grid, final_grid_path)

    # Stage the grid for execution, possible remove other output
    MUST.stage_namelists(grid, clean_namelists=clean_namelists, clean_logs=clean_logs)
end

#=============================================================================#


begin
    TSO.end_timing!()
end
