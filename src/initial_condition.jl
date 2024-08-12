#==============================================================================#
# Initial condition generator
#-> generate initial conditions based on stellar parameters. 
#==============================================================================#

include("modelgrids.jl")        # for creating initial models
include("prepare4dispatch.jl")  # for opacity binning and input namelist

using Dates

#==============================================================================#
# 1. Initial Models
#==============================================================================#

# global grid files that can be used for interpolation
staggergrid = MUST.Atmos1DGrid(abspath(joinpath(dirname(@__FILE__), "..", "initial_grids", "Stagger", "stagger_grid.mgrid")))
marcsgrid = MUST.Atmos1DGrid(abspath(joinpath(dirname(@__FILE__), "..", "initial_grids", "MARCS", "marcs_grid.mgrid")))


opacityguess(name, root) = begin
    all = MUST.glob("$(name)*", root)
    if length(all) >= 1
        first(all)
    else
        error("There is no table with name, root: $(name),$(root)")
    end
end

"""
    initialmodel(paras; grid=staggergrid, eos=:closest, savedir="av_models", kwargs...)

Create initial models with the given parameters (Matrix with paras[:, 1]=teff, 2=logg and 3=feh), 
interpolated from the given grid. 
When `eos=:closest` the EoS that is closest to a grid point is chosen. You can use a different EoS
than those that are used for the grid models by specifying the path to the eos __file__. It is assumed that the opacities 
corresponding to this EoS are located in the same folder.
__Note__: The paths to the av_models in the grid are assumed to be relative paths to the location 
of the grid! In the fault case this is the examples/initial_models folder, but if you use a different 
grid the path of that grid will be used. 
"""
function initialmodel(paras::AbstractMatrix; grid=staggergrid, eos=[:closest for _ in axes(paras, 1)], opa="combined_opacities", sopa="combined_sopacities", savedir="av_models", kwargs...)
    # check if parameters are fine and within the grid
    @assert all(checkparameters(grid, paras))

    # convert to absolute paths in the grid
    grid_local = deepcopy(grid)
    MUST.absolute_path!(grid_local)

    # pick the right EoS
    eospath = [e==:closest ? find_closest_eos(grid_local, teff, feh, from="matching_eos") : abspath(e) for (teff, feh, e) in zip(paras[:, 1], paras[:, 3], eos)]
    
    eosroot = dirname.(eospath)
    eosloaded = [reload(TSO.SqEoS, eospath[i]) for i in axes(paras, 1)]
    opapath = [opacityguess("$(opa)*", eosroot[i]) for i in axes(paras, 1)]
    sopapath = [opacityguess("$(sopa)*", eosroot[i]) for i in axes(paras, 1)]

    # interpolate the initial models and save them 
    interpolatedGrid = interpolate_from_grid(
        grid_local; 
        teff=round.(paras[:, 1], sigdigits=4), 
        logg=round.(paras[:, 2], sigdigits=4), 
        feh=paras[:, 3],
        eos=eosloaded,
        folder=savedir,
        kwargs...
    )

    # make sure that the eos root dir is saved, so that we can later load the opacities
    interpolatedGrid.info[!, "eos_root"] .= eosroot
    interpolatedGrid.info[!, "matching_eos"] .= eospath
    interpolatedGrid.info[!, "matching_opa"] = opapath
    interpolatedGrid.info[!, "matching_sopa"] = sopapath
    interpolatedGrid
end

"""
    initialmodel(grid)

Create a initial model grid from the given grid of average models by filling in the EoS given.
If no `eos` given, the matching EoS is chosen by default. Note that relative paths
(that are needed for the base models, so that they can be used on different systems) are
replaced by absolute paths now.
"""
initialmodel(grid=staggergrid::MUST.Atmos1DGrid; eos=:matching, opa="combined_opacities", sopa="combined_sopacities") = begin
    grid_local = deepcopy(grid)
    MUST.absolute_path!(grid_local)

    eospath = if eos==:matching
        grid_local.info[!, "matching_eos"]
    else
        [abspath(eos) for _ in 1:nrow(grid_local.info)]
    end
    eosroot = dirname.(eospath)
    grid_local.info[!, "eos_root"] .= eosroot

    opapath = [opacityguess("$(opa)*", grid_local.info[i, "eos_root"]) for i in 1:MUST.nrow(grid_local.info)]
    sopapath = [opacityguess("$(sopa)*", grid_local.info[i, "eos_root"]) for i in 1:MUST.nrow(grid_local.info)]

    grid_local.info[!, "matching_eos"] .= eospath
    grid_local.info[!, "matching_opa"] = opapath
    grid_local.info[!, "matching_sopa"] = sopapath
    grid_local
end


#==============================================================================#
# 2. Opacity Binning & Input Namelists
#==============================================================================#

checkpaths!(grid::MUST.Atmos1DGrid) = begin
	# within the EoS root folder, check if the desired EoS and opacities are available
    iwarn = Ref(false)
    for i in 1:nrow(grid.info)
        if !isfile(grid.info[i, "matching_opa"])
            error(
				"For grid entry $(i) there is no Opacitiy table at the given `opa_path` at the given `eos_root`!"
			)
        end
        if !iwarn[]
            if !isfile(grid.info[i, "matching_sopa"])
                @warn "For grid entry $(i) there is no scattering at the given `sopa_path` at the given `eos_root`!"
                iwarn[] = true
            end
        end
    end

    # save the paths for later
    if !("eos_original" in names(grid.info))
        grid.info[!, "eos_original"] = basename.(grid.info[!, "matching_eos"])
    end
    if !("opa_original" in names(grid.info)) 
        grid.info[!, "opa_original"] = basename.(grid.info[!, "matching_opa"])
    end
    if !("sopa_original" in names(grid.info))
        grid.info[!, "sopa_original"] = basename.(grid.info[!, "matching_sopa"])
    end
    if !("eos_root" in names(grid.info))
        grid.info[!, "eos_root"] = dirname.(grid.info[!, "matching_eos"])
    end
    if (iwarn[])
        grid.info[!, "sopa_original"] = ["" for _ in 1:nrow(grid.info)]
    end
end

function prepare(
	g; 
	name_extension="",
	version = "v1.0",
    final_grid_path = name_extension*"_dispatch_$(string(Dates.today())).mgrid",
	use_adiabat=false, 
	use_avnewz=true, 
	patch_size = 14,                 # Points per patch
    τ_up = -6.0,                     # Upper limit of simulation domain
    τ_surf = 0.0,                    # Optical surface of simulation domain
    τ_down = 6.0,                    # Lower limit of simulation domain
    τ_ee0 = -4.5,                    # Newton cooling placement (energy)
    τ_eemin = τ_up,                  # Mininmum energy of initial condition
    τ_zee0 = -2.5,                   # Newton cooling placement (height)
    τ_rho0 =  0.0,                   # Density normaliztion height
    dxdz_max = 3.0,                  # how much bigger is box in x than z (max)
    scale_resolution = 0.75,         # Down or upsampling of simulation domain
	namelist_kwargs=Dict(),
	recompute_ross=false,
	Nbins = 8,
	skip_binning = false,
	skip_formation = false,
	clean = false,
	clean_namelists = false,
    clean_logs = false,
	make_quadrants = (name, eos_root, opa_path) -> begin
		qlim = round(
			quadrantlimit(name, eos_root, opa_path, λ_lim=5.0), 
			sigdigits=3
		)
	   quadrants = [ 
			TSO.Quadrant((0.0, 4.0), (qlim, 4.5), 2, stripes=:κ),
			TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
			TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
			TSO.Quadrant((0.0, 100.0), (-100, qlim), 4, stripes=:λ),
		]
	end
)
	# make sure that the grid exists
	grid = deepcopy(g)
	deleteat!(grid.info, .!isfile.(grid["av_path"]))

	# save the paths for later
	checkpaths!(grid)
    
    # compute the rosseland opacity if wanted
    do_ross = [(i==1) ? recompute_ross : false for i in 1:nrow(grid.info)]
    eos_final = [
        do_ross[i] ? 
        compute_rosseland_opacities(
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

	formation_and_bin(args) = begin
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
		(!skip_formation) && formation_opacities(
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
			bin_opacities(
				name, eos_root, eos_path, opa_path, av_path; 
				logg=logg, method=method, Nbins=Nbins, scattering_path=sopa_path,
				name_extension=name_extension,
				version=version,
				quadrants=quadrants
			)
		end
	
		# remove the formation opacities
		clean && clean(eos_root, name)
	end
    
    # End-to-end binning, with clean-up
    map(formation_and_bin, args)
    
    # Save the eos info
    grid.info[!, "name_extension"]   = [name_extension for _ in 1:nrow(grid.info)]
    grid.info[!, "version"]          = [version for _ in 1:nrow(grid.info)]
    grid.info[!, "binned_tables"]    = [TSO.join_full(name_extension, grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "binned_E_tables"]  = [TSO.join_full(name_extension, "E", grid.info[i, "name"], version, add_start=false) for i in 1:nrow(grid.info)]
    grid.info[!, "rad_bins"]         = [Nbins for _ in 1:nrow(grid.info)]

    # compute the resolution and the rounded size of the box
    # use the EoS that was just created for this
    resolution!(
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

    fromT_toE.(grid.info[!, "binned_tables"], grid.info[!, "binned_E_tables"], upsample=4096)

    # Copy the average model in the same folder so that we can link it all to the right place
    for i in 1:nrow(grid.info)
        if use_adiabat & use_avnewz
            error("You can not use the adiabat and the av model with new z scale as initial condition!")
        elseif use_adiabat
            adiabat(
                joinpath(grid.info[i, "binned_E_tables"], "eos.hdf5"),
                grid.info[i, "av_path"], 
                grid.info[i, "logg"],
                saveat=joinpath(grid.info[i, "binned_E_tables"], "inim.dat"),
                ee_min=grid.info[i, "ee_min"],
                common_size=grid.info[i, "initial_model_size"]
            )
        elseif use_avnewz
            new_zscale(
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
        resolution!(
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

    create_namelist!(grid; namelist_kwargs...)
    for i in 1:nrow(grid.info)
        cp(grid.info[i, "namelist_name"], joinpath(grid.info[i, "binned_E_tables"], "ininml.dat"), force=true)
        cp(joinpath(grid.info[i, "binned_tables"], "bin_assignment.hdf5"), joinpath(grid.info[i, "binned_E_tables"], "bin_assignment.hdf5"), force=true)
    end

    MUST.save(grid, final_grid_path)

    # copy over namelists to dispatch
    for i in 1:nrow(grid.info)
        fold = joinpath(grid.info[i, "binned_E_tables"], "ininml.dat")
        f = dirname(fold)
        name = last(split(f, '/', keepempty=false))
        pnew = MUST.@in_dispatch "$(name).nml"
        @info "copy $(fold) to $(pnew)."
        cp.(fold, pnew, force=true)
    end

    # Stage the grid for execution, possible remove other output
    #MUST.stage_namelists(grid, clean_namelists=clean_namelists, clean_logs=clean_logs)

    grid
end