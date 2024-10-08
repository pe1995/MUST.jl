#= General =#
begin
    # Path of the dispatch installation
    dispatch_location = "../../../dispatch2"

    # Under what name should the binned opacities be saved
    name_extension    = "StaggerExtrapolated/SE1"

    # interpolated models
    initial_grid_path = "StaggerExtrapolated/SE1_magg_m0_a0_vmic1.mgrid"
    initial_cl_path   = "StaggerExtrapolated/SE1_magg_m0_a0_vmic1_avail.mgrid"
    initial_mod_path  = "StaggerExtrapolated/SE1_magg_m0_a0_vmic1_solar.mgrid"
    final_grid_path   = "StaggerExtrapolated/SE1_magg_m0_a0_vmic1_090824.mgrid"

    # clean namelists in dispatch folder (other than new ones)
    clean_namelists = false
    clean_logs = false

    # replace the initial condition with the corresponding adiabat
    use_adiabat  = false

    # add new z scale to initial model based on rosseland opacity (for dispatch only)
    use_avnewz = true
end


#= Dispatch setup =#
begin
    patch_size = 14                 # Points per patch
    τ_up = -5.0                     # Upper limit of simulation domain
    τ_surf = 0.0                    # Optical surface of simulation domain
    τ_down = 6.0                    # Lower limit of simulation domain
    τ_ee0 = -4.5                    # Newton cooling placement (energy)
    τ_eemin = τ_up                  # Mininmum energy of initial condition
    τ_zee0 = -2.5                   # Newton cooling placement (height)
    τ_rho0 =  0.0                   # Density normaliztion height
    dxdz_max = 3.0                  # how much bigger is box in x than z (max)
    scale_resolution = 0.60         # Down or upsampling of simulation domain
    namelist_kwargs = Dict(         # Additional modifications in namelist
        :newton_time=>100.0,        
        :friction_time=>120.0,     
        :newton_decay_scale=>10.0,  
        :friction_decay_scale=>10.0,  
        :courant_target=>0.25,
        :courant_hd=>0.25,
        :courant_rt=>0.5,
        :newton_params=>(           #   Optional: Give namelist field = NamedTuple 
            :on=>true,              #   for direct namelist replacement
            :delay_rt=>true
        ),
        :io_params=>(
            :out_time=>0.25,
            :end_time=>400
        ),
        :aux_params=>(
            :select=>["flux"],
        ),
        :boundary_params=>(
            :upper_bc=>7,
            :radiation_bc=>3,
            :smallr=>1e-9
        ),
        :an_params=>(
            :smallr=>1e-8,
            :dlnr_limit=>0.7
        ),
        :patch_params=>(
            :grace=>0.01,
            :nt=>5,
            :dt_small=>1e-10
        ),
        :sc_rt_params=>(
            :rt_grace=>0.01,
            :rt_freq=>0.0,
            :cdtd=>1.0
            #:timestep_limiter=>1.1 
        ),
        #:dispatcher0_params=>(
        #    :retry_stalled=>60,
        #),
        #:friction_params=>(
        #    :slope_factor=>0.7,
        #    :slope_decay_scale=>55.0
        #)
    )
end

#= Opacities =#
begin
    # Location of the opacity table
    #mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"
    #mother_table_path = "../../../opacity_tables/TSO_MARCS_magg_m0_a0_v1.8"
    extension = "magg_m0_a0_vmic1"
    mother_table_path = "../../../opacity_tables/TSO_M3D_$(extension)_v5.0"
    eos_path = "combined_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_sopacities_"*extension*".hdf5"

    # verion of the binned table
    # v1.0 - 8 Bins
    version = "v1.0"

    # Number of bins
    Nbins = 8

    # Skip binning procedure (assumes it has already been done)
    skip_binning = false

    # Skip formation opacity procedure (assumes it has already been done)
    skip_formation = false

    # recompute the rosseland optical depth for the first model in the grid
    recompute_ross = false

    # The binning quadrants where opacity bins should be placed
    make_quadrants(name, eos_root, opa_path) = begin
        # Find the formation height of the given wavelength
        qlim = round(
            prepare4dispatch.quadrantlimit(name, eos_root, opa_path, λ_lim=5.0), 
            sigdigits=3
        )

        # 4 MURaM bins
        #=quadrants = [ 
            TSO.Quadrant((0.0, 100.0), (-100, 0.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (0.0, 2.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (2.0, 4.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (4.0, 100.0), 1, stripes=:κ)
        ]=#

        # 8 bins
        quadrants = [ 
            TSO.Quadrant((0.0, 4.0), (qlim, 4.5), 2, stripes=:κ),
            TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
            TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (-100, qlim), 4, stripes=:λ),
        ]
    end

    # remove formation opacities after usage to save disk space
    clean = true
end