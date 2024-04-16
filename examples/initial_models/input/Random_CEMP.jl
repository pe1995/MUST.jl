#= General =#
begin
    # Path of the dispatch installation
    dispatch_location = "/home/eitner/shared/model_grid/dispatch2"

    # Under what name should the binned opacities be saved
    name_extension    = "CEMP/R1"

    # PLATO models
    initial_grid_path = "CEMP/random_CEMP_1.mgrid"
    initial_cl_path   = "CEMP/random_CEMP_1_avail.mgrid"
    initial_mod_path  = "CEMP/random_CEMP_1_avail.mgrid"
    final_grid_path   = "CEMP/random_CEMP_1_120424.mgrid"

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
    τ_up = -4.5                     # Upper limit of simulation domain
    τ_surf = 0.0                    # Optical surface of simulation domain
    τ_down = 6.0                    # Lower limit of simulation domain
    τ_ee0 = -0.6                   # Newton cooling placement (energy)
    τ_eemin = τ_up                  # Mininmum energy of initial condition
    τ_zee0 = -0.3                   # Newton cooling placement (height)
    τ_rho0 = -0.1                   # Density normaliztion height
    scale_resolution = 0.7          # Down or upsampling of simulation domain
    namelist_kwargs = Dict(         # Additional modifications in namelist
        :newton_time=>100.0,        #   Optional: Give namelist field = NamedTuple 
        :newton_decay_scale=>20.0,  #   for direct namelist replacement
        :courant_target=>0.3,
        :courant_rt=>0.3,
        :newton_params=>(
            :on=>true,
            :delay_rt=>true
        ),
        :io_params=>(
            :out_time=>1.0,
        ),
        :aux_params=>(
            :select=>["dt_rt", "flux"],
        ),
        :boundary_params=>(
            :upper_bc=>2,
            :smallr=>1e-16
        ),
        :an_params=>(
            :smallr=>1e-8,
            :smallc=>0.1, 
            :dlnr_limit=>0.5
        ),
        :patch_params=>(
            :grace=>0.1,
            :nt=>5
        ),
        :sc_rt_params=>(
            :rt_grace=>0.05,
            :rt_freq=>2.0 
        )
    )
end

#= Opacities =#
begin
    # Location of the opacity table
    mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_magg_m1_a0_c1_v3.0"
    extension = "magg_m1_a0_c0"
    eos_path = "combined_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_sopacities_"*extension*".hdf5"

    # verion of the binned table
    # v1.0 - 8 Bins, Carbon enhanced (+1)
    #   v1.0.1 - 8 Bins, Carbon not enhanced (+1)
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