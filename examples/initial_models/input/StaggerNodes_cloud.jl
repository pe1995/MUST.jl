#= General information =#
begin
    # Output name of the models
    name_extension    = "DIS_MARCS"

    # Location of the dispatch installation
    dispatch_location = "/home/ubuntu/DISPATCH/dispatch2"

    # input and output names of the grid
    initial_grid_path = "stagger_grid_full_o.mgrid"
    initial_cl_path   = "stagger_grid_avail.mgrid"
    initial_mod_path  = "stagger_grid_solar.mgrid"
    final_grid_path   = "dispatch_grid.mgrid"

    # clean namelists in dispatch folder (other than new ones)
    clean_namelists = false
    clean_logs = false

    # replace initial model with adiabat (for dispatch only)
    use_adiabat = false

    # add new z scale to initial model based on rosseland opacity (for dispatch only)
    use_avnewz = true
end

#= Dispatch setup =#
begin
    patch_size = 22                 # Points per patch
    τ_up = -3.8                     # Upper limit of simulation domain
    τ_surf = 0.0                    # Optical surface of simulation domain
    τ_down = 7.0                    # Lower limit of simulation domain
    τ_ee0 = -0.5                    # Newton cooling placement (energy)
    τ_eemin = τ_up                  # Mininmum energy of initial condition
    τ_zee0 = -2.0                   # Newton cooling placement (height)
    τ_rho0 = 0.5                    # Density normaliztion height
    scale_resolution = 0.9          # Down or upsampling of simulation domain
    namelist_kwargs = Dict(         # Additional modifications in namelist
        :newton_time=>100.0,         #   Optional: Give namelist field = NamedTuple 
        :newton_decay_scale=>10.0,   #   for direct namelist replacement
        :courant_target=>0.2,
        :courant_rt=>0.1,
        :newton_params=>(
            :on=>true,
            :delay_rt=>true
        ),
        :io_params=>(
            :out_time=>0.5,
        )            
    )
end

#= Opacities =#
begin
    # recompute the rosseland optical depth for the first model in the grid
    recompute_ross = false

    # Location of the opacity table
    mother_table_path = "/home/ubuntu/DISPATCH/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.7"
    extension = "magg_m0_a0"
    eos_path = "combined_ross_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"

    # opacity table version (output)
    version = "v0.4"

    # Number of bins in the opacity table (output)
    Nbins = 8

    # Skip binning procedure (assumes it has already been done)
    skip_binning = true

    # Skip formation opacity procedure (assumes it has already been done)
    skip_formation = true

    # remove formation opacities after binning (save disk space)
    clean = false

    # The binning quadrants where opacity bins should be placed
    # For each quadrant in list: 1. log10 λ range, 2. log10 τ range, 3. Nbins
    make_quadrants(name, eos_root, opa_path) = begin
        qlim = round(
            prepare4dispatch.quadrantlimit(name, eos_root, opa_path, λ_lim=5.0), 
            sigdigits=3
        )

        # 8 bins
        quadrants = [ 
            TSO.Quadrant((0.0, 4.0), (qlim, 4.5), 2, stripes=:κ),
            TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
            TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (-100, qlim), 4, stripes=:λ),
        ]

        # 4 MURaM bins
        #=quadrants = [ 
            TSO.Quadrant((0.0, 100.0), (-100, 0.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (0.0, 2.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (2.0, 4.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (4.0, 100.0), 1, stripes=:κ)
        ]=#
    end
end