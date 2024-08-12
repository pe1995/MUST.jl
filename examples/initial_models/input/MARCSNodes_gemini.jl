#= General information =#
begin
    # Output name of the models
    name_extension    = "MARCSExtrapolated/ME1"

    # Location of the dispatch installation
    dispatch_location = "../../../dispatch2/"

    # input and output names of the grid
    #initial_grid_path = "stagger_grid_full_o.mgrid"
    #initial_grid_path = "stagger_grid_sun.mgrid"
    #initial_grid_path = "stagger_grid_full_subgiant.mgrid"
    #initial_grid_path = "stagger_grid_full_solar.mgrid"
    #initial_grid_path = "stagger_grid_red.mgrid"
    #initial_grid_path = "stagger_grid_solar.mgrid"
    initial_grid_path = "MARCSExtrapolated/ME1_magg_m0_a0_vmic1.mgrid"

    initial_cl_path   = "MARCSExtrapolated/ME1_magg_m0_a0_vmic1.mgrid"
    initial_mod_path  = "MARCSExtrapolated/ME1_magg_m0_a0_vmic1.mgrid"
    final_grid_path   = "MARCSExtrapolated/dispatch_ME1_magg_m0_a0_vmic1.mgrid"

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
    patch_size = 14                 # Points per patch
    τ_up = -5.0                     # Upper limit of simulation domain
    τ_surf = 0.0                    # Optical surface of simulation domain
    τ_down = 6.0                    # Lower limit of simulation domain
    τ_ee0 = -4.5                    # Newton cooling placement (energy)
    τ_eemin = τ_up                  # Mininmum energy of initial condition
    τ_zee0 = -2.5                   # Newton cooling placement (height)
    τ_rho0 =  0.0                   # Density normaliztion height
    dxdz_max = 3.0                  # how much bigger is box in x than z (max)
    scale_resolution = 0.55         # Down or upsampling of simulation domain
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
            :out_time=>1.0,
        ),
        :aux_params=>(
            :select=>["flux"],
        ),
        :boundary_params=>(
            :upper_bc=>7,
            :radiation_bc=>3,
            :smallr=>1e-8
        ),
        :an_params=>(
            :smallr=>1e-8,
            :smallc=>0.1, 
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
        :dispatcher0_params=>(
            :retry_stalled=>60,
        ),
        #:friction_params=>(
        #    :slope_factor=>0.7,
        #    :slope_decay_scale=>55.0
        #)
    )
end

#= Opacities =#
begin
    # recompute the rosseland optical depth for the first model in the grid
    recompute_ross = false

    # Location of the opacity table
    extension = "magg_m0_a0_vmic1"
    mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/storage/opacity_tables/TSO_M3D_$(extension)_v5.0"
    eos_path = "combined_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_sopacities_"*extension*".hdf5"

    # opacity table version (output)
    # v0.5 -> 8 bins (MARCS)
    #    v0.5+1 -> 8 bins (MARCS) - redo for debugging (logg 4.44 for sun)
    #    v0.5+3 -> default 8 bins + 1 bin in lines, set line limiter to 5.0
    #    v0.5+4 -> 12 bins, set line limiter to 5.5 to test influence on lines
    #    v0.5.1 -> Grey (MARCS)
    #    v0.5.2 -> 4 MURaM bins (MARCS)
    #    v0.5.3 -> 5 bins in paper-setup
    #    v0.5.4 -> 8 bins in paper-setup
    #    v0.5.5 -> 12 bins in paper-setup
    #    v0.5.6 -> 7 bins in paper-setup
    #    v0.5.7 -> 12 bins in Co5bold setup
    # v0.6 -> 8 bins (M3D)
    version = "v0.6"

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

        # 9 bins
        #=quadrants = [ 
            TSO.Quadrant((0.0, 4.0), (qlim, 5.0), 3, stripes=:κ),
            TSO.Quadrant((0.0, 4.0), (4.5, 100), 1, stripes=:κ),
            TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (-100, qlim), 4, stripes=:λ),
        ]=#

        # 12 bins
        #=quadrants = [ 
            TSO.Quadrant((0.0, 4.0), (qlim, 5.5), 5, stripes=:κ),
            TSO.Quadrant((0.0, 4.0), (5.5, 100), 1, stripes=:κ),
            TSO.Quadrant((4.0, 100.0), (qlim, 100), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (-100, qlim), 5, stripes=:λ),
        ]=#

        # 4 MURaM bins
        #=quadrants = [ 
            TSO.Quadrant((0.0, 100.0), (-100, 0.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (0.0, 2.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (2.0, 4.0), 1, stripes=:κ),
            TSO.Quadrant((0.0, 100.0), (4.0, 100.0), 1, stripes=:κ)
        ]=#

        # 12 CO5BOLD bins
        #=quadrants = [ 
            TSO.Quadrant(log10.((0.0,    5500.0)),     reverse(-1 .* (0.15,  99.0)), 1, stripes=:κ)
            TSO.Quadrant(log10.((5500.0, 10000000.0)),  reverse(-1 .* (0.15,  99.0)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    6000.0)),     reverse(-1 .* (0.00,  0.15)), 1, stripes=:κ)
            TSO.Quadrant(log10.((6000.0, 10000000.0)),  reverse(-1 .* (0.00,  0.15)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    6500.0)),     reverse(-1 .* (-0.75, 0.00)), 1, stripes=:κ)
            TSO.Quadrant(log10.((6500.0, 10000000.0)),  reverse(-1 .* (-0.75, 0.00)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-1.50, -0.75)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-2.25, -1.50)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-3.00, -2.25)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-3.75, -3.00)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-4.50, -3.75)), 1, stripes=:κ)
            TSO.Quadrant(log10.((0.0,    10000000.0)),  reverse(-1 .* (-99.0, -4.50)), 1, stripes=:κ)
        ]=#

        # grey
        #=quadrants = [ 
            TSO.Quadrant((0.0, 100.0), (-100, 100), 1)
        ]=#

        # N bins paper style
        #=quadrants=[ 
            TSO.Quadrant((0.0, 4.0),   (-100, 4.5), 4),
            TSO.Quadrant((0.0, 4.0),   (4.5, 100), 1),
            TSO.Quadrant((4.0, 100.0), (-100, 100), 2)
        ]=#
    end
end