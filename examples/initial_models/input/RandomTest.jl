host = "gemini"

if host == "raven"
    name_extension    = "DIS_MARCS"
    dispatch_location = "/u/peitner/DISPATCH/dispatch2/"

    initial_grid_path = "random_grid.mgrid"
    initial_cl_path   = "random_grid_avail.mgrid"
    initial_mod_path  = "random_grid_solar.mgrid"
    final_grid_path   = "random_models.mgrid"

    mother_table_path = "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6"

    extension         = "magg22"
    eos_path = "combined_ross_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"

    version = "v0.1"
    Nbins = 8
    clean = false
    use_adiabat = false

elseif host == "gemini"
    name_extension    = "DIS_MARCS"
    dispatch_location = "/home/eitner/shared/model_grid/dispatch2"

    ## Random test models
    initial_grid_path = "random_setup.mgrid"
    final_grid_path   = "random_grid.mgrid"
    initial_grid_path = "node_setup.mgrid"
    final_grid_path   = "node_grid.mgrid"

    ## This is only used in case no EoS is given in grid. BUT THERE SHOULD ALWAYS BE ONE FROM THE INTERPOLATION
    #mother_table_path = "/home/eitner/shared/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_v1.6"
    #mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_asplund_m0_a0_v1.6"
    mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.6"
    #mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.5"
    
    # Because there are multiple EoS per opacity table possible (e.g. to test influence of EoS with same opacities)
    # you need to provide the name extension under which the EoS will be searched. In the given root_eos, the tables
    # will be loaded from
    # mother_table_path * "combined_ross_eos" * extension * ".hdf5"
    # mother_table_path * "combined_opacities" * extension * ".hdf5"
    extension = "magg_m0_a0"
    eos_path = "combined_ross_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"

    # v0.3 - 8 Bins
    # v0.4 - 10 Bins
    # v0.5 - 4 bins
    version = "v0.3"
    Nbins = 8
    clean = false
    use_adiabat = false
elseif host == "cloud"
    name_extension    = "DIS_MARCS"
    dispatch_location = "/home/ubuntu/DISPATCH/dispatch2"

    #initial_grid_path = "random_setup.mgrid"
    #final_grid_path   = "random_grid.mgrid"
    #initial_grid_path = "node_setup.mgrid"
    #final_grid_path   = "node_grid.mgrid"

    mother_table_path = "/home/ubuntu/DISPATCH/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.7"
    
    extension = "magg_m0_a0"
    eos_path = "combined_ross_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"

    version = "v0.4"
    Nbins = 8
    clean = false
    use_adiabat  = true
end

# The binning quadrants where opacity bins should be placed
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
end