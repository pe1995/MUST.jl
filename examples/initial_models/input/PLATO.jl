#= General =#
begin
    # Path of the dispatch installation
    dispatch_location = "/home/eitner/shared/model_grid/dispatch2"

    # Under what name should the binned opacities be saved
    name_extension    = "PLATO/M"

    # PLATO models
    initial_grid_path = "PLATO/plato_initial_121223.mgrid"
    initial_cl_path   = "PLATO/plato_initial_121223_avail.mgrid"
    initial_mod_path  = "PLATO/plato_initial_121223_solar.mgrid"
    final_grid_path   = "PLATO/plato_121223.mgrid"

    # clean namelists in dispatch folder (other than new ones)
    clean_namelists = false
    clean_logs = false

    # replace the initial condition with the corresponding adiabat
    use_adiabat  = false
end


#= Dispatch setup =#
begin
    patch_size = 22 
    τ_up = -4.0 
    τ_surf = 0.0 
    τ_down = 6.0
    scale_resolution = 0.9
    namelist_kwargs = Dict()
end

#= Opacities =#
begin
    # This is only used in case no EoS is given in grid. BUT THERE SHOULD ALWAYS BE ONE FROM THE INTERPOLATION
    mother_table_path = "/mnt/beegfs/gemini/groups/bergemann/users/eitner/TS_opacity_tables/TSO.jl/examples/converting_tables/TSO_MARCS_magg_m0_a0_v1.6"
        
    # Opacity tables that should be used for the binning (unbinned tables)
    extension = "magg_m0_a0"
    eos_path = "combined_ross_eos_"*extension*".hdf5"
    opa_path = "combined_opacities_"*extension*".hdf5"
    sopa_path = "combined_Sopacities_"*extension*".hdf5"

    # verion of the binned table
    # v0.3 - 8 Bins
    # v0.4 - 10 Bins
    # v0.5 - 4 bins
    version = "v0.3"

    # Number of bins
    Nbins = 8

    # Skip binning procedure (assumes it has already been done)
    skip_binning = false

    # Skip formation opacity procedure (assumes it has already been done)
    skip_formation = false

    # recompute the rosseland optical depth for the first model in the grid
    recompute_ross = true

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
    clean = false
end