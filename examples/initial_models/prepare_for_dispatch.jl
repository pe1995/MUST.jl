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
    addprocs(3)
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
    name_extension    = "DIS_MARCS"
    version           = "0.1"
    dispatch_location = "/u/peitner/DISPATCH/dispatch2/"
    MUST.@import_dispatch dispatch_location
end





#= Additional functionality =#

@everywhere function resolution(av_model, min_x, max_x, min_z, max_z, hres, 
                                patch_size=20; cut_bottom=0.3)
    nz = ceil(length(av_model.z) * (1-cut_bottom)) * 1.2
    dz = round((max_z - min_z) * (1-cut_bottom), sigdigits=2)
    target_res = round(dz/nz, sigdigits=2)
    @show target_res dz nz
    nz = ceil(nz / patch_size) * patch_size

    ## how many patches do we need then
    n_patches = div(nz, patch_size)

    ## So we need to make sure that we have an integer multiple of this number of
    ## patches in the other dimensions
    dx = round((max_x - min_x) * (1-cut_bottom), sigdigits=2)
    
    ## we need this many more patches in x dimension
    dxz_frac = dx / dz
    n_xpatches = dxz_frac * n_patches
    n_xpatches_orig = n_xpatches

    ## we round this to the next quater
    to_y = 0.1
    
   

    n_xpatches_m = find_integer(-, dxz_frac, n_patches, to_y; digits=1)


    dxz_frac = round(round_to_next(dx / dz, to_y), digits=1)
    n_xpatches_p = find_integer(+, dxz_frac, n_patches, to_y; digits=1)


    dxz_frac = round(round_to_next(dx / dz, to_y), digits=1)
    n_xpatches_hr = find_integer(+, dxz_frac, n_patches+1, to_y; digits=1)

    diff_orig(x) = abs(x / n_patches - n_xpatches_orig / n_patches)
    n_xpatches = [n_xpatches_m, n_xpatches_p, n_xpatches_hr][argmin([diff_orig(n_xpatches_m), diff_orig(n_xpatches_p), diff_orig(n_xpatches_hr)])]


    nz = n_patches * patch_size
    nx = n_xpatches * patch_size

    target_res = round(dz/nz, sigdigits=2)
    dz = round(nz * target_res, sigdigits=2)
    res_points = nx / nz
    dx = res_points * dz
    
    dx, nx, dz, nz
end

@everywhere function resolutionHD(av_model, min_x, max_x, min_z, max_z, hres, 
                                    patch_size=30; cut_bottom=0.3, scale_resolution=1.2)
    nx = ceil(abs(max_x - min_x)/abs(hres)) * scale_resolution
    dx = round((max_x - min_x) * (1-cut_bottom), sigdigits=2)
    target_res = round(dx/nx, sigdigits=2)
    nx = ceil(nx / patch_size) * patch_size

    ## how many patches do we need then
    n_patches = div(nx, patch_size)

    ## So we need to make sure that we have an integer multiple of this number of
    ## patches in the other dimensions
    dz = round((max_z - min_z) * (1-cut_bottom), sigdigits=2) 
    
    ## we need this many more patches in x dimension
    dxz_frac = dz / dx
    n_zpatches = dxz_frac * n_patches
    n_zpatches_orig = n_zpatches

    ## we round this to the next to_y
    to_y = 0.1
    n_zpatches_m = find_integer(-, dxz_frac, n_patches, to_y; digits=1)

    dxz_frac = round(round_to_next(dz / dx, to_y), digits=1)
    n_zpatches_p = find_integer(+, dxz_frac, n_patches, to_y; digits=1)

    dxz_frac = round(round_to_next(dz / dx, to_y), digits=1)
    n_zpatches_hr = find_integer(+, dxz_frac, n_patches+1, to_y; digits=1)

    diff_orig(x) = abs(x / n_patches - n_zpatches_orig / n_patches)
    n_zpatches = [n_zpatches_m, n_zpatches_p, n_zpatches_hr][argmin([diff_orig(n_zpatches_m), diff_orig(n_zpatches_p), diff_orig(n_zpatches_hr)])]


    nx = n_patches * patch_size
    nz = n_zpatches * patch_size

    target_res = round(dx/nx, sigdigits=2)
    dx = round(nx * target_res, sigdigits=2)
    res_points = nz / nx
    dz = res_points * dx
    
    @show dx, nx, dz, nz

    dx, nx, dz, nz
end

@everywhere resolution!(grid::MUST.AbstractMUSTGrid; 
                            patch_size=30, cut_bottom=0.3) = begin
    xr, zr = zeros(Int, nrow(grid.info)), zeros(Int, nrow(grid.info))
    xd, zd = zeros(Float64, nrow(grid.info)), zeros(Float64, nrow(grid.info))


    # At this point, it might be usefull to use the EoS
    # We should -> compute optical depth of model
    #           -> shift the box down by so much that roughy tau=-5 is the upper edge
    #           -> also save the density and the z at the optical surface, we can use those as scaling
    eos = [reload(SqEoS, joinpath(grid.info[i, "binned_tables"], "eos.hdf5")) for i in 1:nrow(grid.info)]
    opa = [reload(SqOpacity, joinpath(grid.info[i, "binned_tables"], "binned_opacities.hdf5")) for i in 1:nrow(grid.info)]

    models = [@optical(Average3D(eos[i], grid.info[i, "av_path"], logg=grid.info[i, "logg"]), eos[i], opa[i]) for i in 1:nrow(grid.info)]


    rho_norm = zeros(nrow(grid.info))
    l_norm = zeros(nrow(grid.info))
    z_up = zeros(nrow(grid.info))
    eemin = zeros(nrow(grid.info))

    z_lo = zeros(nrow(grid.info))
    d_lo = zeros(nrow(grid.info))
    T_lo = zeros(nrow(grid.info))


    @show maximum(log10.(models[1].τ)) minimum(log10.(models[1].τ))

    for i in 1:nrow(grid.info)
        ## What should the resolution be
        xd[i], xr[i], zd[i], zr[i] = resolutionHD(models[i], 
                                                grid.info[i, "mi_x"], grid.info[i, "ma_x"], grid.info[i, "mi_z"], grid.info[i, "ma_z"], grid.info[i, "hres"],
                                                patch_size; cut_bottom=cut_bottom)
        
        ## Where is the optical surface
        it0 = argmin(abs.(log10.(models[i].τ) .- 0.0))

        ## where is the upper edge (roughly)
        itup = argmin(abs.(log10.(models[i].τ) .- -5.5))

        ## Where can we start with the adiabat
        itlo = argmin(abs.(log10.(models[i].τ) .- 5))


        rho_norm[i] = models[i].lnρ[it0] |> exp |> log10 |> round |> exp10 
        l_norm[i]   = zd[i] |> abs |> log10 |> floor |> exp10
        z_up[i]     = abs(models[i].z[itup] - models[i].z[it0])
        eemin[i]    = exp.(models[i].lnEi[itup])

        z_lo[i] = round(models[i].z[itlo], sigdigits=5)
        d_lo[i] = round(exp.(models[i].lnρ[itlo]), sigdigits=5)
        T_lo[i] = round(exp.(models[i].lnT[itlo]), sigdigits=5)

    end

    xd, xr, zd, zr

    grid.info[!, "x_resolution"] = xr
    grid.info[!, "z_resolution"] = zr
    grid.info[!, "x_size"] = xd
    grid.info[!, "z_size"] = zd
    grid.info[!, "patch_size"] = [patch_size for _ in 1:nrow(grid.info)]
    grid.info[!, "rt_patch_size"] = [patch_size*2 for _ in 1:nrow(grid.info)]


    grid.info[!, "rho_norm"] = rho_norm
    grid.info[!, "l_norm"] = l_norm
    grid.info[!, "z_up"] = z_up
    grid.info[!, "ee_min"] = eemin
    grid.info[!, "initial_model_size"] = [length(models[i].z) for i in 1:nrow(grid.info)]

    grid.info[!, "z_lo"] = z_lo
    grid.info[!, "d_lo"] = d_lo
    grid.info[!, "T_lo"] = T_lo


    nothing
end

@everywhere round_to_next(x, to_y) = begin
    x_new = x
    x_new = if (x_new % to_y) >= to_y/2
        div(x_new, to_y) * to_y + to_y
    else
        div(x_new, to_y) * to_y
    end
end

@everywhere find_integer(f, dxz_frac, n_patches, to_y = 0.1; digits=1) = begin
    dxz_frac = round(round_to_next(dxz_frac, to_y), digits=digits)

    ## now we decrease the number of patches in x direction until we get an even number
    n_xpatches = dxz_frac * n_patches
    i = 0
    while !isinteger(n_xpatches)
        dxz_frac = f(dxz_frac, to_y)
        dxz_frac = round(dxz_frac, digits=digits)
        n_xpatches = dxz_frac * n_patches

        i += 1
    end

    n_xpatches
end

@everywhere function formation_opacities(mother_table, av_path, logg, name)

    # The EoS has already been smoothed in the running process
    table_folder = mother_table
    opacities    = reload(SqOpacity, joinpath(table_folder, "combined_opacities.hdf5"))
    eos          = reload(SqEoS,     joinpath(table_folder, "combined_eos.hdf5"))
    aos          = @axed eos


    # For any model, compute the internal energy given T, rho by bisecting using the EoS
    model = Average3D(av_path, logg=logg)


    # Recompute the rosseland optical depth scale based on the combined opacities
    # This should be more accurate, because the rosseland optical depth requires 
    # integration over frequency, which is not available in the individual runs
    if !isfile(joinpath(table_folder, "combined_ross_eos.hdf5"))
        @info "computing rosseland ($name)"
        rosseland_opacity!(aos.eos.lnRoss, aos, opacities; weights=ω_midpoint(opacities))
        transfer_rosseland!(aos.eos, opacities)
        save(aos.eos,   joinpath(table_folder, "combined_ross_eos.hdf5"))
        save(opacities, joinpath(table_folder, "combined_opacities.hdf5"))
    end


    # Compute the optical depth scale + the formation height
    #if isfile(joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"))
    #    @warn "skipping formation opacities for $(name)"
    #    return
    #end

    τ_ross, τ_λ = optical_depth(aos, opacities, model)
    d_ross, d_κ = formation_height(model, aos, opacities, τ_ross, τ_λ)


    # Save the results
    formation_opacities = SqOpacity(d_κ, d_ross, opacities.src, opacities.λ, true)
    save(formation_opacities, joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"))

end

@everywhere function bin_opacities(mother_table, av_path, logg, name)
    table_folder   = mother_table


    # TS quantities after post-processing
    eos           = reload(SqEoS,     joinpath(table_folder, "combined_ross_eos.hdf5")) # ross for others
    opacities     = reload(SqOpacity, joinpath(table_folder, "combined_opacities.hdf5"), mmap=true)
    formOpacities = reload(SqOpacity, joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"), mmap=true)

    # λ Integration weights
    weights = ω_midpoint(opacities)

    # Load a model for the transition between optically thin and thick regime
    model = Average3D(av_path, logg=logg)

    #bins_semistagger = StaggerBinning(TSO.SemiStaggerBins, 
    #                                opacities=opacities, 
    #                                formation_opacity=-log10.(formOpacities.κ_ross), κ_bins=3, Nbins=7)

    #bins_tabgen = TabgenBinning(TSO.UniformTabgenBins, opacities=opacities, formation_opacity=-log10.(formOpacities.κ_ross), Nbins=5, line_bins=3)
    TSO.Clustering.Random.seed!(42)    
    bins = ClusterBinning(TSO.KmeansBins;
                            opacities=opacities, 
                            formation_opacity=-log10.(formOpacities.κ_ross), 
                            Nbins=12, maxiter=1000, λ_split=(4.0, 5))

    bins = ClusterBinning(
        TSO.KmeansBins;
        opacities=opacities, 
        formation_opacity=-log10.(formOpacities.κ_ross), 
        Nbins=7, maxiter=1000, 
        quadrants=[ 
            TSO.Quadrant((0.0, 4.0),   (-100, 4.5), 4),
            TSO.Quadrant((0.0, 4.0),   (4.5, 100), 1),
            TSO.Quadrant((4.0, 100.0), (-100, 100), 2),
        ],
        display=:none
    )

    if isdir("$(name_extension)_$(name)_v$(version)")
        @warn "skipping binning for $(name)"
        return
    end

    bin8  = binning(bins, opacities, -log10.(formOpacities.κ_ross)) 


    # Save the binned opacities only
    function save_table(binned_opacities, version)
        eos_table_name = "$(name_extension)_$(name)_v$(version)"

        !isdir(eos_table_name) && mkdir(eos_table_name) 

        save(binned_opacities.opacities, joinpath(eos_table_name, "binned_opacities.hdf5"))

        # Copy the eos for convenience. Usually not a big deal because rather small
        cp(joinpath(table_folder, "combined_ross_eos.hdf5"), joinpath(eos_table_name, "eos.hdf5"), force=true)
    end

    binned_opacities8 = tabulate(bin8, weights, eos, opacities, transition_model=model)
    save_table(binned_opacities8, version)
end

@everywhere bin_opacities(args) = bin_opacities(args...)


@everywhere function fromT_toE(table_folder, folder_new; upsample=2048)
    eos = reload(SqEoS,     joinpath(table_folder, "eos.hdf5"))
    opa = reload(SqOpacity, joinpath(table_folder, "binned_opacities.hdf5"));

    aos = @axed eos
    eosE, opaE = switch_energy(aos, opa, upsample=upsample);
    aosE = @axed eosE

    TSO.fill_nan!(aosE, opaE)
    TSO.set_limits!(aosE, opaE)

    !isdir(folder_new) && mkdir(folder_new) 

    for_dispatch(eosE, opaE.κ, opaE.src, ones(eltype(opaE.src), size(opaE.src)...), folder_new)

    save(opaE, joinpath(folder_new, "binned_opacities.hdf5"))
    save(eosE, joinpath(folder_new, "eos.hdf5"))

    #mv("tabparam.in", joinpath(folder_new, "tabparam.in"), force=true)
    #mv("eostable.dat", joinpath(folder_new, "eostable.dat"), force=true)
    #mv("rhoei_radtab.dat", joinpath(folder_new, "rhoei_radtab.dat"), force=true);

end

@inline patches(res, patch_size) = Int(floor(res / patch_size))

function create_namelist(name, x_resolution, z_resolution, x_size, z_size, 
                        patch_size, rt_patch_size, δz, l_cgs, d_cgs, logg, teff, 
                        eemin, nz, initial_path, n_bin, eos_table, tscale)
    ngrid = nrow(grid.info)
    phase = grid.name

    dummy_nml = MUST.StellarNamelist("stellar_default.nml")

    # name of namelists
    name = "grid_$(name).nml"

    # copy the dummy namelist
    nml = deepcopy(dummy_nml)

    # set the box size + number of patches and points per patch
    # in the rt params, set the upper limit of rt correctly (and lower)
    # shift the box down by the right amount (there needs to be a parameter from the initial model for this)
    # set the EoS and the number of bins correctly 
    # also set the path and size of the initial model (also an info for this needs to be added.)
    dup = round(z_size /2 -δz, sigdigits=2) *1.
    
    courant_rt = if (logg >= 4) & (teff<6500)
        1.0
    elseif ((logg < 4) & (logg >= 3)) | (teff>=6500)
        0.8
    else
        0.7
    end

    courant_hd = if logg >= 4
        0.3
    elseif (logg < 4) & (logg >= 3)
        0.25
    else
        0.2
    end

    larger_than_sun = round(x_size/l_cgs / 4.6, sigdigits=1)
    newton_time = 100 #* larger_than_sun
    friction_time = 100 #* larger_than_sun
    newton_scale = 0.1 #* larger_than_sun

    # experiment with t scale based on size only
    test_tscale = false
    tscale = if test_tscale
        round(x_size / (4.6 * 1e8), sigdigits=1) * 100.0
    else
        larger_than_sun*tscale
    end

    @show tscale
    @show dup δz z_size/2

    MUST.set!(nml, cartesian_params=(:size=>[round(x_size/l_cgs, digits=2), round(x_size/l_cgs, digits=2), round(z_size/l_cgs, digits=2)], 
                                        :dims=>[patches(x_resolution, patch_size), 
                                                patches(x_resolution, patch_size), 
                                                patches(z_resolution, patch_size)],
                                    :position=>[0,0,round(-dup/l_cgs, digits=2)]),
                    patch_params=(:n=>[patch_size, patch_size, patch_size],),
                    scaling_params=(:l_cgs=>l_cgs, :d_cgs=>d_cgs, :t_cgs=>tscale),
                    stellar_params=(:g_cgs=>round(exp10(logg), digits=5), 
                                    :ee_min_cgs=>round(log(eemin), digits=5), 
                                    :nz=>nz-1, 
                                    :initial_path=>initial_path),
                    friction_params=(:end_time=>friction_time,),
                    gravity_params=(:constant=>-round(exp10(logg), digits=5),),
                    newton_params=(:ee0_cgs=>round(log(eemin), digits=5), 
                                    :position=>round((z_size/2 - 1.2*dup)/l_cgs, digits=2), 
                                    :end_time=>newton_time, 
                                    :scale=>newton_scale),
                    sc_rt_params=(  :rt_llc=>[-round(x_size/2/l_cgs, digits=2), -round(x_size/2/l_cgs, digits=2), -round((z_size/2 + dup)/l_cgs, digits=2)], 
                                    :rt_urc=>[ round(x_size/2/l_cgs, digits=2),  round(x_size/2/l_cgs, digits=2),  round((z_size/2 - dup)/l_cgs, digits=2)], 
                                    :n_bin=>n_bin,
                                    :courant=>courant_rt,
                                    :rt_res=>[-1,-1,rt_patch_size]),
                    an_params=(:courant=>courant_hd,),
                    eos_params=(:table_loc=>eos_table,)
                )
    
    # write the namelists
    MUST.write(nml, MUST.@in_dispatch(name))

    # save the summary
    MUST.@in_dispatch(name)
end

create_namelist!(grid::MUST.StaggerGrid) = begin
    g(i, v) = grid.info[i, v]

    names = []
    for i in 1:nrow(grid.info)
        name = create_namelist( g(i, "name"), 
                                g(i, "x_resolution"), 
                                g(i, "z_resolution"), 
                                g(i, "x_size"),
                                g(i, "z_size"),
                                g(i, "patch_size"),
                                g(i, "rt_patch_size"),
                                g(i, "z_up"),
                                g(i, "l_norm"),
                                g(i, "rho_norm"),
                                g(i, "logg"),
                                g(i, "teff"),
                                g(i, "ee_min"),
                                g(i, "initial_model_size"),
                                joinpath("input_data/grd/", 
                                            g(i, "binned_E_tables"), 
                                            "inim.dat"),
                                g(i, "rad_bins"),
                                joinpath("input_data/grd/", 
                                            g(i, "binned_E_tables")),
                                g(i, "tscale"))
        append!(names, [name])
    end

    grid.info[!, "namelist_name"] = names
end








begin
#= Step (A): The Grid =#

    grid = MUST.StaggerGrid("stagger_grid.mgrid")

    ## Check for opacity table field
    mother_table_path = "/u/peitner/DISPATCH/opacity_tables/TSO_MARCS_v1.6"
    if !("eos_root" in names(grid.info))
        @warn "The given grid does not contain information about the root EoS and opacity source. Adding the fallback table."
        grid.info[!, "eos_root"] = [mother_table_path for _ in 1:nrow(grid.info)]
    end



#= Step (B): Fomation opacities =#

    for i in 1:nrow(grid.info)
        formation_opacities(grid.info[i, "eos_root"], grid.info[i, "av_path"], grid.info[i, "logg"], grid.info[i, "name"])
    end

    ## Do the binning in parallel across many workers, this is the most time consuming part
    args = [
        (
            grid.info[i, "eos_root"], 
            grid.info[i, "av_path"], 
            grid.info[i, "logg"], 
            grid.info[i, "name"]
        ) for i in 1:nrow(grid.info)
    ]
    Distributed.pmap(bin_opacities, args)    

    ## Save the eos info
    grid.info[!, "name_extension"]   = [name_extension for _ in 1:nrow(grid.info)]
    grid.info[!, "version"]          = [version for _ in 1:nrow(grid.info)]
    grid.info[!, "binned_tables"]    = ["$(name_extension)_$(grid.info[i, "name"])_v$(version)" for i in 1:nrow(grid.info)]
    grid.info[!, "binned_E_tables"]  = ["$(name_extension)_E_$(grid.info[i, "name"])_v$(version)" for i in 1:nrow(grid.info)]
    grid.info[!, "rad_bins"]         = [7 for _ in 1:nrow(grid.info)]


    ## compute the resolution and the rounded size of the box
    ## use the EoS that was just created for this
    resolution!(grid, patch_size=15, cut_bottom=0.45)
    @show grid.info[!, "x_resolution"]  grid.info[!, "z_resolution"] grid.info[!, "x_size"]  grid.info[!, "z_size"]
  



#= Step (C): Conversion =#
  
    fromT_toE.(grid.info[!, "binned_tables"], grid.info[!, "binned_E_tables"])

    ## Copy the average model in the same folder so that we can link it all to the right place
    for i in 1:nrow(grid.info)
        cp(grid.info[i, "av_path"], joinpath(grid.info[i, "binned_E_tables"], "inim.dat"), force=true)
    end




#= Step (D): Input namelists for dispatch =#
    create_namelist!(grid)

    for i in 1:nrow(grid.info)
        cp(grid.info[i, "namelist_name"], joinpath(grid.info[i, "binned_E_tables"], "ininml.dat"), force=true)
    end

    MUST.save(grid, "dispatch_grid.mgrid")
end