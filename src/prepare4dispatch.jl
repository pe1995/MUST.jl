using MUST
using TSO
using DataFrames
using DelimitedFiles

#= Utilities =#

round_to_next(x, to_y) = begin
    x_new = x
    x_new = if (x_new % to_y) >= to_y/2
        div(x_new, to_y) * to_y + to_y
    else
        div(x_new, to_y) * to_y
    end
end

find_integer(f, dxz_frac, n_patches, to_y = 0.1; digits=1) = begin
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

@inline patches(res, patch_size) = Int(floor(res / patch_size))






#= Resolution =# 

"""
    resolution(av_model, min_x, max_x, min_z, max_z, hres, patch_size=20; cut_bottom=0.3)

Compute the needed resolution based on input parameters from other 3D simulation.
"""
function resolution(av_model, min_x, max_x, min_z, max_z, hres, 
                                patch_size=20; cut_bottom=0.3)
    nz = ceil(length(av_model.z) * (1-cut_bottom)) * 1.2
    dz = round((max_z - min_z) * (1-cut_bottom), sigdigits=2)
    target_res = round(dz/nz, sigdigits=2)
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

"""
    resolution(av_model, min_x, max_x, min_z, max_z, hres, patch_size=20; cut_bottom=0.3)

Compute the needed resolution based on input parameters from other 3D simulation.
Put more emphasis on the horizontal extend.
"""
function resolutionHD(av_model, min_x, max_x, min_z, max_z, hres, 
                                    patch_size=30; cut_bottom=0.3, scale_resolution=1.2)
    nx = ceil(abs(max_x - min_x)/abs(hres)) * scale_resolution
    dx = round((max_x - min_x) * (1-cut_bottom), sigdigits=3)
    target_res = round(dx/nx, sigdigits=3)
    nx = ceil(nx / patch_size) * patch_size

    ## how many patches do we need then
    n_patches = div(nx, patch_size)

    ## So we need to make sure that we have an integer multiple of this number of
    ## patches in the other dimensions
    dz = round((max_z - min_z) * (1-cut_bottom), sigdigits=3) 
    
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

    target_res = round(dx/nx, sigdigits=3)
    dx = round(nx * target_res, sigdigits=3)
    res_points = nz / nx
    dz = res_points * dx
    
    dx, nx, dz, nz
end







#= opacities =#

formation_opacities(args...; kwargs...) = TSO.compute_formation_opacities(args...; kwargs...)
bin_opacities(args...; kwargs...) = TSO.bin_opacity_table(args...; kwargs...)
fromT_toE(args...; kwargs...) = TSO.convert_fromT_toE(args...; kwargs...)






#= Namelist generation =#

function create_namelist(name, x_resolution, z_resolution, x_size, z_size, 
                        patch_size, rt_patch_size, δz, l_cgs, d_cgs, logg, teff, 
                        eemin, nz, initial_path, n_bin, eos_table, 
                        tscale, vmax, vmin)
    #ngrid = nrow(grid.info)
    #phase = grid.name

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

    larger_than_sun = x_size/l_cgs / 4.6
    newton_time = 100 #* larger_than_sun
    friction_time = 150 #* larger_than_sun
    newton_scale = 0.1 #* larger_than_sun

    l_cgs = MUST.roundto(l_cgs * larger_than_sun, 0.1, magnitude=1e3)

    # experiment with t scale based on size only
    test_tscale = true
    Δt(R, u, c=1.0) = c * R / u

    tscale = if test_tscale
        # Chose dt, such that with the given scaling c is a given value at the highest 
        # absolute velocity in the given stagger model. There will be higher 
        # velocities during the simulation, so there has to be room for larger
        # currants. Round to the next quater
        max(100.0, MUST.roundto(Δt(l_cgs, max(abs(vmax), abs(vmin)), 0.9), 0.25, magnitude=1e2))
    else
        #round(larger_than_sun*tscale, sigdigits=2)
        tscale
    end

    x = round(z_size/l_cgs, sigdigits=3) * 2
    MUST.set!(
        nml, 
        cartesian_params=(:size=>[x, x, round(z_size/l_cgs, sigdigits=3)], 
                        :dims=>[patches(x_resolution, patch_size), 
                                patches(x_resolution, patch_size), 
                                patches(z_resolution, patch_size)],
                        :position=>[0,0,round(-dup/l_cgs, sigdigits=3)]),
        patch_params=(:n=>[patch_size, patch_size, patch_size],),
        scaling_params=(:l_cgs=>l_cgs, :d_cgs=>d_cgs, :t_cgs=>tscale),
        stellar_params=(:g_cgs=>round(exp10(logg), digits=5), 
                        :ee_min_cgs=>round(log(eemin), digits=5), 
                        :nz=>nz-1, 
                        :initial_path=>initial_path),
        friction_params=(:end_time=>friction_time, :decay_scale=>10.0),
        gravity_params=(:constant=>-round(exp10(logg), digits=5),),
        newton_params=(:ee0_cgs=>round(log(eemin), digits=5), 
                        :position=>round((z_size/2 - 1.2*dup)/l_cgs, digits=2), 
                        :end_time=>newton_time, 
                        :decay_scale=>25.0,
                        :scale=>newton_scale),
        sc_rt_params=(  :rt_llc=>[-x, -x, -round((z_size/2 + dup)/l_cgs, sigdigits=3)], 
                        :rt_urc=>[ x,  x,  round((z_size/2 - dup)/l_cgs, sigdigits=3)], 
                        :n_bin=>n_bin,
                        :courant=>courant_rt,
                        :start_time=>newton_time,
                        :decay_scale=>25.0,
                        :rt_freq=>0.0,
                        :rt_res=>[-1,-1,rt_patch_size]),
        an_params=(:courant=>courant_hd,),
        eos_params=(:table_loc=>eos_table,)
    )
    
    # write the namelists
    MUST.write(nml, MUST.@in_dispatch(name))

    # save the summary
    MUST.@in_dispatch(name)
end







#= Modification of initial grid (interface) =#

resolution!(grid::MUST.AbstractMUSTGrid; 
                            patch_size=30, cut_bottom=0.3) = begin
    xr, zr = zeros(Int, nrow(grid.info)), zeros(Int, nrow(grid.info))
    xd, zd = zeros(Float64, nrow(grid.info)), zeros(Float64, nrow(grid.info))


    # At this point, it might be usefull to use the EoS
    # We should -> compute optical depth of model
    #           -> shift the box down by so much that roughy tau=-5 is the upper edge
    #           -> also save the density and the z at the optical surface, we can use those as scaling
    eos = [reload(SqEoS, joinpath(grid.info[i, "binned_tables"], "eos.hdf5")) for i in 1:nrow(grid.info)]
    opa = [reload(SqOpacity, joinpath(grid.info[i, "binned_tables"], "binned_opacities.hdf5")) for i in 1:nrow(grid.info)]

    models = []
    for i in 1:nrow(grid.info)
        append!(models, [@optical(Average3D(eos[i], grid.info[i, "av_path"], logg=grid.info[i, "logg"]), eos[i], opa[i])])
    end

    #models = [@optical(Average3D(eos[i], grid.info[i, "av_path"], logg=grid.info[i, "logg"]), eos[i], opa[i]) for i in 1:nrow(grid.info)]

    rho_norm = zeros(nrow(grid.info))
    l_norm = zeros(nrow(grid.info))
    z_up = zeros(nrow(grid.info))
    eemin = zeros(nrow(grid.info))

    z_lo = zeros(nrow(grid.info))
    d_lo = zeros(nrow(grid.info))
    T_lo = zeros(nrow(grid.info))

    τ_up = -4.5
    τ_surf = 0.0
    τ_down = 4.5

    for i in 1:nrow(grid.info)
        ## What should the resolution be
        xd[i], xr[i], zd[i], zr[i] = resolutionHD(models[i], 
                                                grid.info[i, "mi_x"], grid.info[i, "ma_x"], grid.info[i, "mi_z"], grid.info[i, "ma_z"], grid.info[i, "hres"],
                                                patch_size; cut_bottom=cut_bottom)
        ## Where is the optical surface
        it0 = argmin(abs.(log10.(models[i].τ) .- τ_surf))

        ## where is the upper edge (roughly)
        itup = argmin(abs.(log10.(models[i].τ) .- τ_up))

        ## Where can we start with the adiabat
        itlo = argmin(abs.(log10.(models[i].τ) .- τ_down))

        ## interpolator
        ip_r = MUST.linear_interpolation(log10.(models[i].τ), models[i].lnρ)
        ip_z = MUST.linear_interpolation(log10.(models[i].τ), models[i].z)
        ip_E = MUST.linear_interpolation(log10.(models[i].τ), models[i].lnEi)
        ip_T = MUST.linear_interpolation(log10.(models[i].τ), models[i].lnT)

        #rho_norm[i] = models[i].lnρ[it0] |> exp |> log10 |> round |> exp10 
        rho_norm[i] = ip_r(τ_surf) |> exp |> log10 |> round |> exp10 

        l_norm[i] = zd[i] |> abs |> log10 |> floor |> exp10
        
        #z_up[i] = abs(models[i].z[itup] - models[i].z[it0])
        z_up[i] = abs(ip_z(τ_up) - ip_z(τ_surf))

        #eemin[i] = exp.(models[i].lnEi[itup])
        eemin[i] = exp.(ip_E(τ_up))

        #z_lo[i] = round(models[i].z[itlo], sigdigits=5)
        z_lo[i] = round(ip_z(τ_down), sigdigits=5)

        #d_lo[i] = round(exp.(models[i].lnρ[itlo]), sigdigits=5)
        d_lo[i] = round(exp.(ip_r(τ_down)), sigdigits=5)

        #T_lo[i] = round(exp.(models[i].lnT[itlo]), sigdigits=5)
        T_lo[i] = round(exp.(ip_T(τ_down)), sigdigits=5)
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
                                g(i, "tscale"),
                                g(i, "vmax"),
                                g(i, "vmin"))
        append!(names, [name])
    end

    grid.info[!, "namelist_name"] = names
end
