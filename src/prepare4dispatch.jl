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

find_integer_simple(f, dxdz_approx, desired_n_patches; step=0.1) = begin
    actual_x_patches = desired_n_patches * dxdz_approx
    while !isinteger(actual_x_patches)
        dxdz_approx = f(dxdz_approx, step)
        dxdz_approx = round(dxdz_approx, sigdigits=2)
        actual_x_patches = desired_n_patches * dxdz_approx
    end

    actual_x_patches
end

find_integer_simple(dxdz_approx, desired_n_patches; step=0.1) = begin
    p = find_integer_simple(+, dxdz_approx, desired_n_patches)
    m = find_integer_simple(-, dxdz_approx, desired_n_patches)
    a = dxdz_approx * desired_n_patches

    [p, m][argmin([abs(p - a), abs(m - a)])]
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


"""
    resolution(av_model, min_x, max_x, min_z, max_z, hres, patch_size=20; cut_bottom=0.3)

Compute the needed resolution based on input parameters from other 3D simulation.
Make everything easier and avoid rounding.
"""
function resolutionSimple(av_model, min_x, max_x, min_z, max_z, τ_top, τ_surf, τ_bottom, hres,
                                    patch_size=30; scale_resolution=1.0)
    dx = abs(max_x - min_x)
    dz = abs(max_z - min_z)
    dxdz = dx / dz

    # This is the rt resolution we need. In HD we use half the points
    vres = minimum(abs.(diff(av_model.z))) *2.0 

    mask = sortperm(av_model.τ)
    ip_z = MUST.linear_interpolation(
        MUST.Interpolations.deduplicate_knots!(log10.(av_model.τ[mask])), 
        av_model.z[mask], 
        extrapolation_bc=MUST.Flat()
    )

    actual_top = ip_z(τ_top)
    actual_bottom = ip_z(τ_bottom)
    actual_dz = abs(actual_top - actual_bottom)

    # now we need to add so many patches of the given size that is matches at least
    # the desired resolution
    desired_resolution = vres / scale_resolution
    desired_n_points = actual_dz / desired_resolution
    desired_n_patches = ceil(desired_n_points / patch_size)

    # we need to make sure that x patches are integer, so we need to round here
    dxdz_approx = round(dxdz, sigdigits=2)
    step = 0.1

    actual_x_patches = find_integer_simple(dxdz_approx, desired_n_patches, step=step)
    
    # For some reason the fraction needs to be integer as well...
    actual_x_patches = round(dxdz, sigdigits=1) * desired_n_patches
    actual_dx = actual_x_patches / desired_n_patches * actual_dz

    actual_dx, actual_x_patches*patch_size, actual_dz, desired_n_patches*patch_size
end







#= opacities =#

formation_opacities(args...; kwargs...) = TSO.compute_formation_opacities(args...; kwargs...)
bin_opacities(args...; kwargs...) = TSO.bin_opacity_table(args...; kwargs...)
fromT_toE(args...; kwargs...) = TSO.convert_fromT_toE(args...; kwargs...)
compute_rosseland_opacities(args...; kwargs...) = TSO.compute_rosseland_opacities(args...; kwargs...)


quadrantlimit(name, table_folder, opa_path; λ_lim=5.0) = begin
    opa = reload(
        SqOpacity, 
		joinpath(table_folder, opa_path), 
		mmap=true
    )
	
    fopa = reload(
        SqOpacity, 
		joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"), 
		mmap=true
    )

    TSO.median(-log10.(fopa.κ_ross)[log10.(opa.λ) .> λ_lim])
end

clean(table_folder, name) = begin
	rm(joinpath(table_folder, "combined_formation_opacities_$(name).hdf5"))
end





#= Namelist generation =#

function create_namelist(name, x_resolution, z_resolution, x_size, z_size, 
                        patch_size, rt_patch_size, δz, zs, zh, z_lo, l_cgs, d_cgs, logg, teff, 
                        eemin, ee0, zee0, nz, initial_path, n_bin, eos_table, 
                        tscale, vmax, vmin; 
                        courant_rt=0.2, courant_hd=0.2, newton_time=100, friction_time=150,
                        newton_scale=0.1, newton_decay_scale=15.0, friction_decay_scale=10.0,
                        courant_target=0.25, kwargs...)
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
    # this is how much the cube needs to be shifted down so that the amount of star above the surface is equal
    dup = (z_size /2 -δz) *1.0

    # z_lo should be the bottom z value. The same should be true for the simulation domain
    # this means that we need to shift it down until that value is reached
    ddown = abs((z_size/2.0 - abs(z_lo)))
    
    #=courant_rt = if (logg >= 4) & (teff<6500)
        1.0
    elseif ((logg < 4) & (logg >= 3)) | (teff>=6500)
        0.8
    else
        0.7
    end=#
    courant_rt = courant_rt

    #=courant_hd = if logg >= 4
        0.3
    elseif (logg < 4) & (logg >= 3)
        0.25
    else
        0.2
    end=#
    courant_hd = courant_hd


    larger_than_sun = x_size /1.0e8 /6.0
    newton_time = newton_time #* larger_than_sun
    friction_time = friction_time #* larger_than_sun
    newton_scale = newton_scale #* larger_than_sun

    l_cgs_raw = 1e8 * larger_than_sun
    l_cgs = round(1e8 * larger_than_sun, sigdigits=3) #MUST.roundto(l_cgs * larger_than_sun, 0.1, magnitude=1e3)

    # experiment with t scale based on size only
    test_tscale = true
    Δt(R, u, c=1.0) = c * R / u

    velocity_ratio = max(abs(vmax), abs(vmin)) / 1e6
    dynamic_scale_ratio = velocity_ratio / larger_than_sun

    tscale = if test_tscale
        # Chose dt, such that with the given scaling c is a given value at the highest 
        # absolute velocity in the given stagger model. There will be higher 
        # velocities during the simulation, so there has to be room for larger
        # currants. Round to the next quater
        dx = x_size / (patches(x_resolution, patch_size) * patch_size)
        round(Δt(dx, max(abs(vmax), abs(vmin)), courant_target) / 5e-3, sigdigits=3)
    else
        tscale
    end

    stellar_w = round(0.1 * velocity_ratio, sigdigits=3)
    strength = round(0.1 * 100.0 / tscale / velocity_ratio, sigdigits=5)  #round(0.1 / velocity_ratio, sigdigits=3)
    tbot = round(0.01 * 100.0 / tscale / velocity_ratio, sigdigits=5)  #round(0.1 / velocity_ratio, sigdigits=3)

    x = round(z_size/l_cgs_raw, sigdigits=3) * patches(x_resolution, patch_size) / patches(z_resolution, patch_size)
    MUST.set!(
        nml, 
        cartesian_params=(:size=>[x, x, round(z_size/l_cgs_raw, sigdigits=3)], 
                        :dims=>[patches(x_resolution, patch_size), 
                                patches(x_resolution, patch_size), 
                                patches(z_resolution, patch_size)],
                        :position=>[0,0,round(-ddown/l_cgs_raw, sigdigits=3)]),
        patch_params=(:n=>[patch_size, patch_size, patch_size], :grace=>0.2),
        scaling_params=(:l_cgs=>l_cgs, :d_cgs=>d_cgs, :t_cgs=>tscale),
        experiment_params=(:t_bot=>tbot,),
        stellar_params=(:g_cgs=>round(exp10(logg), sigdigits=5), 
                        :ee_min_cgs=>round(log(eemin), sigdigits=7), 
                        :nz=>nz-1, 
                        :w_perturb=>stellar_w,
                        :initial_path=>initial_path),
        friction_params=(:end_time=>friction_time, :decay_scale=>friction_decay_scale, :time=>strength,),
        gravity_params=(:constant=>-round(exp10(logg), sigdigits=5),),
        newton_params=(:ee0_cgs=>round(log(ee0), sigdigits=7), 
                        :position=>round(zee0/l_cgs_raw, sigdigits=3),#round((z_size/2 - 1.2*dup)/l_cgs_raw, sigdigits=3), 
                        :end_time=>newton_time, 
                        :decay_scale=>newton_decay_scale,
                        :time=>strength,
                        :scale=>newton_scale),
        sc_rt_params=(  :rt_llc=>[-x/2, -x/2, -round((z_size/2 + ddown)/l_cgs_raw, sigdigits=3)], 
                        :rt_urc=>[ x/2,  x/2,  round((z_size/2 - ddown)/l_cgs_raw, sigdigits=3)], 
                        :n_bin=>n_bin,
                        :courant=>courant_rt,
                        :rt_freq=>0.0,
                        :rt_grace=>0.1,
                        :rt_res=>[-1,-1,rt_patch_size]),
        an_params=(:courant=>courant_hd,),
        eos_params=(:table_loc=>eos_table, :gamma=>1.666667)
    )

    # set additional kwargs if wanted
    MUST.set!(nml; kwargs...)
    
    # write the namelists
    MUST.write(nml, MUST.@in_dispatch(name))

    # save the summary
    MUST.@in_dispatch(name)
end




#= adiabats =#

#================================================================== Adiabats =#

function adiabat(eospath, av_path, logg; common_size=1000, saveat=nothing, kwargs...)
	eos = reload(
        SqEoS,
        joinpath(eospath)
    )

	data = TSO.flip(Average3D(av_path, logg=logg))
	start_point = TSO.pick_point(data, 1)
	end_point = TSO.pick_point(data, length(data.z))
	
	a = TSO.upsample(TSO.adiabat(start_point, end_point, eos; kwargs...), common_size)
	a = TSO.flip(a, depth=false)
	if !isnothing(saveat)
		open(saveat, "w") do f
			MUST.writedlm(f, [a.z exp.(a.lnT) a.lnρ])
		end
	end

	a
end

function new_zscale(eospath, av_path, logg; common_size=1000, saveat=nothing, kwargs...)
    eos = reload(
        SqEoS,
        joinpath(eospath)
    )
    aos = @axed eos

    # interpolate to equally spaced in z
	# for this we need to construct a new z scale 
	mnew = TSO.flip(@optical(Average3D(aos, av_path, logg=logg), aos), depth=true)

    a = if maximum(log10.(mnew.τ)) < 0.0
        @warn "New Z scale could not be computed for $(basename(av_path))."
        mnew
    else
        # recompute z scale from opacity
        znew = TSO.rosseland_depth(aos, mnew)
        mnew.z .= znew

        # find the optical surface
        TSO.optical_surface!(mnew)
        TSO.flip!(mnew)

        # We now can interpolate everthing to this new z scale
        TSO.flip(
            TSO.interpolate_to(
                mnew, 
                z=collect(range(maximum(mnew.z), minimum(mnew.z), length=common_size))
            ), 
            depth=false
        )
    end

    if !isnothing(saveat)
		open(saveat, "w") do f
			MUST.writedlm(f, [a.z exp.(a.lnT) a.lnρ])
		end
	end

    a
end




#= Modification of initial grid (interface) =#

resolution!(grid::MUST.AbstractMUSTGrid; 
                            patch_size=30, scale_resolution=1.0, 
                            τ_up=-4.5, τ_surf=0.0, τ_down=5.5,
                            τ_ee0=-1.0, τ_eemin=-1.0, τ_zee0=-1.0, τ_rho0=-2.0, use_inim=false) = begin
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
        if !use_inim
            append!(models, [@optical(Average3D(eos[i], grid.info[i, "av_path"], logg=grid.info[i, "logg"]), eos[i], opa[i])])
        else
            append!(models, [@optical(Average3D(eos[i], joinpath(grid.info[i, "binned_E_tables"], "inim.dat"), logg=grid.info[i, "logg"]), eos[i], opa[i])])
        end
    end

    rho_norm = zeros(nrow(grid.info))
    l_norm = zeros(nrow(grid.info))
    z_up = zeros(nrow(grid.info))
    eemin = zeros(nrow(grid.info))
    ee0 = zeros(nrow(grid.info))
    zee0 = zeros(nrow(grid.info))


    z_lo = zeros(nrow(grid.info))
    z_s = zeros(nrow(grid.info))
    z_h = zeros(nrow(grid.info))
    d_lo = zeros(nrow(grid.info))
    T_lo = zeros(nrow(grid.info))

    for i in 1:nrow(grid.info)
        ## What should the resolution be
        #xd[i], xr[i], zd[i], zr[i] = resolutionHD(models[i], 
        #                                        grid.info[i, "mi_x"], grid.info[i, "ma_x"], grid.info[i, "mi_z"], grid.info[i, "ma_z"], grid.info[i, "hres"],
        #                                        patch_size; cut_bottom=cut_bottom)
        ## Where is the optical surface
        #it0 = argmin(abs.(log10.(models[i].τ) .- τ_surf))
        ### where is the upper edge (roughly)
        #itup = argmin(abs.(log10.(models[i].τ) .- τ_up))
        ### Where can we start with the adiabat
        #itlo = argmin(abs.(log10.(models[i].τ) .- τ_down))

        ## interpolator
        m = TSO.flip(models[i])
        mask = sortperm(m.τ)
        τ_down = min(τ_down, maximum(log10.(m.τ)))

        ip_r = MUST.linear_interpolation(
            MUST.Interpolations.deduplicate_knots!(log10.(m.τ[mask])),
            m.lnρ[mask], 
            extrapolation_bc=MUST.Flat()
        )
        ip_z = MUST.linear_interpolation(
            MUST.Interpolations.deduplicate_knots!(log10.(m.τ[mask])),
            m.z[mask], 
            extrapolation_bc=MUST.Flat()
        )
        ip_E = MUST.linear_interpolation(
            MUST.Interpolations.deduplicate_knots!(log10.(m.τ[mask])),
            m.lnEi[mask], 
            extrapolation_bc=MUST.Flat()
        )
        ip_T = MUST.linear_interpolation(
            MUST.Interpolations.deduplicate_knots!(log10.(m.τ[mask])),
            m.lnT[mask], 
            extrapolation_bc=MUST.Flat()
        )

        rho_norm[i] = round(exp.(ip_r(τ_rho0)), sigdigits=5) 

        l_norm[i] = zd[i] |> abs 
        
        z_up[i] = abs(ip_z(τ_up) - ip_z(τ_surf))
        z_s[i] = ip_z(τ_surf)
        z_h[i] = ip_z(τ_up)

        ee0[i] = exp.(ip_E(τ_ee0))
        zee0[i] = ip_z(τ_zee0)
        eemin[i] = exp.(ip_E(τ_eemin)) 

        z_lo[i] = ip_z(τ_down)
        d_lo[i] = exp.(ip_r(τ_down))
        T_lo[i] = exp.(ip_T(τ_down))

        xd[i], xr[i], zd[i], zr[i] = resolutionSimple(
            m, 
            grid.info[i, "mi_x"], 
            grid.info[i, "ma_x"], 
            grid.info[i, "mi_z"], 
            grid.info[i, "ma_z"], 
            τ_up,
            τ_surf,
            τ_down,
            grid.info[i, "hres"],
            patch_size, 
            scale_resolution=scale_resolution
        )
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
    grid.info[!, "z_s"] = z_s
    grid.info[!, "z_h"] = z_h
    grid.info[!, "ee_min"] = eemin
    grid.info[!, "initial_model_size"] = [length(models[i].z) for i in 1:nrow(grid.info)]

    grid.info[!, "ee0"] = ee0
    grid.info[!, "zee0"] = zee0

    grid.info[!, "z_lo"] = z_lo
    grid.info[!, "d_lo"] = d_lo
    grid.info[!, "T_lo"] = T_lo

    nothing
end

create_namelist!(grid::MUST.StaggerGrid; kwargs...) = begin
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
                                g(i, "z_s"),
                                g(i, "z_h"),
                                g(i, "z_lo"),
                                g(i, "l_norm"),
                                g(i, "rho_norm"),
                                g(i, "logg"),
                                g(i, "teff"),
                                g(i, "ee_min"),
                                g(i, "ee0"),
                                g(i, "zee0"),
                                g(i, "initial_model_size"),
                                joinpath("input_data/grd/", 
                                            g(i, "binned_E_tables"), 
                                            "inim.dat"),
                                g(i, "rad_bins"),
                                joinpath("input_data/grd/", 
                                            g(i, "binned_E_tables")),
                                g(i, "tscale"),
                                g(i, "vmax"),
                                g(i, "vmin");
                                kwargs...)
        append!(names, [name])
    end

    grid.info[!, "namelist_name"] = names
end
