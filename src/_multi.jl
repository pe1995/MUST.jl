#= Wrapper for the MULTI output =#

struct M3DISRun{P<:PythonCall.Py}
    run::P
end

M3DISRun(path::String) = begin
    isnothing(multi_location) && error("No Multi module has been loaded.")
    M3DISRun(m3dis.read(@in_m3dis(path)))
end






#= Utility functions =#

Base.getproperty(m::M3DISRun, arg::Symbol) = begin
    if arg in fieldnames(typeof(m))
        getfield(m, arg)
    else
        v = Base.getproperty(m.run, arg)
        try
            pyconvert.(Any, v)
        catch
            v
        end
    end
end





#= running M3DIS =#

"""
	whole_spectrum(model_name; [namelist_kwargs, m3dis_kwargs])

Submit a job to the M3DIS code, which will compute the outgoing flux across the entire wavelength range.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function whole_spectrum(model_name::String; namelist_kwargs=Dict(), m3dis_kwargs=Dict(), slurm=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")

    # Create the default namelist (with adjustments)
    nml = whole_spectrum_namelist(model_name; namelist_kwargs...)
    write(nml, @in_m3dis("$(model_name).nml"))

    # run multi (with waiting)
    if slurm
        srun_m3dis("$(model_name).nml"; wait=true, m3dis_kwargs...)
    else
        run_m3dis("$(model_name).nml"; wait=true, m3dis_kwargs...)
    end

    # read the output
    M3DISRun(joinpath(nml.io_params["datadir"], model_name))
end

"""
	whole_spectrum(model_names...; [namelist_kwargs, m3dis_kwargs])

Submit jobs to the M3DIS code, which will compute the outgoing flux across the entire wavelength range.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function whole_spectrum(model_names::AbstractVector{String}; namelist_kwargs=Dict(), m3dis_kwargs=Dict(), slurm=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")

    data_dir = whole_spectrum_namelist(model_names |> first; namelist_kwargs...).io_params["datadir"]

    if (!slurm) & (length(model_names) > 1)
        @warn "You are submitting $(length(model_names)) jobs without resource management (slurm=false)!"
    end

    # Create the default namelist (with adjustments)
    results = []
    for model_name in model_names
        nml = whole_spectrum_namelist(model_name; namelist_kwargs...)
        write(nml, @in_m3dis("$(model_name).nml"))

        # run multi (with waiting)
        if slurm
            append!(results, [srun_m3dis("$(model_name).nml"; wait=false, m3dis_kwargs...)])
        else
            append!(results, [run_m3dis("$(model_name).nml"; wait=false, m3dis_kwargs...)])
        end
    end

    # wait for success
    for (i, r) in enumerate(results)
        s = success(r)
        @info "$(model_names[i]) finished with success status $(s)."
    end

    # read the output
    [M3DISRun(joinpath(data_dir, model_name)) for model_name in model_names]
end



"""
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit a job to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_name::String; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")

    # Create the default namelist (with adjustments)

    # Check if NLTE and twostep, in that case we follow up with an LTE run
    nml = if NLTE & twostep
        nml_nlte = spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)
        nml_spectrum = spectrum_namelist(model_name; NLTE=false, namelist_kwargs...)
        path_nlte = joinpath(nml_nlte.io_params["datadir"], "$(model_name)_departure")

        set!(
            nml_spectrum; 
            atom_params=(:precomputed_depart=>path_nlte,),
            m3d_params=(:long_scheme=>"lobatto",)
        )
        set!(
            nml_nlte; 
            m3d_params=(:long_scheme=>"disk_center",)
        )

        # run the initial namelist first
        write(nml_nlte, @in_m3dis("$(model_name)_departure.nml"))
        @info "Running M3D (departure)."
        if slurm
            srun_m3dis("$(model_name)_departure.nml"; wait=true, m3dis_kwargs...)
        else
            run_m3dis("$(model_name)_departure.nml"; wait=true, m3dis_kwargs...)
        end
        @info "M3D completed (departure)."

        nml_spectrum
    else
        spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)  
    end

    write(nml, @in_m3dis("$(model_name).nml"))

    # run multi (with waiting)
    @info "Running M3D with NLTE=$(NLTE)."
    if slurm
        srun_m3dis("$(model_name).nml"; wait=true, m3dis_kwargs...)
    else
        run_m3dis("$(model_name).nml"; wait=true, m3dis_kwargs...)
    end
    @info "M3D completed."

    # read the output
    M3DISRun(joinpath(nml.io_params["datadir"], model_name))
end

"""
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit jobs to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_names::AbstractVector{String}; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")

    data_dir = whole_spectrum_namelist(model_names |> first; namelist_kwargs...).io_params["datadir"]

    if (!slurm) & (length(model_names) > 1)
        @warn "You are submitting $(length(model_names)) jobs without resource management (slurm=false)!"
    end
    
    # Create the default namelist (with adjustments)
    results = []
    results_nlte = []
    nml = []
    for model_name in model_names
        # Create the default namelist (with adjustments)
        nml_l = if NLTE & twostep
            path_nlte = joinpath(data_dir, "$(model_name)_departure")
            nml_nlte = spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)
            nml_spectrum = spectrum_namelist(model_name; NLTE=false, namelist_kwargs...)

            set!(
                nml_spectrum; 
                atom_params=(:precomputed_depart=>path_nlte,),
                m3d_params=(:long_scheme=>"lobatto",)
            )
            set!(
                nml_nlte; 
                m3d_params=(:long_scheme=>"disk_center",)
            )
    
            # run the initial namelist first
            write(nml_nlte, @in_m3dis("$(model_name)_departure.nml"))

            @info "Running M3D (departure)."
            r = if slurm
                srun_m3dis("$(model_name)_departure.nml"; wait=false, m3dis_kwargs...)
            else
                run_m3dis("$(model_name)_departure.nml"; wait=false, m3dis_kwargs...)
            end
            append!(results_nlte, [r])
            @info "M3D completed (departure)."
    
            nml_spectrum
        else
            spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)  
        end

        append!(nml, [nml_l])
    end

    # wait for success
    for (i, r) in enumerate(results_nlte)
        s = success(r)
        @info "$(model_names[i]) finished with success status $(s). (departure)"
    end


    for (i, model_name) in enumerate(model_names)
        write(nml[i], @in_m3dis("$(model_name).nml"))

        # run multi (with waiting)
        @info "Running M3D with NLTE=$(NLTE)."
        r = if slurm
            srun_m3dis("$(model_name).nml"; wait=false, m3dis_kwargs...)
        else
            run_m3dis("$(model_name).nml"; wait=false, m3dis_kwargs...)
        end

        append!(results, [r])
    end

    # wait for success
    for (i, r) in enumerate(results)
        s = success(r)
        @info "$(model_names[i]) finished with success status $(s)."
    end

    # read the output
    [M3DISRun(joinpath(data_dir, model_name)) for model_name in model_names]
end



"""
    heating(model_name, args...; kwargs...)

Compute the radiative heating within the given model.
"""
function heating(model_name::String, opacity_table=nothing; slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict())
    isnothing(multi_location) && error("No Multi module has been loaded.")

    binned_ext = isnothing(opacity_table) ? "_unbinned" : "_binned"

    # Create the default namelist (with adjustments)
    nml = heating_namelist(model_name, opacity_table; namelist_kwargs...)
    write(nml, @in_m3dis("$(model_name)$(binned_ext).nml"))

    # run multi (with waiting)
    if slurm
        srun_m3dis("$(model_name)$(binned_ext).nml"; wait=true, m3dis_kwargs...)
    else
        run_m3dis("$(model_name)$(binned_ext).nml"; wait=true, m3dis_kwargs...)
    end

    # read the output
    M3DISRun(joinpath(nml.io_params["datadir"], "$(model_name)$(binned_ext)"))
end





#= Effective temperature computations =#

window(line::PythonCall.Py, args...; kwargs...) = begin
    pyconvert.(Array, line.crop(args...; kwargs...))
end

window(run::M3DISRun, iline, args...; kwargs...) = begin
    line = run.line[iline]
    window(line, args...; kwargs...)
end





flux(run; norm=true) = begin
    norm ? (run.lam, run.flux ./ run.flux_cnt) : (run.lam, run.flux)
end

"""
    flux(run, norm=true)
Return and optionally normalize the flux coming from a m3dis run.
"""
flux


"""
    B(λ, T)
Planck function (cgs, λ).
"""
B(λ, T) = begin
    @. 2*HPlanck*CLight^2/(aa_to_cm*λ)^5 * 1 / 
                (exp(HPlanck*CLight /(aa_to_cm*λ*KBoltzmann*T)) -1 )
end

"""
    wien(T)
Maximum of Planck function (cgs).
"""
wien(T) = 2898e-4 /aa_to_cm / T



Teff(λ, F) = begin
    ν = reverse(CLight ./ (λ*aa_to_cm))
    Fl = reverse(F)

    (integrate(ν, Fl) /σ_S)^(1/4)
end

Teff(run) = Teff(run.lam, run.flux)

"""
    Teff(run)
Effective temperature from the M3D run.

    Teff(λ, F)
Effective temperature from the flux itself.
"""
Teff