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
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit a job to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_name::String; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict())
    isnothing(multi_location) && error("No Multi module has been loaded.")

    # Create the default namelist (with adjustments)
    nml = spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)
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





#= Effective temperature computations =#

flux(run, norm=true) = begin
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