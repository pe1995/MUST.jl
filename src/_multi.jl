#= Wrapper for the MULTI output =#

struct M3DISRun{P<:PythonCall.Py}
    run::P
end

M3DISRun(path::String) = begin
    isnothing(multi_location) && error("No Multi module has been loaded.")
    M3DISRun(m3dis.read(@in_m3dis(path)))
end

M3DISRun(path::String, folder::String) = begin
    M3DISRun(joinpath(folder, path))
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

find_model(runs, kind, abundance) = begin
	models = [split(dirname(join(r.sfolder)), "/")[end-1] for r in runs]
	found = []
	for (i, model) in enumerate(models)
		
		if occursin(kind, model) 
			abund = pyconvert(Any, runs[i].atom.abnd)
			if abund ≈ abundance
				append!(found, [runs[i]])
			end
		end
	end

	if length(found) == 0
		@show kind, abundance
	end
	
	if length(found) == 1
		first(found)
	else
		found
	end
end

extension_results(extension, folder="data") = begin
	all_results = glob("*/", @in_m3dis(folder))
	paths = []
	
	for (i, r) in enumerate(all_results)
		if occursin(extension, r)
			append!(paths, [r])
		end
	end

	joinpath.(folder, last.(split.(dirname.(paths), "/")))
end





getabundances(run) = pyconvert(Float64, run.atom.abnd)


equivalentwidth(line::PythonCall.Py; LTE=false, kwargs...) = pyconvert(Any, line.calc_weq(; LTE=LTE, kwargs...))
equivalentwidth(lines::AbstractArray, abundances::AbstractArray; kwargs...) = begin
    m = sortperm(abundances)
    linear_interpolation(abundances[m], equivalentwidth.(lines; kwargs...)[m])
end
equivalentwidth(runs::AbstractArray; line=1, kwargs...) = begin
    lines = [run.line[line] for run in runs]
    abund = getabundances.(runs)
    equivalentwidth(lines, abund; kwargs...)
end


abundance(lines::AbstractArray, abundances::AbstractArray; kwargs...) = begin
    ew = equivalentwidth.(lines; kwargs...)
    m = sortperm(ew)
    linear_interpolation(ew[m], abundances[m])
end
abundance(runs::AbstractArray; line=1, kwargs...) = begin
    lines = [run.line[line] for run in runs]
    abund = getabundances.(runs)
    abundance(lines, abund; kwargs...)
end  


ΔNLTE(lines::AbstractArray, abundances::AbstractArray; at, reference=:LTE) = begin
    fLTE = equivalentwidth(lines, abundances, LTE=true)
    fNLTE = equivalentwidth(lines, abundances, LTE=false)

    aLTE = abundance(lines, abundances, LTE=true)
    aNLTE = abundance(lines, abundances, LTE=false)

    if reference == :LTE
        aNLTE(at |> fLTE) - at
    elseif reference == :NLTE
        at - aLTE(at |> fNLTE)
    else
        error("Pass either :LTE or :NLTE as reference.")
    end
end

ΔNLTE(runs::AbstractArray; at, reference=:LTE, line=1) = begin
    lines = [run.line[line] for run in runs]
    abund = getabundances.(runs)

    ΔNLTE(lines, abund; at=at, reference=reference)
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

opacityTable(models; folder, linelist, λs, λe, δλ, H_atom="input_multi3d/atoms/atom.h20",
				in_log=true, slurm=false, kwargs...) = begin
    whole_spectrum(
		models, 
		namelist_kwargs=(
			:model_folder=>folder,
			:linelist=>nothing,
			:absmet=>nothing,
			:linelist_params=>(:line_lists=>linelist,),
			:atom_params=>(:atom_file=>H_atom, ),
			:spectrum_params=>(:daa=>δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>in_log),
			:atmos_params=>(
				:dims=>1, 
				:atmos_format=>"Text",
				:use_density=>true, 
				:use_ne=>false
			),
			:m3d_params=>(
				:n_nu=>1, 
				:ilambd=>0,
				:quad_scheme=>"disk_center",
				:long_scheme=>"disk_center",
				:make_opac_table=>true
			),
            kwargs...
		),
		slurm=slurm
	)
end








"""
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit a job to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_name::String; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true, name="")
    isnothing(multi_location) && error("No Multi module has been loaded.")

    # Create the default namelist (with adjustments)
    mn(n) = join([n, name], "_")

    # Check if NLTE and twostep, in that case we follow up with an LTE run
    nml = if NLTE & twostep
        nml_nlte = spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)
        nml_spectrum = spectrum_namelist(model_name; NLTE=false, namelist_kwargs...)
        path_nlte = joinpath(nml_nlte.io_params["datadir"], "$(mn(model_name))_departure")

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
        write(nml_nlte, @in_m3dis("$(mn(model_name))_departure.nml"))
        @info "Running M3D (departure)."
        if slurm
            srun_m3dis("$(mn(model_name))_departure.nml"; wait=true, m3dis_kwargs...)
        else
            run_m3dis("$(mn(model_name))_departure.nml"; wait=true, m3dis_kwargs...)
        end
        @info "M3D completed (departure)."

        nml_spectrum
    else
        spectrum_namelist(model_name; NLTE=NLTE, namelist_kwargs...)  
    end

    write(nml, @in_m3dis("$(mn(model_name)).nml"))

    # run multi (with waiting)
    @info "Running M3D with NLTE=$(NLTE)."
    if slurm
        srun_m3dis("$(mn(model_name)).nml"; wait=true, m3dis_kwargs...)
    else
        run_m3dis("$(mn(model_name)).nml"; wait=true, m3dis_kwargs...)
    end
    @info "M3D completed."

    # read the output
    M3DISRun(joinpath(nml.io_params["datadir"], mn(model_name)))
end

"""
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit jobs to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_names::AbstractVector{String}; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true, wait=true, name="")
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
            path_nlte = joinpath(data_dir, join(["$(model_name)", name, "departure"], "_"))
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
            write(nml_nlte, @in_m3dis(join(["$(model_name)", name, "departure"], "_")*".nml"))

            @info "Running M3D (departure)."
            r = if slurm
                srun_m3dis(join(["$(model_name)", name, "departure"], "_")*".nml"; wait=false, m3dis_kwargs...)
            else
                run_m3dis(join(["$(model_name)", name, "departure"], "_")*".nml"; wait=false, m3dis_kwargs...)
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
        write(nml[i], @in_m3dis(join([model_name, name], "_")*".nml"))

        # run multi (with waiting)
        @info "Running M3D with NLTE=$(NLTE)."
        r = if slurm
            srun_m3dis(join([model_name, name], "_")*".nml"; wait=false, m3dis_kwargs...)
        else
            run_m3dis(join([model_name, name], "_")*".nml"; wait=false, m3dis_kwargs...)
        end

        append!(results, [r])
    end

    # wait for success
    if wait
        for (i, r) in enumerate(results)
            s = success(r)
            @info "$(model_names[i]) finished with success status $(s)."
        end
        # read the output
        [M3DISRun(joinpath(data_dir, model_name)) for model_name in model_names]

    else
        results
    end
end

"""
    spectrum(model_names::AbstractVector{String}, namelist_kwargs::AbstractVector; NLTE=false, slurm=true, m3dis_kwargs=Dict(), twostep=true)

Run `spectrum` for every model in `model_names` with different `namelist_kwargs`.
"""
function spectrum(model_names::AbstractVector{String}, runnames::Dict; kwargs...)
    names = keys(runnames) |> collect
    namelist_kwargs = [runnames[k] for k in names]

    results = []
    for (i, input_params) in enumerate(namelist_kwargs)
        #@show input_params names[i]
        r = spectrum(model_names; namelist_kwargs=input_params, name=names[i], wait=false, kwargs...)
        append!(results, r)
    end

    for (i, r) in enumerate(results)
        s = success(r)
    end

    runs = []
    data_dir = whole_spectrum_namelist(model_names |> first; first(namelist_kwargs)...).io_params["datadir"]

    for (i, model_name) in enumerate(model_names)
        for (j, input_params) in enumerate(namelist_kwargs)
            r = M3DISRun(joinpath(data_dir, join([model_name, names[j]], "_")))
            append!(runs, [r])
        end
    end

    runs
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

"""
    heating(model_name, args...; kwargs...)

Compute the radiative heating within the given models in parallel.
"""
function heating(model_names::AbstractVector{String}, opacity_table=nothing; namelist_kwargs=Dict(), m3dis_kwargs=Dict(), slurm=true, name="")
    isnothing(multi_location) && error("No Multi module has been loaded.")

    binned_ext = isnothing(opacity_table) ? "_unbinned" : "_binned"
    data_dir = heating_namelist(model_names |> first, opacity_table; namelist_kwargs...).io_params["datadir"]

    if (!slurm) & (length(model_names) > 1)
        @warn "You are submitting $(length(model_names)) jobs without resource management (slurm=false)!"
    end

    # Create the default namelist (with adjustments)
    results = []
    for model_name in model_names
        nml = heating_namelist(model_name, opacity_table; namelist_kwargs...)
        write(nml, @in_m3dis("$(model_name)$(name)$(binned_ext).nml"))

        # run multi (with waiting)
        if slurm
            append!(results, [srun_m3dis("$(model_name)$(name)$(binned_ext).nml"; wait=false, m3dis_kwargs...)])
        else
            append!(results, [run_m3dis("$(model_name)$(name)$(binned_ext).nml"; wait=false, m3dis_kwargs...)])
        end
    end

    # wait for success
    for (i, r) in enumerate(results)
        s = success(r)
        @info "$(model_names[i])$(name)$(binned_ext) finished with success status $(s)."
    end

    # read the output
    [M3DISRun(joinpath(data_dir, "$(model_name)$(name)$(binned_ext)")) for model_name in model_names]
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



Teff(λ::AbstractVector{T}, F::AbstractVector{T}) where {T<:AbstractFloat} = begin
    ν = reverse(CLight ./ (λ*aa_to_cm))
    Fl = reverse(F)

    (integrate(ν, Fl) /σ_S)^(1/4)
end

Teff(F::AbstractVector{T}) where {T<:AbstractFloat} = (sum(F) /σ_S)^(1/4)

Teff(run::M3DISRun) = Teff(run.lam, run.flux)

"""
    Teff(run)
Effective temperature from the M3D run.

    Teff(λ, F)
Effective temperature from the flux itself.

    Teff(F)
Effective temperature from the binned flux itself.
"""
Teff