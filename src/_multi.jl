#= Wrapper for the MULTI output =#

struct TUMULTRun{P<:PythonCall.Py}
    run::P
end

TUMULTRun(path::String; kwargs...) = begin
    isnothing(multi_location) && error("No Multi module has been loaded.")
    p = isdir(path) ? path : @in_m3dis(path)
    @assert isdir(p)
    TUMULTRun(tumult.read(p; kwargs...))
end

TUMULTRun(path::String, folder::String; kwargs...) = begin
    TUMULTRun(joinpath(folder, path); kwargs...)
end

# alias for old Multi package
M3DISRun = TUMULTRun








#= Utility functions =#

Base.getproperty(m::TUMULTRun, arg::Symbol) = begin
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




#= NLTE/LTE correction =#

_spread_NLTE_corrections!(correction_χ, correction_S, χ_fraction, S_fraction, lnρ, lnT, idistR, dist, σ) = begin
    for j in eachindex(lnρ)
        for i in eachindex(lnT)
            correction_χ[i, j, :] .= (χ_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
            correction_S[i, j, :] .= (S_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
        end
    end
end

"""
    NLTE_grid_correction(result1, result2; λ_new, lnρ_new, lnT_new)

Compute NLTE/LTE corrections for opacity and source function and interpolate the
result to the given λ, lnρ and lnT grids. This can be used to create a 1D NLTE
opacity table before binning.
"""
function NLTE_grid_correction(result_NLTE, result_LTE; λ_new, lnρ_new, lnT_new)
    λ = pyconvert(Array, result_NLTE.lam)
    χ_NLTE = pyconvert(Array, 
        result_NLTE.run.read_patch_save("chi_01", concat=true, fdim=0, lazy=false)[0]
    )
    S_NLTE = pyconvert(Array, 
        result_NLTE.run.read_patch_save("Snu_01", concat=true, fdim=0, lazy=false)[0]
    )
    χ_LTE = pyconvert(Array, 
		result_LTE.run.read_patch_save("chi_01", concat=true, fdim=0, lazy=false)[0]
	)
	S_LTE = pyconvert(Array, 
		result_LTE.run.read_patch_save("Snu_01", concat=true, fdim=0, lazy=false)[0]
	)

    lnρ = log.(pyconvert(Array, result_NLTE.run.rho))
	lnT = log.(pyconvert(Array, result_NLTE.run.temp))

    # Interpolate in wavelength first
    χ_fraction = zeros(eltype(lnT_new), length(lnT), length(λ_new)) 
    S_fraction = zeros(eltype(lnT_new), length(lnT), length(λ_new)) 
    frac = similar(λ)

    for i in eachindex(lnρ)
        # interpolate in wavelength
        frac .= χ_NLTE[:, 1, 1, i] ./ χ_LTE[:, 1, 1, i]
        f_χ_NLTE = linear_interpolation(λ, frac, extrapolation_bc=NaN)
        χ_fraction[i, :] .= f_χ_NLTE.(λ_new)
        χ_fraction[i, isnan.(@view(χ_fraction[i, :]))] .= 1.0

        frac .= S_NLTE[:, 1, 1, i] ./ S_LTE[:, 1, 1, i]
        f_S_NLTE = linear_interpolation(λ, frac, extrapolation_bc=NaN)
        S_fraction[i, :] .= f_S_NLTE.(λ_new)
        S_fraction[i, isnan.(@view(S_fraction[i, :]))] .= 1.0
    end

    # Interpolate in rho-T second. 
    # For now we just pick the closest point in rho, and apply a gaussian filter
    # based on the 2D distance from that point.
    correction_χ = zeros(eltype(lnT_new), length(lnT_new), length(lnρ_new), length(λ_new))
    correction_S = zeros(eltype(lnT_new), length(lnT_new), length(lnρ_new), length(λ_new))
    dist = zeros(eltype(lnT_new), length(lnT_new), length(lnρ_new))
    idist = zeros(Int, size(dist)...)
    idistR = zeros(Int, size(dist)...)
    
    for j in eachindex(lnρ_new)
        for i in eachindex(lnT_new)
            idist[i, j] = argmin(
                (lnρ_new[j]/4 .- lnρ/4).^2 .+ (lnT_new[i] .- lnT).^2 
            )
            idistR[i, j] = argmin(
                (lnρ_new[j] .- lnρ).^2 
            )
            dist[i, j] = (lnρ_new[j]/4 .- lnρ[idist[i, j]]/4).^2 .+ (lnT_new[i] .- lnT[idist[i, j] ]).^2 
        end
    end
    
    σ_rho = (abs(maximum(lnρ_new)/4 - minimum(lnρ_new)/4) / 16) ^2
    σ_T = (abs(maximum(lnT_new) - minimum(lnT_new)) / 16) ^2
    σ = σ_rho + σ_T
    
    _spread_NLTE_corrections!(
        correction_χ, correction_S, 
        χ_fraction, S_fraction, 
        lnρ_new, lnT_new, 
        idistR, dist, σ
    )

    correction_χ, correction_S
end








#= running M3DIS =#

"""
    abund_abundances(;α=0.0, default="./input_multi3d/abund_magg", eles...)

Create an abund file with the given abundance ratios.
"""
function abund_abundances(;α=0.0, default="./input_multi3d/abund/abund_magg", eles...)
	# read the default abundances
	abund_default = readdlm(@in_m3dis(default))
	abund_new = deepcopy(abund_default)
    ele_names = lowercase.(abund_new[:, 1])
	new_name = default*"_mod"
	if α != 0.0
		for ele in ["o", "ne", "mg", "si", "s", "ar", "ca"]
			iele = findfirst(lowercase(ele) .== ele_names)
			if !isnothing(iele)
				abund_new[iele, 2] = abund_default[iele, 2] + α
			else
				@warn "element $(ele) not found in absdat $(default)."
			end
		end
		new_name *= "_a$(α)"
	end

	for (eleS, val) in eles
		ele = lowercase(string(eleS))
		iele = findfirst(ele .== ele_names)
		if !isnothing(iele)
			abund_new[iele, 2] = abund_default[iele, 2] + val
			new_name *= "_$(ele)$(val)"
		else
			@warn "element $(string(eleS)) not found in absdat $(default)."
		end
	end

	if new_name != default
        abs_new_name_dir = dirname(@in_m3dis(default))
        random_file_in_dir = tempname(abs_new_name_dir)
		open(random_file_in_dir, "w") do f
			for i in axes(abund_new, 1)
				line = @sprintf "%-4s %-.3f\n" abund_new[i, 1] abund_new[i, 2]
				write(f, line)
			end
		end
        joinpath(dirname(default), basename(random_file_in_dir))
	else
		default
	end
end







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
    TUMULTRun(joinpath(nml.io_params["datadir"], model_name), read_atmos=false)
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
    [TUMULTRun(joinpath(data_dir, model_name), read_atmos=false) for model_name in model_names]
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
				:use_rho=>true, 
				:use_ne=>false
			),
			:m3d_params=>(
				:n_nu=>1, 
				:ilambd=>0,
				:short_scheme=>"disk_center",
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
function spectrum(model_name::String; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true, name="", cleanup=true)
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


    if cleanup
        nmls_created = glob("$(mn(model_name)).nml*", @in_m3dis(""))
        rm.(nmls_created)
    end


    # read the output
    TUMULTRun(joinpath(nml.io_params["datadir"], mn(model_name)))
end

"""
	spectrum(model_name; NLTE=false, slurm=true, [namelist_kwargs, m3dis_kwargs])

Submit jobs to the M3DIS code, which will compute the outgoing flux across in LTE or NLTE.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function spectrum(model_names::AbstractVector{String}; NLTE=false, slurm=true, namelist_kwargs=Dict(), m3dis_kwargs=Dict(), twostep=true, wait=true, name="", cleanup=true)
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

        if cleanup
            nmls_created = glob(join([model_name, name], "_")*".nml*", @in_m3dis(""))
            rm.(nmls_created)
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
        [TUMULTRun(joinpath(data_dir, model_name)) for model_name in model_names]

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
            r = TUMULTRun(joinpath(data_dir, join([model_name, names[j]], "_")))
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
    TUMULTRun(joinpath(nml.io_params["datadir"], "$(model_name)$(binned_ext)"))
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
    [TUMULTRun(joinpath(data_dir, "$(model_name)$(name)$(binned_ext)")) for model_name in model_names]
end







"""
    multimodel(model_name::String; namelist_kwargs=Dict(), m3dis_kwargs=Dict(), name="", cleanup=true)

Submit a job to the M3D code, which will compute the chi500 and store the atmosphere.
IMPORTANT: Make sure you have loaded m3dis in advance using the @import_m3dis macro.
"""
function multimodel(model_name::String; namelist_kwargs=Dict(), m3dis_kwargs=Dict(), name="", cleanup=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")

    # Create the default namelist (with adjustments)
    mn(n) = join([n, name], "_")

    # Check if NLTE and twostep, in that case we follow up with an LTE run
    nml = tau500_namelist(model_name; namelist_kwargs...)  
    write(nml, @in_m3dis("$(mn(model_name)).nml"))

    # run multi (with waiting)
    @info "Running M3D."
    run_m3dis("$(mn(model_name)).nml"; wait=true, m3dis_kwargs...)
    @info "M3D completed."

    if cleanup
        nmls_created = glob("$(mn(model_name)).nml*", @in_m3dis(""))
        rm.(nmls_created)
    end

    # read the output
    TUMULTRun(joinpath(nml.io_params["datadir"], mn(model_name)))
end















#= Effective temperature computations =#

window(line::PythonCall.Py, args...; kwargs...) = begin
    pyconvert.(Array, line.crop(args...; kwargs...))
end

window(run::TUMULTRun, iline, args...; kwargs...) = begin
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

Teff(run::TUMULTRun) = Teff(run.lam, run.flux)

"""
    Teff(run)
Effective temperature from the M3D run.

    Teff(λ, F)
Effective temperature from the flux itself.

    Teff(F)
Effective temperature from the binned flux itself.
"""
Teff





#= Create a 3D cube from a run =#

_is_1D(x) = ndims(x) == 1
_set_size_1D_to_3D(x) = _is_1D(x) ? reshape(x, 1, 1, length(x)) : x

function Box(result::TUMULTRun)
    x = Base.convert.(Float32, pyconvert(Array, result.run.xx) .* 1e8)
    y = Base.convert.(Float32, pyconvert(Array, result.run.yy) .* 1e8)
    z = Base.convert.(Float32, pyconvert(Array, result.run.zz) .* 1e8)

    data = Dict{Symbol, Any}(
        :T => _set_size_1D_to_3D(pyconvert(Array, result.run.temp)),
        :d => _set_size_1D_to_3D(pyconvert(Array, result.run.rho)),
        :ux => _set_size_1D_to_3D(pyconvert(Array, result.run.vx)),
        :uy => _set_size_1D_to_3D(pyconvert(Array, result.run.vy)),
        :uz => _set_size_1D_to_3D(pyconvert(Array, result.run.vz)),
        :τ500 => exp10.(pyconvert(Array, result.run.ltau)),
    )
    xx, yy, zz = meshgrid(x, y, z)
    Box(xx, yy, zz, data, AtmosphericParameters())
end