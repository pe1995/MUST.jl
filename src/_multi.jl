# ============================================================================
# Wrapper for the MULTI output
# ============================================================================
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

M3DISRun = TUMULTRun

# ============================================================================
# Utility functions
# ============================================================================

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

run_folder(m::TUMULTRun) = pyconvert(String, pystr(m.run.folder))

# ============================================================================
# Running OpenMULTI3D
# ============================================================================

"""
    synthesize(name, settings=Dict(); presets=(), threads=nothing, overwrite=true)

Run OpenMULTI3D with `tumult.launch` in the installation loaded with `@import_tumult`.
The namelist is the tumult `presets` (e.g. `("1d",)` or `("3d", "nlte")`) updated
with `settings`, a `Dict` of namelist group => `Dict` of key => value. The run is
`data/<name>` of the installation.
"""
function synthesize(name, settings=Dict(); presets=(), threads=nothing, overwrite=true)
    isnothing(multi_location) && error("No Multi module has been loaded.")
    preset_dir = joinpath(dirname(pyconvert(String, tumult.__file__)), "presets")
    defaults = pylist([joinpath(preset_dir, "$(p).nml") for p in presets])
    out = try
        tumult.launch(name; nml_settings=_pynamelist(settings), default=defaults,
            overwrite=overwrite, threads=threads, cwd=multi_location, verbose=1)[0]
    finally
        rm(joinpath(multi_location, "$(name).nml"), force=true)
        d = dirname(joinpath(multi_location, name))
        while (normpath(d) != normpath(multi_location)) && isdir(d) && isempty(readdir(d))
            rm(d)
            d = dirname(d)
        end
    end
    TUMULTRun(pyconvert(String, pystr(out)))
end

"""
    spectrum(atmos_file; format, name, composition=Dict(), settings=Dict(), presets=("1d",), NLTE=false, threads=nothing, overwrite=true)

Spectrum of the model atmosphere `atmos_file` in the OpenMULTI3D `atmos_format`
`format` (e.g. "text", "marcs", "multi"), for `composition` as in
`AtmosphericParameters` (see `composition_settings`). `settings` update the
namelist of the `presets`, with the "nlte" preset added for `NLTE`. In NLTE, the
lines of the model atoms in the line list get its departure coefficients
(`match_nlte_levels`), as line lists like VALD have no level indices.
"""
function spectrum(atmos_file::AbstractString; format, name, composition=Dict(), settings=Dict(), presets=("1d",), NLTE=false, threads=nothing, overwrite=true)
    s = merge_settings(
        composition_settings(composition),
        Dict("atmos_params" => Dict("atmos_file" => abspath(atmos_file), "atmos_format" => format)),
        NLTE ? Dict("linelist_params" => Dict("match_nlte_levels" => true)) : Dict(),
        settings
    )
    synthesize(name, s; presets=(presets..., (NLTE ? ("nlte",) : ())...), threads=threads, overwrite=overwrite)
end

const _ALPHA_ELEMENTS = ("O", "Ne", "Mg", "Si", "S", "Ar", "Ca", "Ti")

"""
    composition_settings(composition)

Namelist settings for a composition as in `AtmosphericParameters`: `[Fe/H]`,
`[alpha/Fe]` and `[X/Fe]` of individual elements. The [X/Fe] are the final values:
OpenMULTI3D adds the element offsets on top of its [Fe/H] scaling (which skips He and
Li) and of [alpha/Fe] (O, Ne, Mg, Si, S, Ar, Ca, Ti), so those are compensated.
Entries that are missing are left to OpenMULTI3D (e.g. the [Fe/H] of a MARCS model).
"""
function composition_settings(composition)
    c = Dict(replace(String(k), "α" => "alpha") => Float64(v) for (k, v) in composition)
    feh = get(c, "[Fe/H]", 0.0)
    alpha = get(c, "[alpha/Fe]", 0.0)
    elements = sort([
        (uppercasefirst(lowercase(m[1])), v)
        for (k, v) in c for m in (Base.match(r"^\[(\w+)/Fe\]$", k),)
        if !isnothing(m) && (m[1] != "alpha")
    ])
    offset(e, v) = v + (e in ("He", "Li") ? feh : 0.0) - (e in _ALPHA_ELEMENTS ? alpha : 0.0)
    s = Dict{String, Dict{String, Any}}("composition_params" => Dict{String, Any}(
        "abundances" => join(["$(e)=$(offset(e, v))" for (e, v) in elements], " "),
        "relative_abundances" => true
    ))
    haskey(c, "[Fe/H]") && (s["atmos_params"] = Dict{String, Any}("FeH" => feh))
    haskey(c, "[alpha/Fe]") && (s["composition_params"]["alpha_enhancement"] = c["[alpha/Fe]"])
    s
end

"""
    merge_settings(settings...)

Merge namelist settings group by group, later ones take precedence.
"""
function merge_settings(settings...)
    m = Dict{String, Dict{String, Any}}()
    for s in settings, (g, group) in s
        merge!(get!(m, String(g), Dict{String, Any}()), Dict(String(k) => v for (k, v) in Dict(group)))
    end
    m
end

_pyvalue(v::AbstractVector) = pylist(v)
_pyvalue(v) = v
_pynamelist(settings) = pydict(Dict(
    String(g) => pydict(Dict(String(k) => _pyvalue(v) for (k, v) in Dict(group)))
    for (g, group) in settings
))

# ============================================================================
# Equivalent width & abundance calculations from Tumult
# ============================================================================

_line(run::TUMULTRun, i) = run.run.atoms[0].bbtrans[i-1]

getabundances(run::TUMULTRun) = pyconvert(Float64, run.run.atoms[0].abnd)

equivalentwidth(line::PythonCall.Py; LTE=false, kwargs...) = pyconvert(
    Float64, line.crop(; LTE=LTE).get_weq(; target_unit=pyimport("astropy.units").AA, kwargs...).value
)
equivalentwidth(lines::AbstractArray, abundances::AbstractArray; kwargs...) = begin
    m = sortperm(abundances)
    linear_interpolation(abundances[m], equivalentwidth.(lines; kwargs...)[m])
end
equivalentwidth(runs::AbstractArray; line=1, kwargs...) = begin
    lines = [_line(run, line) for run in runs]
    abund = getabundances.(runs)
    equivalentwidth(lines, abund; kwargs...)
end

abundance(lines::AbstractArray, abundances::AbstractArray; kwargs...) = begin
    ew = equivalentwidth.(lines; kwargs...)
    m = sortperm(ew)
    linear_interpolation(ew[m], abundances[m])
end
abundance(runs::AbstractArray; line=1, kwargs...) = begin
    lines = [_line(run, line) for run in runs]
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
    lines = [_line(run, line) for run in runs]
    abund = getabundances.(runs)

    ΔNLTE(lines, abund; at=at, reference=reference)
end

# ============================================================================
# NLTE / LTE corrections
# ============================================================================

_spread_NLTE_corrections!(correction_χ, correction_S, χ_fraction, S_fraction, lnρ, lnT, idistR, dist, σ) = begin
    for j in eachindex(lnρ)
        for i in eachindex(lnT)
            correction_χ[i, j, :] .= (χ_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
            correction_S[i, j, :] .= (S_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
        end
    end
end

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

    lnρ = log.(pyconvert(Array, result_NLTE.run.atmos.rho))
	lnT = log.(pyconvert(Array, result_NLTE.run.atmos.temp))

    χ_fraction = zeros(eltype(lnT_new), length(lnT), length(λ_new))
    S_fraction = zeros(eltype(lnT_new), length(lnT), length(λ_new))
    frac = similar(λ)

    for i in eachindex(lnρ)
        frac .= χ_NLTE[:, 1, 1, i] ./ χ_LTE[:, 1, 1, i]
        f_χ_NLTE = linear_interpolation(λ, frac, extrapolation_bc=NaN)
        χ_fraction[i, :] .= f_χ_NLTE.(λ_new)
        χ_fraction[i, isnan.(@view(χ_fraction[i, :]))] .= 1.0

        frac .= S_NLTE[:, 1, 1, i] ./ S_LTE[:, 1, 1, i]
        f_S_NLTE = linear_interpolation(λ, frac, extrapolation_bc=NaN)
        S_fraction[i, :] .= f_S_NLTE.(λ_new)
        S_fraction[i, isnan.(@view(S_fraction[i, :]))] .= 1.0
    end

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

# ============================================================================
# Effective temperature computations
# ============================================================================

window(line::PythonCall.Py, args...; kwargs...) = begin
    pyconvert.(Array, line.crop(args...; kwargs...))
end

window(run::TUMULTRun, iline, args...; kwargs...) = window(_line(run, iline), args...; kwargs...)

# ============================================================================
# Helper functions
# ============================================================================

flux(run; norm=true) = begin
    norm ? (run.lam, run.flux ./ run.flux_cnt) : (run.lam, run.flux)
end

"""
    flux(run, norm=true)
Return and optionally normalize the flux coming from a m3dis run.
"""
flux

B(λ, T) = begin
    @. 2*HPlanck*CLight^2/(aa_to_cm*λ)^5 * 1 /
                (exp(HPlanck*CLight /(aa_to_cm*λ*KBoltzmann*T)) -1 )
end

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

# ============================================================================
# Create a 3D cube from a run
# ============================================================================

_is_1D(x) = ndims(x) == 1
_set_size_1D_to_3D(x) = _is_1D(x) ? reshape(x, 1, 1, length(x)) : x

"""
    Box(result::TUMULTRun)

The atmosphere of the run on its own mesh (temperature, density, velocities, τ500).
"""
function Box(result::TUMULTRun)
    a = result.run.atmos
    x = Base.convert.(Float32, pyconvert(Array, a.xx) .* 1e8)
    y = Base.convert.(Float32, pyconvert(Array, a.yy) .* 1e8)
    z = Base.convert.(Float32, pyconvert(Array, a.zz) .* 1e8)

    data = Dict{Symbol, Any}(
        :T => _set_size_1D_to_3D(pyconvert(Array, a.temp)),
        :d => _set_size_1D_to_3D(pyconvert(Array, a.rho)),
        :ux => _set_size_1D_to_3D(pyconvert(Array, a.vx)),
        :uy => _set_size_1D_to_3D(pyconvert(Array, a.vy)),
        :uz => _set_size_1D_to_3D(pyconvert(Array, a.vz)),
        :τ500 => _set_size_1D_to_3D(exp10.(pyconvert(Array, a.ltau))),
    )
    xx, yy, zz = meshgrid(x, y, z)
    Box(xx, yy, zz, data, AtmosphericParameters())
end
