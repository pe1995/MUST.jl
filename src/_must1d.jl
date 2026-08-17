# =============================================================================
# The MUST 1D model format
# =============================================================================

struct MUST1DColumn
    key         ::Union{Symbol,Nothing}
    name        ::String
    unit        ::String
    description ::String
    log         ::Bool
    logaverage  ::Bool
end

MUST1DColumn(key, name, unit, description; log=false, logaverage=log) = MUST1DColumn(key, name, unit, description, log, logaverage)

# =============================================================================
# Constant MUST1D columns and name mapping
# =============================================================================

const must1d_fixed_columns = MUST1DColumn[
    MUST1DColumn(:z,      "z",                       "cm",       "geometrical depth, increasing downwards"),
    MUST1DColumn(:τ_ross, "tau_ross",                "1",        "Rosseland optical depth",              logaverage=true),
    MUST1DColumn(:τ500,   "tau_500",                 "1",        "optical depth at 500 nm",              logaverage=true),
    MUST1DColumn(:T,      "T",                       "K",        "gas temperature"),
    MUST1DColumn(:d,      "log10_density",           "g cm-3",   "mass density",                         log=true),
    MUST1DColumn(:Pg,     "log10_gas_pressure",      "dyn cm-2", "gas pressure",                         log=true),
    MUST1DColumn(:Ne,     "log10_electron_density",  "cm-3",     "electron number density",              log=true),
    MUST1DColumn(:kr,     "log10_rosseland_opacity", "cm2 g-1",  "Rosseland opacity per unit mass",      log=true),
    MUST1DColumn(:κ500,   "log10_opacity_at_500nm",  "cm2 g-1",  "opacity per unit mass at 500 nm",      log=true),
    MUST1DColumn(:ee,     "internal_energy",         "erg g-1",  "internal energy per unit mass"),
    MUST1DColumn(nothing, "vmic",                    "km s-1",   "microturbulence, 1/3*sqrt(<u2> - (<ux>^2+<uy>^2+<uz>^2))")
]

const must1d_name_map = Dict{Symbol,MUST1DColumn}(
    :ux             => MUST1DColumn(:ux,            "velocity_x",              "cm s-1",       "velocity, x component"),
    :uy             => MUST1DColumn(:uy,            "velocity_y",              "cm s-1",       "velocity, y component"),
    :uz             => MUST1DColumn(:uz,            "velocity_z",              "cm s-1",       "velocity, z component (vertical)"),
    :e              => MUST1DColumn(:e,             "internal_energy_density", "erg cm-3",     "internal energy per unit volume"),
    :gamma          => MUST1DColumn(:gamma,         "gamma",                   "1",            "adiabatic index"),
    :P_turb         => MUST1DColumn(:P_turb,        "turbulent_pressure",      "dyn cm-2",     "turbulent pressure"),
    :bx             => MUST1DColumn(:bx,            "magnetic_field_x",        "G",            "magnetic field, x component"),
    :by             => MUST1DColumn(:by,            "magnetic_field_y",        "G",            "magnetic field, y component"),
    :bz             => MUST1DColumn(:bz,            "magnetic_field_z",        "G",            "magnetic field, z component (vertical)"),
    :flux           => MUST1DColumn(:flux,          "radiative_flux",          "erg cm-2 s-1", "vertical radiative energy flux"),
    :F_rad          => MUST1DColumn(:F_rad,         "radiative_flux",          "erg cm-2 s-1", "vertical radiative energy flux"),
    :F_conv         => MUST1DColumn(:F_conv,        "convective_flux",         "erg cm-2 s-1", "vertical convective energy flux"),
    :F_err_rel      => MUST1DColumn(:F_err_rel,     "flux_error_relative",     "1",            "relative deviation of the total flux from sigma*Teff^4"),
    :qr             => MUST1DColumn(:qr,            "radiative_heating",       "erg g-1 s-1",  "radiative heating rate"),
    :Q_rad          => MUST1DColumn(:Q_rad,         "radiative_heating",       "erg g-1 s-1",  "radiative heating rate"),
    :J              => MUST1DColumn(:J,             "mean_intensity",          "erg cm-2 s-1 sr-1", "angle averaged mean intensity"),
    :P_rad          => MUST1DColumn(:P_rad,         "radiation_pressure",      "dyn cm-2",     "radiation pressure"),
    :g_rad          => MUST1DColumn(:g_rad,         "radiative_acceleration",  "cm s-2",       "radiative acceleration"),
    :sourcefunction => MUST1DColumn(:sourcefunction,"source_function",         "erg cm-2 s-1 sr-1", "radiative source function"),
    :kapparho       => MUST1DColumn(:kapparho,      "opacity_times_density",   "cm-1",         "absorption coefficient"),
    :chi_min        => MUST1DColumn(:chi_min,       "chi_min",                 "cm-1",         "minimum absorption coefficient over all opacity bins"),
    :chi_max        => MUST1DColumn(:chi_max,       "chi_max",                 "cm-1",         "maximum absorption coefficient over all opacity bins"),
    :dT             => MUST1DColumn(:dT,            "temperature_correction",  "K",            "temperature correction of the last 1D iteration"),
    :dFconv_dT      => MUST1DColumn(:dFconv_dT,     "dFconv_dT",               "erg cm-2 s-1 K-1", "derivative of the convective flux w.r.t. temperature"),
    :dt_rt          => MUST1DColumn(:dt_rt,         "dt_rt",                   "s",            "radiative time step")
)

const must1d_ignored_keys = Symbol[:x, :y, :i_patch, :u2, :newscale, :extrapolated]
const must1d_suffixes = ("must1d_z", "must1d_tauross", "must1d_tau500")

# =============================================================================
# Composition and sanitization
# =============================================================================

"""
    must1d_composition(c)

`[X/Fe]` composition as a `Dict`, from a `Dict` or a `"Mn_-0.1,Co_0.1".
"""
must1d_composition(c::AbstractString) = Dict{String,Float64}(
    String(first(split(a, "_", keepempty=false))) => parse(Float64, last(split(a, "_", keepempty=false)))
    for a in split(strip(c), ",", keepempty=false)
)
must1d_composition(c::AbstractDict) = Dict{String,Float64}(string(k) => Float64(v) for (k, v) in c)
must1d_composition(::Nothing) = Dict{String,Float64}()

must1d_sanitized_name(key) = begin
    s = String(key)
    for (from, to) in ("τ"=>"tau", "κ"=>"kappa", "ρ"=>"rho", "σ"=>"sigma", "λ"=>"lambda",
                       "μ"=>"mu", "ν"=>"nu", "α"=>"alpha", "β"=>"beta", "γ"=>"gamma",
                       "δ"=>"delta", "ε"=>"epsilon", "χ"=>"chi", "Δ"=>"delta")
        s = replace(s, from=>to)
    end
    replace(replace(s, r"\s+"=>"_"), r"[^A-Za-z0-9_]"=>"_")
end

# =============================================================================
# Data computation
# =============================================================================

must1d_column(key::Symbol) = get(
    must1d_name_map, key,
    MUST1DColumn(key, must1d_sanitized_name(key), "?", "quantity '$(key)' of the MUST.Box object")
)

function must1d_columns(b::Box; only_default=true)
    columns = copy(must1d_fixed_columns)
    only_default && return columns

    used_keys  = Set{Symbol}(c.key for c in columns if !isnothing(c.key))
    used_names = Set{String}(c.name for c in columns)

    rest = MUST1DColumn[]
    for key in keys(b.data)
        (key in used_keys) && continue
        (key in must1d_ignored_keys) && continue
        (typeof(b.data[key]) <: AbstractArray) || continue
        (size(b.data[key]) == size(b.z)) || continue

        c = must1d_column(key)
        if c.name in used_names
            c = MUST1DColumn(key, must1d_sanitized_name(key), c.unit, c.description, c.log, c.logaverage)
        end
        push!(used_names, c.name)
        push!(rest, c)
    end

    sort!(rest, by=c->c.name)
    append!(columns, rest)
end

_must1d_available(b::Box, key::Symbol) = haskey(b.data, key) || (key in (:x, :y, :z))

function must1d_profile(b::Box, c::MUST1DColumn)
    key = c.key
    (isnothing(key) || !_must1d_available(b, key)) && return fill(NaN, size(b.z, 3))

    logaverage = (c.log || c.logaverage) && all(x -> x > 0, b[key])
    _, y = if logaverage
        profile(mean, b, :z, Symbol("log10", key), filter_out=:extrapolated)
    else
        profile(mean, b, :z, key, filter_out=:extrapolated)
    end

    y = Base.convert(Vector{Float64}, y)

    if c.log && !logaverage
        y = [yi > 0 ? log10(yi) : NaN for yi in y]
    elseif !c.log && logaverage
        y = exp10.(y)
    end

    y
end

function must1d_vmic(b::Box)
    all(haskey(b.data, k) for k in (:ux, :uy, :uz)) || return fill(NaN, size(b.z, 3))

    had_u2 = haskey(b.data, :u2)
    had_u2 || (b.data[:u2] = b[:ux] .^2 .+ b[:uy] .^2 .+ b[:uz] .^2)

    _, vx = profile(mean, b, :z, :ux, filter_out=:extrapolated)
    _, vy = profile(mean, b, :z, :uy, filter_out=:extrapolated)
    _, vz = profile(mean, b, :z, :uz, filter_out=:extrapolated)
    _, v2 = profile(mean, b, :z, :u2, filter_out=:extrapolated)

    had_u2 || delete!(b.data, :u2)

    1/3 .* sqrt.(max.(v2 .- (vx .^2 .+ vy .^2 .+ vz .^2), 0)) ./ 1e5
end

_must1d_field_width(c::MUST1DColumn) = max(length(c.name) + 2, 22)

# =============================================================================
# Header
# =============================================================================

_must1d_sorting(data) = begin
    d = data["log10_density"]
    ref = all(isnan, d) ? data["T"] : d
    (first(ref) < last(ref)) ? collect(1:length(ref)) : collect(length(ref):-1:1)
end

function _write_must1d_header(f, b::Box, columns, header, information, snapshot, feh, alpha, composition, abund)
    p = b.parameter
    version = try string(pkgversion(@__MODULE__)) catch; "unknown" end

    info = ["cube = $(size(b, 1)) x $(size(b, 2)) x $(size(b, 3))"]
    p.teff > 0    && push!(info, "Teff [K] = $(round(p.teff, digits=2))")
    p.logg > -90  && push!(info, "log(g) = $(round(p.logg, digits=3))")
    p.time > -90  && push!(info, "time [s] = $(p.time)")
    isnothing(snapshot) || push!(info, "snapshot = $(snapshot)")

    write(f, "# MUST 1D model: $(isnothing(header) ? "unnamed" : header)\n")
    write(f, "# $(join(info, ", "))\n")

    comp = String[]
    isnothing(feh)   || push!(comp, "[Fe/H] = $(feh)")
    isnothing(alpha) || push!(comp, "[alpha/Fe] = $(alpha)")
    cdict = must1d_composition(composition)
    for (k, v) in sort(collect(cdict), by=first)
        push!(comp, "[$(k)/Fe] = $(v)")
    end
    # an empty composition is not an unknown one: it means the default abundances
    # of the abund file were used. Only a composition that was never passed at
    # all (nothing) is unknown.
    (!isnothing(composition) && isempty(cdict)) && push!(comp, "default composition")
    (isnothing(abund) || isempty(strip(string(abund)))) || push!(comp, "abund file = $(abund)")
    write(f, "# $(isempty(comp) ? "composition unknown" : join(comp, ", "))\n")
    write(f, "# written $(Libc.strftime("%Y-%m-%d %H:%M:%S", time())) by MUST.jl $(version)\n")
    if !isnothing(information) && length(strip(information)) > 0
        for l in split(information, "\n")
            write(f, "# $(replace(l, r"^\s*\*\s*"=>""))\n")
        end
    end
    write(f, "#\n")
    for (i, c) in enumerate(columns)
        write(f, @sprintf("# %3d  %-24s %-20s %-15s %s\n", i, c.name, c.unit,
                          isnothing(c.key) ? "(derived)" : String(c.key), c.description))
    end
end

# =============================================================================
# Interface
# =============================================================================

"""
    save_must1d(b::Box, f_new; kwargs...)

Write the MUST 1D representation of the `Box` `b` to `f_new`. The averages are
taken at constant height scale of `b`, which follows from the cube you pass
(`box_sn*`, `box_tau_sn*`, `box_tau500_sn*`).

- `only_default`      : write only `must1d_fixed_columns`, `false` appends every
                        other quantity of the cube
- `header`            : name of the model
- `information`       : free form text added to the header
- `recompute_ne_func` : `(rho, T) -> ne`, recomputes the electron density from
                        the averaged rho and T instead of averaging it
- `vmic`              : microturbulence [km s-1] for models without a velocity field
- `feh`, `alpha`      : [Fe/H] and [alpha/Fe], for the header
- `composition`       : [X/Fe] of individual elements, `Dict` or `"Mn_-0.1,Co_0.1"`
- `snapshot`          : snapshot number, for the header
"""
function save_must1d(b::Box, f_new;
        only_default=true, header=nothing, information=nothing,
        recompute_ne_func=nothing, vmic=nothing, snapshot=nothing,
        feh=nothing, alpha=nothing, composition=nothing, abund=nothing)
    columns = must1d_columns(b, only_default=only_default)

    data = Dict{String,Vector{Float64}}()
    for c in columns
        data[c.name] = (c.name == "vmic") ?
            Base.convert(Vector{Float64}, must1d_vmic(b)) :
            must1d_profile(b, c)
    end

    if !isnothing(vmic)
        data["vmic"] = Base.convert(Vector{Float64}, (vmic isa Number) ? fill(vmic, size(b.z, 3)) : vmic)
    end
    if !isnothing(recompute_ne_func) && haskey(data, "log10_electron_density")
        ne = recompute_ne_func(exp10.(data["log10_density"]), data["T"])
        data["log10_electron_density"] = [n > 0 ? log10(n) : NaN for n in ne]
    end

    order = _must1d_sorting(data)
    for k in keys(data)
        data[k] = data[k][order]
    end
    if !all(isnan, data["z"]) && (first(data["z"]) > last(data["z"]))
        data["z"] = -1 .* data["z"]
    end

    open(f_new, "w") do f
        _write_must1d_header(f, b, columns, header, information, snapshot, feh, alpha, composition, abund)

        for c in columns
            Printf.format(f, Printf.Format("%$(_must1d_field_width(c))s"), c.name)
        end
        write(f, "\n")

        for i in 1:size(b.z, 3)
            for c in columns
                Printf.format(f, Printf.Format("%$(_must1d_field_width(c)).10E"), data[c.name][i])
            end
            write(f, "\n")
        end
    end

    f_new
end

"""
    save_must1d_all(b::Box, path; boxes=(b, b, b), kwargs...)

Write `\$(path)_must1d_z.txt`, `\$(path)_must1d_tauross.txt` and
`\$(path)_must1d_tau500.txt`. `boxes` are the geometrical, Rosseland and 500 nm
versions of the model, i.e. normally `box_sn*`, `box_tau_sn*`, `box_tau500_sn*`.
For a genuine 1D model the default of using `b` for all three is what you want.
"""
save_must1d_all(b::Box, path; boxes=(b, b, b), kwargs...) = [
    save_must1d(bi, "$(path)_$(suffix).txt"; kwargs...)
    for (bi, suffix) in zip(boxes, must1d_suffixes)
]

"""
    read_must1d(path)

Read a MUST 1D model file. Returns the columns as a `Dict{String,Vector{Float64}}`
and the commented header.
"""
function read_must1d(path)
    header = String[]
    names  = String[]
    rows   = Vector{Float64}[]

    open(path, "r") do f
        for l in eachline(f)
            sl = strip(l)
            isempty(sl) && continue
            if startswith(sl, "#")
                push!(header, l)
            elseif isempty(names)
                names = String.(split(sl))
            else
                push!(rows, parse.(Float64, split(sl)))
            end
        end
    end

    isempty(names) && error("$(path) does not contain a column name line.")

    Dict{String,Vector{Float64}}(n => [r[i] for r in rows] for (i, n) in enumerate(names)),
    join(header, "\n")
end
