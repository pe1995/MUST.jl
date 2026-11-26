# ============================================================================
# Spectra of a model
# ============================================================================

"""
    Spectrum

One spectrum of a model, e.g. of one OpenMULTI3D run: the wavelength `λ` [Å, air],
the ray directions `μ` (n_rays × 3) with their integration `weights`, and in `data`
the horizontally averaged `:intensity` and `:continuum_intensity` (λ × n_rays), the
`:flux` and `:continuum_flux`, for NLTE runs also `:intensity_lte` and `:flux_lte`,
and optionally the `:resolved_intensity` and `:resolved_continuum`
(λ × nx × ny × n_rays) on the columns `:resolved_x` and `:resolved_y`.
`composition` holds `[Fe/H]`, `[alpha/Fe]` and `[X/Fe]`, `meta` e.g. `:nlte`,
`:atom`, `:snapshot`, `:teff`, `:logg`.

In a file, a spectrum is the group `spectra/<name>`, so one model can hold the
spectra of any number of setups (`spectra(path)` lists them).
"""
struct Spectrum
    name::String
    λ::Vector{Float64}
    μ::Matrix{Float64}
    weights::Vector{Float64}
    data::Dict{Symbol, Any}
    composition::Dict{Symbol, Float64}
    meta::Dict{Symbol, Any}
end

Base.getindex(s::Spectrum, k) = s.data[Symbol(k)]
Base.haskey(s::Spectrum, k) = haskey(s.data, Symbol(k))
Base.keys(s::Spectrum) = keys(s.data)
Base.show(io::IO, s::Spectrum) = print(io,
    "Spectrum \"$(s.name)\": $(length(s.λ)) wavelengths $(round(minimum(s.λ), digits=2))-$(round(maximum(s.λ), digits=2)) Å, ",
    isempty(s.composition) ? "composition unknown" : composition_string(s.composition),
    get(s.meta, :nlte, false) ? ", NLTE" : ""
)

"""
    composition_string(composition)

`[Fe/H]=..,[alpha/Fe]=..,[X/Fe]=..` of a composition `Dict`.
"""
function composition_string(composition)
    first_keys = (Symbol("[Fe/H]"), Symbol("[alpha/Fe]"))
    ks = vcat([k for k in first_keys if haskey(composition, k)], sort([k for k in keys(composition) if !(k in first_keys)]))
    join(["$(k)=$(composition[k])" for k in ks], ",")
end

"""
    flux(s::Spectrum; norm=true, lte=false)

Wavelength and flux, normalized by the continuum flux with `norm`.
"""
flux(s::Spectrum; norm=true, lte=false) = (s.λ, norm ? s[lte ? :flux_lte : :flux] ./ s[:continuum_flux] : s[lte ? :flux_lte : :flux])

"""
    μ_z(s::Spectrum)

The distinct μ of the rays, from the disk center outwards.
"""
μ_z(s::Spectrum) = sort(unique(s.μ[:, 3]), rev=true)

"""
    intensity(s::Spectrum, μ; norm=true, lte=false)

Wavelength and intensity averaged over the rays with that μ, normalized by the
continuum intensity with `norm`.
"""
intensity(s::Spectrum, μ; norm=true, lte=false) = begin
    m = s.μ[:, 3] .≈ μ
    i = vec(mean(s[lte ? :intensity_lte : :intensity][:, m], dims=2))
    (s.λ, norm ? i ./ vec(mean(s[:continuum_intensity][:, m], dims=2)) : i)
end

# column positions [cm] of the run's mesh, from its first slab (all slabs share x and y)
function _columns(folder)
    l = readlines(joinpath(folder, "atmos_mesh.txt"))
    nx, ny = parse.(Int, split(l[1])[1:2])
    dx, dy = parse.(Float64, split(l[2])[1:2])
    collect(0:nx-1) .* dx .* 1e8, collect(0:ny-1) .* dy .* 1e8
end

_horizontal_mean(I) = reshape(mean(I, dims=(2, 3)), size(I, 1), size(I, 4))

"""
    Spectrum(run::TUMULTRun, name; composition=Dict(), meta=Dict(), resolved=false)

The spectrum of an OpenMULTI3D run. Runs with resolved intensities (3D) are averaged
horizontally, and keep them with `resolved`.
"""
function Spectrum(run::TUMULTRun, name; composition=Dict(), meta=Dict(), resolved=false)
    r = run.run
    w = pyconvert(Vector{Float64}, r.wts)
    data = Dict{Symbol, Any}()
    if pyhasattr(r, "i3")
        I, C = pyconvert(Array, r.i3), pyconvert(Array, r.c3)
        data[:intensity], data[:continuum_intensity] = _horizontal_mean(I), _horizontal_mean(C)
        pyhasattr(r, "i3_lte") && (data[:intensity_lte] = _horizontal_mean(pyconvert(Array, r.i3_lte)))
        if resolved
            data[:resolved_intensity], data[:resolved_continuum] = I, C
            data[:resolved_x], data[:resolved_y] = _columns(run_folder(run))
        end
    else
        data[:intensity], data[:continuum_intensity] = pyconvert(Matrix{Float64}, r.ie), pyconvert(Matrix{Float64}, r.ie_cnt)
        pyconvert(String, r.mode) == "NLTE" && (data[:intensity_lte] = pyconvert(Matrix{Float64}, r.ie_lte))
    end
    data[:flux] = data[:intensity] * w
    data[:continuum_flux] = data[:continuum_intensity] * w
    haskey(data, :intensity_lte) && (data[:flux_lte] = data[:intensity_lte] * w)

    m = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in meta)
    m[:nlte] = haskey(data, :intensity_lte)
    Spectrum(
        replace(String(name), "/" => "-"),
        pyconvert(Vector{Float64}, r.lam), pyconvert(Matrix{Float64}, r.vec), w,
        data, Dict{Symbol, Float64}(Symbol(k) => v for (k, v) in composition), m
    )
end

# ============================================================================
# Reading and writing
# ============================================================================

"""
    save(s::Spectrum, path)

Store `s` as the group `spectra/<name>` of the HDF5 file `path`, replacing a
spectrum of the same name.
"""
function save(s::Spectrum, path::String)
    HDF5.h5open(path, isfile(path) ? "r+" : "w") do fid
        g = haskey(fid, "spectra") ? fid["spectra"] : HDF5.create_group(fid, "spectra")
        haskey(g, s.name) && HDF5.delete_object(g, s.name)
        h = HDF5.create_group(g, s.name)
        h["wavelength"] = s.λ
        h["mu"] = s.μ
        h["weights"] = s.weights
        for (k, v) in s.data
            h[String(k)] = v
        end
        _write_composition!(h, s.composition)
        a = HDF5.attributes(h)
        for (k, v) in s.meta
            a[String(k)] = v
        end
    end
    path
end

"""
    spectra(path)

Names of the spectra in the HDF5 file `path`.
"""
spectra(path::String) = HDF5.h5open(path, "r") do fid
    names = String[]
    if haskey(fid, "spectra")
        g = fid["spectra"]
        append!(names, [k for k in keys(g) if g[k] isa HDF5.Group])
        append!(names, _legacy_tags(g))
    end
    unique(append!(names, _legacy_tags(fid)))
end

"""
    Spectrum(path, name; resolved=false)

Load the spectrum `name` from the HDF5 file `path`, with `resolved` also the
horizontally resolved intensities.
"""
Spectrum(path::String, name; resolved=false) = HDF5.h5open(path, "r") do fid
    n = String(name)
    g = haskey(fid, "spectra") ? fid["spectra"] : nothing
    if !isnothing(g) && haskey(g, n) && (g[n] isa HDF5.Group)
        _read_spectrum(g[n], n, resolved)
    elseif !isnothing(g) && haskey(g, "$(n)_wavelength")
        _read_legacy(g, n, resolved)
    elseif haskey(fid, "$(n)_wavelength")
        _read_legacy(fid, n, resolved)
    else
        error("No spectrum $(n) in $(path). Available: $(spectra(path))")
    end
end

"""
    load_spectra(path; resolved=false)

All spectra in the HDF5 file `path`.
"""
load_spectra(path::String; resolved=false) = [Spectrum(path, n; resolved=resolved) for n in spectra(path)]

function _read_spectrum(h, name, resolved)
    skip = ("wavelength", "mu", "weights", "composition")
    data = Dict{Symbol, Any}(
        Symbol(k) => read(h[k]) for k in keys(h)
        if !(k in skip) && (resolved || !startswith(k, "resolved_"))
    )
    a = HDF5.attributes(h)
    Spectrum(
        name, read(h["wavelength"]), read(h["mu"]), read(h["weights"]), data,
        _read_composition(h, Float64), Dict{Symbol, Any}(Symbol(k) => read(a[k]) for k in keys(a))
    )
end

# spectra written before the group layout: datasets `<tag>_<key>` with the
# composition as a string `[Fe/H]=..,[α/Fe]=..,[X/Fe]=..`
_legacy_tags(g) = [k[1:end-11] for k in keys(g) if endswith(k, "_wavelength") && haskey(g, k[1:end-11] * "_composition")]

function _read_legacy(g, tag, resolved)
    r(k) = haskey(g, "$(tag)_$(k)") ? read(g["$(tag)_$(k)"]) : nothing
    I, C = r("intensity"), r("continuum")
    data = Dict{Symbol, Any}()
    data[:intensity] = something(r("meanIntensity"), isnothing(I) ? nothing : _horizontal_mean(I))
    data[:continuum_intensity] = isnothing(C) ? data[:intensity] ./ r("meanIntensityNorm") : _horizontal_mean(C)
    data[:flux] = r("meanFlux")
    data[:continuum_flux] = r("meanFlux") ./ r("meanFluxNorm")
    if !isnothing(r("meanFluxLTE"))
        data[:flux_lte] = r("meanFluxLTE")
        data[:intensity_lte] = r("meanIntensityLTE")
    end
    if resolved && !isnothing(I)
        data[:resolved_intensity], data[:resolved_continuum] = I, C
    end
    composition = Dict{Symbol, Float64}(
        Symbol(replace(strip(first(kv)), "α" => "alpha")) => parse(Float64, last(kv))
        for kv in split.(split(r("composition"), ',', keepempty=false), '=') if length(kv) == 2
    )
    Spectrum(tag, r("wavelength"), r("mu"), r("weights"), data, composition, Dict{Symbol, Any}(:nlte => haskey(data, :flux_lte)))
end

"""
    snapshot_spectra(folder; name=nothing)

The spectra of the converted snapshots in `folder` as `Dict(snapshot => Vector{Spectrum})`,
only those called `name` if given.
"""
function snapshot_spectra(folder; name=nothing)
    found = Dict{Any, Vector{Spectrum}}()
    for file in glob("*.hdf5", normpath(folder))
        m = try first(_layout(file)) catch; continue end
        haskey(m, :snapshot) || continue
        names = try spectra(file) catch; String[] end
        isnothing(name) || (names = filter(==(String(name)), names))
        for n in names
            s = Spectrum(file, n)
            s.meta[:snapshot] = m[:snapshot]
            push!(get!(found, m[:snapshot], Spectrum[]), s)
        end
    end
    found
end

"""
    save_text(s::Spectrum, path; norm=false, lte=false, info=String[])

Write the flux and the intensity of each μ as columns of a text file, with the
lines of `info`, the composition and Teff/logg (if known) in the header.
"""
function save_text(s::Spectrum, path; norm=false, lte=false, info=String[])
    μs = μ_z(s)
    open(path, "w") do f
        println(f, norm ? "# Spectrum: Normalized Flux & Intensity" : "# Spectrum: Flux & Intensity")
        for l in info
            println(f, "# ", l)
        end
        haskey(s.meta, :teff) && println(f, "# Teff = $(round(Int, s.meta[:teff])) K")
        haskey(s.meta, :logg) && println(f, "# logg = $(round(s.meta[:logg], digits=2))")
        println(f, "# Composition: $(composition_string(s.composition))")
        println(f, get(s.meta, :nlte, false) && !lte ? "# NLTE atom: $(get(s.meta, :atom, ""))" : "# LTE")
        println(f, "# Column 1: Wavelength [Å]")
        println(f, "# Column 2: $(norm ? "Normalized Flux" : "Flux [erg/s/cm^2/Å]")")
        println(f, "# Columns 3-$(length(μs) + 2): $(norm ? "Normalized intensity" : "Intensity [erg/s/cm^2/Å/sr]")")
        println(f, "# lambda flux" * join([" mu=$(round(m, digits=4))" for m in μs]))
        F = last(flux(s; norm=norm, lte=lte))
        I = [last(intensity(s, m; norm=norm, lte=lte)) for m in μs]
        for i in eachindex(s.λ)
            println(f, @sprintf("%10.5f %10.5E", s.λ[i], F[i]) * join([@sprintf(" %10.5E", c[i]) for c in I]))
        end
    end
    path
end

# ============================================================================
# Combining spectra
# ============================================================================

"""
    average(spectra; name)

Mean of spectra of the same setup, e.g. of the same spectrum for several snapshots.
They need the same wavelengths, rays and composition. Resolved intensities are left out.
"""
function average(spectra::AbstractVector{Spectrum}; name=first(spectra).name)
    s1 = first(spectra)
    for s in spectra
        (length(s.λ) == length(s1.λ)) && (s.λ ≈ s1.λ) || error("$(s.name) has other wavelengths than $(s1.name).")
        (s.μ ≈ s1.μ) || error("$(s.name) has other rays than $(s1.name).")
        (s.composition == s1.composition) || error("$(s.name) has another composition than $(s1.name).")
    end
    ks = filter(k -> !startswith(String(k), "resolved_") && all(s -> haskey(s, k), spectra), collect(keys(s1)))
    meta = copy(s1.meta)
    meta[:snapshots] = join(string.(get.(getfield.(spectra, :meta), :snapshot, "")), ",")
    delete!(meta, :snapshot)
    Spectrum(name, s1.λ, s1.μ, s1.weights, Dict{Symbol, Any}(k => mean(s[k] for s in spectra) for k in ks), s1.composition, meta)
end

"""
    average_spectra(folder, name; snapshots=nothing)

Average the spectrum `name` over the snapshots in `folder` that have it (all, or
those in `snapshots`).
"""
function average_spectra(folder, name; snapshots=nothing)
    found = snapshot_spectra(folder; name=name)
    snaps = sort([s for s in keys(found) if isnothing(snapshots) || (s in snapshots)])
    isempty(snaps) && error("No snapshot in $(folder) has the spectrum $(name).")
    average(vcat([found[s] for s in snaps]...))
end

"""
    interpolate_spectra(spectra, element => value; lte=false)

Normalized flux at the [X/Fe] `value` of `element`, interpolated linearly between
spectra that differ in it (a spectrum without that element has [X/Fe] = 0).
"""
function interpolate_spectra(spectra::AbstractVector{Spectrum}, abundance; lte=false)
    key = Symbol("[$(first(abundance))/Fe]")
    ab = [get(s.composition, key, 0.0) for s in spectra]
    length(unique(ab)) > 1 || error("The spectra need at least two different $(key).")
    m = sortperm(ab)
    spectra, ab = spectra[m], ab[m]

    λ_common = range(maximum(minimum(s.λ) for s in spectra), minimum(maximum(s.λ) for s in spectra), length=length(first(spectra).λ)) |> collect
    vals = reduce(hcat, [linear_interpolation(s.λ, last(flux(s; lte=lte)), extrapolation_bc=Line()).(λ_common) for s in spectra])
    λ_common, [linear_interpolation(ab, vals[i, :], extrapolation_bc=Line())(last(abundance)) for i in eachindex(λ_common)]
end

"""
    equivalentwidth(s::Spectrum; λs=nothing, λe=nothing, lte=false)

Equivalent width [Å] of the normalized flux between `λs` and `λe`.
"""
equivalentwidth(s::Spectrum; λs=nothing, λe=nothing, lte=false) = begin
    λ, F = flux(s; norm=true, lte=lte)
    m = something(λs, minimum(λ)) .<= λ .<= something(λe, maximum(λ))
    integrate(λ[m], 1 .- F[m])
end

"""
    curve_of_growth(element, reference, models; λs=nothing, λe=nothing, use_log=true, lte=false)

Abundance correction of `element` from the curve of growth of the spectra `models`,
which differ in its [X/Fe], at the equivalent width of the `reference` spectrum.
Returns the reference [X/Fe] and equivalent width, the [X/Fe] and equivalent widths of
the models, and the correction.
"""
function curve_of_growth(element, reference::Spectrum, models::AbstractVector{Spectrum}; λs=nothing, λe=nothing, use_log=true, lte=false)
    key = Symbol("[$(element)/Fe]")
    λs = something(λs, max(minimum(reference.λ), (minimum(s.λ) for s in models)...))
    λe = something(λe, min(maximum(reference.λ), (maximum(s.λ) for s in models)...))

    ew_ref = equivalentwidth(reference; λs=λs, λe=λe, lte=lte)
    ab_ref = get(reference.composition, key, 0.0)

    abundances = [get(s.composition, key, 0.0) for s in models]
    equivalent_widths = [equivalentwidth(s; λs=λs, λe=λe, lte=lte) for s in models]

    fl = use_log ? log10 : identity
    m = sortperm(equivalent_widths)
    f_ab = linear_interpolation(fl.(equivalent_widths[m]), abundances[m], extrapolation_bc=Line())
    a = sortperm(abundances)

    ab_ref, ew_ref, abundances[a], equivalent_widths[a], f_ab(fl(ew_ref)) - ab_ref
end
