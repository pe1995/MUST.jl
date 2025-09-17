"""
    MeanSpectrum(λ,spectrum,composition,kind,model_info)

Mean spectrum type for any kind of average spectrum.
"""
struct MeanSpectrum
    λ
    spectrum
    composition
    kind
    model_info
end

"""
    MeanSpectra(path; kind=nothing)

Read mean spectra from file. The format is assumed to follow (see `average_spectra.jl` script in dispatch):
# Average spectra for model_name
# Teff_min[K],Teff_max[K],logg,feh,alpha,n_snapshots,dt[h]
# para1,para2,...,paraN
# wavelength,abundance1,abundance2...,abundanceM
... data in csv format ...

Note that abundances are expected as [X/Fe]=xxx_[Y/Fe]=yyy, etc.
"""
MeanSpectra(path; kind=nothing) = begin
    data = CSV.read(path, DataFrame, header=4)
    header = open(path, "r") do f
        readlines(f)[1:4]
    end
    modelname = strip(split(header[1], ' ', keepempty=false)[end])
    params = split(strip(header[2][2:end]), ',', keepempty=false)
    params_vals = split(strip(header[3][2:end]), ',', keepempty=false)
    model_info = Dict{String}{Any}()
    for (p, v) in zip(params, params_vals)
        model_info[p] = try
            parse(Float64, v) 
        catch
            v
        end
    end
    model_info["modelname"] = modelname

    kind = if isnothing(kind)
        try 
            #s = String.(split(basename(path), '_', keepempty=false)[2:end])
            #s[end] = split(s[end], ".csv", keepempty=false)[1]

            istart = findfirst("mean_", basename(path))
            basename(path)[first(istart):end-4]
        catch
            basename(path)
        end
    else
        kind
    end

    tag = try
        String(split(basename(path), '_')[1])
    catch
        ""
    end
    model_info["tag"] = tag

    comps = [c for c in [n for n in names(data) if !(n=="# wavelength")]]
    sort_on_composition([MeanSpectrum(data[!,"# wavelength"][:], data[!, c][:], abundance_from_composition_string(c, split_on='_'), kind, model_info) for c in comps])
end


"""
    sort_on_composition(spectra, element=nothing)

Sort a list of spectra based on the composition of the given element ([X/Fe]).
Uses the elements in order of number of appearances in all spectra entries of the list, if no element is given.
"""
sort_on_composition(spectra, element=nothing) = begin
    c = [s.composition for s in spectra]
    all_elements = unique(vcat([keys(comp)|>collect for comp in c]...))
    how_many_times = Dict(
        e => count([e in keys(comp) for comp in c])
        for e in all_elements
    )
    # create a sorting tuple in the order of number of appearances
    n_appear = [how_many_times[e] for e in all_elements]
    appearance_sorted_elements = all_elements[sortperm(n_appear, rev=true)]

    mask = if isnothing(element)
        # sort based on this
        st = [Tuple(e in keys(comp) ? comp[e] : 99.0 for e in appearance_sorted_elements) for comp in c]
        sortperm(st)
    else
        # check if element is available in all, if not put those to the back
        st = [element in keys(comp) ? Tuple([comp[element], [comp[e] for e in appearance_sorted_elements if !(e==element)]...]) : Tuple([99.0, [comp[e] for e in appearance_sorted_elements if !(e==element)]...]) for comp in c]
        sortperm(st)
    end

    spectra[mask]
end

sortperm_on_composition_strings(spectra, element=nothing; kwargs...) = begin
    c = [abundance_from_composition_string(s; kwargs...) for s in spectra]
    all_elements = unique(vcat([keys(comp)|>collect for comp in c]...))
    how_many_times = Dict(
        e => count([e in keys(comp) for comp in c])
        for e in all_elements
    )
    # create a sorting tuple in the order of number of appearances
    n_appear = [how_many_times[e] for e in all_elements]
    appearance_sorted_elements = all_elements[sortperm(n_appear, rev=true)]

    mask = if isnothing(element)
        # sort based on this
        st = [Tuple(e in keys(comp) ? comp[e] : 99.0 for e in appearance_sorted_elements) for comp in c]
        sortperm(st)
    else
        # check if element is available in all, if not put those to the back
        st = [element in keys(comp) ? Tuple([comp[element], [comp[e] for e in appearance_sorted_elements if !(e==element)]...]) : Tuple([99.0, [comp[e] for e in appearance_sorted_elements if !(e==element)]...]) for comp in c]
        sortperm(st)
    end
end

"""
    interpolate(spectra, abundance)

Interpolate in `MeanSpectra` to the given abundance. Please specify e.g. ("[C/Fe]"=>val)
"""
function interpolate(spectra::Vector{S}, abundance) where {S<:MeanSpectrum}
    el = Symbol(first(abundance))
    target = last(abundance)
    spectra = [s for s in spectra if el in keys(s.composition)]

    if length(spectra) == 0
        error("No spectrum for the given composition found.")
    end

    ab = [s.composition[el] for s in spectra]
    sortmask = sortperm(ab)
    ab = ab[sortmask]
    spectra = spectra[sortmask]

    λ_min = max([minimum(s.λ) for s in spectra]...)
    λ_max = min([maximum(s.λ) for s in spectra]...)
    λ_common = range(λ_min, λ_max, length=length(spectra[1].λ)) |> collect

    vals = zeros(eltype(λ_common), length(spectra), length(λ_common))
    for (i, s) in enumerate(spectra)
        vals[i, :] = linear_interpolation(s.λ, s.spectrum, extrapolation_bc=Line()).(λ_common)
    end

    # interpolate in abundance wavelength by wavelength
    #ip = [linear_interpolation(ab, vals[:, i]).(target) for i in eachindex(λ_common)]
    ip = zeros(eltype(vals), length(λ_common))
    @inbounds for i in eachindex(λ_common)
        ip[i] = linear_interpolation(ab, vals[:, i], extrapolation_bc=Line())(target)
    end

    MeanSpectrum(λ_common, ip, Dict(el=>target), first(spectra).kind, first(spectra).model_info)
end


"""
    is_mean_spectrum(path, tag; kind="flux", norm=false) = begin

Check if the given path is a MeanSpectrum compatible file.
"""
is_mean_spectrum(path, tag; kind="flux", norm=false) = begin
	p = basename(path)
	t = String(tag)
	norm ? (p[1:length(t)] == t) && occursin(kind, p) && (occursin("norm", p)) : (p[1:length(t)] == t) && occursin(kind, p) && (!occursin("norm", p))
end
