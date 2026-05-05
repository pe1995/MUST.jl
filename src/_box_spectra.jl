# =============================================================================
# Detecting spectra stored in cubes
# =============================================================================

"""
    spectra_key_from_tag(key, tag)

Name tag of the corresponding key in the data dict of Box.
"""
spectra_key_from_tag(key, tag) = Symbol(join([String(tag), String(key)], '_'))

"""
    spectra_tags(box::Box)

List all the spectra tags of this snapshot.
"""
spectra_tags(box::Box; check_for=["wavelength","composition"]) = begin
    spectra_tags = []
    keys_box = keys(box.data)
    for (para, val) in box.data
        if occursin(check_for[1], String(para))
            # check for first key, if found use the string in the back as tag
            key = join(split(String(para), "_", keepempty=false)[1:end-1], '_')

            # check all other check_fors for this key
            if all([spectra_key_from_tag(check, key) in keys_box for check in check_for])
                append!(spectra_tags, [Symbol(key)])
            end
        end
    end

    spectra_tags
end

"""
    spectra_tags(path; kwargs...)

Find all snapshots at a given path that contain a spectrum.
"""
spectra_tags(path::String; check_for=["wavelength","composition"], kwargs...) = begin
    spectra_tags(converted_snapshots(p; kwargs...); check_for=check_for)
end

"""
    spectra_tags(converted_snapshots)

Find all snapshots at a given path that contain a spectrum.
"""
function spectra_tags(cs::Dict; check_for=["wavelength","composition"])
    st = Dict()

    for snapname in list_snapshots(cs)
        try
            b, _ = pick_snapshot(cs, snapname)
            st[snapname] = spectra_tags(b; check_for=check_for)
        catch
            st[snapname] = []
        end
    end

    st
end

spectra_keys(b::Box, tag::Symbol; kwargs...) = begin 
    ks = [String(k) for k in keys(b.data) if is_idential_tag(tag, k)]
    [Symbol(split(k, "_", keepempty=false) |> last) for k in ks]
end

# =============================================================================
# Check tags
# =============================================================================

is_idential_tag(tag, to_check) = begin
    # the key needs to be idential
    # we thus need to check if the beginning of the key is the tag, and if the rest is really not belonging to the tag
    main_tag = String(tag)
    to_check_key = String(to_check)
    nl = length(String(tag))

    # fist check if it is even the main tag
    im = is_main_tag(tag, to_check)

    if im
        # it is the main tag, so we need to check if the rest of the key is just the extension
        # that means that it has to be followed by exactly one underscore and the rest of the key
        rest_of_key = to_check_key[length(main_tag)+1:end]
        stars_with_underscore = rest_of_key[1] == '_'
        only_one_underscore = count(==('_'), rest_of_key) == 1
        stars_with_underscore && only_one_underscore
    else
        # it is not even the main tag
        false
    end
end

"""
    is_main_tag(tag, to_check)
    
Check if `tag` is the main tag of `to_check`. Example: Different abundance 
runs may be called `CHGBand+0` and `CHGBand+0.2`, etc. In this case, 
`is_main_tag("CHGBand", "CHGBand+0")` would return true. The result is always true,
if `tag` is exactly what is writtin at the beginning of `to_check`! Be careful,
something like `CHGBandLTE` and `CHGBandNLTE` both have the same main tag `CHGBand`!
So looking for this main tag will group them together! Be specific about the main tags.
"""
is_main_tag(tag, to_check) = begin
    #to_check_maintag = split(String(to_check), '_')[1:end-1]
    to_check_maintag = String(to_check)
    # the main tag has to be identical.
    # there are possibly extensions added behind the main tag, but not before,
    # otherwise this counts as a differernt tag
    nl = length(String(tag))
    to_check_maintag[1:min(nl, length(to_check_maintag))] == String(tag)
end

# =============================================================================
# Add spectra from one box to another
# =============================================================================

add_spectra!(to::Box, from::Box) = begin
    tags = spectra_tags(from)

    for (key, val) in from.data
        for tag in tags
            if occursin(String(tag)*"_", String(key))
                # this needs to be added
                to[key] = deepcopy(val)
            end
        end
    end
end

# =============================================================================
# Interface
# =============================================================================

"""
    spectra(b::Box, tag::Symbol, key::Symbol) 

Read the spectrum entry with tag and key from the box.
"""
spectra(b::Box, tag::Symbol, key::Symbol) = b[spectra_key_from_tag(key, tag)]

"""
    spectra(b::Box, tag::Symbol) 

Read the spectrum (all keys) at the entry with tag from the box.
"""
spectra(b::Box, tag::Symbol) = Dict(
   k => b[spectra_key_from_tag(k, tag)] for k in spectra_keys(b, tag)
)

"""
    MeanSpectrum(b::Box, tag::Symbol, kind="meanFluxNorm")

Create a MeanSpectrum from the data that is stored in `b` under `tag`.
"""
MeanSpectrum(b::Box, tag::Symbol, kind="meanFluxNorm") = begin
    λ, F = try
        b[spectra_key_from_tag("wavelength", tag)], b[spectra_key_from_tag("meanFluxNorm", tag)]
    catch
        error("$(kind) not found in box. Compute it and add before calling this function.")
    end

    MeanSpectrum(
        λ, 
        F, 
        abundance_from_composition_string(b[spectra_key_from_tag("composition", tag)]), 
        kind, 
        Dict(p=>getfield(b.parameter, p) for p in fieldnames(AtmosphericParameters))
    )
end

# =============================================================================
# Spectra related quantities
# =============================================================================

get_μs(box, tag) = unique(box[MUST.spectra_key_from_tag("mu", tag)][:, 3])

"""
    angular_intensity(box::Box, tag, μ; norm=true)

Compute the mean of the intensity in direction μ.
"""
angular_intensity(box::Box, tag, μ; norm=true, intensity="intensity") = begin
    λ = box[spectra_key_from_tag("wavelength", tag)]
    I = box[spectra_key_from_tag("intensity", tag)]
	c = box[spectra_key_from_tag("continuum", tag)]
    mu = box[spectra_key_from_tag("mu", tag)]
	muz = mu[:, 3]
    mu_mask = muz .≈ μ

    z = norm ? I[:,:,:,mu_mask] ./ c[:,:,:,mu_mask] : I[:,:,:,mu_mask]
	(λ, reshape(mean(z, dims=4), size(I)[1:3]...))
end

"""
    mean_angular_intensity(box::Box, tag, μ; norm=true)

Compute the horizontal mean of the intensity in direction μ.
"""
mean_angular_intensity(box::Box, tag, μ; norm=true, intensity="intensity") = begin
    λ = box[spectra_key_from_tag("wavelength", tag)]
    I = box[spectra_key_from_tag(intensity, tag)]
	c = box[spectra_key_from_tag("continuum", tag)]
    mu = box[spectra_key_from_tag("mu", tag)]
	muz = mu[:, 3]
    mu_mask = muz .≈ μ

    z = norm ? I[:,:,:,mu_mask] ./ c[:,:,:,mu_mask] : I[:,:,:,mu_mask]
	(λ, reshape(mean(mean(z, dims=4), dims=(2, 3)), :))
end

"""
    mean_angular_intensity(box::Box, tag, μ; norm=true)

Compute the horizontal mean of the intensity in direction μ.
"""
mean_intensity(box::Box, tag; norm=true, intensity="intensity") = begin
    λ = box[spectra_key_from_tag("wavelength", tag)]
    I = box[spectra_key_from_tag(intensity, tag)]
	c = box[spectra_key_from_tag("continuum", tag)]

    z = norm ? I ./ c : I
	(λ, reshape(mean(z, dims=(2, 3)), length(λ), size(I, 4)))
end
 
"""
    mean_integrated_flux(box::Box, tag::Symbol; norm=true)

Integrate angular intensity using the weights stored in box. Then average horizontally.
"""
mean_integrated_flux(box::Box, tag; norm=true, intensity="intensity") = begin
    λ = box[spectra_key_from_tag("wavelength", tag)]
	I = box[spectra_key_from_tag(intensity, tag)]
	c = box[spectra_key_from_tag("continuum", tag)]
	w = box[spectra_key_from_tag("weights", tag)]

	F = zeros(eltype(I), size(I)[1])
	C = zeros(eltype(I), size(I)[1])
	
    @inbounds for i in axes(I, 1)
        F[i] = sum(reshape(mean(@view(I[i, :, :, :]), dims=(1, 2)), :) .* w)
        C[i] = sum(reshape(mean(@view(c[i, :, :, :]), dims=(1, 2)), :) .* w)
    end
		
    (λ, norm ? F ./ C : F)
end

"""
    integrate_flux(box::Box, tag::Symbol; norm=true)

Integrate angular intensity using the weights stored in box.
"""
function integrated_flux(box::Box, tag; norm=true, intensity="intensity")
    λ = box[spectra_key_from_tag("wavelength", tag)]
	I = box[spectra_key_from_tag(intensity, tag)]
	c = box[spectra_key_from_tag("continuum", tag)]
	w = box[spectra_key_from_tag("weights", tag)]

	F = zeros(eltype(I), size(I)[1:3]...)
	C = zeros(eltype(I), size(I)[1:3]...)
	@inbounds for k in axes(I, 3)
		@inbounds for j in axes(I, 2)
			@inbounds for i in axes(I, 1)
				F[i, j, k] = sum(I[i, j, k, :] .* w)
				C[i, j, k] = sum(c[i, j, k, :] .* w)
			end
		end
	end

    (λ, norm ? F ./ C : F)
end

# =============================================================================
# Compute average spectra and store them to file
# =============================================================================

"""
    average_spectra(available_run, nsnaps, nbatches; modelname1D="", datafolder="data")

Computes and saves average spectra.

- For 3D runs (modelname1D is empty), it processes the run in `nbatches` batches,
  each containing `nsnaps` snapshots, starting from the end and moving backward.
- For 1D runs, it processes the single specified model.
"""
function average_spectra(available_run, selectedSpecTagGrid, nsnaps, nbatches; datafolder="data", modelname1D="", selectedSnaps=[], elements=nothing, skip_from_end=0)
    #================================================================================#
    # Helper Function: Contains the core logic to process a given set of snapshots. #
    # This avoids code duplication and ensures 1D and 3D cases are treated identically. #
    #================================================================================#
    function _process_and_write_batch(selectedSnaps_raw, header_info, filename_prefix, convspec, specTags, o; elements=nothing)
        if isempty(selectedSnaps_raw)
            @warn "No snapshots to process for prefix $(filename_prefix). Skipping."
            return false # Indicate that processing was skipped
        end

        selectedSpecTagGridSym = Symbol(selectedSpecTagGrid)

        matchingTags = unique(vcat([
            [t for t in specTags[s] if is_main_tag(selectedSpecTagGridSym, t)] 
            for s in selectedSnaps_raw if haskey(specTags, s)
        ]...))

        matchingTagsPerSnap = Dict(
            s=>[t for t in matchingTags if (t in specTags[s])] for s in selectedSnaps_raw 
        )

        @info "Including spectra with tags $(String.(matchingTags)) in averaging." 
        selectedSnapsGrid = [s for s in selectedSnaps_raw if length(matchingTagsPerSnap[s]) > 0]
        
        if isempty(selectedSnapsGrid)
            @warn "There are no snapshots that contain the tag $(selectedSpecTagGrid) for this batch."
            return false
        else
            @info "$(header_info.batch_text) Averaging spectra from $(length(selectedSnapsGrid)) cubes with the tag $(selectedSpecTagGrid)."
        end
        
        specboxesGrid = Dict(s => pick_snapshot(convspec, s) |> first for s in selectedSnapsGrid)
        spectraGrid = Dict(s => Dict(t=>spectra(b, t) for t in matchingTagsPerSnap[s]) for (s, b) in specboxesGrid)

        is1D_local = header_info.is1D
        timesGrid = is1D_local ? [-99.0] : [first(pick_snapshot(o, s)).parameter.time for s in selectedSnapsGrid] ./ (60*60)
        teffsGrid = is1D_local ? [-99.0] : [first(pick_snapshot(o, s)).parameter.teff for s in selectedSnapsGrid]

        ncG(name, tag) = spectra_key_from_tag(name, tag)
        for (s, sp) in spectraGrid
            for t in matchingTagsPerSnap[s]
                specSnap = sp[t]
                if !haskey(specSnap, Symbol("meanFlux"))
                    if !is1D_local
                        @warn "There is no mean flux available for snapshot $s. Computing."
                    end
                    b = specboxesGrid[s]
                    b.data[ncG("meanFlux", t)] = mean_integrated_flux(b, t, norm=false) |> last
                    b.data[ncG("meanFluxNorm", t)] = mean_integrated_flux(b, t, norm=true) |> last
                    b.data[ncG("meanIntensity", t)] = mean_intensity(b, t, norm=false) |> last
                    b.data[ncG("meanIntensityNorm", t)] = mean_intensity(b, t, norm=true) |> last
                    spectraGrid[s][t] = spectra(b, t)
                end
            end
        end

        clean_comp_string(cs) = begin
            s = composition_string_from_abundance(; skip_eles=["[Fe/H]","[alpha/Fe]"], abundance_from_composition_string(cs)...)
            if s == ""
                "[X/Fe]=0.0"
            else
                s
            end
        end
        params_from_string(cs) = (abundance_from_composition_string(cs, "[Fe/H]"), abundance_from_composition_string(cs, "[alpha/Fe]"))
        contains_elements(cs) = isnothing(elements) ? true : all([occursin(e, cs) for e in elements])  

        abundanceGrid = sort(unique(vcat([[clean_comp_string(replace(sp[t][:composition], "α"=>"alpha")) for t in matchingTagsPerSnap[s]] for (s, sp) in spectraGrid]...)), rev=true)
        abundanceGrid = abundanceGrid[contains_elements.(abundanceGrid)]
        abundanceGrid = abundanceGrid[sortperm_on_composition_strings(abundanceGrid, split_on=',')]
        chemParamGrid = unique(vcat([[params_from_string(replace(sp[t][:composition], "α"=>"alpha")) for t in matchingTagsPerSnap[s]] for (s, sp) in spectraGrid]...)) |> first

        spectraGridAbundance = Dict(a=>[] for a in abundanceGrid)
        for (s, sp) in spectraGrid
            for t in matchingTagsPerSnap[s]
                c = clean_comp_string(replace(sp[t][:composition], "α"=>"alpha"))
                if c in abundanceGrid
                    append!(spectraGridAbundance[c], [sp[t]])
                end
            end
        end

        dt = length(timesGrid) > 1 ? diff(timesGrid)|> first : 0
        time1 = first(timesGrid)
        time2 = last(timesGrid)
        alpha_ref = chemParamGrid[2]
        feh_ref = chemParamGrid[1]
        
        header_line1 = "# Average spectra for $(available_run)\n"
        header_line11(k, n, a) = "# $(k),$(selectedSpecTagGrid),$(n),$(a)\n"
        header_line2 = "# Teff_min[K],Teff_max[K],logg,feh,alpha,n_snapshots,dt[h],firstSnap,lastSnap,firstSnapTime,lastSnapTime\n# "
        header_line3 = if is1D_local
            @sprintf "%.2f,%.2f,%.3f,%.3f,%.3f,%i,%.7E,%s,%s,%.7E,%.7E" minimum(teffsGrid) maximum(teffsGrid) -99.0 feh_ref alpha_ref length(selectedSnapsGrid) dt first(selectedSnapsGrid) last(selectedSnapsGrid) time1 time2
        else
            @sprintf "%.2f,%.2f,%.3f,%.3f,%.3f,%i,%.7E,%i,%i,%.7E,%.7E" minimum(teffsGrid) maximum(teffsGrid) ModelInformation(available_run).logg feh_ref alpha_ref length(selectedSnapsGrid) dt first(selectedSnapsGrid) last(selectedSnapsGrid) time1 time2
        end
        
        write_spectra(fname, specDict; header) = begin
            open(fname, "w") do f
                write(f, header)
                λ = specDict["wavelength"]
                for i in eachindex(λ)
                    line_str = @sprintf "%.7f," λ[i]
                    for ab in abundanceGrid
                        line_str *= @sprintf "%.7E," specDict[ab][i]
                    end
                    write(f, line_str[1:end-1] * "\n")
                end
            end
            @info "File $(basename(fname)) written."
        end
        
        getfilename(p) = joinpath(datafolder, available_run, "$(filename_prefix)_$(p).csv")
        
        nanmean(x) = begin
            nanx = [all(.!isnan.(x[i])) for i in eachindex(x)]
            count(nanx) > 0 ? mean(x[nanx]) : [NaN for _ in eachindex(first(x))]
        end

        wvl_dims = Dict(a=>median([length(s[:wavelength]) for s in sp]) for (a, sp) in spectraGridAbundance)
        wvl_dims = floor(Int, median([w for (_,w) in wvl_dims]))
        # filter those that have a different wavelength dimension
        for (a, sp) in spectraGridAbundance
            s_new = []
            for (i, s) in enumerate(sp)
                if length(s[:wavelength]) != wvl_dims
                    @warn "Entry for abundance $(a) has a different wavelength grid: $(length(s[:wavelength])) ≠ $(wvl_dims)"
                else
                    append!(s_new, [s])
                end
            end
            spectraGridAbundance[a] = s_new
        end

        fluxDictGrid = Dict(a=>nanmean([s[:meanFlux] for s in sp]) for (a, sp) in spectraGridAbundance)
        fluxDictGrid["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
        for (k, v) in  fluxDictGrid
            @assert length(v) == length(fluxDictGrid["wavelength"])
        end
        full_header = header_line1 * header_line11("flux", "absolute", "integrated") *header_line2 * header_line3 * "\n"
        ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
        write_spectra(getfilename("mean_flux"), fluxDictGrid, header=ab_header)

        fluxDictGridNorm = Dict(a=>nanmean([s[:meanFluxNorm] for s in sp]) for (a, sp) in spectraGridAbundance)
        fluxDictGridNorm["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
        full_header = header_line1 * header_line11("flux", "normalized", "integrated") *header_line2 * header_line3 * "\n"
        ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
        write_spectra(getfilename("mean_flux_norm"), fluxDictGridNorm, header=ab_header)

        # See if there is LTE as well
        LTEAvail = try
            Dict(a=>nanmean([s[:meanFluxLTE] for s in sp]) for (a, sp) in spectraGridAbundance)
            true
        catch
            false
        end
        if LTEAvail
            fluxDictGrid = Dict(a=>nanmean([s[:meanFluxLTE] for s in sp]) for (a, sp) in spectraGridAbundance)
            fluxDictGrid["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
            for (k, v) in  fluxDictGrid
                @assert length(v) == length(fluxDictGrid["wavelength"])
            end
            full_header = header_line1 * header_line11("flux_LTE", "absolute", "integrated") *header_line2 * header_line3 * "\n"
            ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
            write_spectra(getfilename("mean_flux_LTE"), fluxDictGrid, header=ab_header)


            fluxDictGridNorm = Dict(a=>nanmean([s[:meanFluxNormLTE] for s in sp]) for (a, sp) in spectraGridAbundance)
            fluxDictGridNorm["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
            full_header = header_line1 * header_line11("flux_LTE", "normalized", "integrated") *header_line2 * header_line3 * "\n"
            ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
            write_spectra(getfilename("mean_flux_LTE_norm"), fluxDictGridNorm, header=ab_header)
        end

        mus = first(spectraGridAbundance[first(abundanceGrid)])[:mu]
        muzs = sort(unique(mus[:, 3]), rev=true)
        for muz in muzs
            mu_mask = mus[:, 3] .≈ muz
            mustring = @sprintf "%.3f" muz
            
            fg = Dict(a=>nanmean([reshape(mean(s[:meanIntensity][:, mu_mask], dims=2), :) for s in sp]) for (a, sp) in spectraGridAbundance)
            fg["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
            full_header = header_line1 * header_line11("intensity", "absolute", "$(mustring)") *header_line2 * header_line3 * "\n"
            ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
            write_spectra(getfilename("mean_intensity_mu$(mustring)"), fg, header=ab_header)

            fg_norm = Dict(a=>nanmean([reshape(mean(s[:meanIntensityNorm][:, mu_mask], dims=2), :) for s in sp]) for (a, sp) in spectraGridAbundance)
            fg_norm["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
            full_header = header_line1 * header_line11("intensity", "normalized", "$(mustring)") *header_line2 * header_line3 * "\n"
            ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
            write_spectra(getfilename("mean_intensity_mu$(mustring)_norm"), fg_norm, header=ab_header)

            # See if there is LTE as well
            LTEAvail = try
                Dict(a=>nanmean([reshape(mean(s[:meanIntensityLTE][:, mu_mask], dims=2), :) for s in sp]) for (a, sp) in spectraGridAbundance)
                true
            catch
                false
            end
            if LTEAvail
                fg = Dict(a=>nanmean([reshape(mean(s[:meanIntensityLTE][:, mu_mask], dims=2), :) for s in sp]) for (a, sp) in spectraGridAbundance)
                fg["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
                full_header = header_line1 * header_line11("intensity_LTE", "absolute", "$(mustring)") *header_line2 * header_line3 * "\n"
                ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
                write_spectra(getfilename("mean_intensity_mu$(mustring)_LTE"), fg, header=ab_header)

                fg_norm = Dict(a=>nanmean([reshape(mean(s[:meanIntensityNormLTE][:, mu_mask], dims=2), :) for s in sp]) for (a, sp) in spectraGridAbundance)
                fg_norm["wavelength"] = first(spectraGridAbundance[first(abundanceGrid)])[:wavelength]
                full_header = header_line1 * header_line11("intensity_LTE", "normalized", "$(mustring)") *header_line2 * header_line3 * "\n"
                ab_header = full_header * "# wavelength," * join(replace.(abundanceGrid, ','=>'_'), ',') * "\n"
                write_spectra(getfilename("mean_intensity_mu$(mustring)_LTE_norm"), fg_norm, header=ab_header)
            end
        end
        return true # Indicate success
    end

    #=====================================================#
    # Main Function Body: Setup and execution control.    #
    #=====================================================#
    o = joinpath(datafolder, available_run)
    is1D = length(modelname1D) > 0

    convspec = converted_snapshots(o, box_name="box_m3dis")
    all_available_snaps = if length(selectedSnaps)==0 
        list_snapshots(convspec, numbered_only=!is1D)
    else
        selectedSnaps
    end
    specTags = spectra_tags(convspec)
    specTags = Dict(isnap => haskey(specTags, isnap) ? specTags[isnap] : [] for isnap in all_available_snaps)

    if is1D
        # --- 1D Case: Process a single model, no looping ---
        @info "Processing 1D model: $(modelname1D)"
        selectedSnaps_raw = [modelname1D]
        filename_prefix = "$(String(selectedSpecTagGrid))_$(modelname1D)"
        header_info = (is1D=true, batch_text="")
        
        _process_and_write_batch(selectedSnaps_raw, header_info, filename_prefix, convspec, specTags, o; elements=elements)
    else
        # --- 3D Case: Loop over batches ---
        snaps = sort([s for s in all_available_snaps if (s > 0 && haskey(specTags, s) && !isempty(specTags[s]))])
        
        end_idx = length(snaps) - skip_from_end
        batch_num = 1

        @info "Starting 3D batch processing for run $(available_run). Total batches to run: $(nbatches)."

        while batch_num <= nbatches && end_idx > 0
            @info "--- Computing Batch $(batch_num) of $(nbatches) ---"
            
            start_idx = max(1, end_idx - nsnaps + 1)
            batch_snaps = snaps[start_idx:end_idx]

            filename_prefix = "$(String(selectedSpecTagGrid))_$(snaps[start_idx])-$(snaps[end_idx])"
            header_info = (is1D=false, batch_text=" - Batch $(batch_num)")

            _process_and_write_batch(batch_snaps, header_info, filename_prefix, convspec, specTags, o, elements=elements)

            # Update loop variables for the next iteration
            end_idx -= nsnaps
            batch_num += 1
            
            if end_idx < 1 && batch_num <= nbatches
                @info "Reached the beginning of the snapshot list before completing all requested batches. Stopping."
            end
        end
        @info "Batch processing finished."
    end
end

# =============================================================================
# Extract composition
# =============================================================================

abundance_from_composition_string(cs, element) = begin
    if occursin(element, cs)
        cs_split = split(cs, ',', keepempty=false)
        comp = cs_split[findfirst(occursin.(element, cs_split))]
        parse(Float64, split(comp, '=', keepempty=false) |> last)
    else
        nothing
    end
end

abundance_from_composition_string(cs; split_on=',') = begin
	ss = split(cs, split_on, keepempty=false)
	eles = [first(split(s, '=')) for s in ss]
	abund = [last(split(s, '=')) for s in ss]
	Dict(Symbol(e)=>parse(Float64, a) for (e, a) in zip(eles, abund))
end

abundance_from_tag(model, tag, element) = begin
    comp_ref = model[
        spectra_key_from_tag(
            "composition", tag
        )
    ]
    abundance_from_composition_string(comp_ref, element)
end

composition_string_from_abundance(;skip_eles=[], ab...) = begin
	s = []
	for (ele, abund) in ab
		if !(String(ele) in skip_eles)
			eleabund = join([String(ele), "$(abund)"], '=')
			append!(s, [eleabund])
		end
	end
    join(s, ',')
end

# =============================================================================
# Equivalent width
# =============================================================================

equivalentwidth(model::Box, tag; λs=nothing, λe=nothing) = begin
    # get reference wavelength
    λ_ref, F_ref = try
        model[spectra_key_from_tag("wavelength", tag)],model[spectra_key_from_tag("meanFluxNorm", tag)]
    catch
        @warn "No mean flux available. Computing..."
        ll, ff = mean_integrated_flux(model, tag, norm=true)
        model[spectra_key_from_tag("meanFluxNorm", tag)] = ff
        ll, ff
    end

    # masking region
    λs = isnothing(λs) ? minimum(λ_ref) : λs
    λe = isnothing(λe) ? maximum(λ_ref) : λe
    λmask = λs .<= λ_ref .<= λe
    λ_ref = λ_ref[λmask]
    F_ref = F_ref[λmask]
    
    # get reference abundance and equivalent width
    integrate(λ_ref, 1 .- F_ref) 
end

equivalentwidth(spectrum::MeanSpectrum; λs=nothing, λe=nothing) = begin
    # get reference wavelength
    λ_ref, F_ref = spectrum.λ, spectrum.spectrum

    # masking region
    λs = isnothing(λs) ? minimum(λ_ref) : λs
    λe = isnothing(λe) ? maximum(λ_ref) : λe
    λmask = λs .<= λ_ref .<= λe
    λ_ref = λ_ref[λmask]
    F_ref = F_ref[λmask]
    
    # get reference abundance and equivalent width
    integrate(λ_ref, 1 .- F_ref) 
end

# =============================================================================
# Curve of growth
# =============================================================================

"""
    curve_of_growth(element, reference_model, reference_tag, model, model_tags; λs=nothing, λe=nothing)

Compute the CoG (+ abundance corrections) for a given element (e.g. `[C/Fe]`) from a reference model to the target model.
The spectra obtained via `reference_tag` is used to compute the difference in equivalent width to all spectra in `model_tags`.
Abundances are obtained from the `composition` entry of those spectra. The correction is `abund_model(ew_ref) - abund_ref`.
Returns: ab_ref, ew_ref, ab, ew, Δab
"""
function curve_of_growth(element, reference_model::Box, reference_tag, model::Box, model_tags; λs=nothing, λe=nothing, interpolate_reference_to=nothing, use_log=true)
	ew_ref, ab_ref = if isnothing(interpolate_reference_to)
        ew_ref = equivalentwidth(reference_model, reference_tag; λs=λs, λe=λe)
        ab_ref = abundance_from_tag(reference_model, reference_tag, element)

        ew_ref, ab_ref
    else
        # interpolate reference tags to the given abundance
        ew_ref = []
        ab_ref = []
        for i in eachindex(reference_tag)
            try
                ew_l = equivalentwidth(reference_model, reference_tag[i]; λs=λs, λe=λe)
                ab_l = abundance_from_tag(reference_model, reference_tag[i], element)
                append!(ew_ref, [ew_l])
                append!(ab_ref, [ab_l])
            catch
                @warn "Error in tag $(reference_tag[i]). Skipping."
            end
        end

        # interpolation
        f_ab = linear_interpolation(
            ab_ref[sortperm(ab_ref)], ew_ref[sortperm(ab_ref)], extrapolation_bc=Line()
        )

        target_ab = interpolate_reference_to

        f_ab(target_ab), target_ab
    end

	# Compute the same thing for all other models
	abundances = []
	equivalent_widths = []
	for tag in model_tags
        try
            ew = equivalentwidth(model, tag; λs=λs, λe=λe)
            ab = abundance_from_tag(model, tag, element)

            append!(equivalent_widths, [ew])
            append!(abundances, [ab])
        catch
            @warn "Error in tag $(tag). Skipping."
        end
	end

	# sorting
	abmask = sortperm(abundances)
	ewmask = sortperm(equivalent_widths)

    fl = use_log ? log10 : identity

	# interpolation
	f_ab = linear_interpolation(
		fl.(equivalent_widths[ewmask]), abundances[ewmask], extrapolation_bc=Line()
	)
	#f_ew = linear_interpolation(abundances[abmask], equivalent_widths[abmask], extrapolation_bc=Line())
	Δab = f_ab(fl(ew_ref)) - ab_ref

	ab_ref, ew_ref, abundances[abmask], equivalent_widths[abmask], Δab
end

"""
    curve_of_growth(element, reference_model::MeanSpectrum, models::Vector{MeanSpectrum}; λs=nothing, λe=nothing)


Compute the CoG (+ abundance corrections) for a given element (e.g. `[C/Fe]`) from a reference model to the target model.
The spectra is used to compute the difference in equivalent width to all spectra in `models`.
Abundances are obtained from the `composition` entry of those spectra. The correction is `abund_model(ew_ref) - abund_ref`.
Returns: ab_ref, ew_ref, ab, ew, Δab
"""
function curve_of_growth(element, reference_model::S, models::Vector{S}; λs=nothing, λe=nothing, use_log=true) where {S<:MeanSpectrum}
    element = Symbol(element)

    λ_min = max(minimum(reference_model.λ),[minimum(s.λ) for s in models]...)
    λ_max = min(maximum(reference_model.λ),[maximum(s.λ) for s in models]...)

    λs = isnothing(λs) ? λ_min : λs
    λe = isnothing(λe) ? λ_max : λe

    ew_ref = equivalentwidth(reference_model; λs=λs, λe=λe)
    ab_ref = reference_model.composition[element]

	# Compute the same thing for all other models
	abundances = []
	equivalent_widths = []
	for model in models
        try
            ew = equivalentwidth(model; λs=λs, λe=λe)
            ab =  model.composition[element]

            append!(equivalent_widths, [ew])
            append!(abundances, [ab])
        catch
            @warn "Error in composition $(model.composition[element]). Skipping."
        end
	end

	# sorting
	abmask = sortperm(abundances)
	ewmask = sortperm(equivalent_widths)

    fl = use_log ? log10 : identity
	# interpolation
	f_ab = linear_interpolation(
		fl.(equivalent_widths[ewmask]), abundances[ewmask], extrapolation_bc=Line()
	)
	#f_ew = linear_interpolation(abundances[abmask], equivalent_widths[abmask], extrapolation_bc=Line())
	Δab = f_ab(fl(ew_ref)) - ab_ref

	ab_ref, ew_ref, abundances[abmask], equivalent_widths[abmask], Δab
end