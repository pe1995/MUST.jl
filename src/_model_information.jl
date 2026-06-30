"""
	ModelInformation

Information from model headers in grid naming convention.
"""
Base.@kwdef mutable struct ModelInformation
	category = ""
	teff = ""
	logg = ""
	feh = ""
	alpha = ""
	vmic = ""
	extension = ""
	wrong_format = false
	original_name = ""
end

# ============================================================================= 
# Information of grid models from Name
# =============================================================================

_find_datastring(parts) = begin
	i_part = Ref(-1)
	for (i, x) in enumerate(parts)
		if occursin("t", x)&occursin("g", x)&occursin("m", x)
			p = try 
				_parameters_from_string(x)
			catch
				nothing
			end
			if isnothing(p)
				continue
			else
				i_part[] = i
				break
			end
		else
			continue
		end
	end

	i_part[]
end

_parameters_from_string(para) = begin
	t1 = findfirst(x->x=='t', para) + 1
	t2 = findfirst(x->x=='g', para) - 1
	t3 = findfirst(x->x=='m', para) - 1

	t = parse(Float64, para[t1:t2]) .* 100
	g = parse(Float64, para[t2+2:t3]) ./10
	m = parse(Float64, para[t3+2:end])

	d = Dict("teff"=>t, "logg"=>g, "feh"=>m)
end

# ============================================================================= 
# Comparison of model information
# =============================================================================

"""
	same_parameters(a, b)

Check if a and b have the same stellar parameters.
"""
same_parameters(a::ModelInformation, b::ModelInformation) = (a.teff≈b.teff) & (a.logg≈b.logg) & (a.feh≈b.feh)

"""
	same_id(a, b, field)

Check if a and b have the same `field` entry. `False` if entries are empty.
"""
same_id(a::ModelInformation, b::ModelInformation, field) = if (length(getfield(a, field)) > 0) & (length(getfield(b, field)) > 0)
	getfield(a, field) == getfield(b, field)
else
	false
end

"""
	same_id(a, b, field)

Check if a and b have the same `category` and `extension` entries.
"""
same_class(a::ModelInformation, b::ModelInformation) = same_id(a, b, :category) & same_id(a, b, :extension)

# ============================================================================= 
# Building new name string from information
# =============================================================================

"""
	name_string(minfo::ModelInformation; add_E=false)

Build a new name string from a model information.
"""
name_string(minfo::ModelInformation; add_E=false) = begin
	tname = MUST.@sprintf "%.2f" minfo.teff/100
	gname = MUST.@sprintf "%.2f" minfo.logg*10
	fname = MUST.@sprintf "%.3f" minfo.feh
	name = "t$(tname)g$(gname)m$(fname)"
	
	components = if add_E
		[minfo.category, "E", name, minfo.extension]
	else
		[minfo.category, name, minfo.extension]
	end

	components = [c for c in components if length(c)>0]
	join(components, '_')
end

pretty_name(minfo::ModelInformation, snapshot=nothing; format="new") = begin
	# Convert from internal name to pretty name
    a_str = (minfo.alpha == "" || minfo.alpha === nothing) ? "" : MUST.@sprintf("_a%+.2f", minfo.alpha)
    v_str = (minfo.vmic == "" || minfo.vmic === nothing) ? "" : MUST.@sprintf("_vmic%+.2f", minfo.vmic)
    t_val = (minfo.teff == "" || minfo.teff === nothing) ? 0 : round(Int, minfo.teff)
    g_val = (minfo.logg == "" || minfo.logg === nothing) ? 0.0 : minfo.logg
    z_val = (minfo.feh == "" || minfo.feh === nothing) ? 0.0 : minfo.feh

    if format == "new"
        if isnothing(snapshot)
            MUST.@sprintf "%s_t%05d_g%+.2f_z%+.2f%s%s" minfo.category t_val g_val z_val a_str v_str
        else
            MUST.@sprintf "%s_t%05d_g%+.2f_z%+.2f%s%s_sn%i" minfo.category t_val g_val z_val a_str v_str snapshot
        end
    else
        if isnothing(snapshot)
            MUST.@sprintf "%s_t%i_g%.3f_z%.3f" minfo.category minfo.teff minfo.logg minfo.feh
        else
            MUST.@sprintf "%s_t%i_g%.3f_z%.3f_sn%i" minfo.category minfo.teff minfo.logg minfo.feh snapshot
        end
    end
end

pretty_name(name::String, args...; format="new", kwargs...) = pretty_name(ModelInformation(name; format=(format=="new" ? "auto" : format)), args...; format=format, kwargs...)

"""
	ModelInformation(name; format="auto")

Extract category, teff, logg, feh and version from standard m3dis format names.
The format is `category_teff/1000_logg*10_±feh_extension` or `category_txxxxx_g+-x.xx_z+-x.xx_a+-x.xx_vmicx.xx`.
"""
ModelInformation(name; format="auto") = begin
	name_parts = split(name |> strip, "_")

	# if there is an E in the parameter name we remove it
	name_parts = [p for p in name_parts if !(p=="E")]

    if format == "auto" || format == "new"
        # Check for new format: e.g. M3DIS_t04000_g-0.50_z-2.00_a+0.40_vmic+2.00
        t_idx = findfirst(x -> startswith(x, "t") && length(x) > 1 && tryparse(Float64, x[2:end]) !== nothing, name_parts)
        g_idx = findfirst(x -> startswith(x, "g") && length(x) > 1 && tryparse(Float64, x[2:end]) !== nothing, name_parts)
        z_idx = findfirst(x -> startswith(x, "z") && length(x) > 1 && tryparse(Float64, x[2:end]) !== nothing, name_parts)
        
        if !isnothing(t_idx) && !isnothing(g_idx) && !isnothing(z_idx)
            teff = parse(Float64, name_parts[t_idx][2:end])
            logg = parse(Float64, name_parts[g_idx][2:end])
            feh = parse(Float64, name_parts[z_idx][2:end])
            
            prefix = join(name_parts[1:t_idx-1], '_')
            
            last_param_idx = maximum([t_idx, g_idx, z_idx])
            
            a_idx = findfirst(x -> startswith(x, "a") && length(x) > 1 && tryparse(Float64, x[2:end]) !== nothing, name_parts)
            alpha = ""
            if !isnothing(a_idx)
                alpha = parse(Float64, name_parts[a_idx][2:end])
                last_param_idx = max(last_param_idx, a_idx)
            end
            
            v_idx = findfirst(x -> startswith(x, "vmic") && length(x) > 4 && tryparse(Float64, x[5:end]) !== nothing, name_parts)
            vmic = ""
            if !isnothing(v_idx)
                vmic = parse(Float64, name_parts[v_idx][5:end])
                last_param_idx = max(last_param_idx, v_idx)
            end
            
            extension = join(name_parts[last_param_idx+1:end], '_')
            
            return ModelInformation(
                teff=teff,
                logg=logg,
                feh=feh,
                alpha=alpha,
                vmic=vmic,
                category=prefix,
                extension=extension,
                original_name=name
            )
        end
    end

    if format == "auto" || format == "legacy"
	    try
		    # find the data string
		    idatastring = _find_datastring(name_parts)

		    if idatastring == -1
			    # this has the wrong format
			    ModelInformation(wrong_format=true, original_name=name)
		    else
			    # everyting before the datastring is prefix, the rest is extension
			    prefix = try
				    join(name_parts[1:idatastring-1], '_')
			    catch
				    ""
			    end

			    extension = try
				    join(name_parts[idatastring+1:end], '_')
			    catch
				    ""
			    end
	
			    p = _parameters_from_string(name_parts[idatastring])
			    ModelInformation(
				    teff=p["teff"],
				    logg=p["logg"],
				    feh=p["feh"],
				    category=prefix,
				    extension=extension,
				    original_name=name
			    )
		    end
	    catch
		    ModelInformation(wrong_format=true, original_name=name)
	    end
    end
end
