# ============================================================================= 
# Namelist Types
# =============================================================================

abstract type AbstractNamelist end

mutable struct StellarNamelist <: AbstractNamelist
    snapshot_params      ::Dict{String,Any}      
    cartesian_params     ::Dict{String,Any}             
    patch_params         ::Dict{String,Any}  
    experiment_params    ::Dict{String,Any}  
    scaling_params       ::Dict{String,Any}  
    stellar_params       ::Dict{String,Any}               
    newton_params        ::Dict{String,Any}  
    friction_params      ::Dict{String,Any}  
    boundary_params      ::Dict{String,Any}                   
    gravity_params       ::Dict{String,Any}  
    pressure_node_params ::Dict{String,Any}  
    rt_integral_params   ::Dict{String,Any}  
    sc_rt_params         ::Dict{String,Any}           
    restart_params       ::Dict{String,Any}  
    io_params            ::Dict{String,Any}                  
    aux_params           ::Dict{String,Any}  
    eos_params           ::Dict{String,Any}  
    bifrost_params       ::Dict{String,Any}                
    ramses_params        ::Dict{String,Any}  
    stagger_params       ::Dict{String,Any}           
    divb_clean_params    ::Dict{String,Any}  
    halo_params          ::Dict{String,Any}  
    timer_params         ::Dict{String,Any}  
    lock_params          ::Dict{String,Any} 
    task_list_params     ::Dict{String,Any} 
    an_params            ::Dict{String,Any} 
    dispatcher0_params   ::Dict{String,Any} 
end

mutable struct SnapshotNamelist <: AbstractNamelist
    io_nml       ::Dict{String,Any}  
    snapshot_nml ::Dict{String,Any}  
    idx_nml      ::Dict{String,Any}  
end

"""
    M3DNamelist()

Default namelist for M3D computations.
"""
mutable struct M3DNamelist <: MUST.AbstractNamelist
    io_params          ::Dict{String,Any}  
    timer_params       ::Dict{String,Any}  
    atmos_params       ::Dict{String,Any}  
    stagger_params     ::Dict{String,Any}  
    atom_params        ::Dict{String,Any}  
    m3d_params         ::Dict{String,Any}  
    linelist_params    ::Dict{String,Any}  
    spectrum_params    ::Dict{String,Any}  
    dispatcher0_params ::Dict{String,Any}  
    task_list_params   ::Dict{String,Any}   
    composition_params ::Dict{String,Any}   
    line_mask          ::Dict{String,Any}   
end

"""
    FreeNamelist()

Namelist with no limitations on the fields it contains.
"""
mutable struct FreeNamelist <: MUST.AbstractNamelist
    data ::Dict{Symbol,Any}   
end

# ============================================================================= 
# Interface
# =============================================================================

nmlFields(nml::AbstractNamelist) = fieldnames(typeof(nml))
nmlFields(nml::FreeNamelist) = keys(nml.data)
nmlField(nml::AbstractNamelist, key) = getfield(nml, key)
nmlField(nml::FreeNamelist, key) = getindex(nml.data, key)

nmlSetField!(nml::AbstractNamelist, key, value) = setfield!(nml, key, value)
nmlSetField!(nml::FreeNamelist, key, value) = nml.data[key] = value

nmlParameter(nml::AbstractNamelist, para) = begin
    matches = Dict()

    for f in nmlFields(nml)
        d = nmlField(nml, f)
        if para in keys(d)
            matches[f] = d[para]
        end
    end

    if length(keys(matches)) == 0
        @warn "No paramter $(para) found."
    elseif length(matches) > 1
        #@warn "More than one paramter $(para) found. Consider using `nmlField` with the correct field."
    end

    matches
end

nmlValue(nml::AbstractNamelist, para) = begin
    v = values(nmlParameter(nml, para))
    if length(v) > 0
        first(v)
    else
        nothing
    end
end

# ============================================================================= 
# Constructors
# =============================================================================

M3DNamelist() = M3DNamelist([Dict{String,Any}() for f in fieldnames(M3DNamelist)]...)
M3DNamelist(path::String) = begin
    s = M3DNamelist()
    read!(s, path)
    s
end

StellarNamelist() = StellarNamelist([Dict{String,Any}() for f in fieldnames(StellarNamelist)]...)
StellarNamelist(path::String) = begin 
    s = StellarNamelist()
    read!(s, path)
    s
end

SnapshotNamelist() = SnapshotNamelist([Dict{String,Any}() for f in fieldnames(SnapshotNamelist)]...)
SnapshotNamelist(path::String) = begin 
    s = SnapshotNamelist()
    read!(s, path)
    s
end

FreeNamelist() = FreeNamelist(Dict{Symbol, Any}())
FreeNamelist(path::String; buzzword=nothing) = begin 
    s = FreeNamelist()
    read!(s, path, buzzword=buzzword)
    s
end
FreeNamelist(path::String, stack) = begin 
    i = 1
    nmls = []
    while true
        bw = (stack, i)
        s = FreeNamelist(path, buzzword=bw)

        if length(nmlFields(s)) == 0
            break
        end

        append!(nmls, [s])
        i+=1
    end

    nmls
end

# ============================================================================= 
# Utility Functions for Namelists
# =============================================================================

function read!(nml::AbstractNamelist, path::String; buzzword=nothing)
    content = Dict{Symbol, Dict{String, Any}}()
    raw_text = read(path, String)
    
    text_no_comments = Char[]
    masked_text = Char[]
    
    in_sq, in_dq = false, false
    in_comment = false
    
    for c in raw_text
        if in_comment
            if c == '\n'
                in_comment = false
                push!(text_no_comments, c)
                push!(masked_text, c)
            end
            continue
        end
        
        # Toggle string states
        if c == '\'' && !in_dq; in_sq = !in_sq; end
        if c == '"'  && !in_sq; in_dq = !in_dq; end
        
        # Skip comments (unless the '!' is inside a string)
        if c == '!' && !in_sq && !in_dq
            in_comment = true
            continue
        end
        
        push!(text_no_comments, c)
        
        if (in_sq || in_dq) && (c == '/' || c == '=' || c == '&')
            push!(masked_text, ' ')
        else
            push!(masked_text, c)
        end
    end
    
    clean_text = String(text_no_comments)
    safe_text = String(masked_text)
    
    buzzcount = 0
    block_matches = collect(eachmatch(r"&([a-zA-Z_][a-zA-Z0-9_]*)", safe_text))
    
    for b_match in block_matches
        para = Symbol(lowercase(b_match.captures[1]))
        
        if !isnothing(buzzword)
            if para == buzzword[1]
                buzzcount += 1
            end
            if para != buzzword[1] || buzzcount != buzzword[2]
                continue
            end
        end
        
        start_search_idx = b_match.offset + length(b_match.match)
        end_idx = findnext('/', safe_text, start_search_idx)
        end_idx = isnothing(end_idx) ? length(safe_text) : end_idx
        
        block_safe = safe_text[start_search_idx:end_idx-1]
        block_clean = clean_text[start_search_idx:end_idx-1]
        
        content[para] = Dict{String, Any}()
        
        key_regex = r"([a-zA-Z_][a-zA-Z0-9_]*)\s*="
        k_matches = collect(eachmatch(key_regex, block_safe))
        
        for (k_idx, k_match) in enumerate(k_matches)
            key = lowercase(k_match.captures[1])
            
            val_start = k_match.offset + length(k_match.match)
            val_end = (k_idx < length(k_matches)) ? k_matches[k_idx+1].offset - 1 : length(block_clean)
            
            raw_val = block_clean[val_start:val_end]
            clean_val = strip(raw_val, [' ', '\t', '\n', '\r', ','])
            
            content[para][key] = parse_from_namelist(clean_val)
        end
        
        if !isnothing(buzzword)
            break
        end
    end

    for key in keys(content)
        if !(key in nmlFields(nml))
            if _nmlWarn(nml)
                continue
            end
        end
        nmlSetField!(nml, key, content[key])
    end
end

function parse_from_namelist(value_string)
    if length(strip(value_string)) == 0
        return ""
    end
    
    v_str = ("$(strip(value_string)[end])" == "/") ? strip(strip(value_string)[1:end-1]) : strip(value_string)
    
    val = if occursin(",", v_str) || occursin("*", v_str)
        v = parse_multiple(v_str) 
        length(v) == 1 ? first(v) : v
    elseif occursin(".", v_str) || occursin("E", uppercase(v_str)) || occursin("D", uppercase(v_str))
        try
            safe_float_str = replace(uppercase(v_str), "D" => "E")
            parse(Float64, safe_float_str)
        catch
            String(v_str)
        end
    else
        try
            parse(Int64, v_str)
        catch
            String(v_str)
        end
    end

    val = if typeof(val) <: AbstractString
        v_lower = lowercase(strip(val))
        if v_lower in ["t", ".true."]
            true
        elseif v_lower in ["f", ".false."]
            false
        else
            # Strip delimiters, remove line breaks, and strip the resulting Fortran space padding
            cleaned = strip(val, ['\'', '"'])
            cleaned = replace(cleaned, r"[\r\n]+" => "")
            strip(cleaned)
        end
    elseif (typeof(val) <: AbstractArray) && !isempty(val) && (typeof(first(val)) <: AbstractString)
        v = []
        for iv in val
            iv_lower = lowercase(strip(iv))
            vi = if iv_lower in ["t", ".true."]
                true
            elseif iv_lower in ["f", ".false."]
                false
            else
                # Do the same for string arrays
                cleaned = strip(iv, ['\'', '"'])
                cleaned = replace(cleaned, r"[\r\n]+" => "")
                strip(cleaned)
            end
            push!(v, vi)
        end
        v
    else
        val
    end

    return val
end

function parse_multiple(s)
    sComponents = []
    
    sComma = String[]
    in_sq, in_dq = false, false
    current_chunk = Char[]
    
    for c in s
        if c == '\'' && !in_dq; in_sq = !in_sq; end
        if c == '"'  && !in_sq; in_dq = !in_dq; end
        
        if c == ',' && !in_sq && !in_dq
            push!(sComma, strip(String(current_chunk)))
            empty!(current_chunk)
        else
            push!(current_chunk, c)
        end
    end
    push!(sComma, strip(String(current_chunk)))
    filter!(!isempty, sComma)

    sMulti = []
    for si in sComma
        in_sq, in_dq = false, false
        found_star = false
        left_chunk = Char[]
        right_chunk = Char[]
        
        for c in si
            if c == '\'' && !in_dq; in_sq = !in_sq; end
            if c == '"'  && !in_sq; in_dq = !in_dq; end
            
            if c == '*' && !in_sq && !in_dq && !found_star
                found_star = true
                continue
            end
            
            if found_star
                push!(right_chunk, c)
            else
                push!(left_chunk, c)
            end
        end
        
        if found_star
            push!(sMulti, [strip(String(left_chunk)), strip(String(right_chunk))])
        else
            push!(sMulti, [strip(String(left_chunk))])
        end
    end

    for component in sMulti
        c = if length(component) > 1
            repeat([last(component)], parse(Int, first(component)))
        else
            component
        end
        append!(sComponents, c)
    end

    sComponentsParsed = []
    for component in sComponents
        c_new = if occursin(".", component) || occursin("E", uppercase(component)) || occursin("D", uppercase(component))
            try 
                safe_float_str = replace(uppercase(component), "D" => "E")
                parse(Float64, safe_float_str)
            catch
                component
            end
        else
            try 
                parse(Int, component)
            catch
                component
            end
        end
        push!(sComponentsParsed, c_new)
    end
        
    return sComponentsParsed
end

function reverse_parse(value)
    val_str = ""
    if typeof(value) <: AbstractArray
        for val in value
            val_str = val_str * "$(reverse_parse(val)),"
        end
        val_str = val_str[1:end-1]
    elseif typeof(value) <: String
        if (value in ["t", "f", ".true.", ".false."]) || (occursin('*', value))
            v = strip(value)
            v = (v[1] == '\'') ? v[2:end] : v
            v = (v[end] == '\'') ? v[1:end-1] : v
            val_str = "$(strip(v))"
        elseif (occursin("'", value)) || (occursin('"', value))
            v = strip(value)
            v = (v[1] == '\'') ? v[2:end] : v
            v = (v[end] == '\'') ? v[1:end-1] : v
            val_str = "$(strip(v))"
        else
            can_be_float = try 
                parse(Float64, value)
                true
            catch
                false
            end

            if can_be_float
                val_str = "$(value)"
            else
                val_str = "'$(value)'"
            end
        end
    elseif typeof(value) <: Bool
        val_str = ".$(value)."
    else
        val_str = "$(value)"
    end
    return val_str
end

# ============================================================================= 
# Modified warning
# =============================================================================

_nmlWarn(nml::AbstractNamelist) = begin
    @warn "$(key) is not a valid fieldname"
    true
end
_nmlWarn(nml::FreeNamelist) = false

"""
Write a namelist to path.
"""
function write(nml::AbstractNamelist, path::String; soft_line_limit=85)
    para_names   = String.(nmlFields(nml))
    longest_name = maximum(length.(para_names))

    f = open(path, "w")

    for (i,p) in enumerate(para_names)
        spaces = longest_name - length(p) +1
        line   = "&$(p)" * repeat(" ", spaces) 
        d = nmlField(nml, Symbol(p))
        
        lline = longest_name+2
        for word in keys(d)
            add_to_line = "$(word)=" * reverse_parse(d[word]) * " "
            
            lline += length(add_to_line)
            if (lline >= soft_line_limit) 
                line  = line*"\n" * repeat(" ", longest_name+2) 
                lline = longest_name+2
            end
            
            line = line * add_to_line
        end 

        line = line * "/\n"

        Base.write(f, line)
    end

    close(f)
end

function get_restart_snap_nml(nml::StellarNamelist)
    # the restart snapshot
    i_snap    = nml.restart_params["snapshot"]
    folder    = @in_dispatch "data/$(strip(nml.restart_params["run"], [''', ' ']))"
    s = SnapshotNamelist(joinpath(_snapshot_folder(i_snap, glob("*/", folder)), "snapshot.nml"))
    s
end

"""
Set fields of the Namelist. Enter field=(parameter=>value,...)
"""
function set!(nml::AbstractNamelist; kwargs...)
    for (field,paras) in kwargs
        if !(typeof(paras) <: Tuple) 
            @warn "Argument to $(field) is not a Tuple. Add trailing comma if there is only one entry."
            continue
        end

        for (p,v) in paras
            nmlField(nml, field)[String(p)] = v
        end
    end
end

read_eos_params(path::String) = begin 
    open(path, "r") do f
        lines = readlines(f)
        lines = [strip.(split(strip(l), "=")) for l in lines]
        lines = Dict(kv[1]=>kv[2] for kv in lines if length(kv)==2)
    end
end

# ============================================================================= 
# Default namelists
# =============================================================================

_whole_spectrum_namelist!(nml::M3DNamelist; 
	io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
	timer_params=(:sec_per_report=>120,),
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_rho=>true),
    atom_params=(#:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true),
	m3d_params=(:verbose=>1, :fcheck=>1, :pcheck=>[1,1,1], :linecheck=>1, 
				:lvlcheck=>1,
                :n_nu=>32, :maxiter=>0, :decouple_continuum=>true),
	#linelist_params=(:dlam=>1.0,),
    composition_params=(:abund_file=>"./input_multi3d/abund_magg",),#:absdat_file=>"./input_multi3d/TS_absdat.dat"),
	#spectrum_params=(:daa=>2., :aa_blue=>2000, :aa_red=>40000)
    ) = begin

	set!(
		nml; 
		io_params=io_params, 
		timer_params=timer_params,
		atmos_params=atmos_params,
		m3d_params=m3d_params,
		#linelist_params=linelist_params,
		#spectrum_params=spectrum_params,
        atom_params=atom_params,
        composition_params=composition_params
	)
end

"""
    whole_spectrum_namelist(model_name; kwargs...)

Create a default namelist for the given model, that
can be used to compute the whole spectrum emerging this model.
This is handy when you want to compute the effective temperature.
"""
function whole_spectrum_namelist(model_name::String; 
					model_folder="./input_multi3d/MUST",
					linelist="./input_multi3d/vald_2490-25540.list",
					absmet=nothing,
					kwargs...)	
	# Create an empty Namelist
	nml = M3DNamelist()

	# Fill in the defaults for this, they can be modified
	_whole_spectrum_namelist!(nml)

	# path of the model
	model_path = joinpath(model_folder, model_name)

	# additionally apply the model specific fields
	set!(nml; atmos_params=(:atmos_file=>model_path,))

    if (isnothing(absmet)) & (!isnothing(linelist))
        set!(nml; linelist_params=(:linelist_file=>linelist,))
    elseif (!isnothing(absmet)) & (isnothing(linelist))
        set!(nml, composition_params=(:absmet_file=>absmet,))
        set!(nml, m3d_params=(:absmet_in_output=>true, ))
    elseif (isnothing(absmet)) & (isnothing(linelist))
        nothing
    else
        error("Need to specify either linelist or absmet (not both).")
    end

    # apply custom input parameters
    set!(nml; kwargs...)

	nml
end

_spectrum_namelist_lte!(nml::M3DNamelist; 
	io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
	timer_params=(:sec_per_report=>120,),
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_rho=>true),
    atom_params=(#:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true),
	m3d_params=(:verbose=>1, :fcheck=>1, :linecheck=>1, 
				:lvlcheck=>1, :short_scheme=>"set_a2", :long_scheme=>"lobatto",
                :n_nu=>32, :maxiter=>0, :ilambd=>0)) = begin
    # :decouple_continuum=>true
	set!(
		nml; 
		io_params=io_params, 
		timer_params=timer_params,
		atmos_params=atmos_params,
        m3d_params=m3d_params,
        atom_params=atom_params
	)
end

_spectrum_namelist_nlte!(nml::M3DNamelist; 
	io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
	timer_params=(:sec_per_report=>120,),
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_rho=>true),
    atom_params=(#:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true,
                :convlim=>1e-2),
	m3d_params=(:verbose=>1, :fcheck=>1, :linecheck=>1, 
				:lvlcheck=>1, :short_scheme=>"set_a2", :long_scheme=>"lobatto",
                :n_nu=>32, :maxiter=>100)) = begin
    # :decouple_continuum=>true
	set!(
		nml; 
		io_params=io_params, 
		timer_params=timer_params,
		atmos_params=atmos_params,
        m3d_params=m3d_params,
        atom_params=atom_params
	)
end

"""
    spectrum_namelist(model_name; kwargs...)

Create a default namelist for the given model, that
can be used to compute the spectrum emerging this model.
You can specify NLTE=true in order to load the default NLTE namelist.
"""
function spectrum_namelist(model_name::String; NLTE=false,
					model_folder="./input_multi3d/MUST",
					linelist=nothing,#"./input_multi3d/vald_2490-25540.list",
                    absmet=nothing,
					kwargs...)	
	# Create an empty Namelist
	nml = M3DNamelist()

	# Fill in the defaults for this, they can be modified
    NLTE ? _spectrum_namelist_nlte!(nml) : _spectrum_namelist_lte!(nml)

	# path of the model
	model_path = joinpath(model_folder, model_name)

	# additionally apply the model specific fields
	set!(nml; atmos_params=(:atmos_file=>model_path,))

    if (isnothing(absmet)) & (!isnothing(linelist))
        set!(nml; linelist_params=(:linelist_file=>linelist,))
    elseif (!isnothing(absmet)) & (isnothing(linelist))
        set!(nml, composition_params=(:absmet_file=>absmet,))
        set!(nml, m3d_params=(:absmet_in_output=>true, ))
    elseif (isnothing(absmet)) & (isnothing(linelist))
        nothing
    else
        error("Need to specify either linelist or absmet (not both).")
    end

    # apply custom input parameters
    set!(nml; kwargs...)

	nml
end

"""
    heating_namelist(model_name, args...; kwargs...)

Compute the radiative heating namelist within the given model.
"""
heating_namelist(model_name::String, opacity_table=nothing, args...; kwargs...) = begin
    nml = spectrum_namelist(model_name, args...; kwargs...)

    if !isnothing(opacity_table)
        set!(nml, composition_params=(:opac_table=>opacity_table,))
    end

    #set!(nml, spectrum_params=(:daa=>1., :aa_blue=>1500, :aa_red=>9000))    
    #set!(nml; kwargs...)

    set!(nml, m3d_params=(:save_patch_kind=>"mean_rad", :save_Qrad=>true, :n_nu=>1))
    set!(nml, atmos_params=(:dims=>1,))

    nml
end

# ============================================================================= 
# Namelist for simple tau500 computation
# =============================================================================

_tau500_namelist!(nml::M3DNamelist,
    io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
    timer_params=(:sec_per_report=>1e8, :verbose=>2),
    atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_rho=>true),
    m3d_params=(:verbose=>2, :long_scheme=>"none")) = begin
    set!(
        nml; 
        io_params=io_params, 
        timer_params=timer_params,
        atmos_params=atmos_params,
        m3d_params=m3d_params,
    )
end

function tau500_namelist(model_name::String; model_folder="./input_multi3d/MUST", kwargs...)	
    # Create an empty Namelist
    nml = M3DNamelist()

    # Fill in the defaults for this, they can be modified
    _tau500_namelist!(nml)

    # path of the model
    model_path = joinpath(model_folder, model_name)

    # additionally apply the model specific fields
    set!(nml; atmos_params=(:atmos_file=>model_path,))

    # apply custom input parameters
    set!(nml; kwargs...)

    nml
end