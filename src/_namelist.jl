#= Namelist Types =#

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
end





#= Constructors =#

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
    for field in fieldnames(typeof(s))
        d = getfield(s, field)
        for key in keys(d)
            !(typeof(d[key])<:AbstractArray) ? continue : nothing
            d[key] = parse_from_namelist.(d[key][1:end-1])
        end
    end
    s
end






#= Utility Functions for Namelists =#

"""
Read a namelist from path.
"""
function read!(nml::AbstractNamelist, path::String)
    content = Dict{Symbol, Dict{String, Any}}()

    # Read the data line by line and save the parameter
    lines = readlines(path)
    para  = :nothing
    for line in lines
        occursin("!", line) ? continue : nothing
        
        # If there is a & it is a new parameter, then remove this from the line
        if occursin("&", line) 
            para          = Symbol(lowercase(split(line)[1][2:end]))
            content[para] = Dict{String, Any}()
            line          = line[length(String(para))+2:end]
        end

        # Split the lines at the = sign. Every second entry is value then
        line_s = split_namelist_line(line)
        for i in 1:2:length(line_s)
            (length(line_s) >= 2) ? key = lowercase(strip(line_s[i])) : continue
            content[para][key] = parse_from_namelist(strip(line_s[i+1])) 
        end
    end

    for key in keys(content)
        if !(key in fieldnames(typeof(nml))) 
            @warn "$(key) is not a valid fieldname"
            continue
        end
        setfield!(nml, key, content[key])
    end
end

function parse_from_namelist(value_string)
    if length(value_string) == 0
        return ""
    end
    value_string = ("$(strip(value_string)[end])" == "/") ? strip(strip(value_string)[1:end-1]) : value_string
    val = nothing 
    if occursin(",", value_string) 
        val_spl = split(value_string, ",")
        
        if any(occursin.(".", val_spl))
            val = try
                parse.(Float64, val_spl)
            catch
                String.(val_spl)
            end
        else
            val = try
                parse.(Int64, val_spl)
            catch
                String.(val_spl)
            end
        end

    elseif occursin(".", value_string)
        val = try
            parse(Float64, value_string)
        catch
            String(value_string)
        end
    else
        val = try
            parse(Int64, value_string)
        catch
            String(value_string)
        end
    end
    return val
end

function reverse_parse(value)
    val_str = ""
    if typeof(value) <: AbstractArray
        for val in value
            val_str = val_str * "$(val),"
        end
        val_str = val_str[1:end-1]
    elseif typeof(value) <:String
        if (occursin('"', value)) | (value in ["t", "f", ".true.", ".false."]) | (occursin('*', value))
            val_str = "$(value)"
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
    elseif typeof(value) <:Bool
        val_str = ".$(value)."
    else
        val_str = "$(value)"
    end
    return val_str
end

function split_namelist_line(line)
    # First split at =
    line_s = split(strip(line), "=")
    out    = String[]

    for i in 1:length(line_s)
        if (i>1) & (i<length(line_s))
            line_s_strip = strip(line_s[i])
            i_last = first(findlast(" ", line_s_strip))

            append!(out, [strip(line_s_strip[1:i_last])])
            append!(out, [strip(line_s_strip[i_last:end])])
        else
            append!(out, [strip(line_s[i])])
        end
    end

    out
end

"""
Write a namelist to path.
"""
function write(nml::AbstractNamelist, path::String; soft_line_limit=85)
    para_names   = String.(fieldnames(typeof(nml)))
    longest_name = maximum(length.(para_names))

    f = open(path, "w")

    for (i,p) in enumerate(para_names)
        spaces = longest_name - length(p) +1
        line   = "&$(p)" * repeat(" ", spaces) 
        d      = getfield(nml, Symbol(p))
        
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

"""Set fields of the Namelist. Enter field=(parameter=>value,...)"""
function set!(nml::AbstractNamelist; kwargs...)
    for (field,paras) in kwargs
        if !(typeof(paras) <: Tuple) 
            @warn "Argument to $(field) is not a Tuple. Add trailing comma if there is only one entry."
            continue
        end

        for (p,v) in paras
            getfield(nml, field)[String(p)] = v
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





#= Default namelists =#

_whole_spectrum_namelist!(nml::M3DNamelist; 
	io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
	timer_params=(:sec_per_report=>120,),
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_density=>true),
    atom_params=(:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true),
	m3d_params=(:verbose=>1, :fcheck=>1, :pcheck=>[1,1,1], :linecheck=>1, 
				:lvlcheck=>1,
                :n_nu=>32, :maxiter=>0, :decouple_continuum=>true),
	#linelist_params=(:dlam=>1.0,),
    composition_params=(:abund_file=>"./input_multi3d/abund_magg",:absdat_file=>"./input_multi3d/TS_absdat.dat"),
	spectrum_params=(:daa=>2., :aa_blue=>2000, :aa_red=>40000)) = begin

	set!(
		nml; 
		io_params=io_params, 
		timer_params=timer_params,
		atmos_params=atmos_params
		,m3d_params=m3d_params,
		#linelist_params=linelist_params,
		spectrum_params=spectrum_params,
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
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_density=>true),
    atom_params=(:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true),
	m3d_params=(:verbose=>1, :fcheck=>1, :pcheck=>[1,1,1], :linecheck=>1, 
				:lvlcheck=>1,
                :n_nu=>32, :maxiter=>0)) = begin
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
	atmos_params=(:dims=>23, :atmos_format=>"MUST", :use_density=>true),
    atom_params=(:atom_file=>"./input_multi3d/atoms/atom.h20", 
                :exclude_trace_cont=>true, :exclude_from_line_list=>true,
                :convlim=>1e-2),
	m3d_params=(:verbose=>1, :fcheck=>1, :pcheck=>[1,1,1], :linecheck=>1, 
				:lvlcheck=>1,
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

    set!(nml, m3d_params=(:save_patch_kind=>"mean_rad", :save_Qrad=>true, :ilambd=>1, :n_nu=>1))
    set!(nml, atmos_params=(:dims=>1,))

    nml
end