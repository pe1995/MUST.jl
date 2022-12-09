abstract type AbstractEOS end

# Type for squaregas EOS tables
SqGTable = Union{Nothing, Array{F,N} where {N,F<:AbstractFloat}}
RangeOrVector = Union{R where {R<:AbstractRange}, Vector{F} where {F<:Integer}}

"""
Wrapper for python EOS.
"""
struct PythonEOS <: AbstractEOS
    eos          :: PyObject
    variables    :: Vector{Symbol}
    log_variable :: Vector{Bool}
end

"""Squaregas EOS in the Bifrost format."""
struct SquareGasEOS{NE<:SqGTable, RO<:SqGTable, 
                    EO<:SqGTable, EX<:SqGTable, 
                    LN<:SqGTable, TH<:SqGTable,
                    AX<:SqGTable, BX<:SqGTable, 
                    T<:AbstractFloat} <: AbstractEOS

    params           ::Dict{String, Any}
    variables        ::Vector{Symbol}
    lnRho_axis       ::Vector{T}
    lnEi_axis        ::Vector{T}
    lnTg_axis        ::AX
    lnNe_axis        ::BX
    ne_epstable      ::NE
    ne_temtable      ::NE
    ne_opatable      ::NE
    epstable         ::RO
    temtable         ::RO
    opatable         ::RO
    eostable         ::EO
    expieos_err      ::EX  
    lndlnT_table     ::LN  
    theta_rho_table  ::TH
end





#=== Init routines ===#

"""
Create a SquareGasEOS from the Tabgen input file.
"""
function SquareGasEOS(para_path)
    d         = para_path
    para_path = joinpath(d, "tabparam.in")
    p = read_squaregas_params(para_path)

    # Record length
    p["RhoEi_recl"]          = p["nEiBin"] *p["nRhoBin"] *4
    p["RhoEiRadTable_recl"]  = p["nEiBin"] *p["nRhoBin"] * p["nRadBins"]

    if haskey(p, "nNeBin") && haskey(p, "nTBin")
        p["NeTgRadTable_recl"] = p["nNeBin"] *p["nTBin"] * p["nRadBins"] * 2
    else
        p["NeTgRadTable_recl"] = -1
    end

    # Read the EOS tables
    # Ne - Tg Rad tables
    if p["NeTgRadTable_recl"] != -1
        path = joinpath(d, p["NeTgRadTableFile"])
        try
            ne_epstable = read_squaregas_table(path, record = 1, recl = p["NeTgRadTable_recl"], dim = (p["nTBin"],p["nNeBin"],p["nRadBins"],2))  
            ne_temtable = read_squaregas_table(path, record = 2, recl = p["NeTgRadTable_recl"], dim = (p["nTBin"],p["nNeBin"],p["nRadBins"],2))  
            ne_opatable = read_squaregas_table(path, record = 3, recl = p["NeTgRadTable_recl"], dim = (p["nTBin"],p["nNeBin"],p["nRadBins"],2)) 
        catch
            ne_epstable = nothing
            ne_temtable = nothing
            ne_opatable = nothing
        end
    else
        ne_epstable = nothing
        ne_temtable = nothing
        ne_opatable = nothing
    end
    
    # Rho - Ei - rad tables
    path = joinpath(d, p["RhoEiRadTableFile"])
    epstable = read_squaregas_table(path, record = 1, recl = p["RhoEiRadTable_recl"], dim = (p["nEiBin"],p["nRhoBin"],p["nRadBins"]))  
    temtable = read_squaregas_table(path, record = 2, recl = p["RhoEiRadTable_recl"], dim = (p["nEiBin"],p["nRhoBin"],p["nRadBins"]))  
    opatable = read_squaregas_table(path, record = 3, recl = p["RhoEiRadTable_recl"], dim = (p["nEiBin"],p["nRhoBin"],p["nRadBins"]))  

    # EOS tables
    path = joinpath(d, p["EOSTableFile"])
    eostable = read_squaregas_table(path, record = 1, recl = p["RhoEi_recl"], dim = (p["nEiBin"],p["nRhoBin"],4))  

    # Type of tables
    TT = eltype(eostable)

    # The axes of the tables
    lnRho_axis = Vector(range(TT(log(p["RhoMin"])), TT(log(p["RhoMax"])); length=p["nRhoBin"]))
    lnEi_axis  = Vector(range(TT(log(p["EiMin"])),  TT(log(p["EiMax"]));  length=p["nEiBin"]))

    # Add Tg and Ne axes
    if haskey(p,"TMin") && haskey(p,"TMax")
        lnTg_axis  = Vector(range(TT(log(p["EiMin"])),  TT(log(p["EiMax"]));  length=p["nEiBin"]))
    else
        lnTg_axis = nothing
    end

    if haskey(p,"NeMin") && haskey(p,"NeMax")
        try
            lnNe_axis  = Vector(range(TT(log(p["NeMin"])),  TT(log(p["NeMax"]));  length=p["nNeBin"]))
        catch
            lnNe_axis = nothing
        end
    else
        lnNe_axis = nothing
    end

    SquareGasEOS(p, [:d, :ee], lnRho_axis, lnEi_axis, lnTg_axis, lnEi_axis, ne_epstable, ne_temtable, ne_opatable, epstable, temtable, opatable, eostable, nothing, nothing, nothing)
end

function read_squaregas_params(file_name::String)
    f = open(file_name,"r")
    l = readlines(f);
    l = [ strip(i) for i in l if !isempty(strip(i)) ]; # remove empty str
    l = [ i for i in l if !occursin("/", i) ]; # remove header
    l = [ i for i in l if !(i[1]=='&') ]; # remove header
    l = [ strip(i) for i in l if (i[1] != ';')]
    l = [ i[1]=='!' ? strip(i[2:end]) : i for i in l ]; # remove empty str
    l = [ replace(i,"'" => "\"") for i in l ];
    l = [ split(i,'=') for i in l ]; # remove lines which starts with ';'
    params = Dict{String,Any}()
    for x in l
        key = strip(x[1])
        val = eval(Meta.parse(x[2]))
        params[key] = val
    end
    return params
end

function read_squaregas_table(path; record, recl, dim, out_type=Float32)
    ispath(path) ? nothing : return nothing
    f = FortranFile(path,"r",access="direct",recl=recl*4);
    read(f,rec=record,(out_type,dim));
end





#=== Lookup routines ===#

"""
Lookup function. Log is done automaticall if requested in EOS
"""
lookup(eos::PythonEOS, parameter, variables::AbstractVector...) = begin
    in_var = [eos.log_variable[i] ? log.(v) : v for (i,v) in enumerate(variables)]
    eos.eos.lookup("$(parameter)", in_var...)
end

lookup(eos::PythonEOS, parameter, variables::T...) where {T<:AbstractFloat} = lookup(eos, parameter, [T[v] for v in variables]...)

"""
Lookup function for SquareGas EOS.
    The input is assumed to be CGS. 

    Example:
        lookup(sqg_eos, :T, [1e-7, 2e-7, 3e-7], [4e-12, 5e-12, 6e-12])

    The EOS tabe is listed in log units, so you can also pass
        lookup(sqg_eos, :T, log.([1e-7, 2e-7, 3e-7]), log.([4e-12, 5e-12, 6e-12]), to_log=false)
"""
function lookup(eos::SquareGasEOS, parameter::Symbol, d::AbstractVector{T}, ee::AbstractVector{T}; to_log=true) where {T<:AbstractFloat}
    if to_log
        d  = log.(d)
        ee = log.(ee)
    end
    
    # The EOS table is searched based on the regular grid imposed by the corresponding axis
    result = if ((parameter == :P) || (parameter == :Pg))
        exp.(interpolate_in_table(eos, eos.eostable, d, ee, 1))
    elseif ((parameter == :T) || (parameter == :Temp))
             interpolate_in_table(eos, eos.eostable, d, ee, 2)
    elseif parameter == :Ne
        exp.(interpolate_in_table(eos, eos.eostable, d, ee, 3))
    elseif ((parameter == :kr) || (parameter == :ross))
        exp.(interpolate_in_table(eos, eos.eostable, d, ee, 4))
    elseif ((parameter == :src) || (parameter == :Src))
        exp.(interpolate_in_table(eos, eos.temtable, d, ee, 1:eos.params["nRadBins"]))
    elseif parameter == :rk
        exp.(interpolate_in_table(eos, eos.opatable, d, ee, 1:eos.params["nRadBins"]))
    else
        error("The given parameter $(parameter) is not implemented for the SquareGasEOS.")
    end

    return result
end

"""
Interpolate in the given table using the axis in the EOS. 
"""
function interpolate_in_table(eos::SquareGasEOS, table::SqGTable, d::AbstractVector{T}, ee::AbstractVector{T}, index::RangeOrVector) where {T<:AbstractFloat}
    result = zeros(eltype(d), length(d), length(index))

    ei_min  = eos.lnEi_axis[1]
    rho_min = eos.lnRho_axis[1]

    dlnRho = eos.lnRho_axis[2] - eos.lnRho_axis[1]
    dlnEi  = eos.lnEi_axis[2] - eos.lnEi_axis[1]

    il=0; jl=0; ifr=0.0; jfr=0.0; tt=0.0; u=0.0

    @inbounds for i in 1:length(d)
        ifr = (ee[i] - ei_min) / dlnEi + 1.0
        il  = floor(Int, ifr)
        tt  = Base.mod(ifr, 1.0) 
        m1t = 1.0 - tt
    
        jfr = (d[i] - rho_min) / dlnRho + 1.0
        jl  = floor(Int, jfr)
        u   = Base.mod(jfr, 1.0) 
        m1u = 1.0 - u

        if ((ifr<0) || (jfr<0))
            result[i,:] .= Base.convert(eltype(d),NaN) 
        else
            result[i,:] .= m1t .* m1u .* table[il,  jl  ,index] +
                          tt  .* m1u .* table[il+1,jl  ,index] +
                          tt  .*   u .* table[il+1,jl+1,index] +
                          m1t .*   u .* table[il  ,jl+1,index] 
        end
    end

    result
end

interpolate_in_table(eos::SquareGasEOS, table::SqGTable, d::AbstractVector{T}, ee::AbstractVector{T}, index::I) where {T<:AbstractFloat, I<:Integer} = interpolate_in_table(eos, table, d, ee, [index])[:,1]

lookup(eos::SquareGasEOS, parameter::Symbol, d::T, ee::T; to_log=true) where {T<:AbstractFloat} = first(lookup(eos, parameter, [d], [ee]; to_log=to_log))
lookup(eos::SquareGasEOS, parameter::String, args...; kwargs...) = lookup(eos, Symbol(parameter), args...; kwargs...)
lookup(eos::AbstractEOS, parameter::Symbol, d::AbstractArray{T,2}, ee::AbstractArray{T,2}, args...; kwargs...) where {T<:AbstractFloat} = begin
    dp     = dimension(eos, parameter)
    result = rand(eltype(d), size(d)..., dp)

    @inbounds for i in 1:size(d, 2)
        di  = view(d,  :, i)
        eei = view(ee, :, i)

        dp > 1 ? 
            result[:,i,:] .= lookup(eos, parameter, di, eei, args...; kwargs...) : 
            result[:,i,1] .= lookup(eos, parameter, di, eei, args...; kwargs...)
            
    end

    dp > 1 ? result : result[:,:,1]
end

lookup(eos::AbstractEOS, parameter::Symbol, d::AbstractArray{T,3}, ee::AbstractArray{T,3}, args...; kwargs...) where {T<:AbstractFloat} = begin
    dp     = dimension(eos, parameter)
    result = rand(eltype(d), size(d)..., dp)

    @inbounds for k in 1:size(d, 3)
        @inbounds for j in 1:size(d, 2)
            di  = view(d, :, j, k)
            eei = view(ee, :, j, k)

            dp > 1 ? 
                result[:,j,k,:] .= lookup(eos, parameter, di, eei, args...; kwargs...) : 
                result[:,j,k,1] .= lookup(eos, parameter, di, eei, args...; kwargs...)
        end   
    end

    dp > 1 ? result : result[:,:,:,1]
end

dimension(eos::SquareGasEOS, parameter::Symbol) = begin
    d = if parameter == :rk
        eos.params["nRadBins"]
    elseif ((parameter == :src) || (parameter == :Src))
        eos.params["nRadBins"]
    else
        1
    end

    d
end

dimension(eos::PythonEOS, parameter::Symbol) = 1


"""
Inverse lookup parameter in the EOS tables. 
"""
function inverse_lookup(eos::SquareGasEOS, iterations=50, antilinear=false; kwargs...)
end





#=== Python EOS aliases ===#

"""
Generate a legacy Stagger EOS (Python).
"""
macro legacyPythonEOS(path="", table_name="table.dat") 
    p = esc(path)
    t = esc(table_name)
    d = esc(:EOS)
    
    quote
        if $(p) == ""
            PythonEOS($(d).stagger(top=dispatch_location, file=$(t)), [:d,:ee], [:true,:false])
        else
            PythonEOS($(d).stagger(top=$(p), file=$(t)), [:d,:ee], [:true,:false])
        end
    end
end

"""
Generate a Squaregas EOS (Python).
"""
macro squaregasPythonEOS(path="", table_name="eostable.dat", dir="/experiments/stellar_atmospheres/input_data/square_gas/") 
    p = esc(path)
    t = esc(table_name)
    f = esc(dir)
    d = esc(:EOS)
    
    quote
        if $(p) == ""
            PythonEOS($(d).bifrost(top=dispatch_location, file=$(t), dir=$(f)), [:d,:ee], [:true,:true])
        else
            PythonEOS($(d).bifrost(top=$(p), file=$(t), dir=$(f)), [:d,:ee], [:true,:true])
        end
    end
end

"""
Bisect the EOS to get a value of either of the variables 
on the basis of the other one an one parameter. Enter the names and values as kwargs.
The one that shall be found needs to be given a (start, end).
"""
function bisect(eos::AbstractEOS, iterations=50, antilinear=false; kwargs...)
    bisect_var = Symbol("")
    normal_var = Symbol("")
    para       = Symbol("")
    for (key,value) in kwargs
        if length(value) == 2
            bisect_var = key
        elseif !(key in eos.variables)
            para = key
        else
            normal_var = key
        end
    end

    # now loop and find the requested values
    bisect_res = 0
    b0, b1     = kwargs[bisect_var]
    
    bisect_var_pos = first(findfirst(bisect_var .== eos.variables))
    normal_var_pos = first(findfirst(normal_var .== eos.variables))
   
    args = Float64[0.0,0.0]
    args[normal_var_pos] = kwargs[normal_var]

    for i in 1:iterations
        bisect_res           = (b0+b1)/2
        args[bisect_var_pos] = bisect_res

        para_res = lookup(eos, String(para), args...)
        para_res = size(para_res) == 0 ? para_res : para_res[1]
     
        if !(antilinear)
            if para_res > kwargs[para]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        else
            if para_res < kwargs[para]
                b1 = bisect_res
            else
                b0 = bisect_res
            end
        end
    end

    return bisect_res
end


function bisect_cube(eos::AbstractEOS; iterations=50, kwargs...)
    key_values = [[],[],[]]
    i = 1
    for (key,value) in kwargs
        if ndims(value) == 4
            key_values[3] = [key, value]
        else
            key_values[i] = [key, value]
            i += 1
        end
    end
    s = size(key_values[1][2])
    r = zeros(eltype(key_values[1][2]), s...)

    p = zeros(2)
    @inbounds for k in 1:s[3]
        @inbounds for j in 1:s[2]
            @inbounds for i in 1:s[1]
                p   .= key_values[3][2][i, j, k, :]
                args = (key_values[1][1]=>key_values[1][2][i,j,k], key_values[2][1]=>key_values[2][2][i,j,k], key_values[3][1]=>p)

                r[i,j,k] = bisect(eos, iterations; args...)
            end
        end
    end

    r
end

bisect_energy(eos::SquareGasEOS, s, args...; kwargs...) = begin
    _, _, eemin, eemax = limits(eos)
    e_cube_lims = zeros(s..., 2)
    e_cube_lims[:,:,:,1] .= eemin
    e_cube_lims[:,:,:,2] .= eemax
    bisect_cube(eos, args...; ee=e_cube_lims, kwargs...)
end

limits(eos::PythonEOS) = begin
    d0 = eos.eos.scale0.min
    d1 = eos.eos.scale0.max
    e0 = eos.eos.scale1.min
    e1 = eos.eos.scale1.max
    return d0,d1,e0,e1
end

limits(eos::SquareGasEOS) = begin
    emin = exp(eos.lnEi_axis[1])
    emax = exp(eos.lnEi_axis[end])
    dmin = exp(eos.lnRho_axis[1])
    dmax = exp(eos.lnRho_axis[end])
    return dmin, dmax, emin, emax
end

# Avoid Broadcast over EoS
broadcastable(eos::AbstractEOS) = Ref(eos)