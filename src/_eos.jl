abstract type AbstractEOS end

"""
Wrapper for python EOS.
"""
struct PythonEOS <: AbstractEOS
    eos          :: PyObject
    variables    :: Vector{Symbol}
    log_variable :: Vector{Bool}
end

"""
Lookup function. Log is done automaticall if requested in EOS
"""
lookup(eos::PythonEOS, parameter, variables...) = begin
    in_var = [eos.log_variable[i] ? log.(v) : v for (i,v) in enumerate(variables)]
    eos.eos.lookup(parameter, in_var...)
end

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