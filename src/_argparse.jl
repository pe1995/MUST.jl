function MUSTGridArgs()
    s = ArgParse.ArgParseSettings()  # The parser

    ArgParse.@add_arg_table s begin
        "--info_path", "-p"
            help = "path where the info/summary csv will be stored (leave out .csv).\n"*
                   "The summary can be given to a MUSTGrid in the info field."
            arg_type = String
            default = "summary"

        "--index_range", "-i"
            help     =  "The Namelists of the grid will be named according to this range.\n"*
                        "The length of this range is the number of atmospheres to create."
            arg_type = Vector{Int}
            default  = [1,10]
    end

    return ArgParse.parse_args(s;as_symbols=true)
end

function ArgParse.parse_item(::Type{Vector{T}},x::AbstractString) where {T}
    splitted_string = strip.(split(x[2:end-1],","))
    return [ArgParse.parse_item(T,s) for s in splitted_string]
end

function ArgParse.parse_item(::Type{Dict{S,T}},x::AbstractString) where {S,T}
    splitted_string = strip.(split(x[2:end-1],","))
    splitted_string = split.(splitted_string,"=>")
    return Dict(ArgParse.parse_item(S,s[1])=>ArgParse.parse_item(T,s[2]) for s in splitted_string)
end