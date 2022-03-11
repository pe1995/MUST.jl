macro in_MUST(source_of_help)
    source_of_help_l = esc(source_of_help)
    folder_of_file   = dirname(@__FILE__)
    joinpath_l       = esc(Base.joinpath)
    :($joinpath_l($(folder_of_file),String($(source_of_help_l))*".jl"))
end

macro get_help(source_of_help)
    #source_of_help_l = esc(source_of_help)
    source_of_help_s = :($source_of_help)
    folder_of_file   = dirname(@__FILE__)
    #joinpath_l       = esc(Base.joinpath)
    include_l        = esc(Base.include)
    main             = esc(Main)
    quote
        path = joinpath($(folder_of_file),String($(QuoteNode(source_of_help_s)))*".jl")
        if !isdefined($main, $(QuoteNode(source_of_help_s)))
            $include_l($main,path)
        end
    end
end

function _get_help_py(mod ,dir=dirname(@__FILE__))
	sys   = pyimport("sys")
	dir in sys."path" ? nothing : append!(sys."path",[dir])
    pyimport(String(mod))
end

macro get_help_py(mod)
    mod_e = esc(mod)
    mod_s = :($mod)
    path  = dirname(@__FILE__)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path))
end

macro get_help_py(mod, dir)
    mod_e    = esc(mod)
    path_esc = esc(dir)
    mod_s = :($mod)
    :($(mod_e) = _get_help_py($(QuoteNode(mod_s)), $path_esc))
end