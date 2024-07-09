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
	Py(dir) in sys."path" ? nothing : sys."path".append(Py(dir))
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



include_helper(name) = joinpath(dirname(@__FILE__), name)


"""
    ingredients(path::String)

Include the julia script given into the current namespace.
The file be looked for in the src folder where ` MUST`` is located.
"""
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	path = include_helper(path)
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end



"""
Split array arr in nsplits roughly equal junks.
If mask=true the indices will be returned.
"""
function split_similar(arr, nsplits; mask=false)
    nsplits = length(arr) < nsplits ? length(arr) : nsplits
    splits  = div(length(arr), nsplits)
    nrest   = length(arr) % nsplits
    split_sizes = Int[splits for _ in 1:nsplits]
    
    for i in 1:nrest
        split_sizes[i] = split_sizes[i] + 1
    end

    split_masks = []
    last_idx    = 0
    for i in 1:nsplits
        append!(split_masks, [[ (last_idx+1:last_idx+split_sizes[i])... ]])
        last_idx = last_idx+split_sizes[i]
    end

    if mask
        return split_masks
    else
        return [arr[mask] for mask in split_masks]
    end
end

"""
create random number between a and b
"""
randrange(a,b,args...) = begin
    xmin,xmax = min(a,b),max(a,b)
    rand(args...) .* (xmax-xmin) .+ xmin
end

function meshgrid(ax...)
    grids = []
    space = Iterators.product(ax...)
    for i in eachindex(ax)
        append!(grids, [getindex.(space, i)])
    end

    grids
end

uniqueidx(v) = unique(i -> v[i], eachindex(v))

"""
    roundto(x, y)

Round the given value `x` to multiples of `y`.

# Examples
```jldoctest
julia> roundto(π, 0.01)
3.14
julia> roundto(π, 0.25)
3.25
julia> roundto(2.6086956525e8, 0.25)
2.5e8
```
"""
roundto(x, y; magnitude=nothing) = begin
    order_of_magnitude = if isnothing(magnitude) 
        exp10(round(log10(x)))
    else
        exp10(round(log10(x))) / magnitude
    end

    x = x / order_of_magnitude
    x = if (x % y) >= y/2
        div(x, y) * y + y
    else
        div(x, y) * y
    end

    x = x * order_of_magnitude
end




#= Integration functions =#

"""
    integrate(x, y; [method])

Integrate the values in the y array defined on the grid of the x array.
"""
function integrate(x, y; method=QuadGKJL())
	mask = sortperm(x)
	xs, ys = x[mask], y[mask]

	ip = Interpolations.linear_interpolation(xs, ys)
	func(xi, p) = ip(xi)

	prob = IntegralProblem(func, first(xs), last(xs))
	solve(prob, method).u
end





#= stuff for timing =#

"""
    activate_timing!(t)

Activate timing of the given timer. See `timers` for a list of available timers.
"""
activate_timing!(t) = begin
    reset_timer!(generalTimer)
    t[] = true
end

"""
    deactivate_timing!(t)

Deactivate timing of the given timer. See `timers` for a list of available timers.
"""
deactivate_timing!(t) = begin
    reset_timer!(generalTimer)
    t[] = false
end

"""
    start_timing!(t=generalTimer)

Reset the timer.
"""
start_timing!(t=generalTimer) = reset_timer!(t) 

"""
    end_timing!(t=generalTimer)

Reset the timer.
"""
end_timing!(t=generalTimer) = begin
    println("")
    show(t)
    println("")
    reset_timer!(t) 
end

macro optionalTiming(name, exp)
    name_e = esc(name)
    ex = esc(exp)
    name_string = "$(name)"
    quote
        if $(name_e)[]
            @timeit generalTimer $(name_string) begin
                $(ex)
            end
        else
            $(ex)
        end
    end
end



"""
    parametersFromName(name)

Try to extract stellar parameters from the name. Parameters need to be given as
tXXXXgXXXmXXXX. Please make sure to include the minus sign in m.
"""
parametersFromName(name) = begin
	names = split(name, "_", keepempty=false)
	i_name = 0
	for (i, split) in enumerate(names)
		if occursin("t", split) & occursin("g", split) & occursin("m", split)
			i_name = i
		end
	end

	t, g, m = if i_name == 0
		@warn "could not guess atmospheric parameters from name."
		-1, -1, 0.0
	else
        relevantPart = names[i_name]
        t1 = findfirst(x->x=='t', relevantPart) + 1
        t2 = findfirst(x->x=='g', relevantPart) - 1
        t3 = findfirst(x->x=='m', relevantPart) - 1

        t = parse(Float64, relevantPart[t1:t2])
        g = parse(Float64, relevantPart[t2+2:t3]) ./10
        m = parse(Float64, relevantPart[t3+2:end])

		t, g, m
	end

	t, g, m
end