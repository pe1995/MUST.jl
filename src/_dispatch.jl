dispatch_location   = nothing
experiment_location = "experiments/stellar_atmospheres/"

"""
Load the dispatch module. Either the global location is used or 
a specific path is given to the function. If submodules are listed, those are
imported only.
"""
macro import_dispatch(location=nothing,submodules...)
	location_l    = esc(location)
	submodules_l  = esc(submodules)
	dispatch_name = esc(:dispatch)
	submodules_e  =[esc(m) for m in submodules]
	expressions   = []
	if length(submodules) >0
		append!(expressions,[:(all_modules = _import_dispatch($(location_l),$(submodules_l)))])
		for i in 1:length(submodules)
			append!(expressions,[:($(submodules_e[i]) = all_modules[$i])])
		end
		return quote
			$(expressions...)
		end
	else
		return quote
			$(dispatch_name) = _import_dispatch($(location_l))
		end
	end
end

"""
Load the dispatch module. Either the global location is used or 
a specific path is given to the function
"""
function _import_dispatch(location::Union{Nothing,String}=nothing, submodules=[])
    isnothing(location) ? location = joinpath(pwd(),"../dispatch2") : nothing
    
    # Import sys and os, add to path and import dispatch
    os  = pyimport("os")
	sys = pyimport("sys")
	global dispatch_location = abspath(os.getenv("D", location))
	dispatch_python_path = os."path".join(abspath(os.getenv("D", location)),
											"utilities/python/")
	isdir(dispatch_python_path) ? nothing : throw(error("$(dispatch_python_path) is not a dispatch installation."))
	dispatch_python_path in sys."path" ? nothing :
										append!(sys."path",[dispatch_python_path])
	
	if length(submodules) == 0
		return pyimport("dispatch")
	else
		return [pyimport("dispatch.$(String(m))") for m in submodules]
	end
end

"""
The given path is converted from relative to the expermiment to absolute.
"""
macro in_dispatch(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_dispatch($(relative_path_l)))
end
_in_dispatch(relative_path) = abspath(joinpath(dispatch_location,experiment_location,relative_path))

function list_of_snapshots(list_of_files::Vector{String})
	los = []
	for i in 1:length(list_of_files)
		if (isdir(list_of_files[i]))
			try
				append!(los, [parse(Int, basename(dirname(list_of_files[i])))])
			catch
				continue
			end
		end
	end
	los
end
