dispatch_location   = nothing
experiment_location = "experiments/stellar_atmospheres/"

multi_location   = nothing
multi_experiment_location = ""

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
	global dispatch_location = abspath(pyconvert(String, os.getenv("D", location)))
	dispatch_python_path = joinpath(dispatch_location, "utilities/python/")

	dispatch_python_path = if isdir(dispatch_python_path)
		os."path".abspath(os."path".join(os.getenv("D", location), "utilities/python/"))
	else
		throw(error("$(dispatch_location) is not a dispatch installation."))
	end
		
	if !(dispatch_python_path in sys."path")
		sys."path".append(dispatch_python_path)
		sys."path".append(os."path".join(dispatch_python_path, Py("dispatch")))
	end
	
	
	if length(submodules) == 0
		return pyimport("dispatch")
	else
		return [pyimport("$(String(m))") for m in submodules]
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

function _snapshot_folder(i_snap::Int, list_of_files::Vector{String})
	snap_int = -1
	snap_str = ""
	for i in 1:length(list_of_files)
		if (isdir(list_of_files[i]))
			try
				snap_int = parse(Int, basename(dirname(list_of_files[i])))
			catch
				continue
			end
			
			if (snap_int == i_snap) 
				snap_str=list_of_files[i] 
				break 
			end
		end
	end
	
	(snap_str == "") ? error("Snapshot $(i_snap) not found.") : nothing
	snap_str
end





#= MULTI related functionality =#

"""
Load the dispatch module. Either the global location is used or 
a specific path is given to the function. If submodules are listed, those are
imported only.
"""
macro import_m3dis(location=nothing,submodules...)
	location_l    = esc(location)
	submodules_l  = esc(submodules)
	dispatch_name = esc(:m3dis)
	submodules_e  =[esc(m) for m in submodules]
	expressions   = []
	if length(submodules) >0
		append!(expressions,[:(all_modules = _import_m3dis($(location_l),$(submodules_l)))])
		for i in 1:length(submodules)
			append!(expressions,[:($(submodules_e[i]) = all_modules[$i])])
		end
		return quote
			$(expressions...)
		end
	else
		return quote
			$(dispatch_name) = _import_m3dis($(location_l))
		end
	end
end

"""
Load the dispatch module. Either the global location is used or 
a specific path is given to the function
"""
function _import_m3dis(location::Union{Nothing,String}=nothing, submodules=[])
    # Import sys and os, add to path and import dispatch
    os  = pyimport("os")
	sys = pyimport("sys")
	global multi_location = abspath(pyconvert(String, os.getenv("D", location)))

	multi_python_path = if isdir(joinpath(multi_location, "experiments/Multi3D/m3dis"))
		global multi_experiment_location = joinpath(multi_location, "experiments/Multi3D")
		os."path".abspath(os."path".join(os.getenv("D", location), "experiments/Multi3D"))
	elseif isdir(joinpath(multi_location, "m3dis"))
		os."path".abspath(os.getenv("D", location))
	else
		throw(error("$(multi_location) is not a m3dis installation."))
	end
		
	if !(multi_python_path in sys."path")
		sys."path".append(multi_python_path)
		#append!(sys."path",[joinpath(multi_python_path, "m3dis")])
	end

	if length(submodules) == 0
		return pyimport("m3dis")
	else
		return [pyimport("$(String(m))") for m in submodules]
	end
end

"""
The given path is converted from relative to the expermiment to absolute.
"""
macro in_m3dis(relative_path)
	relative_path_l = esc(relative_path)
	:(_in_m3dis($(relative_path_l)))
end
_in_m3dis(relative_path) = abspath(joinpath(multi_location,multi_experiment_location,relative_path))
