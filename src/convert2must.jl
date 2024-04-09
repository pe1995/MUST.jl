using MUST
using TSO

import MUST.Box
dispatch = MUST.dispatch

MUST.Box(simulation, snapshot; data=MUST.@in_dispatch("data"), legacy=false, kwargs...) = begin
	if legacy
		_Box_py(simulation, snapshot; data=data, kwargs...)
	else
		_Box_jl(simulation, snapshot; data=data, kwargs...)
	end
end

"""
	_Box_py(simulation, snapshot; eos, opacity, add_selected=false, save=true)

Python + Julia version:
Convert a snapshot with number `snapshot` from the simulation output `simulation` to a ```MUST.jl``` ```Box```. EoS and opacity tables in the ```TSO.jl``` format are used to add T and κ-ross.
"""
function _Box_py(simulation, snapshot; eos, opacity, data, add_selected=false, save=true)
	data = joinpath(data, simulation)
	snap = dispatch.snapshot(snapshot, data=data)
    units = MUST.StandardUnits(snap)

	s = if add_selected
		si = MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg)
		MUST.convert!(
			si, units; 
			d=:d, ee=:ee, e=:e, 
			qr=:qr, pg=:p,
			ux=:u, uy=:u, uz=:u,
			x=:l, y=:l, z=:l, time=:t
		)

		si
	else
		si = MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)
		MUST.convert!(
			si, units; 
			d=:d, ee=:ee, e=:e, 
			ux=:u, uy=:u, uz=:u,
			x=:l, y=:l, z=:l, time=:t
		)
		
		si
	end

	
    bs = MUST.Box(s)
	s = nothing

	# Get rosseland opacity and T from Opacity table
	κross = lookup(eos, opacity, :κ_ross, log.(bs[:d]), log.(bs[:ee]))
	T = exp.(lookup(eos, :lnT, log.(bs[:d]), log.(bs[:ee])))
	MUST.add!(bs, κross, :kr)
	MUST.add!(bs, T, :T)

	# Compute an vertical optical depth
	τ = MUST.optical_depth(bs, opacity=:kr, density=:d)
    MUST.add!(bs, τ, :τ_ross)

    # Save the box in the same folder
    if save
	    MUST.save(bs; name="box_sn$(snapshot)", folder=simulation)
	end

    bs_τ = MUST.height_scale_fast(bs, :τ_ross)
    if save
        MUST.save(bs_τ; name="box_tau_sn$(snapshot)", folder=simulation)
    end

	bs, bs_τ
end


"""
	_Box_jl(simulation, snapshot; eos, opacity, add_selected=false, save=true)

Julia version:
Convert a snapshot with number `snapshot` from the simulation output `simulation` to a ```MUST.jl``` ```Box```. EoS and opacity tables in the ```TSO.jl``` format are used to add T and κ-ross.
"""
function _Box_jl(simulation, snapshot; data, mmap=false)
    bs = MUST.Box(snapshot, run=simulation, data=data, eos_reader=sqeos_reader, lookup_generator=lookup_function_generator, mmap=mmap)

	# Compute an vertical optical depth
	τ = MUST.optical_depth(bs, opacity=:kr, density=:d)
    MUST.add!(bs, τ, :τ_ross)

    bs_τ = MUST.height_scale_fast(bs, :τ_ross)
 	MUST.save(bs_τ; name="box_tau_sn$(snapshot)", folder=joinpath(data, simulation))
	bs_τ = MUST.Box("box_tau_sn$(snapshot)", folder=joinpath(data, simulation))

	bs, bs_τ
end


"""
	sqeos_reader(run)

Read the EoS for a given run from the name of the run.
Return the names that should be saved to the cube. 
Either give e.g. :T or :lnT, in agreement with the 
`lookup_function_generator`.
"""
function sqeos_reader(run)
	inputNamelist = MUST.@in_dispatch run
	eos_path = if isfile(inputNamelist*".nml")
    	inputNamelist = MUST.StellarNamelist(inputNamelist*".nml")    
		p = replace(
			MUST.nmlField(inputNamelist, :eos_params)["table_loc"], 
			"'"=>""
		)
		joinpath(p, "eos.hdf5")
	else
		nothing
	end
	if (!isnothing(eos_path)) && isfile(MUST.@in_dispatch(eos_path))
		TSO.reload(SqEoS, MUST.@in_dispatch(eos_path)), [:T, :kr, :Pg, :Ne]
	else
		@warn "No EoS found!"
		nothing, []
	end
end

"""
	lookup_function_generator(eos, q)

Create a generator for the TSO lookup functions.
Given the EoS and quantity, return the interpolated function
and the respective function to call for the final lookup.
"""
function lookup_function_generator(eos, q)
	q = if q == :T
		:lnT
	elseif q == :Pg
		:lnPg
	elseif q == :Ne
		:lnNe
	elseif q == :kr
		:lnRoss
	else
		q
	end
		
	TSO.lookup_function(TSO.@axed(eos), q), (args...) -> exp.(TSO.lookup(args...))
end

