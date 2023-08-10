using MUST
using TSO

import MUST.Box
dispatch = MUST.dispatch

"""
	dispatch2box(simulation, snapshot; eos, opacity, add_selected=false)

Convert a snapshot with number `snapshot` from the simulation output `simulation` to a ```MUST.jl``` ```Box```. EoS and opacity tables in the ```TSO.jl``` format are used to add T and κ-ross.
"""
function Box(simulation, snapshot; eos, opacity, add_selected=false, save=true)
	snap = dispatch.snapshot(snapshot, data=simulation)
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