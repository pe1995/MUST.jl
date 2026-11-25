using MUST
using TSO

import MUST.Box
dispatch = MUST.dispatch

MUST.Box(simulation, snapshot; data=MUST.@in_dispatch("data"), legacy=false, eos_path=nothing, kwargs...) = begin
	if legacy
		_Box_py(simulation, snapshot; data=data, kwargs...)
	else
		_Box_jl(simulation, snapshot; data=data, eos_path=eos_path, kwargs...)
	end
end

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

	κross = lookup(eos, opacity, :κ_ross, log.(bs[:d]), log.(bs[:ee]))
	T = exp.(lookup(eos, :lnT, log.(bs[:d]), log.(bs[:ee])))
	MUST.add!(bs, κross, :kr)
	MUST.add!(bs, T, :T)

	τ = MUST.optical_depth(bs, opacity=:kr, density=:d)
    MUST.add!(bs, τ, :τ_ross)

    merge!(bs.meta, Dict(:snapshot=>snapshot, :kind=>"geo", :run=>last(splitpath(abspath(simulation)))))
    box_path = save ? MUST.save(bs; name="box_sn$(snapshot)", folder=simulation) : nothing

    bs_τ = MUST.height_scale_fast(bs, :τ_ross)
    if save
        MUST.save!(box_path; tau=bs_τ)
    end

	bs, bs_τ
end

function _Box_jl(simulation, snapshot; data, mmap=false, eos_path=nothing)
	MUST.snapshotBox(
		snapshot,
		eos_path=eos_path,
		folder=joinpath(data, simulation),
		eos_reader=(x)->sqeos_reader(x, eos_path=eos_path, inputNamelist=joinpath(folder, "params")),
		lookup_generator=lookup_function_generator,
		use_mmap=mmap,
		legacy=false
	)
end

function sqeos_reader(run; inputNamelist=@in_dispatch(run), eos_path=nothing)
	eos_path = if isnothing(eos_path)
		if isfile(inputNamelist*".nml")
			inputNamelist = MUST.StellarNamelist(inputNamelist*".nml")
			p = replace(
				MUST.nmlField(inputNamelist, :eos_params)["table_loc"],
				"'"=>""
			)
			MUST.@in_dispatch(joinpath(p, "eos.hdf5"))
		else
			nothing
		end
	else
		eos_path
	end
	if (!isnothing(eos_path)) && isfile(eos_path)
		TSO.reload(SqEoS, eos_path), [:T, :kr, :Pg, :Ne]
	else
		@warn "No EoS found!"
		nothing, []
	end
end

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

	lf = TSO.lookup_function(TSO.@axed(eos), q)
	lkp!(res, args...) = begin
		res .= exp.(TSO.lookup(args...))
	end

	lf, lkp!
end
