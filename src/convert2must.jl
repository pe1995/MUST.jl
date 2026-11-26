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

function _Box_py(simulation, snapshot; eos, data, add_selected=false, save=true)
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

	lnRoss, lnT = TSO.sample(eos, :lnRoss, :lnT; lnRho=log.(bs[:d]), lnEi=log.(bs[:ee]))
	T = exp.(lnT)
	MUST.add!(bs, exp.(lnRoss), :kr)
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
			MUST.@in_dispatch(replace(MUST.nmlField(inputNamelist, :eos_params)["table_loc"], "'"=>""))
		else
			nothing
		end
	else
		eos_path
	end
	if !isnothing(eos_path) && ispath(eos_path)
		TSO.dispatch_eos_reader(eos_path)(run)
	else
		@warn "No EoS found!"
		nothing, []
	end
end

lookup_function_generator(eos, q) = TSO.dispatch_lookup_generator(eos, q)
