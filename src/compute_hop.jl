using MUST

MUST.@import_hop

# ============================================================================
# Compute EoS + Opacities from the standalone HOP
# ============================================================================

function compute_eos_HOP(
	T, rho, lam;
	HOP_kwargs=Dict{Symbol,Any}(:feh => 0.0),
	mmap=true,
	linelist=nothing,
	line_margin=100.0, include_absmet=false, line_scattering=true, vmic=nothing,
	nthreads=1, nprocs=nothing,
	t_block=length(T), rho_block=length(rho),
	Teltype=Float32,
	verbose=false,
)
	chi, scat, E2d, pg2d, ne2d, chi500_2d = HOP.compute_eos(
		T, rho, lam;
		HOP_kwargs=HOP_kwargs,
		mmap=mmap,
		linelist=linelist,
		line_margin=line_margin, include_absmet=include_absmet, line_scattering=line_scattering, vmic=vmic,
		nthreads=nthreads, nprocs=nprocs,
		t_block=t_block, rho_block=rho_block,
		Teltype=Teltype,
		verbose=verbose
	)

	_nan_out_of_range!(chi)
	_nan_out_of_range!(scat)

	_nan_out_of_range!(E2d)
	_nan_out_of_range!(pg2d)
	_nan_out_of_range!(ne2d)
	_nan_out_of_range!(chi500_2d)

	rho_axis  = Teltype.(log.(rho))
	temp_axis = Teltype.(log.(T))
	Elog      = Teltype.(log.(E2d))
	pglog     = Teltype.(log.(pg2d))
	nelog     = Teltype.(log.(ne2d))
	chi500log = Teltype.(log.(chi500_2d))
	l_axis    = Teltype.(lam)

	rho_axis, temp_axis, Elog, pglog, nelog, chi500log, l_axis, chi, scat
end

function _nan_out_of_range!(A)
	@inbounds for i in eachindex(A)
		v = A[i]
		if v < 1e-30 || v > 1e30 || isnan(v)
			A[i] = NaN
		end
	end
	A
end

# ============================================================================
# Constants
# ============================================================================

const HOP_TILE_ELEMS_AUTO = 1 << 28
const HOP_TILE_ELEMS_CEIL = 1 << 30
