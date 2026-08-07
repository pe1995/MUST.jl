using MUST

MUST.@import_hop

# ============================================================================
# Compute EoS + Opacities from the standalone HOP
# ============================================================================

"""
    compute_eos_HOP(T, rho, lam; HOP_kwargs=Dict(:feh=>0.0), mmap=true, kwargs...)

Run the standalone HOP opacity engine over a (T, ρ) grid and return the raw EoS +
opacity arrays, in the exact layout `TSO.assemble_opacity` expects.

Returns the 9-tuple `(lnRho, lnT, lnEi, lnPg, lnNe, lnχ500, λ, κ, σ)`: log density
axis (nρ), log temperature axis (nT), log internal energy / gas pressure / electron
density / reference opacity (each nT×nρ), the wavelength axis (nλ), and the
absorption and scattering cubes `κ`/`σ` (nT×nρ×nλ, per-gram cm²/g).

Driven through the Julia port `<OpenMULTI3D>/tumult/HOP.jl` (loaded from the tumult
install set by `@import_tumult`). One HOP instance is created and the line list loaded
once. The (T, ρ) grid is walked in `t_block × rho_block` tiles **sequentially**; within
each tile the EOS is solved once and the **line opacity is threaded over the tile's T-ρ
cells** (`HOP.opacity_threaded`, `julia -t N`) — continuum computed once, then each
cell-chunk adds its lines on the one shared instance, then a serial combine. No line
list is copied and no EOS re-solved per thread, and there is no OpenMP in Fortran.
`nchunks` sets the cell-chunks per tile (default: the thread count); `nprocs` specifies
the number of worker processes to spawn for multiprocessing over wavelength chunks.

Requires a loaded tumult install (`@import_tumult "<OpenMULTI3D>"`). The (T, ρ) grid
is given by the 1D Julia arrays `T` [K] and `rho` [g/cm³]; `lam` is vacuum-Å.

The (nT, nρ) grid is walked **sequentially** in tiles of `t_block × rho_block` points
(so a large table streams through bounded RAM); each tile is one HOP `load_state` per
worker. `mmap` only decides where the opacity cubes live: with `mmap=true` (default)
each tile is streamed into a Julia-side `Mmap`ed file so the full `(nT, nρ, nλ)` cube
never sits in RAM; with `mmap=false` the cubes are plain in-RAM arrays.

This-function keyword arguments
- `mmap`       : opacity cubes as memory-mapped files (`true`) or in-RAM arrays (`false`)
- `t_block`, `rho_block` : sequential tile sizes in grid points (default: whole grid as one tile; shrink `rho_block` to stream a large table in bounded RAM)
- `energy_offset` : absolute-energy shift added to HOP's internal energy ([`HOP_ENERGY_OFFSET`])
- `Teltype`    : element type of the opacity/EoS arrays (default `Float32`)
- `verbose`    : log each tile as it finishes
- `hop`        : override the HOP module (default: the Julia port from the tumult install)

`HOP_kwargs` — `Dict{Symbol,Any}` forwarded to `HOP.HopEOS`:
- `:absdat_file`  : continuum absorption data (`input/TS_absdat.dat`)
- `:abund_file`   : abundance table (`input/abund/abund_magg22`)
- `:feh`          : metallicity [Fe/H] (`0.0`)
- `:isotope_file` : isotope data (`input/isotopes_blp`)
- `:molec_file`   : molecular data (`input/molec_data_BC.txt`)
- `:tmolim`       : molecular temperature limit [K] (`-99.0`, i.e. unset)
- `:alpha_enhancement`     : [α/Fe] enhancement (`0.0`)
- `:r_process_enhancement` : r-process enhancement (`0.0`)
- `:verbose`      : `false` = silent; `0` = line-count, `1` = progress, `2` = per-solve

Opacity keywords — forwarded to `HOP.opacity(eos, lam; ...)` (defaults in parentheses):
- `linelist`        : line-list file / list / folder (`nothing` → continuum only)
- `dlam`            : line-wing search radius [Å] (`100.0`)
- `include_absmet`  : add the sampled metal-line background (`false`)
- `line_scattering` : fold the (1−ε)-weighted line opacity into the scattering part (`true`)
- `vmic`            : microturbulence [km/s]
- `nthreads`        : cell-chunks the line opacity is threaded over per tile (`1`)

Run with `julia -t N` for threaded cell computation, and set `nprocs > 1` for
multiprocessing over wavelengths. `include_absmet=true` is not supported by the threaded path.
"""
function compute_eos_HOP(
	T, rho, lam;
	HOP_kwargs=Dict{Symbol,Any}(:feh => 0.0),
	mmap=true,
	linelist=nothing,
	dlam=100.0, include_absmet=false, line_scattering=true, vmic=nothing,
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
		dlam=dlam, include_absmet=include_absmet, line_scattering=line_scattering, vmic=vmic,
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