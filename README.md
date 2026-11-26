<img src="must_logo2.png" height=100% width=100%/> 
Julia package for reading, converting and analysing 3D stellar atmosphere models computed with `DISPATCH` (M3DIS). It converts the native `DISPATCH` output to data cubes (`MUST.Box`), computes statistics while a simulation is running, and prepares models for the MULTI3D radiative transfer code. It is complementary to the `TSO.jl` package, which handles the EoS and opacity tables.

The complete pipeline (tables, binning, initial conditions, running, shipping and spectra) is run through `m3dis.jl` in the `stellar_atmospheres` experiment; see the M3DIS README for that. This page gives an overview of what you can do with `MUST.jl` directly.

-------------------------

# Shortcuts

1. [Installation](#installation)
2. [Boxes](#boxes)
    1. [File Layout](#file-layout)
    2. [Finding Snapshots](#finding-snapshots)
    3. [Converting DISPATCH Snapshots](#converting-dispatch-snapshots)
    4. [Restarting from a Box](#restarting-from-a-box)
3. [Atmosphere Analysis](#atmosphere-analysis)
4. [Monitoring a Run](#monitoring-a-run)
5. [MULTI3D Models](#multi3d-models)
6. [Spectra](#spectra)
    1. [Reading Spectra](#reading-spectra)
    2. [Averaging and Combining](#averaging-and-combining)
    3. [Abundance Corrections](#abundance-corrections)

-------------------------
# Installation

`MUST.jl` is included in the M3DIS repository and installed together with the `stellar_atmospheres` environment. It can also be added on its own through the package manager

```
pkg> add https://github.com/pe1995/MUST.jl.git
```

To get access please contact eitner@mpia.de. After adding the package it is available as

```julia
using MUST
```

-------------------------
# Boxes

All models are represented as a `MUST.Box`, which holds the data as 3D arrays:

```julia
# x, y, z meshgrids
b.x, b.y, b.z

# data arrays are stored in the data dictionary, e.g. temperature
b.data[:T]

# or access it with indexing (this copies the array)
b[:T]

# all variables in data
keys(b.data)

# time, Teff, logg and composition of the snapshot
b.parameter
```

## File Layout

Every snapshot is stored in one `HDF5` file that can hold the same snapshot on several scales, or anything else that belongs to it. Each of them is a group of the file with its own `x`, `y`, `z` and data arrays, so a group can also be on a different grid (e.g. the τ scale); the spectra of the snapshot are groups of their own:

```
box_sn198.hdf5
├── geo/       geometrical height scale
├── tau/       Rosseland optical depth scale
├── tau500/    500nm optical depth scale (added when shipping)
├── spectra/   one group per spectrum, spectra/<name> (added by spectrum.jl)
└── restart/   namelists needed to restart DISPATCH from this snapshot
```

The top level holds what is the same for all groups: the snapshot number and run name as attributes, and the time, Teff, logg and chemical composition. These are the same for every group, so each group uses them. The composition is the group `composition` with one attribute per entry, `[Fe/H]`, `[alpha/Fe]` and `[X/Fe]` of individual elements, taken from the opacity table the model was computed with (`b.parameter.composition`). Files that store it as the datasets `composition_e` and `composition_v` are read as well.

Loading a file only needs its path. By default you get the geometrical cube (or the only one in the file):

```julia
b = MUST.Box("path/to/box_sn198.hdf5")

availableData(b)                 # [:geo, :tau, ...]
b_τ = MUST.Box(path; data=:tau)  # or: switch(b; data=:tau)
switch!(b; data=:tau)            # in place: b now holds the optical depth cube
```

`switch` returns a new `Box` and leaves the original untouched, while `switch!` replaces the data of `b` itself. The data arrays are loaded as memory maps automatically (`mmap=false` to read them into memory).

Information that is not data, like the snapshot number, the run, the scale and the file the box was loaded from, is in the `meta` dictionary of the box:

```julia
MUST.meta(b)               # Dict(:snapshot => 198, :kind => "geo", :path => ..., ...)
MUST.meta(b, :snapshot)    # 198
```

You can save a box to a new file, add any box as a new group, or write a modified box back to the group it was loaded from:

```julia
path = MUST.save(b; name="box_sn198", folder=folder)  # new file, group geo
MUST.save!(path; tau=b_τ, mytest=b_test)              # replaces groups of that name
MUST.save(b; mode="r+")                               # write b back where it came from
```

## Finding Snapshots

To browse a folder for converted snapshots use

```julia
b, b_τ = pick_snapshot(folder, :recent)
```

which returns the geometrical and the optical depth cube of a snapshot (the latter is `nothing` if it is not there). `:recent` returns the last available snapshot, but you can pass any snapshot number that is converted, and e.g. `tau_name="tau500"` to get the 500nm scale instead. A list of converted snapshots is given by

```julia
converted_snapshots(folder)
```

which returns a dictionary containing the folder, and for every snapshot number the corresponding boxes, as `BoxRef`s that can be loaded with `Box(ref)`. `MUST.snapshot_path(folder, 198)` returns the file of a given snapshot.

## Converting DISPATCH Snapshots

Snapshots in the native `DISPATCH` format are converted with `snapshotBox`. The location of the `DISPATCH` installation has to be known first, so that the run folders can be found:

```julia
using MUST
@import_dispatch "path/to/M3DIS"
folder = @in_dispatch "data/yourmodel"

# list of all snapshots in that folder (dispatch format)
snapshots = list_of_snapshots(folder)

# convert a snapshot to a Box, also on the optical depth scale
b, b_τ = snapshotBox(snapshots[end]; folder=folder)
```

With `save_snapshot=true` (the default), both cubes are stored in `box_sn<number>.hdf5` in the folder of the run, together with everything that is needed to restart `DISPATCH` from this snapshot later. In practice this is done by the watchdog (`--save_box`) or the `ship` mode of `m3dis.jl`, which also use the binned `TSO.jl` EoS of the run for the conversion.

## Restarting from a Box

The cubes contain everything that is needed to restart `DISPATCH`, so the native `DISPATCH` output does not have to be kept. To restart, the watchdog can write a restart folder from a single box and exit:

```bash
julia monitor.jl your_model --snapshot 198 --restart
```

This creates `<run>/m3dis_198.restart` (converting the snapshot first, if needed) and prints the `restart_params` to add to the namelist, e.g.

```
&restart_params data='/path/to/data/your_model/', run='m3dis_198.restart', snapshot=198 /
```

From within Julia the same can be done with `MUST.write_restart(MUST.snapshot_path(folder, 198), "path/to/m3dis_198.restart")`.

-------------------------
# Atmosphere Analysis

The most common use case is the average stratification of a cube. `plane_statistic` applies a statistic to every horizontal plane:

```julia
# mean temperature of the snapshot box
T = plane_statistic(mean, b, :T)

# plane-wise std of temperature
stdT = plane_statistic(std, b, :T)

# the corresponding axis
z = axis(b, :z)
```

The shortcut `profile` returns the axis and the statistic together, where `log10` or `log` can be put in front of any variable:

```julia
τ, T = profile(mean, b_τ, :log10τ_ross, :T)
plot(τ, T)
```

Note that this is only correct if the `Box` is indeed on the scale you ask for. Averages on the geometrical scale plotted against the average optical depth are **not** the same as averages on planes of constant optical depth, so use the `tau` group (see [File Layout](#file-layout)) for profiles in optical depth.

To get a horizontal slice at a given optical depth, interpolate every column:

```julia
plane = MUST.interpolate_to(b, :T, τ_ross=1, logspace=true)
```

-------------------------
# Monitoring a Run

The watchdog follows a running `DISPATCH` simulation, converts every new snapshot and computes statistics (averages, RMS, Teff, ...) from it, which are stored in the `monitoring` folder of the run. It is started next to the simulation with

```bash
nohup julia monitor.jl your_model -k 10 --save_box > your_model.nml.monitor &
```

where `-k 10` keeps only the last 10 native snapshots and `--save_box` keeps the converted cubes. See `--help` for more options, e.g. the automatic rebinning of the opacities during the run. The results can be loaded with

```julia
w = MUST.WatchDog("your_model", folder=@in_dispatch("data/"))
monitoring = MUST.reload!(w)

# e.g. the effective temperature over time
teff = MUST.timeevolution(monitoring, "atmosphericParameters", "teff")
```

The `analysis.jl` Pluto notebook in the `stellar_atmospheres/m3dis` folder gives a convenient interface to all of this, including the converted 3D cubes.

-------------------------
# MULTI3D Models

A `Box` can be converted to the MULTI3D format. MULTI3D needs the electron density, which can be added from the EoS with `TSO.jl`:

```julia
using TSO

eos = TSO.EquationOfState("path/to/table"; parts=(:eos,))
b.data[:ne] = exp.(TSO.sample(eos, :lnNe; rho=b.data[:d], T=b.data[:T]))

MUST.multiBox(b, "output_name", downsample_xy=2)
```

where `downsample_xy` picks only every nth point in the horizontal. If `:ne` is not available, set `use_ne=.false.` in MULTI3D. The model can be read back as a `Box` with `MUST.multiBox("output_name")`. Spectra for converted snapshots are computed with `m3dis.jl spectrum` (see the M3DIS README), which stores each of them as `spectra/<name>` of the snapshot file: `spectra(path)` lists them, `Spectrum(path, name)` loads one, with its wavelength, flux, intensities, composition and setup.

You can also run OpenMULTI3D for a converted model yourself from Julia:

```julia
MUST.@import_tumult "path/to/OpenMULTI3D"
run = MUST.spectrum("m3dis_198"; 
    format="multi", 
    name="my_spectrum", 
    presets=("3d",),
    composition=b.parameter.composition,
    settings=Dict(
        "spectrum_params"=>Dict("aa_blue"=>5160.0, "aa_red"=>5170.0, "daa"=>0.01)
    )
)
s = Spectrum(run, "my_spectrum"; composition=b.parameter.composition)
```

`settings` can hold any entries of the OpenMULTI3D namelist.

-------------------------
# Spectra

Spectra are computed with OpenMULTI3D, usually through `m3dis.jl spectrum` (see the M3DIS README), or with `MUST.spectrum` as above. Each spectrum is stored under its name as the group `spectra/<name>` of the model's `HDF5` file (for 3D models the file of the snapshot), so one model can hold the spectra of any number of setups. A spectrum (`MUST.Spectrum`) contains

- the wavelengths `s.λ` [Å, air], the ray directions `s.μ` and their integration weights `s.weights`,
- the flux and the intensities averaged over the surface, with their continua, and for NLTE runs also the LTE spectrum of the same run,
- optionally the horizontally resolved intensities of 3D models,
- the composition it was computed with, `s.composition` ([Fe/H], [α/Fe] and [X/Fe]), and its setup, `s.meta` (model, snapshot, Teff, logg, window, NLTE, atom).

## Reading Spectra

```julia
path = "my_run/box_sn481.hdf5"

spectra(path)                       # names of the spectra in the file
s = Spectrum(path, "GridMn1")       # one spectrum; resolved=true also loads the 3D intensities
ss = load_spectra(path)             # all of them

λ, F = MUST.flux(s)                 # flux normalized by the continuum
λ, F = MUST.flux(s; norm=false)     # absolute flux [erg/s/cm²/Å]
λ, F = MUST.flux(s; lte=true)       # LTE spectrum of an NLTE run
μ_z(s)                              # the μ of the rays
λ, I = intensity(s, 1.0)            # intensity at μ = 1 (norm and lte as for the flux)
s[:resolved_intensity]              # λ × x × y × rays, if loaded with resolved=true
```

Spectra written by older versions (`box_m3dis_snN.hdf5` files and `<tag>_<quantity>` entries) are read the same way.

## Averaging and Combining

The name of a spectrum ties together the runs of the same setup. Compute it with the same name for several snapshots and average it over the snapshots:

```julia
snapshot_spectra("my_run")                        # Dict(snapshot => spectra of that snapshot)
s = average_spectra("my_run", "GridMn1")          # mean over all snapshots with that spectrum
s = average_spectra("my_run", "GridMn1"; snapshots=470:481)
s = MUST.average([s1, s2, s3]; name="mean")       # or any spectra with the same wavelengths, rays and composition
```

The average is a spectrum like any other:

```julia
MUST.save(s, "averages.hdf5")                     # stored as spectra/GridMn1 of that file
Spectrum("averages.hdf5", "GridMn1")
MUST.save_text(s, "GridMn1.txt")                  # flux and intensity per μ as text; norm=true, lte=true as above
```

## Abundance Corrections

Spectra with different compositions are compared through the composition they carry, never through their names. Compute the same window for a few values of an element and a reference spectrum, and get the abundance correction from the curve of growth:

```julia
grid = [average_spectra("my_run", "GridMn$(i)") for i in 1:3]
ref = Spectrum("my_1d_model/my_1d_model.hdf5", "Mn1D")

MUST.equivalentwidth(ref; λs=6012.0, λe=6014.0)   # equivalent width [Å] of the normalized flux
ab_ref, ew_ref, ab, ew, correction = curve_of_growth(:Mn, ref, grid; λs=6012.0, λe=6014.0)
λ, F = MUST.interpolate_spectra(grid, :Mn => -0.1) # normalized flux interpolated to [Mn/Fe] = -0.1
```

`curve_of_growth` returns the [Mn/Fe] and equivalent width of the reference and of the grid, and the change of [Mn/Fe] the grid needs to reproduce the equivalent width of the reference. `lte=true` uses the LTE spectra of NLTE runs.

-------------------------

*Claude Opus 5.5 has been used to improve text quality.*
