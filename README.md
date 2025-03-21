<a href="url"><img src="must_logo2.png" height=100% width=100%/></a> 
Julia package for creating 3D stellar atmosphere models using the DISPATCH framework.
This package is capable of executing simulations via a SLURM distriubution system. It furthermore handles the reading, conversion and post-processing/visualisation of those simulations. It prepares the input for the Multi3D radiative transfer code. This package is complementary to the `TSO.jl` package, which is responsible for opacity related questions.
In the following the main functionality is presented. Examples for most use-cases can be found in the `examples` folder.

-------------------------

# Shortcuts
You can jump through the documentation by using the follwing topic shortcuts. Note that less relevant sections of the documentation may not appear here.

1. [Installation](#installation)
2. [3D Model Atmospheres](#3d-model-atmospheres)
    1. [Reading Dispatch Models](#reading-dispatch-models)
    2. [MULTI3D Models](#multi3d-models)
    3. [MURaM Models](#muram-models)
    4. [Stagger Models](#stagger-models)
3. [Atmosphere Analysis](#atmosphere-analysis)
4. [Running Dispatch](#running-dispatch)
    1. [Running MULTI3D](#running-multi3d)
    2. [Running Stellar Atmospheres](#running-stellar-atmospheres)
4. [Workflow](#workflow) 

-------------------------
# Installation

To use `MUST.jl` a julia installation is requiered. Please avoid using versions older than `julia 1.6`. After downloading julia, the `MUST.jl` package can be added, after permissions to the repository have been granted, by simply adding it though the package manager.

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.2 (2023-07-05)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(@v1.9) pkg> add https://github.com/pe1995/MUST.jl.git
```

To get access please contact eitner@mpia.de. After adding the Package it should be available in your code as 

```julia
using MUST
```

If you want to use `MUST.jl` in python this is straight forward as well. We refer the reader to [juliacall.py and Pythoncall.jl](https://juliapy.github.io/PythonCall.jl/stable/juliacall/).
If you want to access one of the examples, there are local julia environments included for all of them. If you get an error message when using `MUST.jl`, you may try to re-adding it. So `cd` into the `examples` folder and run

```
$ julia --project
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.2 (2023-07-05)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

(initial_models) pkg> rm MUST
(initial_models) pkg> add ../..
```

Which should work in most cases. 

-------------------------
# 3D Model Atmospheres

The main functionality of `MUST.jl` is to read different stellar atmospheres in different formats and convert them to one comman format that can be used for post-processing purposes. For this, there are two different model representations available, the `MUST.Space` and `MUST.Box`. The `MUST.Space` is designed as a loose array of points associated with different coordinates, that can have any shape or orientation. This is important if models with unspecified, ungirdded orientation are present. In most cases, this is not relevant unless there is mesh refinement active in Dispatch. In any case, for any post-processing purposes the model should always be converted in a `MUST.Box` object, which tabulates any data present as 3D arrays.

A convenient interface for reading, converting and analysing is given in the `examples/initial_models/analyze.jl` notebook. It is a Pluto notebook and can be used by starting a `Pluto.jl` session.

```
(@v1.9) pkg> add Pluto
julia> using Pluto; Pluto.run()
```
which will open a Pluto notebook in your default browser. You can then open the `analyze.jl` notebook from here. If you already added Pluto, you can also just run it from the command line, e.g. 
```console
$ julia -e 'using Pluto; Pluto.run()'
```


## Reading Dispatch Models

*Relevant examples: `examples/dispatch2space/convert.jl, prepare_restart.jl, snapshot2box.jl`*

For converting from Dispatch to `MUST.jl` the Dispatch python routines are used. For tnis, the location of the `dispatch2` installation has to be known.

```julia
using MUST
MUST.@import_dispatch "path/to/dispatch2/"
```
This macro will import the dispatch location into the right possition, such that the python modules can be found. Next, the list of snapshots that should be converted has to be loaded. 

```julia
folder = MUST.@in_dispatch "data/yourmodelhere/"
snapshots = sort(MUST.list_of_snapshots(folder))
```

where `MUST.@in_dispatch` will append the correct experiment folder to your dispatch location, so that you can give it relative to the `stellar_atmospheres` experiment. The snapshots themselves can now be read using the python rourtine, and then assambled to a actual 3D object, to avoid intransparent patch-like data structures.

```julia
for i_s in eachindex(snapshots)
    @info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
    try
        # The dispatch snapshot object (Python)
        snap = dispatch.snapshot(snapshots[i_s], data=folder)

        # Units for conversion to CGS
        units = MUST.StandardUnits(snap)

        # Convert its content to pure Julia
        s = MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)

        # Apply the conversion to CGS
        MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                ux=:u, uy=:u, uz=:u,
                                x=:l, y=:l, z=:l, time=:t)
    catch
        @warn "snapshot $(snapshots[i_s]) could not be loaded."
    end
end
```

One can furthermore load the namelist that was used for creating the snapshot and load the EoS. This package has a EoS module identical to the one in Dispatch to be consistent, however in the `TSO.jl` module there is a more julia-like and flexible `EoS` module that is similar to the `Box` type an intermediate representation that fits into the package API. When the `TSO.jl` EoS is converted to be used in Dispatch, it can be loaded in `MUST.jl` in the same way, to allow debugging etc. 

```julia
# Initial namelist
nml_name = MUST.@in_dispatch splitpath(folder)[end]
nml = MUST.StellarNamelist(nml_name*".nml")

# Use the new Squaregas EOS 
eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
eos_sq = MUST.SquareGasEOS(MUST.@in_dispatch(eos_path))
```
With this information one can now add the Temperature and opacity to the data cube, because they are relevant later when analysing the output for example.

```julia
 # Add additional columns are already in CGS after converting
MUST.add_from_EOS!(s, eos, :T)
MUST.add_from_EOS!(s, eos, :kr)
```

This `Space` can now either be saved, if wanted, or should rather be converted to a `Box`. It will assume that the `uniform` function of x, y and z will return the proper 3D grid. If this is not the case, one may also choose a new axis the cube should be interpolated to. This is significantly slower thouch and can not be recommended. This is not optimized in any way because at the moment there is no mesh refinement in Dispatch 3D model atmospheres.

```julia
b_s = MUST.Box(s)

# or specify N points to interpolate to
b_s = MUST.Box(s, 200)

# or give axes directly
b_s = MUST.Box(s, x, y, z)

MUST.save(b_s; name="box_sn$(snapshots[i_s])", folder=folder)
```

If you want to pick snapshots later using the rest of the API, it is important to name snapshots like specified in the `name` kwarg above. This will allow the code later to specify between geometrical and $\rm \tau$. If you want to convert the cube to the optical depth scale (recommended for visualisation) you can do this simply by specifying what fields of the `Box` relate to opacity and desity.

```julia
τ = MUST.optical_depth(b_s, opacity=:kr, density=:d)
MUST.add!(b_s, τ, :τ_ross)

# Convert the height scale from cm to optical depth
b_τ = MUST.height_scale_fast(b_s, :τ_ross)
MUST.save(b_τ; name="box_tau_sn$(snapshots[i_s])", folder=folder)
```

A handy shortcut that contains the above procedure is provided as well. After loading the dispatch2 module through the usage of the macro mentioned above it is possible to use the following function, that incorporates above steps

```julia
using MUST

# Load the dispatch2 module
@import_dispatch "path/to/dispatch2"
folder = @in_dispatch "data/yourdispatchmodel"

# list of all snapshots in that folder (dispatch format)
snapshots = list_of_snapshots(folder)

# convert the given snapshot to a Box, also on optical depth scale
b, bτ = snapshotBox(snapshots[end-1], folder=folder)
```
Note that `snapshotBox` has the kwarg `is_box=true`, which assumes that the given snapshot already is box-like, meaning that it is gridded. If set to false (e.g. because of possible mesh refinement) the cube will first be converted to a `Space` object, and then be interpolated to a regular grid in a new `Box` object.
The naming convention of the saved snapshots will make it possible to browse the folder for converted models on different scales later, such that we can pick specific snapshots using e.g.

```julia
pick_snapshot(folder, :recent)
```

which will return geometrical and optical depth snapshots, if available with above names. `:recent` will return the last available snapshot, but you can pass any snapshot number that is converted. You can obtain a list of converted snapshots from

```julia
converted_snapshots(folder)
```
which will return a dictionary, containing the folder, and for every found snapshot number the corresponding box names, if they were found at all.
You can however pick them any way you like. You can load them later using a simple

```julia
b = MUST.Box(name, folder=nothing)
```

that will check if name is a file, if not it will join

```julia
path = isnothing(folder) ? @in_dispatch(name*".hdf5") : joinpath(folder, name*".hdf5")
```

It will load the data arrays as memory maps automatically. If you want to enforce the naming convention you can also try the automatic detection of the vertical scale using

```julia
MUST.save(b_τ, snapshots[i_s], folder=folder)
```

which will check if all z columns are identical and assume that the scale is cartesian if they are.
From this internal `Box` model we can convert to other formats, like e.g. the `MULTI3D` format.
Note that you can always ensure a consistent orientation (from bottom to top, negative z values at the bottom, increasing upwards) by using the `flip!` function on this `Box` model.

```julia
flip!(b_s, density=:d, uz=:uz)
```

where you need to specify what field density and vetical velocity are. The density will be used to detect the ordering of arrays and z-velocity will be flipped in sign if the axis z needs to be flipped.

## MULTI3D Models

*Relevant examples: `examples/dispatch2multi/box2multi.jl`*

A conversion from the `Box` model to the MULTI3D format can be done conveniently. For this only the path to the hdf5 file of the model needs to be given to the correspoding function. Note that an electron density is needed, however in MULTI3D `use_ne` can be set to `false`. In that case the elecctron density can be filled with any value, it won't be used. If the EoS is loaded in the `MUST.jl` format, it can be passed to the function to direclty compute `Ne`.

```julia
MUST.multiBox(box, eos, output_name, downsample_xy=downsample)
```
where `downsample_xy` can be given to pick only every nth point in the horizontal as a downsampling technique. In case you want to use the more transparend `TSO.jl` interface, you can use it like this

```julia
aos = @axed reload(SqEoS, eos_name)
ne = if !TSO.is_internal_energy(aos)
    lookup(aos, :lnNe, log.(model_box[:d]), log.(model_box[:T]))
else
    lookup(aos, :lnNe, log.(model_box[:d]), log.(model_box[:ee]))
end
model_box.data[:ne] = exp.(ne)
MUST.multiBox(model_box, output_name, downsample_xy=downsample)
```
where you can either use an EoS on the energy or temperature grid. 
To load up the model again you may use the same function

```julia
MUST.multiBox(name)
```

which returns it again as a `Box` model.

There is a convenient collection which uses these models, that is included in this repository. However, because it uses `TSO.jl`, it is not part of the module. It can be loaded however, once you installed `TSO.jl` to your local environment.

```julia
must2multi = MUST.ingredients("convert2multi.jl")
```

which returns the script, which is located within the `src` folder. To use it one may call e.g.

```julia
must2multi.snaps2multi(
	models..., 
	eos=eos, 
	label=["$(l)_$(xres)x$(xres)x$(zres)" for (i, l) in enumerate(labels)],
	name="muram_m2",
	n_horizontal=xres,
	n_vertical=zres,
	outfolder="nameoffolder/"
)
```

Where the EoS needs to be a `TSO.SqEoS`. `Models` is a list of `MUST.Box` models, `labels` are the corresponding labels. Specifying a resolution here will lead to a resampling by linear interpolation, unless an other `method` is specified. See `gresample` function for more info. You can leave out the EoS, in which case it is assumed that the electron density is available as `:ne`. If it is not, it will be set to 1.0. Make sure to use `use_ne=.false.` in M3D in that case especially.

## MURaM Models

*Relevant examples: `examples/Muram/convert2m3dis.jl, convert_muram.jl`*

MURaM models are saved in a convenient format, which is why the conversion is straight forward. 

```julia
# Load a MURaM model
MUST.MURaMBox(model_path)

# interpolate to optical depth scale
MUST.MURaMBox(model_path, opacity)
```

where you can get the opacity e.g. from an `TSO.SqOpacity` table.

## Stagger Models

*Relevant examples: `examples/stagger2box.ipynb`*

Functions to read and convert Stagger snapshots are available. For this reader `msh, dat and aux` files are required within the same folder with the same filename (except the ending). Loading and converting is then very convenient.

```julia
s = MUST.StaggerSnap(filename, folder)
```

If additionally and EoS in the `MUST.jl` format is given, the gas pressure and rosseland opacity are recomputed from this EoS. Otherwise the values from the Stagger snapshot are used.

```julia
b = MUST.Box(s, eos=sqEOS)
```

Models in the `MUST.Box` format can then be converted to MULTI3D format as seen above. Please note that in M3D@DISPATCH the format needs to be specified as `"must"`, because the velocity in the vertical is not flipped in these models here.

----------------

# Atmosphere Analysis

*Relevant examples: `examples/Paper_I/paper_v2.jl, examples/initial_models/analyse.jl`*

The module contains usefull functionality that allow easy and transparent analysis of `MUST.Box` models. The most common usecase is the investigation of the average stratification of various quantities. In general, you can access all data cubes stored within this `Box` by using

```julia
# x, y, z meshgrids
b.x, b.y, b.z

# data arrays are stored in the data dictionary, e.g. temperature
b.data[:T]

# or access it with indexing
b[:T]

# all fields can be shown with
@show fieldnames(typeof(b))

# all variables in data can be shown with 
@show collect(keys(b.data))
```

## Plane Statistics

Plane statistics are usefull to apply a given statistic plane wise to the entire data cube. The most convenient interface is the `plane_statistic` function, that is designed for exactly this usecase

```julia
# mean temperature of the snapshot box
T = plane_statistic(mean, box, :T)

# mean log density 
lnRho = log.(plane_statistic(mean, box, :d))

# plane-wise std of temperature
stdT = plane_statistic(std, box, :T)
```

The corresponding axis, most of the time the z-axis, can be obtained as 1D array.

```julia
z = axis(box, :z)
```

If the axis is not z, but e.g. $\rm \tau$, one may either use 

```julia
# get the 3rd dimension of τ_ross array
axis(box, :τ_ross, 3)
```

or simply use the shortcut

```julia
τ, T = profile(mean, box, :log10τ_ross, :T)
```

where log10 or log can be added in front of any variable to apply the corresponding log directly. This is convenient when you want to past this in the plotting function directly, for example

```julia
plot(profile(mean, box, :log10τ_ross, :T)...)
```

will plot the average temperature as a function of optical depth. Note that this will only be correct, if the given `Box` is indeed on the $\rm \tau$ scale. Computing the average quantites on the geometrical scale and plot against average optical depth is **not** equal to average planes in a cube interpolated to a uniform optical depth scale. One should always compute averages on planes of constant optical depth if one intents to compare on that scale.

## Surfaces

You can either get surfaces which are closest to a given axis value, or interpolate column wise.

```julia
# closest plane to optical surface
MUST.closest(log10.(axis(snap, :τ_ross, 3)), 0)

# interpolate to optical surface
fh = MUST.height_where(τ_ross=1.0)
plane = MUST.reduce_by_column(fh, box)
```

where `reduce_by_column` will create a new Box with size 1 in z-direction. `fh` is a function that takes in a rosseland opacity column and returns the interpolated column as a function of z, such that a call to `fh` with all fields of `box` as kwargs will return a function that can be evaluated at any point in optical depth and will return then corresponding height. This can be used to interpolate the cube to this height. An easier interpolation can be achived by using the `interpolate_to` function

```julia
plane = MUST.interpolate_to(box, :T, τ_ross=0, logspace=true)
```
which is more convenient and will also interpolate by column.

## Time Statistics

*Relevant examples: `examples/initial_models/progress.jl, monitor.jl`*

Time dependent statistics can be computed on the fly by `MUST.WatchDog` capabilities. Such a watchdog can be run in parallel to the running DISPATCH simulation. It will detect new snapshots and perform the statistics presented above directly. Per default it will not save the converted cubes to save disk space, also because many of the intermediate snapshots are not relevant for the end-product and will be ignored anyways. Saving their cubes is hence mostly irrelevant. You can create a `WatchDog` by passing the name of the simulation it should track and the functions it should evaluate on the cubes.

```julia
using MUST
@import_dispatch "path/to/dispatch2"

# default watchdog contains many usefull statistics
w = MUST.defaultWatchDog("grid_t5777g44m00", folder=@in_dispatch("data/"))

# you can also add functions with the following call structure e.g.
mystatistic(watchdog, box, box_t) = begin
    # e.g. average profile, can be aything
    # 'box' is the geometrical, 'box_t' the optical depth 3D cube.
    x, y = MUST.profile(MUST.mean, box, :z, :T)

    # return your result as Dict, so that it will have proper names when loaded
    Dict(
        "x" => x,
	"y" => y
    )
end

w = MUST.defaultWatchDog("grid_t5777g44m00", folder=@in_dispatch("data/"), mystatistic=mystatistic)

# run the monitoring, stops when there is now new snapshot after timeout seconds
MUST.monitor(w; timeout=2*60*60)
```

The default watchdog can also be started by e.g.

```
$ julia --project monitor.jl grid_t5777g44m00
```
or send to the background with `nohup`. You can look at the results by calling `monitoring = MUST.reload(MUST.WatchDog, "grid_t5777g44m00", folder=datafolder)`, or by simply using the `examples/initial_models/progress.jl` notebook, which provides a nice interface.

--------------

# Running Dispatch

It is furthermore possible to run Dispatch (MULTI3D or stellar atmosperes) directly through `MUST.jl`. This includes automatic creation of namelists, submitting the job and fetching the output. This functionality can even be executed within a SLURM allocation, and hence provides easy access to the functional and repeated execution of tasks when considering grids. 

Note that the grid computation functionality of `MUST.jl` is outdated, and needs to be updated according to the new philosophy, as can be seen in `examples/initial_models`. The desired interface should be more similar to the one developed for MULTI3D. 

## Running MULTI3D

*Relevant examples: `examples/Paper_I/running_multi/opacity_tables.jl, effective_temperature.jl`*

Running M3D@DISPATCH from within julia is straight forward. If you have a working, compiled version of MULTI3D available, you can either use the existing high-level functions to automatically create namelists, or run whatever namelist you can think of. The execution task can either be piped to slurm, or just as a normal background task. You may either wait for the execution to finish or submit it to the background and fetch it later. The following gives an example for running M3D for a couple of different snapshots of a MURaM model in full 3D.

```julia
using MUST

# location of the dispatch2 M3D installation
@import_m3dis "path/to/Multi3D"

# folder where snapshots are saved
modelatmosfolder = "./input_multi3d/muram_m2"

# names of the snapshots within this folder
snapshots = [
	"muram_m2_HDm_50x50",
	"muram_m2_SSDm_50x50",
	"muram_m2_HDl_50x50",
	"muram_m2_SSDl_50x50",
	"muram_m2_HDh_50x50",
	"muram_m2_SSDh_50x50"
]

# linelist or absmet file path
linelist = "./input_multi3d/vald_2490-25540.list"
absmet = "./input_multi3d/absmet"

# input parameters that you want to specify (optional)
input_parameters = Dict(
    :model_folder=>modelatmosfolder,
    :linelist=>linelist,
    :absmet=>nothing,
    :atom_params=>(
    	:atom_file=>"./input_multi3d/atoms/atom.li12_col", 
    	:use_atom_abnd=>false,
    	:abundance=>2.1,
    	:exclude_trace_cont=>true,
    	:exclude_from_line_list=>true,
    	:hydrogen_BPO=>false
    ),
    :spectrum_params=>(
    	:daa=>0.01, :aa_blue=>6706, :aa_red=>6710
    ),
    :composition_params=>(
    	:absdat_file=>"./input_multi3d/TS_absdat.dat",
    	:abund_file=>"./input_multi3d/abund_asplund07"
    ),
    :atmos_params=>(
    	:atmos_format=>"MUST", 
    	:use_density=>true, 
    	:use_ne=>false,
    	:FeH=>-2.0,
    	:dims=>16
    ),
    :m3d_params=>(
    	:n_nu=>1, 
    	:quad_scheme=>"set_a2"
    )
)
```

Now from this general setup you can produce a input dict with different run-names and variations of `input_parameters`, and run those in parallel. The run-name will be appended to the output name (name of the atmosphere), such that every model is run for every input setup.

```julia
# modify input parameters as you want
params = Dict("test1" => input_parameters1, "test2" => input_parameters2)

# run all of it in parallel using slurm
MUST.spectrum(
    snapshots, 	    # run these snapshots
    params, 	    # and those paramters
    NLTE=true, 	    # in NLTE
    slurm=true,     # using Slurm (+ wait)
    twostep=false   # first dep. => LTE + dep.
)
```

Note that you can if course also run M3D with individual snapshots and only one setup. For this you only specify the snapshot(s) and `namelist_kwargs` in agreement with the M3D input parameters. Please read the M3D documentation for more info. There are a couple of different namelist templates and higher level functions available (e.g. `spectrum(), whole_spectrum()`, etc.). Take a look at `src/_multi.jl` and `src/namelist.jl`.
In case all of this is to high-level for you, take a look at the direct execution functions `srun_m3dis()` and `run_m3dis()` for slurm, and non-slurm, respectively. These just handle the execution of pre-created namelists. So if you have a different way of creating your namelists you can just call those functions to run them.

The result will be saved in the normal M3D format. There is a thin wrapper available, which uses the excellent `m3dis.py` package to read the output into julia. Note that, if not done automatically, you may need to call `MUST.pyconvert()` onto the result to bring it to julia types.

```julia
# read results stored in resultspaths
m3druns = MUST.M3DISRun.(resultpaths) 

# handle the output for e.g. different abundances
abund = MUST.getabundances.(m3druns)         # get abundances
MUST.equivalentwidth(m3druns, line=1)        # interpolate EW
MUST.abundance(m3druns, line=1)              # interpolate abundance
MUST.ΔNLTE(m3druns, line=1, reference=:LTE)  # compute NLTE corrections

# or execute directly for individual lines
lines = [run.line[10] for run in m3druns]    # Base.getproperty is redirected to m3dis.py
MUST.abundance(lines, abund)                 # interpolate abundance for line 10
```

The `m3dis.run` object is stored in the `m3drun.run` field and can hence be used as a normal python object if needed. `Base.getproperty()`, so e.g. `m3drun.flux` should automatically access the corresponding field of the python object.


## Running Stellar Atmospheres
*Relevant examples: `examples/initial_models/prepare.jl, random_models.jl, exoplanets.jl, run.jl`*

Running MUST@DISPATCH is similar to running M3D@DISPATCH. Namelists can either be run via a slurm job submission or by submission to the background. If slurm submissions are executed within a slurm allocation, they will be executed as job steps automatically. `MUST.jl` also offers the possibility to create namelists automatically. 

You can either load an existing namelist and then modify the fields, or just set the fields manually (Note that currently the namelist type is not very flexible, because the name of the fields are hardcoded. If there are additional fields needed, you need to modify the struct StellarNamelist. This is a bit inconvenient and may need change in the future). The same interface is used in the MULTI case under the hood.

```julia
# load an empty namelist
nml = MUST.StellarNamelist()

# or load an example namelist
nml = MUST.StellarNamelist("stellar_default.nml")

# Set whatever fields you want
MUST.set!(
    nml, 
    patch_params=(:n=>[patch_size, patch_size, patch_size],),
    scaling_params=(:l_cgs=>l_cgs, :d_cgs=>d_cgs, :t_cgs=>tscale)
)

# Save the namelist in dispatch
MUST.write(nml, MUST.@in_dispatch(name))
```

There is a script `prepare.jl` in the mentioned example folder, that is made for creating a large number of namelists for a grid of models. There are multiple steps involved, like the creation of initial 1D models and the corresponding binning of opacites, which involve the `TSO.jl` package. This is why the corresponding source file is available as `ingredient`.

```julia
prepare4dispatch = MUST.ingredients("prepare4dispatch.jl")
```

In order to run the script `prepare.jl`, which will create everything you need, including opactiy tables, you need to have a `MUST.StaggerGrid` available, which contains information about the Stagger grid which will be used to determine initial conditions. There are other initial conditions one can think of, but at the current time this is the preferred one. The steps are as follows.

1. Create the `MUST.StaggerGrid`. For this you need the Stagger grid, which is at the moment not public. A default `stagger_grid.mgrid` is available, which consists of a couple of models from the grid. But you can enhance this table with whatever model you can think of. Have a look at the file, models don't need to be Stagger at all. You only need the lists quantities (like e.g. size of the box, the average model on the geometrical scale, etc.).
2. This grid can be used to create random initial conditions. At the moment, this is done by interpolating in the grid in terms of every quantity, **including** the average model. This may be replaced by a adiabatic initial condition. Also the resolution of the models is interpolated in the grid. This means that the more models you have in the initial grid the better the interpolation will be. Note that the interpolation between models and the construction of a new z scale has to reply on opacities! This means that it is highly recommended to include the EoS already at this stage, when interpolating and constructing the initial condition. It is not a requirement though, you can also perform a point-by-point interpolation by omitting the EoS. However this is strongly discouraged and may lead to unstable intial conditions that are far from being adiabatic. You can provide an EoS for every new model, in case you want them to have different compositions.
3. Also: please make sure to include a column with `matching_eos` in the original grid! This should contain the path the a suitable EoS for each model in the grid, because the rosseland opacity needs to be consistent. If there is a column names `avo_path`, then the model averaged on the optical depth scale stored at the given path are used. The rosseland opacity is recomputed regardless.
```julia
grid = MUST.StaggerGrid("stagger_grid.mgrid")
modelgrids = MUST.ingredients("modelgrids.jl")

# add a column with a suitable EoS for every grid model
grid.info[!, "matching_eos"] = [
	joinpath(eos_root, "ross_combined_eos_magg_m0_a0.hdf5") 
	for _ in 1:nrow(grid.info)
]

# new grid object with the info about the models + average model
ig = modelgrids.interpolate_from_grid(
    grid, 
    teff=[5000.0, 6000.0], 
    logg=[4.0, 4.5], 
    feh=[0.0, 0.0],
    eos=[eos1, eos2]
)

# or add the models to the existing grid
ig = modelgrids.interpolate_from_grid!(...)
```
3. Run `prepare.jl` with this grid. You don't need this particular script, but following the steps in the script is encouraged. This will create a `dispatch_grid.mgrid`, input namelists and binned opacities. You can see input files for this script in the `input` folder. Those can be copied and modified to your needs. 
5. Run the final grid using the `run!()` function.
```julia
@import_dispatch "path/to/dispatch2"

# load the prepared dispatch grid
grid = MUST.StaggerGrid("dispatch_grid.mgrid")

# run it. Threads, mem and timeout are passed to each job step
MUST.run!(grid; threads=40, memMB=90000, timeout="24:00:00", slurm=true)
```

The success of the results will be recorded upon end in the success column, prepended with the name of the grid. Success will be judged based on the availability of all snapshots that should be there. One could also use the return status of the application -- which is available -- however sometimes dispatch may either crash without a false success status or end with a error status, even though it reached the end. You can switch between the two by specifing the `use_status` kwarg. 
Within the grid running procedure, the name of the namelist is used to run dispatch the same way as MULTI.

```julia
# run in slurm
srun_dispatch(nml; threads=40, memMB=1024, timeout="24:00:00")

# or run in shell
run_dispatch(nml; threads=70, wait=true, ddir=@in_dispatch(""))
```

where the first line can be run within a sbatch allocation. You can therefore scrip previous steps, if you already have a namelist, initial model (or adiabat) and opacity table ready. The results of all computations, successfull or not, can be retreived after converting your favourite snapshot, as described [in previous sections](#reading-dispatch-models).

___________

# Workflow

There are many useful tools and programs that help navigating through the usecases of `MUST.jl`. All of these tools are located within the `stellar_atmospheres` folder of `DISPATCH`. To illustrate the main workflow, consider the following example.
Suppose you observed 3 FGK-type stars and took high-resolution optical spectra. Your goal is to analyse those spectra with respect to elemental abundances and to make sure you rely on the physically most realistic models possible you want to perform the spectrum synthesis in 3D, possibly even Non-LTE. Suppose the stars you observed are located on the main-sequence and sub-giant branch and you have estimates of their stellar parameters. To simulate the model atmosphere in 3D you need to generete an initial model to start your simulation from, create the binned opacity table for this model, estimate the size of the simulation domain, and start `DISPATCH`. This first step can be done using the `interpolated_initial.jl` program. 
To make this work, you need to make sure you have monochromatic opacity tables available. Assuming you do not have access to them, you can create them yourself using the `Multi3D` code. The corresponding script is `create_eos.jl`. It can be used like e.g.

```bash
julia --threads=24 create_eos.jl -n table_name  --feh=0.0  --alpha=0.0 --vmic=1 --version=v1.0 --n_lambda=250000 --multi_threads=32 --linelist_dir=input_multi3d/LINE-LISTS/
```

Which will use `Multi3D` to create a table with 250000 wavelength points from 1000Å to 100000Å. Have a look at `--help` to gain more information about the capabilities, like e.g. a specific chemical composition with `--composition`. The new table will be given the name "table_name". In `interpolated_initial.jl` a grid of average models is used to interpolate the new model. The interpolator needs the equation of state, so you need to add your new table to the grid. This can be done with `replace_opacities_in_grid.jl`.  

```bash
julia replace_eos_in_grid.jl --grid=Stagger --name my_new_grid --opacity_tables abs_path/of/all/your/opacity/tables --metallicity_assignment=0=table1,m1=table2...,m5=table6
```

Which will create a new grid file with the opacity tables you specify. Checkout the new file if you want to further replace EoS tables by hand. The new grid `my_new_grid` can now be used for interpolation of your new model. Note that the previous steps only need to be once. As soon as the new grid file has been created it can be used to interpolate in. For this we use the script `interpolated_initial.jl`.

```bash
julia --threads=10 interpolated_initial.jl --teff=5500,6800 --logg=3.7,4.6 --feh=0.0,0.0 --random=10 --adiabatic_extrapolation --tau_bottom=4 --tau_up=-6.0 --tau_down=7.0 --modelFolder=InterpolatedModels/M1 --grid=my_new_grid --scale_resolution=1.0 --patch_size=20 --n_patches=10,10,5 --staggered
```

With this command, we create 10 random initial models for DISPATCH between 5500 and 6000K, 3.7 and 4.6 log(g), with solar metallicity. The physical domain of the box will cover the log optical depth from -6 to 7, and will extrapolate tbe model adiabatically below optical depth of 4, to ensure the bottom part of the model atmosphere is indeed adiabatic. The model will be run with 10x10x5 patches, with 20 points each, such that the frinal resolution will be 200x200x100 in this example. The option `--staggered` makes sure that DISPATCH uses the staggered, conservative radiative transfer solver (strongly recommended) that gerneally performs extraordinary well at low to intermediate resolution. Again see `--help` for more options. The code above also performs the opacity binning (8 bins by default) and created the input namelist for DISPATCH. From here one all that remains is to run DISPATCH either with SLURM or interactively, e.g.

```bash
./dispatch.x your_interpolated_model.nml
```

Note that the watchdog capabilities are very handy, so it is a good idea to run it at the same time. You can either take inspiration from one of the available `.slurm` scripts, or simply write your own script that makes sure resources are spent on it. It can be run as e.g.

```bash
nohup julia monitor.jl your_interpolated_model -k 10 --reverse > your_interpolated_model.nml.monitor &
```

So it will compute statistics of your run while it is running in the background. `-k 10` makes sure that only the last 10 snapshots are saved, the rest is deleted. This is very useful, because snapshots use a lot of disc space while you generally are not interested in all of them but just the last few. Of course the snapshots that get deleted are still included in the monitoring before they are removed.

After your DISPATCH run is finished, it may be a wise idea to collect your output together in one location, so that you can use it any time and don't lose what you have computed. For this you can use the script `ship_models.jl`, which will convert snapshots to `hdf5` datacubes, additionally interpolate them to the optical depth scale, and store them together in one place, including the namelist and EoS tables used to create the model in the first place.

```bash
julia ship_models.jl -r your_interpolated_model --move --n_snaps=10 -n 'shipped' --out_dir='shipped_models'
```

Which will create a new folder in `shipped_models` and store the last 10 snapshots of your model there. It will also add the τ500 opacity to your datacube, if a corresponding table is available in the folder of the EoS table. Note that in the final folder there will be three cubes per snapshot `isnap`, the `box_snisnap.hdf5` on the geometrical scale, `box_tau_snisnap.hdf5` on the rosseland optical depth scale, and `box_tau500_snisnap.hdf5` on the 500nm scale. The latter box can then simply be read back by using `pick_snapshot(folder, isnap, tau_name="tau500")` for . See `--help` for more options. You can e.g. run the shipping for entire grids easily by creating a folder with all the namelists of the models you want to ship, and specify that folder with `-d`. 

Now that your model has been secured, you can compute spectra. There is the very convenient and flexible script `spectrum.jl`, which does the work for you! All you need to do is run e.g. for the CH-G Band:

```bash
julia spectrum.jl -r your_interpolated_model --datafolder=shipped_models --add_m3dis=GBand --remove -n -1 -x 50 -z 180 -s 5880.0 -e 5900.0 -d 0.001 --feh=0.0 --alpha=0.0 --composition= C_3.0,N_0.4,O_0.5 --linelists=input_multi3d/LINE-LISTS/
```

To compute the CH G Band for the snapshot `-1`, which is the most recent one, at a horizontal resolution of `50` points, vertical of `180` points, from `5880 Å` to `5900 Å`, with a spacing of `0.001 Å`, and composition you specified. There are many more options available, so please check out `--help` again. The script will run M3D and store the output in a new data cube, containing the downsampled atmosphere as well as the spectrum. It will be stored as `box_m3dis_snsnapid.hdf5` cube, and can be read back simply by using `pick_snapshot(folder, isnap, box_name="box_m3dis")`.


