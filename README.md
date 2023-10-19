# `MUST.jl`
Julia package for creating 3D stellar atmosphere models using the DISPATCH framework.
This package is capable of executing simulations via a SLURM distriubution system. It furthermore handles the reading, conversion and post-processing/visualisation of those simulations. It prepares the input for the Multi3D radiative transfer code. This package is complementary to the `TSO.jl` package, which is responsible for opacity realted questions.
In the following the main functionality is presented. Examples for most use-cases can be found in the `examples` folder.

-------------------------
## Installation

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
## 3D Model Atmospheres

The main functionality of `MUST.jl` is to read different stellar atmospheres in different formats and convert them to one comman format that can be used for post-processing purposes. For this, there are two different model representations available, the `MUST.Space` and `MUST.Box`. The `MUST.Space` is designed as a loose array of points associated with different coordinates, that can have any shape or orientation. This is important if models with unspecified, ungirdded orientation are present. In most cases, this is not relevant unless there is mesh refinement active in Dispatch. In any case, for any post-processing purposes the model should always be converted in a `MUST.Box` object, which tabulates any data present as 3D arrays.

### Reading Dispatch Models

*Relevant examples: `examples/dispatch2space/convert.jl, prepare_restart.jl`*

For converting from Dispatch to `MUST.jl` the Dispatch python routines are used. For tnis, the location of the `dispatch2` installation has to be known.

```julia
using MUST
MUST.@import_dispatch "path/to/dispatch2/"
```
This macro will import the dispatch location into the right possition, such that the python modules can be found. Next, the list of snapshots that should be converted has to be loaded. 
For this the `Glob` module may be used.

```julia
using Glob
folder = MUST.@in_dispatch "data/yourmodelhere/"
content_of_folder = glob("*/", folder)
snapshots = sort(MUST.list_of_snapshots(content_of_folder))
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

This naming convention will make it possible to browse the folder for converted models on different scales later, such that we can pick specific snapshots using e.g.

```julia
pick_snapshot(out_folder, :recent)
```

which will return geometrical and optical depth snapshots, if available with above names. You can however pick them any way you like. You can load them later using a simple

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
### MULTI3D Models

*Relevant examples: `examples/dispatch2multi/box2multi.jl`*

A conversion from the `Box` model to the MULTI3D format can be done conveniently. For this only the path to the hdf5 file of the model needs to be given to the correspoding function. Note that an electron density is needed, however in MULTI3D `use_ne` can be set to `false`. In that case the elecctron density can be filled with any value, it won't be used. If the EoS is loaded in the `MUST.jl` format, it can be passed to the function to direclty compute `Ne`.

```julia
MUST.multiBox(box, eos, output_name, downsample_xy=downsample)
```
where `downsample_xy` can be given to pick only every nth point in the horizontal as a downsampling technique. In case you want to the more transparend `TSO.jl` interface, you can use it like this

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
To load up the model again you may use the same function.

```julia
MUST.multiBox(name)
```

which returns it again as a `Box` model.

### Reading MURaM Models
### Converting `MUST.Box` Models

