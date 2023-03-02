using Pkg; Pkg.activate("."); 
using MUST
using PyPlot
using DelimitedFiles
using TSO

MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" ;

input_name  = ARGS[1]
output_name = split(last(split(input_name, "/")),".hdf5") |> first
output_name = "m$(output_name)"
eos_name    = ARGS[2] #MUST.@in_dispatch("input_data/$(ARGS[2])")

model_box = MUST.Box(input_name);

#eos       = MUST.SquareGasEOS(eos_name);
#MUST.multiBox(model_box, eos, output_name, downsample_xy=15)

eos = reload(SqEoS, eos_name)
aos = @axed eos
@assert !TSO.is_internal_energy(aos)

ne  = lookup(aos, :lnNe, log.(model_box[:d]), log.(model_box[:T]))

model_box.data[:ne] = exp.(ne)
MUST.multiBox(model_box, output_name, downsample_xy=5)