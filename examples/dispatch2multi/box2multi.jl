using Pkg; Pkg.activate("."); 
using MUST
using PyPlot
using DelimitedFiles
using TSO

MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" ;

input_name  = ARGS[1]
eos_name    = ARGS[2] #MUST.@in_dispatch("input_data/$(ARGS[2])")
eos_dir_name = split(dirname(eos_name), "/") |> last |> String

#output_name = split(dirname(joinpath(input_name, "test")), '/') |> last 
output_name = split(last(split(input_name, "/")),".hdf5") |> first
output_name = "m$(eos_dir_name)_$(output_name)"
@show output_name

downsample  = length(ARGS) > 2 ? parse(typeof(1), ARGS[3]) : 5

model_box = MUST.Box(input_name);

@info "old size: $(size(model_box.x))"

#model_boxes = MUST.Boxes(input_name)
#model_box   = MUST.time_statistic(MUST.mean, model_boxes[1][end-5:end-1]);

#eos       = MUST.SquareGasEOS(eos_name);
#MUST.multiBox(model_box, eos, output_name, downsample_xy=15)

eos = reload(SqEoS, eos_name)
aos = @axed eos
@assert !TSO.is_internal_energy(aos)

ne  = lookup(aos, :lnNe, log.(model_box[:d]), log.(model_box[:T]))

model_box.data[:ne] = exp.(ne)
MUST.multiBox(model_box, output_name, downsample_xy=downsample)

@info "New size: $(size(@view(ne[1:downsample:end, 1:downsample:end, :])))"