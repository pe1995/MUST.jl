using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using MUST

phase = ARGS[1]

MUST.@import_dispatch "../../../dispatch2/"
summary = DataFrame(CSV.File("summary.csv"))
success = rand(Bool, nrow(summary))

for row in 1:nrow(summary)
    nml       = MUST.StellarNamelist(MUST.@in_dispatch(summary.phase1_name[row]))
    name      = summary[row,Symbol("phase$(phase)_name")]
    name      = String(name[1:first(findlast(".", name))-1])
    output_at = MUST.@in_dispatch "data/$(name)"

    success[row] = MUST.check_success(nml, output_at)
    break
end

#summary["phase$(phase)_success"] = success
#CSV.write("summary.csv", DataFrame(summary)) 
