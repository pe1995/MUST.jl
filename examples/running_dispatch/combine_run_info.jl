using Pkg; Pkg.activate(".")
using CSV
using DataFrames

df = DataFrame(CSV.File(ARGS[1]))

for f in ARGS[2:end]
    append!(df, DataFrame(CSV.File(f)))
end

CSV.write("summary_total.csv", df)