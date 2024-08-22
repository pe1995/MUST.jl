using Pkg; Pkg.activate(".")
using MUST
using TSO

marcs_path = ARGS[1]
new_dir = ARGS[2]

function simpleMarcs(path, new_dir)
	name = basename(path)
    new_path = joinpath(new_dir, name)
	
	modelmarcs = MUST.MARCSModel(path)
	model = TSO.Model1D(
		z=modelmarcs.structure["Depth"], 
		lnρ=log.(modelmarcs.structure["Density"]), 
		lnT=log.(modelmarcs.structure["T"]), 
		lnEi=fill!(similar(modelmarcs.structure["Depth"]), 0.0),
		logg=-99.0
	)

    TSO.flip!(model, depth=true)
    z = model.z
    T  = exp.(model.lnT)
    lnr = model.lnρ
    vmic  = parse(Float64, split(modelmarcs.info[4], " ", keepempty=false)[1])
    open(new_path, "w") do f
        write(f, name*"\n")
        write(f, "$(length(z))\n")
        for i in eachindex(z)
            write(f, "$(z[i]) $(T[i]) 1.0 $(exp(lnr[i])) $(vmic) \n")
        end
    end
	
    @info "saved at $(new_path)"
	model
end

marcs_paths = MUST.glob("*.mod", marcs_path)

for mp in marcs_paths
    simpleMarcs(mp, new_dir)
end