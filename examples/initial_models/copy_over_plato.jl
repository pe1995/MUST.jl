using Glob

path_dispatch = "../../../stellar_atmospheres"
folder_to_copy = glob("PLATO/M_E*/")

nml = joinpath.(folder_to_copy, "ininml.dat")
name = [f[end] == '/' ? f[1:end-1] : f for f in folder_to_copy]
name = last.(split.(name, '/', keepempty=false))

for (i, n) in enumerate(name)
    cp.(nml[i], joinpath(path_dispatch, "$(n).nml"), force=true)
end