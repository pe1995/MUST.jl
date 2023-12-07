using Glob

path_dispatch = "../../../stellar_atmospheres"
folder_to_copy = glob("DIS_MARCS_E*/")

nml = joinpath.(folder_to_copy, "ininml.dat")
name = [f[end] == '/' ? f[1:end-1] : f for f in folder_to_copy]
name = [first(split(n, "DIS_MARCS_E_", keepempty=false)) for n in name]
name = [first(split(n, "_", keepempty=false)) for n in name]

for (i, n) in enumerate(name)
    cp.(nml[i], joinpath(path_dispatch, "grid_$(n).nml"), force=true)
end