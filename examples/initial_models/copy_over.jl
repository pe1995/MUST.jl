folder_to_copy = ARGS[1]
nml = joinpath(folder_to_copy, "ininml.dat")

name = folder_to_copy[end] == '/' ? folder_to_copy[1:end-1] : folder_to_copy

cp(nml, joinpath("/u/peitner/DISPATCH/stellar_atmospheres/", "$(name).nml"), force=true)