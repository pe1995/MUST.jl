using Glob

path_dispatch = "../../../stellar_atmospheres"
folder_to_copy = glob(ARGS[2], ARGS[1])

nml = joinpath.(folder_to_copy, "ininml.dat")
name = [f[end] == '/' ? f[1:end-1] : f for f in folder_to_copy]
name = last.(split.(name, '/', keepempty=false))

for (i, n) in enumerate(name)
    @info "copy from $(nml[i]) to $(n).nml"
    cp.(nml[i], joinpath(path_dispatch, "$(n).nml"), force=true)

    if length(ARGS) > 2
        cp.(nml[i], joinpath(path_dispatch, ARGS[3], "$(n).nml"), force=true)
    end
end