using Pkg; Pkg.activate(".")
using MUST

@import_dispatch "../../../dispatch2"

# create a default watchdog for the given simulation
w = MUST.defaultWatchDog(ARGS[1])

# should snapshots be saved
save_box = if length(ARGS) > 1
    ARGS[2] == "--save_box"
else
    false
end

# run the watchdog until it timeouts
MUST.monitor(w, save_box=save_box)