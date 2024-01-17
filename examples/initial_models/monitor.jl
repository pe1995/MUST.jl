using MUST

@import_dispatch "../../../dispatch2"

# create a default watchdog for the given simulation
w = MUST.defaultWatchDog(ARGS[1])

# run the watchdog until it timeouts
MUST.monitor(w)