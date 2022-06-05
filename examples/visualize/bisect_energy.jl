using Pkg; Pkg.activate(".");# Pkg.update()
using MUST

MUST.@import_dispatch "../../../dispatch2" EOS

units = MUST.StaggerCGS()

#leEos = MUST.@legacyPythonEOS MUST.dispatch_location
sqEos = MUST.SquareGasEOS(MUST.@in_dispatch("input_data/square_gas_1/"))

#le_ee_range = [5., 30.] 
sq_ee_range = [Float32(sqEos.params["EiMin"]), Float32(sqEos.params["EiMax"])] 

# MARCS model, ini: Tau=1, cool: Tau=-5
t_ini  = 9298.04
d_ini  = 10^(-6.6383)
t_cool = 4232.14
d_cool = 10^(-9.1195)


# MARCS model 2, ini: Tau=1, cool: Tau=0, min: Tau=-5
t_ini  = 8906.303711
d_ini  = 10^(-6.538473)
t_cool = 6422.082031
d_cool = 10^(-6.530671)
t_min  = 4064.781006
d_min  = 10^(-9.007814)

# MARCS model 2, ini: Tau=1, cool: Tau=-3, min: Tau=-5
t_ini  = 8906.303711
d_ini  = 10^(-6.538473)
t_cool = 4731.717285
d_cool = 10^(-7.940815)
t_min  = 4064.781006
d_min  = 10^(-9.007814)

# MARCS model 2, ini: Tau=1.5, cool: Tau=-6, min: Tau=-6
t_ini  = 9574.467773
d_ini  = 10^(-6.516477)
t_cool = 3875.567139
d_cool = 10^(-9.598085)
t_min  = 3875.567139
d_min  = 10^(-9.598085)


@info "SquareG: Internal Energy (code) at initial point ($(t_ini),$(d_ini)): $(MUST.bisect(  sqEos; ee=sq_ee_range, d=d_ini,  T=t_ini) ./units.ee)"
@info "SquareG: Internal Energy (code) at cooling point ($(t_cool),$(d_cool)): $(MUST.bisect(sqEos; ee=sq_ee_range, d=d_cool, T=t_cool) ./units.ee)"
@info "SquareG: Internal Energy (code) at minimum point ($(t_min),$(d_min)): $(MUST.bisect(sqEos; ee=sq_ee_range, d=d_min, T=t_min) ./units.ee)"