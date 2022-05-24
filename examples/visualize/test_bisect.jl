using Pkg; Pkg.activate(".")
using MUST

MUST.@import_dispatch "../../../dispatch2" EOS

units = MUST.StaggerCGS()

leEos = MUST.@legacyPythonEOS MUST.dispatch_location
sqEos = MUST.@squaregasPythonEOS MUST.dispatch_location

le_ee_range = [5., 30.] 
sq_ee_range = [exp(sqEos.eos.scale1.min), exp(sqEos.eos.scale1.max)] 

@info "Legacy : Internal Energy (code) $(MUST.bisect(leEos; ee=le_ee_range, d=10^(-6.6383) / units.d, T=9298.04))"
@info "SquareG: Internal Energy (code) $(MUST.bisect(sqEos; ee=sq_ee_range, d=10^(-6.6383),           T=9298.04) ./units.ee)"
@info "SquareG: Internal Energy (code) $(MUST.bisect(sqEos; ee=sq_ee_range, d=10^(-9.1195),           T=4232.14) ./units.ee)"