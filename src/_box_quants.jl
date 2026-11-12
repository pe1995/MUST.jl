struct Quantity
    name ::Symbol
    recipe
    conversion
    derived ::Bool
end

standardConversion(name) = begin
    (units) -> getfield(units, name)
end

Quantity(name) = Quantity(name, nothing, standardConversion(name), false)
Quantity(name, recipe) = Quantity(name, recipe, standardConversion(name), true)
Quantity(name, recipe, conversion) = Quantity(name, recipe, conversion, true)

derived(q) = [qi for qi in q if qi.derived]
nonderived(q) = [qi for qi in q if !qi.derived]
varnames(q) = [qi.name for qi in q]

dQ = Quantity(:d)
eQ = Quantity(:e)
BxQ = Quantity(:bx, nothing, standardConversion(:B), false)
ByQ = Quantity(:by, nothing, standardConversion(:B), false)
BzQ = Quantity(:bz, nothing, standardConversion(:B), false)
uxQ = Quantity(:ux, (; d, px, kwargs...)->(px ./ d), standardConversion(:u))
uyQ = Quantity(:uy, (; d, py, kwargs...)->(py ./ d), standardConversion(:u))
uzQ = Quantity(:uz, (; d, pz, kwargs...)->(pz ./ d), standardConversion(:u))
eeQ = Quantity(:ee, (; d, e, kwargs...)->(e ./ d), standardConversion(:ee))
dtQ = Quantity(:dt_rt, nothing, standardConversion(:t))
fluxQ = Quantity(:flux, nothing, standardConversion(:flux))
gammaQ = Quantity(:gamma, nothing, standardConversion(:gamma))
heatingQ = Quantity(:qr, nothing, standardConversion(:qr))
opacityQ = Quantity(:kapparho, nothing, standardConversion(:rk))
srcQ = Quantity(:sourcefunction, nothing, standardConversion(:src))

defaultQuantities = [
    dQ,
    eQ,
    BxQ,
    ByQ,
    BzQ,
    uxQ,
    uyQ,
    uzQ,
    eeQ,
    dtQ,
    fluxQ,
    heatingQ,
    opacityQ,
    srcQ,
    gammaQ
]
