"""
    Quantity

A quantity to read from a DISPATCH patch and include in the `Box` object.
A Quantity always needs a name, a recipe and a conversion.
By default, the conversion can be omitted, which will result in the default conversion, i.e.
a division by the unit with the same `name` (See `MUST.StandardUnits`). Optionally,
a other name can be given as `MUST.standardConversion(name)`. If a combination of 
different units is needed, you can specify your own function as e.g.
```julia
myConversion() = begin
    (units) -> getfield(units, :l) / getfield(units, :t)
end
```
which will then be used as a unit when passed to the respective Quantity.
A `recipe` is a similar thing, which tells the code how to get the quantity in the 
first place. If to recipe is given, the quantity is assumed to be a field of the
patch data anyway. For velocities this e.g. would look like:
```julia
uxQ = Quantity(:ux, (; d, px, kwargs...)->(px ./ d), standardConversion(:u))
```
Where you have to add `kwargs...` because this function will be called with 
all data that is available in the patch.
Collect all quantities you wish (except EoS) in an array and you are good to go.
"""
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
uxQ = Quantity(:ux, (; d, px, kwargs...)->(px ./ d), standardConversion(:u))
uyQ = Quantity(:uy, (; d, py, kwargs...)->(py ./ d), standardConversion(:u))
uzQ = Quantity(:uz, (; d, pz, kwargs...)->(pz ./ d), standardConversion(:u))
eeQ = Quantity(:ee, (; d, e, kwargs...)->(e ./ d), standardConversion(:ee))
dtQ = Quantity(:dt_rt, nothing, standardConversion(:t))
fluxQ = Quantity(:flux, nothing, standardConversion(:flux))
heatingQ = Quantity(:qr, nothing, standardConversion(:qr))
opacityQ = Quantity(:kapparho, nothing, standardConversion(:rk))

defaultQuantities = [
    dQ,
    eQ,
    uxQ,
    uyQ,
    uzQ,
    eeQ,
    dtQ,
    fluxQ,
    heatingQ,
    opacityQ
]
