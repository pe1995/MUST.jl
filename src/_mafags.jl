function read_mafags(file)
    if !isfile(file)
        error("File not found: $file")
    end

    f = FortranFile(file, "r")
    ndepth, nel, nmol, nbreak = read(f, (Int32, 4))
    rewind(f)

    _..., nion, idx_el = read(f, Int32, Int32, Int32, Int32, (Int32, nel), (Int32, nel))
    kion = sum(nion)
    rewind(f)

    _..., nscale, tau, t, kappa_ref, pe, pg, pk, rho = read(
        f, Int32, Int32, Int32, Int32, 
        (Int32, nel), (Int32, nel),   
        (Int32, nbreak),              
        (Float32, ndepth), (Float32, ndepth), (Float32, ndepth), 
        (Float32, ndepth), (Float32, ndepth), (Float32, ndepth), (Float32, ndepth)
    )       

    close(f)

    return (
        exp10.(Base.convert.(Float64, tau)) |> reverse,
        exp10.(Base.convert.(Float64, t)) |> reverse,
        exp10.(Base.convert.(Float64, kappa_ref)) |> reverse,
        exp10.(Base.convert.(Float64, pe)) |> reverse,
        exp10.(Base.convert.(Float64, pg)) |> reverse,
        exp10.(Base.convert.(Float64, pk)) |> reverse,
        exp10.(Base.convert.(Float64, rho)) |> reverse
    )
end

function mafagsBox(path::String)
    tau, t, kappa_ref, pe, pg, pk, rho = read_mafags(path)
    z = geometrical_from_optical(tau, kappa_ref ./ rho, rho)
    optical_surface!(z, tau)

    xx = zeros(1, 1, length(z))
    yy = zeros(1, 1, length(z))
    zz = reshape(z, 1, 1, :)

    d = reshape(rho, 1, 1, :) 
    T = reshape(t, 1, 1, :)    
    k5 = reshape(kappa_ref ./ rho, 1, 1, :)   
    t5 = reshape(tau, 1, 1, :)   
    pg = reshape(pg, 1, 1, :)    
    ne = reshape(pe, 1, 1, :) ./ (KBoltzmann .* T)

    data = Dict{Symbol, Any}(
        :Ne=>ne, :T=>T, :d=>d, :τ500=>t5, :Pg=>pg, :κ500=>k5
    )

    Box(xx, yy, zz, data, AtmosphericParameters())
end