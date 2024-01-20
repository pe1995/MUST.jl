"""
    pressurescaleheight(f, b::Box; kwargs...)

Reduce box to a plane given in kwargs. Then compute the pressure scale height in that plane. 
Apply function f.
"""
function pressurescaleheight(f, b::Box; Pg=:Pg, ρ=:d, logg=b.parameter.logg, kwargs...)
    add!(b, lnPg=log.(b[Pg]))
    add!(b, lnd=log.(b[ρ]))
    planePg = interpolate_to(b, :lnPg; kwargs...)[:lnPg]
    planeRho = interpolate_to(b, :lnd; kwargs...)[:lnd]

    f(exp.(planePg) ./ (exp.(planeRho) .* exp10.(logg)))
end

"""
    pressurescaleheight(b::Box; kwargs...)

Reduce box to a plane given in kwargs. Then compute the pressure scale height mean in that plane.
"""
pressurescaleheight(b::Box; kwargs...) = pressurescaleheight(mean, b; kwargs...)




"""
    convectiveturnovertime(b::Box; v=:uz, kwargs...)

Compute convective turnover time from mean pressure scale height and rms velocity.
"""
function convectiveturnovertime(b::Box; v=:uz, kwargs...)
    hp = pressurescaleheight(b; kwargs...)
    vel = interpolate_to(b, v; kwargs...)[v]

    rms(x) = sqrt(mean(x .^2))
    abs.(hp ./ rms(vel))
end
