"""
    geometrical_from_optical(τ, κ, ρ)

Compute geometrical height scale from optical depth.
"""
function geometrical_from_optical(τ, κ, ρ)
    z = zeros(length(τ))
    ρκ = ρ .* κ 
    for j in eachindex(z)
        if j==1 
            z[1] = 0 + (τ[2] - τ[1]) * 0.5 * (1.0/ρκ[j])
        else
            z[j] = z[j-1] + (τ[j] - τ[j-1]) * 0.5 * (1.0/ρκ[j] + 1.0/ρκ[j-1])
        end
    end

    return z
end


"""
    optical_surface(model)

Return the z coordinate of the optical surface
"""
function optical_surface!(z, τ)
    mask = sortperm(log10.(τ))
    o = linear_interpolation(
        Interpolations.deduplicate_knots!(log10.(τ[mask]), move_knots=true), 
        z[mask],
        extrapolation_bc=Line()
    )(0.0)

    z .= z .- o
end