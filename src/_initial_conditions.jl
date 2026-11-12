# ============================================================================
# Initial conditons (Legacy, should be unused!)
# ============================================================================

struct InitialConditions
    T
    ρ
    g
    hp
end

# ============================================================================
# Initial models
# ============================================================================

struct MarcsInitialModel{T<:AbstractFloat} <:AbstractInitialModel
    Teff ::T
    logg ::T
    μ    ::T
    lgTauR  ::Vector{T}
    lgTau5  ::Vector{T}
    Depth   ::Vector{T}
    T       ::Vector{T}
    Pe      ::Vector{T}
    Pg      ::Vector{T}
    Prad    ::Vector{T}
    Pturb   ::Vector{T}
end

struct SimpleInitialModel{T<:AbstractFloat} <:AbstractInitialModel
    Teff ::T
    logg ::T
    z   ::Vector{T}
    T   ::Vector{T}
    lnρ ::Vector{T}
end

MarcsInitialModel(path; Teff, logg, μ=mean_molecular_weight()) = begin
    model = readdlm(path, skipstart=1)
    MarcsInitialModel(Teff, logg, μ, (model[:, i] for i in axes(model, 2) if i>1)...)
end

SimpleInitialModel(box::Box; Teff, logg) = SimpleInitialModel(Base.convert(eltype(box.z), Teff),
                                                              Base.convert(eltype(box.z), logg),
                                                              collect(axis(box, :z)),
                                                                 plane_statistic(mean, box, :T),
                                                            log.(plane_statistic(mean, box, :d)))

SimpleInitialModel(path::String; Teff, logg) = begin
    model = readdlm(path)
    SimpleInitialModel(Teff, logg, (model[:, i] for i in axes(model, 2))...)
end

SimpleInitialModel(model::MarcsInitialModel) = begin
    ρ = eltype(model.T)[ρg_ideal(model.T[i], model.Pg[i], model.μ) for i in eachindex(model.T)]
    SimpleInitialModel(model.Teff, model.logg, -reverse(model.Depth), reverse(model.T), log.(reverse(ρ)))
end

# ============================================================================
# General Methods
# ============================================================================

mean_molecular_weight(X=0.73, Y=0.26, Z=0.01) = 1.0 /(2.0*X + 3.0/4.0*Y + 1.0/2.0*Z)
Hp(T, g, μ=mean_molecular_weight())           = GeneralUnits.k_B * T / (μ *GeneralUnits.m_H *g)
ρg_ideal(T, Pg, μ=mean_molecular_weight())    = Pg / (GeneralUnits.k_B *T) *μ *GeneralUnits.m_H
granular_size(args...; kwargs...)             = 1.3 * Hp(args...; kwargs...)
box_dimensions(args...;kwargs...)             = Dict(:x=>granular_size(args...; kwargs...)*4 /1e8, :z=>Hp(args...; kwargs...)*14 /1e8)

function cut(model::SimpleInitialModel; kwargs...)
    z, T, ρ = deepcopy(model.z), deepcopy(model.T), deepcopy(model.lnρ)
    mask = trues(length(z))
    for p in keys(kwargs)
        lo, hi = kwargs[p]
        mask .= (mask) .& ((getfield(model, p) .> lo) .& (getfield(model, p) .< hi))
    end

    SimpleInitialModel(model.Teff, model.logg, z[mask], T[mask], ρ[mask])
end

# ============================================================================
# Methods for initial conditions
# ============================================================================

function initial_adiabat(eos, t_ini=1e4, d_ini=3e-7, g_ini=2.75e4;
												  scaling=MUST.StaggerCGS(),
												  dlnd = 0.05,
												  ee_min = 4.0,
												  w_perturb = 0.1,
												  a_perturb = 0.1,
												  nz = 200,
												  i0 = 120,
                                                  z0_position=0.0,
												  n_iter = 30 )
	g0  = g_ini * scaling.t^2/scaling.l
    d0  = d_ini
	tt0 = t_ini
	p0  = d_ini / 1.6726219e-24 * 1.380658e-16 * t_ini
    _, _, ee0, ee1 = limits(eos) ./ scaling.ee
    for iter in 1:20
		ee2 = (ee0+ee1)/2.
		tt0 = MUST.lookup(eos, :T,  d0, ee2*scaling.ee)
		p0  = MUST.lookup(eos, :Pg, d0, ee2*scaling.ee)
		if (tt0 > t_ini)
			ee1 = ee2
		else
			ee0 = ee2
		end
    end

	z  = zeros(nz)
	d  = zeros(nz)
	e  = zeros(nz)
	ee = zeros(nz)
	t  = zeros(nz)
	p  = zeros(nz)

	z[i0]  = z0_position
    d[i0]  = d0 / scaling.d
    ee[i0] = ee0
    t[i0]  = MUST.lookup(eos, :T, d[i0] * scaling.d, ee[i0] * scaling.ee)
	p[i0]  = MUST.lookup(eos, :P, d[i0] * scaling.d, ee[i0] * scaling.ee) / scaling.p

    dee = 0.0
    for i in i0-1:-1:1
		ee[i] = ee[i+1]+dee
		for iter in 1:n_iter
			# ============================================================================
			# Look up new values for temperature and gas pressure, and do averages
			# ============================================================================
			d[i] = d[i+1] *exp(dlnd)
			t[i] = MUST.lookup(eos, :T, d[i]*scaling.d, ee[i] * scaling.ee)
			p[i] = MUST.lookup(eos, :P, d[i]*scaling.d, ee[i] * scaling.ee) / scaling.p
			# ============================================================================
			# Use the step in pressure and the given step in ln(d) to compute dz and
			# ============================================================================
			dp    = p[i] - p[i+1]
			da    = (d[i+1]+d[i]) / 2.
			dz    = dp / (g0*da)
			z[i]  = z[i+1] - dz
			pd    = (p[i+1]+p[i]) / (d[i+1]+d[i])
			dee   = pd*dlnd
			ee[i] = ee[i+1] + dee
		end
	end

	# ============================================================================
	# Integrate from the starting point upwards
	# ============================================================================
    dee = 0.0
    for i in i0+1:nz
      ee[i] = max(ee[i-1]-dee,ee_min)
      for iter in 1:n_iter
        # ============================================================================
        # Look up new values for temperature and gas pressure, and do averages
        # ============================================================================
        d[i] = d[i-1]*exp(-dlnd)
        t[i] = MUST.lookup(eos, :T, d[i]*scaling.d, ee[i] * scaling.ee)
		p[i] = MUST.lookup(eos, :P, d[i]*scaling.d, ee[i] * scaling.ee) / scaling.p
        # ============================================================================
        # Use the step in pressure and the given step in ln(d) to compute dz and
        # ============================================================================
        dp = p[i-1]-p[i]
        da = (d[i-1]+d[i])/2.
        dz = dp/(g0*da)
        z[i] = z[i-1]+dz
        pd = (p[i-1]+p[i])/(d[i-1]+d[i])
        dee = pd*dlnd
        ee[i] = max(ee[i-1]-dee,ee_min)
      end
    end

    z,d*scaling.d,ee,t,p
end

function initial_conditions(model::MarcsInitialModel; kwargs...)
    @assert length(kwargs)==1

    ref_field, ref_val = keys(kwargs)[1], values(kwargs)[1]
    reference_point = argmin(abs.(getfield(model, ref_field) .- ref_val))

    ref_T = model.T[reference_point]

    ref_ρ = ρg_ideal(model.T[reference_point], model.Pg[reference_point], model.μ)

    hp = Hp(model.Teff, exp10(model.logg), model.μ)

    InitialConditions(ref_T, ref_ρ, exp10(model.logg), hp)
end

geometrical_depth(a_goal, a_ini, z) = begin
    smask = sortperm(a_goal)
    linear_interpolation(a_goal[smask], z[smask], extrapolation_bc=Line())(a_ini)
end

match(T_ini, ρ_ini, logg, z_goal, T_goal, d_goal, eos; ee_min=3.5, nz=400, i0=150, hotspot=z_goal[end-1]/1e8, which=:T) = begin
    z0         = which==:T ? geometrical_depth(T_goal, T_ini, z_goal) / 1e8 : geometrical_depth(d_goal, ρ_ini, z_goal) / 1e8
    z,d,ee,t,p = initial_adiabat(eos, T_ini, ρ_ini, exp(logg), ee_min=ee_min, nz=nz, i0=i0, z0_position=z0);
    z_last     = hotspot
    z_low      = argmin(abs.(z .- (z_last)))
    z0        += which==:T ? geometrical_depth(T_goal, t[z_low], z_goal) / 1e8 - (z_last) :
                             geometrical_depth(d_goal, d[z_low], z_goal) / 1e8 - (z_last) ;

    z0
end

function regrid(model::SimpleInitialModel, z_new::AbstractVector{T}) where {T}
    iT = linear_interpolation(model.z, model.T,   extrapolation_bc=Line())
    iR = linear_interpolation(model.z, model.lnρ, extrapolation_bc=Line())

    z_scale = Base.convert.(eltype(model.z), z_new)

    SimpleInitialModel(model.Teff, model.logg, z_scale, iT.(z_scale), iR.(z_scale))
end

function for_dispatch(model::SimpleInitialModel, path)
    open(path, "w") do f
        writedlm(f, hcat(model.z, model.T, model.lnρ))
    end;
end
