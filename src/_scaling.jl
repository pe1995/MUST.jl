struct AtmosUnits
    system ::String
    l      ::Float64
    d      ::Float64
    t      ::Float64
    n      ::Float64
    u      ::Float64
    m      ::Float64
    p      ::Float64
    pm     ::Float64
    ee     ::Float64
    e      ::Float64
    flux   ::Float64
    src    ::Float64
    qr     ::Float64
    k      ::Float64
    rk     ::Float64
    mu     ::Float64
    b      ::Float64
    k_B    ::Float64
    m_H    ::Float64
    m_He   ::Float64
    B      ::Float64
end

"""
StandardUnits units (dispatch code units). 
    code -> cgs: var * Stagger.cgs.var
    code <- cgs: var / Stagger.cgs.var
"""
function StandardUnits(;system = "cgs",                                   
                    l  = 1e8,                                        
                    d  = 1e-7,                                        
                    t  = 1e2,   
                    n  = 1/l^3,                                    
                    u  = l/t,                               
                    m  = d*l^3,                            
                    p  = d*u^2, 
                    pm = d*u,                          
                    ee = u^2,    
                    e  = d*ee,
                    flux = d*u^2/t*l,
                    src = m/(t^3),
                    qr = e/t,
                    k  = l^2/m,
                    rk = 1/l,
                    mu = 1.3,
                    b = u*sqrt(π*4.0*d),
                    k_B  = 1.380658E-16,
                    m_H  = 1.6726219E-24, 
                    m_He = 6.65e-24,
                    B = sqrt(4.0*π*d)*u)
    AtmosUnits(system,l,d,t,n,u,m,p,pm,ee,e,flux,src,qr,k,rk,mu,b,k_B,m_H,m_He,B)
end

"""
Read Units from Dispatch snapshot and save the conversion to CGS.
"""
function StandardUnits(snap::T) where {T<:Py}
    ps = pyconvert(Any, snap.params_list["scaling_params"])
    pnames = keys(ps)
    l_name = "l_cgs" in pnames ? "l_cgs" : "l"
    d_name = "d_cgs" in pnames ? "d_cgs" : "d"
    v_name = "v_cgs" in pnames ? "v_cgs" : "v"
    t_name = "t_cgs" in pnames ? "t_cgs" : "t"

    #pnames = keys(snap.params_list["scaling_nml"])
    #l_name = "l_cgs" in pnames ? "l_cgs" : "l"
    #d_name = "d_cgs" in pnames ? "d_cgs" : "d"
    #v_name = "v_cgs" in pnames ? "v_cgs" : "v"
    #t_name = "t_cgs" in pnames ? "t_cgs" : "t"

    l = ps[l_name]
    d = ps[d_name]

    if v_name in pnames
        v = ps[v_name]
        t = l / v
    end
    if t_name in pnames
        t = ps[t_name]
        v = l / t
    end

    StaggerCGS(l=l,d=d,t=t,u=v)
end

"""
    StandardUnits(name::String; folder=@in_dispatch("data"))

Read Units from Dispatch snapshot with name `name` and save the conversion to CGS.
"""
function StandardUnits(name::String, folder::String)
	rundir = joinpath(folder, name)
	StandardUnits(rundir)
end

StandardUnits(folder::String) = begin
    file = joinpath(folder, "params.nml")
    ps = FreeNamelist(file)

    scalingp, pnames, l_name, d_name, v_name, t_name = try 
        scalingp = nmlField(ps, Symbol("scaling_nml"))
        pnames = keys(scalingp) |> collect

        scalingp, pnames, "l", "d", "v", "t"
    catch
        @warn "SCALING_NML not found in params.nml!"
        scalingp = nmlField(ps, Symbol("scaling_params"))
        pnames = keys(scalingp) |> collect
    
        scalingp, pnames, "l_cgs", "d_cgs", "v_cgs", "t_cgs"
    end

    l = scalingp[l_name]
    d = scalingp[d_name]

	t, v = if (t_name in pnames) & (v_name in pnames) 
        t = scalingp[t_name]
        v = scalingp[v_name]
        t, v
    elseif t_name in pnames
        @warn "Only $(t_name) found."
        t = scalingp[t_name]
        v = l / t
        t, v
	elseif v_name in pnames
        @warn "Only $(v_name) found."
        v = scalingp[v_name]
        t = l / v
        t, v
	else
		error("Either t or v scaling needs to be provided")
    end

    StandardUnits(l=l, d=d, t=t, u=v)
end


DispatchCGS(snap::T) where {T<:Py} = StandardUnits(snap)
StaggerCGS(snap::T) where {T<:Py}  = StandardUnits(snap)
StaggerCGS(args...; kwargs...)     = StandardUnits(args...; kwargs...)

GeneralUnits = StandardUnits()