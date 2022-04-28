struct AtmosUnits
    system ::String
    l      ::Float64
    d      ::Float64
    t      ::Float64
    u      ::Float64
    m      ::Float64
    p      ::Float64
    pm     ::Float64
    ee     ::Float64
    e      ::Float64
    k      ::Float64
    rk     ::Float64
    mu     ::Float64
    b      ::Float64
    k_B    ::Float64
    m_H    ::Float64
    m_He   ::Float64
end

function StaggerCGS(system = "cgs",                                   
                    l  = 1e8,                                        
                    d  = 1e-7,                                        
                    t  = 1e2,                                       
                    u  = l/t,                               
                    m  = d*l^3,                            
                    p  = d*u^2, 
                    pm = d*u,                          
                    ee = u^2,    
                    e  = d*ee,
                    k  = l^2/m,
                    rk = 1/l,
                    mu = 1.3,
                    b = u*sqrt(π*4.0*d),
                    k_B  = 1.380658E-16,
                    m_H  = 1.6726219E-24,
                    m_He = 6.65e-24)
    AtmosUnits(system,l,d,t,u,m,p,pm,ee,e,k,rk,mu,b,k_B,m_H,m_He)
end

function StaggerCGS(snap::T) where {T<:PyCall.PyObject}
    pnames = keys(snap.params_list["scaling_params"])
    l_name = "l_cgs" in pnames ? "l_cgs" : "l"
    d_name = "d_cgs" in pnames ? "d_cgs" : "d"
    v_name = "v_cgs" in pnames ? "v_cgs" : "v"

    l = snap.params_list["scaling_params"][l_name]
    d = snap.params_list["scaling_params"][d_name]
    v = snap.params_list["scaling_params"][v_name]
    t = snap.scaling.t
    StaggerCGS("cgs",l,d,t,v)
end