### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ 17fc2f8e-9a0a-11ec-35dc-117e3cb2fd68
begin
	using Pkg
	Pkg.activate(".")
	using PyCall
	using Plots
end

# ╔═╡ 989addae-bfdc-4df2-8268-ba8523719381
md"
# Visualizing the adiabat production procedure
"

# ╔═╡ b6c7b72a-fbb8-4afa-b756-5938664ba1eb
md"## Imports"

# ╔═╡ 478cade0-1500-46f2-9276-c28c989815a6
begin
	os  = pyimport("os")
	sys = pyimport("sys")
	dispatch_path        = "/lustre/astro/eitner/model_grid/dispatch2-1/"
	dispatch_python_path = os."path".join(os.getenv("D", dispatch_path),
											"utilities/python/")
	
	
	dispatch_python_path in sys."path" ? nothing :
										append!(sys."path",[dispatch_python_path])
	dispatch = pyimport("dispatch")
	EOS      = pyimport("dispatch.EOS")
end;

# ╔═╡ cb84f04c-d899-40f4-a3b7-d86ad4fb4c43
md"## The adiabat procedure"

# ╔═╡ 37225809-dbc5-4d1a-8922-f1b3b397cbbd
md"Choose the EOS"

# ╔═╡ 667f88d3-c296-406b-bf71-da2da2533082
begin 
	eos = EOS.bifrost(top=dispatch_path,verbose=2)
	@show eos.scale0.name
	@show eos.scale1.name
end;

# ╔═╡ 52dc53a9-ad80-459e-9d15-2cda09d12607
begin
		eosStagger = EOS.stagger(top=dispatch_path, file="table.dat").E_to_lnE()
		@show eosStagger.scale0.name
		@show eosStagger.scale1.name
		@show eosStagger.variables
end

# ╔═╡ a01990fb-cc37-4022-92b7-fa787219528f
md"The atmosphere class"

# ╔═╡ d8ca950e-b0c1-4776-94e4-f42bd14738dc
begin
	abstract type AbstractAtmosphere
	end
	
	struct Atmosphere1D{T<:AbstractFloat} <: AbstractAtmosphere
		z ::Vector{T}
		d ::Vector{T}
		e ::Vector{T}
		t ::Vector{T}
		p ::Vector{T}
	end
	
	struct Atmosphere3D{T<:AbstractFloat} <: AbstractAtmosphere
	end
end

# ╔═╡ 764a1169-2681-4aed-8458-3dad8381c586
begin
	struct Scaling
		system  # pass on to Python
		l       # length = 1Mm = 1e8 cm
		d       # photospheric density
		u       # speed = 1 km/s
		t       # time = 1000 sec
		m       # mass
		p       # pressure = energy density
		e       # energy
		ee      # specific energy
		kr      # rosseland opacity [cm2/g]
		mu
		b       # magnetic flux density - this one 
				# is a factor of sqrt(4pi) smaller
		temp    # temperature for given mu
		grav    # constant of gravity
		m_u     
		m_e    
		m_h    
		m_he   
		amu    
		k_b                                                
		h_p    
		c      
		el      
		q      
		stefan 
		G    
		AU     
		pc     
		yr     
		kms    
		mum    
		Myr     
		m_sun   
		r_sun   
		m_earth 
		r_earth 
		mu_0    
		eps_0   
	end
	function Scaling(  ;system = "cgs",                               
						l  = 1e8,                                   
						d  = 1e-7,                                   
						u  = 1e5,                                  
						t  = l/u,                             
						m  = d*l^3,                       
						p  = d*u^2,                       
						e  = d*u^2,                        
						ee = u^2,                              
						kr = 1.0/(d*l),                    
						mu = 1.3,
						b = u*sqrt(π*4.0*d),
						temp = 1.0,
						grav = 6.6743e-8*(m/l^3)*t^2,
						m_u     = 1.6726219e-24,
						m_e     = 9.10938356e-28,
						m_h     = m_u,
						m_he    = 6.65e-24,
						amu     = 1.66054e-24,
						k_b     = 1.380658e-16,                                       
						h_p     = 6.6260755e-27,
						c       = 2.999792e10,
						el      = 4.80320451e-10,
						q       = 4.80320451e-10,
						stefan  = 5.67e-5,
						G       = 6.6743e-8,
						AU      = 1.496e13,
						pc      = 3.086e18,
						yr      = 3.15542e7,
						kms     = 1e5,
						mum     = 1e-4,
						Myr     = 3.15542e13,
						m_sun   = 1.989e33,
						r_sun   = 6.9598e10,
						m_earth = 5.972e27,
						r_earth = 6.371e8,
						mu_0    = 1.0,
						eps_0   = 1.0)
	
		Scaling(system, l, d, u, t, m, p, e, ee, kr, mu, b, temp, grav,
					m_u ,m_e, m_h, m_he, amu, k_b, h_p, c, el, q,      
					stefan, grav, AU, pc, yr, kms, mum, Myr, m_sun,  
					r_sun, m_earth, r_earth, mu_0, eps_0 )    
	end
	staggerScaling = Scaling(u=1e6);
end

# ╔═╡ 536f8500-c8bd-4880-871a-4ddfdf919615
staggerScaling.p

# ╔═╡ d72d9f50-6508-4574-935c-a47aa2341376
"""
Initialize the 1D atmosphere in agreement with the given adiabat and EOS.
"""
function init_adiabat(eos, t_ini=1e4, d_ini=3e-7, g_ini=2.75e4;
												  scaling=staggerScaling,
												  dlnd = 0.05,
												  ee_min = 4.0,
												  w_perturb = 0.1,
												  a_perturb = 0.1,
												  nz = 200,
												  i0 = 120,
												  n_iter = 3 )
	g0  = g_ini
    d0  = d_ini
	tt0 = t_ini
	p0  = d_ini / scaling.m_h * scaling.k_b * t_ini
	ee0 = 3.  
    ee1 = 30.
    for iter in 1:20
		ee2 = (ee0+ee1)/2.
		#@show  log(d0) ee2
		tt0 = eos.lookup("T", log(d0), log(ee2*scaling.ee))
		p0  = eos.lookup("P", log(d0), log(ee2*scaling.ee))
		#@show p0 tt0
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
	
	z[i0]  = 0.0
    d[i0]  = d0 / scaling.d
    ee[i0] = ee0 
    t[i0]  = eos.lookup("T", log(d[i0] * scaling.d), log(ee[i0]) + log(scaling.ee))
	p[i0]  = eos.lookup("P", log(d[i0] * scaling.d), log(ee[i0]) + log(scaling.ee))/
																						scaling.p

	test_x, test_y = zeros(2,n_iter+1),zeros(2,n_iter+1)
    dee = 0.0
    for i in i0-1:-1:1
		ee[i] = ee[i+1]+dee
		if i ==100
			test_x[2,1] = d[i]
			test_y[2,1] = ee[i]
		end
		if i == 99
			test_x[1,1] = d[i]
			test_y[1,1] = ee[i]
		end
		for iter in 1:n_iter
			#-------------------------------------------------------------------------
			# Look up new values for temperature and gas pressure, and do averages
			#-------------------------------------------------------------------------
			d[i] = d[i+1] *exp(dlnd)
			t[i] = eos.lookup("T", log(d[i])  + log(scaling.d), 
								   log(ee[i]) + log(scaling.ee))
			p[i] = eos.lookup("P", log(d[i])  + log(scaling.d), 
								   log(ee[i]) + log(scaling.ee)) / scaling.p
			#@show i p[i]
			#-------------------------------------------------------------------------
			# Use the step in pressure and the given step in ln(d) to compute dz and
			# the step in internal energy per unit mass, dee
			#-------------------------------------------------------------------------
			dp    = p[i] - p[i+1]
			da    = (d[i+1]+d[i]) / 2.
			dz    = dp / (g0*da)
			z[i]  = z[i+1] - dz
			pd    = (p[i+1]+p[i]) / (d[i+1]+d[i])
			dee   = pd*dlnd
			ee[i] = ee[i+1] + dee
			#@show iter p[i] ee[i] d[i]
			#@show scaling.p
			#iter ==3 && @show i,z[i],p[i],dee,ee[i],t[i]
			if i ==100
				test_x[2,iter+1] = d[i]
				test_y[2,iter+1] = pd
			end
			if i == 99
				test_x[1,iter+1] = d[i]
				test_y[1,iter+1] = pd
			end
		end 
	end 
	test_x,test_y
end

# ╔═╡ c7eaf2da-04ee-4ffd-9779-cb1a821d7904
x,y = init_adiabat(eos)

# ╔═╡ a1331084-71d7-42cd-9055-fe310b2c2c43
begin
	plot([1:3...],y[1,2:end]; marker=".")
	#plot!(x[2,:],y[2,:]; marker=".")
end

# ╔═╡ 83375c1e-953f-402d-8d1a-d39b5b100108
x[1,:]

# ╔═╡ Cell order:
# ╟─989addae-bfdc-4df2-8268-ba8523719381
# ╟─b6c7b72a-fbb8-4afa-b756-5938664ba1eb
# ╠═17fc2f8e-9a0a-11ec-35dc-117e3cb2fd68
# ╠═478cade0-1500-46f2-9276-c28c989815a6
# ╟─cb84f04c-d899-40f4-a3b7-d86ad4fb4c43
# ╟─37225809-dbc5-4d1a-8922-f1b3b397cbbd
# ╠═667f88d3-c296-406b-bf71-da2da2533082
# ╠═52dc53a9-ad80-459e-9d15-2cda09d12607
# ╟─a01990fb-cc37-4022-92b7-fa787219528f
# ╠═d8ca950e-b0c1-4776-94e4-f42bd14738dc
# ╠═764a1169-2681-4aed-8458-3dad8381c586
# ╠═536f8500-c8bd-4880-871a-4ddfdf919615
# ╠═d72d9f50-6508-4574-935c-a47aa2341376
# ╠═c7eaf2da-04ee-4ffd-9779-cb1a821d7904
# ╠═a1331084-71d7-42cd-9055-fe310b2c2c43
# ╠═83375c1e-953f-402d-8d1a-d39b5b100108
