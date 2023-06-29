### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 49a4e636-166c-11ee-03ed-c129e8cdc7fc
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	import MUST
	using Plots
end

# ╔═╡ bb45e6a0-2b8b-4731-bd53-190c44a5ea23
md"# Prepare & Run M3DIS"

# ╔═╡ 21951d47-3aef-4297-9edd-cf7afcb4544f
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 4bb87ab7-94fa-460d-9c7b-4b7fe06c8500
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 90818cad-34c7-4ad0-8a8a-530391a3e176
md"## The Namelists
Similar to the automatic execution of the stellar atmospheres code, the Multi code can be executed within a slurm environment. The overall code is the same, so first a namelist has to be created. One can choose a native default namelist, or just construct it from scratch. All parameters that are needed are available as fields. Note that this strictly speaking is not needed, but can easier to fit into the rest of the existing code."

# ╔═╡ 2da3ff54-8b7a-4e7e-8b15-c56e34594875
begin
	mutable struct M3DNamelist <: MUST.AbstractNamelist
		io_params          ::Dict{String,Any}  
		timer_params       ::Dict{String,Any}  
		atmos_params       ::Dict{String,Any}  
		stagger_params     ::Dict{String,Any}  
		atom_params        ::Dict{String,Any}  
		m3d_params         ::Dict{String,Any}  
		linelist_params    ::Dict{String,Any}  
		spectrum_params    ::Dict{String,Any}  
		dispatcher0_params ::Dict{String,Any}  
		task_list_params   ::Dict{String,Any}   
	end

	M3DNamelist() = M3DNamelist(
		[Dict{String,Any}() for f in fieldnames(M3DNamelist)]...
	)

	M3DNamelist(path::String) = begin
		s = M3DNamelist()
		read!(s, path)
		s
	end
end

# ╔═╡ b7d1105d-3a57-4425-83dd-20fa03cd4194
md"A couple of defaults however are usefull when it comes to compute the effective temperature of a model."

# ╔═╡ c490885a-588f-4277-be47-8a567babda27
whole_spectrum_namelist!(nml; 
	io_params=(:datadir=>"data", :gb_step=>10.0, :do_trace=>false),
	timer_params=(:sec_per_report=>120),
	atmos_params=(dims=>12, atmos_format=>"MUST"),
	m3d_params=(:verbose=>1, :fcheck=>1, :pcheck=>[1,1,1], :linecheck=>5, 
				:lvlcheck=>2,
                :n_nu=>100, :maxiter=>0, :decouple_continuum=>true,
                :long_scheme=>"lobatto", :quad_scheme=>"set_a2"),
	linelist_params=(:dlam=>0.5),
	spectrum_params=(:daa=>0.01, :aa_blue=>4840, :aa_red=>4880)) = begin

	set!(
		nml; 
		io_params=io_params, 
		timer_params=timer_params,
		atmos_params=atmos_params
		,m3d_params=m3d_params,
		linelist_params=linelist_params,
		spectrum_params=spectrum_params
	)
end

# ╔═╡ 7fc240a9-4f81-4a08-8734-ccd07b89e0bc
"""
	whole_spectrum(model_path; kwargs...)

Submit a job to the M3DIS code, which will compute the outgoing flux across the entire wavelength range.
"""
function whole_spectrum(model_path::String; 
					linelist="./input_multi3d/nlte_ges_linelist_jmg25jan2023_I_II",
					kwargs...)
	
	# Create an empty Namelist
	nml = M3DNamelist()

	# Fill in the defaults for this, they can be modified
	whole_spectrum_namelist!(nml)

	# additionally apply the model specific fields
	set(nml, linelist_params=(linelist=linelist,))
end

# ╔═╡ Cell order:
# ╟─bb45e6a0-2b8b-4731-bd53-190c44a5ea23
# ╠═49a4e636-166c-11ee-03ed-c129e8cdc7fc
# ╠═21951d47-3aef-4297-9edd-cf7afcb4544f
# ╠═4bb87ab7-94fa-460d-9c7b-4b7fe06c8500
# ╟─90818cad-34c7-4ad0-8a8a-530391a3e176
# ╠═2da3ff54-8b7a-4e7e-8b15-c56e34594875
# ╟─b7d1105d-3a57-4425-83dd-20fa03cd4194
# ╠═c490885a-588f-4277-be47-8a567babda27
# ╠═7fc240a9-4f81-4a08-8734-ccd07b89e0bc
