### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ dc738eee-74a5-11ee-01ab-1b3f09b09f68
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
    using MUST
	using Mmap
    using Glob
    using Distributed
    using Interpolations
	using PythonPlot
end

# ╔═╡ 78859bf7-e996-4e2c-8c2a-a312d29bcfe5
np = MUST.pyimport("numpy")

# ╔═╡ 8644e9cf-6c28-4389-93f0-457a8decce01
plt = matplotlib.pyplot

# ╔═╡ e760784c-b8f2-492a-a51f-89e76ad375e1
begin
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
    #MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" EOS select
end

# ╔═╡ 99bae142-19fa-4a9d-8fa9-1377a34655b6
begin
	folder = MUST.@in_dispatch "data/grid_t50g40m00"
	content_of_folder = glob("*/", folder)
	snapshots = sort(MUST.list_of_snapshots(content_of_folder))[100]
end

# ╔═╡ 9ac522f9-6bc8-4021-8f05-a3352cff6aeb
begin	
	# Name of the namelist of the current folder
	nml_name = MUST.@in_dispatch splitpath(folder)[end]
	
	# Init namelist
	nml = MUST.StellarNamelist(nml_name*".nml")
	
	# Use the new Squaregas EOS 
	eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
	eos_sq = MUST.SquareGasEOS(MUST.@in_dispatch(eos_path))
end

# ╔═╡ e7e2a3e9-13c9-4c42-8695-fd5ddce0b61a
begin 
	b = []
	for i_s in eachindex(snapshots)
		@info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
	
		# The dispatch snapshot object (Python)
		snap = dispatch.snapshot(snapshots[i_s], data=folder)
		
		#=patch = snap.patches[4]
		@show patch.filename 
		@show patch.idx.__dict__.keys() patch.mesh_type snap.io.__dict__.keys()
		
		patch = snap.patches[4]
		@show patch.filename
		@show patch.idx.__dict__.keys() patch.mesh_type snap.io.__dict__.keys()
		

		fname = MUST.pyconvert(Any, patch.filename)
		shp = Tuple(MUST.pyconvert.(Any, patch.ncell))
		off = MUST.pyconvert(Vector{Int}, patch.offset)
		varidx = MUST.pyconvert(Any, patch.idx.__dict__["d"]) + 1
		li = MUST.pyconvert(Vector{Int}, patch.li)
		ui = MUST.pyconvert(Vector{Int}, patch.ui)
		

		@show fname shp off varidx
		v = mmap(fname, Array{Float32, length(shp)}, shp, off[varidx])
		@show v[li[1]:ui[1], li[2]:ui[2], li[3]:ui[3]]

		break
		#v = MUST.pyconvert(Array{Float32, 3}, np.array(patch.var(String(:d))))

		@show patch.var(MUST.Py("d"))=#
		
		# Units for conversion to CGS
		units = MUST.StandardUnits(snap)

		# Convert its content to pure Julia
		s = MUST.Space(snap, :d, :ux, :uy, :uz, :e, :ee)
		
		# Apply the conversion
		MUST.convert!(s, units; d=:d, e=:e, ee=:ee,
								ux=:u, uy=:u, uz=:u,
								x=:l, y=:l, z=:l, time=:t)

		# Add additional columns already in CGS after converting
		MUST.add_from_EOS!(s, eos_sq, :T)
		MUST.add_from_EOS!(s, eos_sq, :kr)

		# Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
		b_s = MUST.Box(s)
		s = nothing
		append!(b, [b_s])
	end
end

# ╔═╡ d93548fb-8dec-4b49-b37b-2479e4e91b5a
b[1][:T][:, : , end]

# ╔═╡ 107dc754-0dbc-4d16-9bb2-abcb73313535
begin
	plt.close()
	f, ax = plt.subplots(1, 1)

	rms5(x) = sqrt.(MUST.mean(x .^2)) ./1e5
	ax.plot(profile(MUST.mean, b[1], :z, :T)...)
	
	gcf()
end

# ╔═╡ Cell order:
# ╠═dc738eee-74a5-11ee-01ab-1b3f09b09f68
# ╠═78859bf7-e996-4e2c-8c2a-a312d29bcfe5
# ╠═8644e9cf-6c28-4389-93f0-457a8decce01
# ╠═e760784c-b8f2-492a-a51f-89e76ad375e1
# ╠═99bae142-19fa-4a9d-8fa9-1377a34655b6
# ╠═9ac522f9-6bc8-4021-8f05-a3352cff6aeb
# ╠═e7e2a3e9-13c9-4c42-8695-fd5ddce0b61a
# ╠═d93548fb-8dec-4b49-b37b-2479e4e91b5a
# ╠═107dc754-0dbc-4d16-9bb2-abcb73313535
