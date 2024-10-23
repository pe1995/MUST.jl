### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 8068b412-910a-11ef-2c08-e7a811075511
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ 82b0b107-7de4-4038-83d5-c4f603e59a86
plt = matplotlib.pyplot

# ╔═╡ 7751c1d4-15a0-4a6c-bb51-4a9e35e56f5c
@import_dispatch "../../../dispatch2"

# ╔═╡ 538710f8-9788-4a43-b838-ead75f2f3be8
datafolder = @in_dispatch "data" 

# ╔═╡ 87d5bd1b-6a69-4b58-91d9-d3b022cd189d


# ╔═╡ 47a6d82b-8a91-4c86-bd3f-4c6d35d5badd
begin
	runname = "ST8_E_t57.77g44.40m0.000_v1.0_mancha_s4"
	#runname = "ST11_E_t57.77g44.40m0.000_v1.0"
	quantities = MUST.defaultQuantities
	
	run, rundir, datadir, params_list, files = MUST._check_files(28, datafolder, runname, nothing)

	# load the units and prepare the converters
	units = MUST.StandardUnits(rundir)
	l_conv = MUST.standardConversion(:l)
	converters = Dict(
		q.name=>q.conversion(units)
		for q in quantities
	)

	# add parameter groups from data/run/NNNNN/*snapshot.nml
	snapshot_nml, variablesSym, idxVars, idxs = MUST._snapshot_variables(files)

	# rank namelists
    rank_files = [f for f in MUST.glob("*", datadir) if occursin("_patches.nml", f)]
	rank_nmls = MUST.FreeNamelist.(rank_files, :patch_nml)

	# Use the Squaregas EOS 
	eos_reader = (x)->MUST._squaregaseos(x, eos_path=nothing)
    eos_sq, eos_quantities = eos_reader(run)

	# First we read the meta data so that we can build the data cube
	patchMeta, patchDataFiles, patchAuxFiles, auxnames = MUST._patchmeta(
		datadir, rank_nmls, snapshot_nml
	)

	# decide if there are quantities to skip because there is no aux data for them
	# and they can not be derived otherwise
	variableMaskOk = [ (!q.derived)|(q.name in Symbol.(auxnames))|(!isnothing(q.recipe)) for q in quantities]
	quantities = quantities[variableMaskOk]

	# prepare the HDF5 file for the final box
	variablesAndAux = [MUST.varnames(quantities)..., eos_quantities...]
	variablesDerived = MUST.derived(quantities)
	
end

# ╔═╡ 247ef560-5e97-4f00-aa0f-7c5025c8caf3
patchDataFiles

# ╔═╡ 708834f3-e81e-4b98-bda8-3a30352b22a6
begin
	for (i, patch) in enumerate(patchMeta)
		#@show patch.li patch.ui patch.li_with_guards patch.ui_with_guards
		#@show patch.li.-patch.ng patch.ui.+patch.ng
		#@show patch.x[patch.li.-patch.ng] patch.x[patch.li.+patch.ng]
	end
end

# ╔═╡ 7f856856-5f24-42cc-a271-0fcb51eba184
function read_patch(patchDataFile, patchMeta, off, ng)
	buffer = Array{Float32, length(Tuple(patchMeta.ncell))}(undef, Tuple(patchMeta.ncell))

	data = Dict()
	open(patchDataFile, "r") do f
		for (j, var) in enumerate(variablesSym)
			varidx = idxs[String(var)] + 1
		    offset = off[varidx]

		    # seek to the correct offset
		    seek(f, offset)
		
		    # Read the full array data
		    read!(f, buffer)

			data[var] = reshape(deepcopy(buffer), ng...)
		end
	end

	data
end

# ╔═╡ 4a621b8b-33ee-48e8-b353-6bbe0b91c32c
begin
	plt.close()
	f, ax = plt.subplots(1, 1)
	patchdata = []
	patchdata_bnd = []

	u = 1
	ax.set_title("z = ui + $(u)")

	zu = maximum([maximum(p.zi) for p in patchMeta])
	for (i, patch) in enumerate(patchMeta)
		if maximum(patch.zi) >= zu
			li = patch.li
			ui = patch.ui
			xx, yy = MUST.meshgrid(patch.xi, patch.yi)
			data = read_patch(patchDataFiles[i], patchMeta[i], patch.offset, patch.n .+ 2 .*patch.ng)

			auxname = patchAuxFiles[i]
			aux = MUST.readaux(auxname)
			for a in aux
				asym = Symbol(a.name)
				if asym in variablesAndAux
					data[asym] = deepcopy(a.data)
				end
			end
			
			#d = data[:e][li[1]:ui[1],li[2]:ui[2],ui[3]+2] ./ data[:d][li[1]:ui[1],li[2]:ui[2],ui[3]+2]
			d = data[:qr][li[1]:ui[1],li[2]:ui[2],ui[3]+u]
			
			append!(patchdata, [(xx, yy, d)])

			xx, yy = MUST.meshgrid(patch.x, patch.y)
			#d = data[:e][:,:,ui[3]+2] ./ data[:d][:,:,ui[3]+2]
			d = data[:qr][:,:,ui[3]+u]
			
			d[li[1]:ui[1],li[2]:ui[2]] .= NaN
			append!(patchdata_bnd, [(xx, yy, d)])
		end
	end

	vmin = minimum([minimum(d[3]) for d in patchdata])
	vmax = maximum([maximum(d[3]) for d in patchdata])
	for i in eachindex(patchdata)
		xx, yy, d = patchdata[i]
		im = ax.pcolormesh(
			xx, 
			yy, 
			d, 
			vmin=vmin, vmax=vmax, rasterized=true, cmap="rainbow"
		)
		xx, yy, d = patchdata_bnd[i]
		im = ax.pcolormesh(
			xx, 
			yy, 
			d, 
			vmin=vmin, vmax=vmax, 
			rasterized=true, cmap="rainbow"
		)
	end
	
	f.colorbar(im, ax=ax)
	f
end

# ╔═╡ b3b2b560-90b3-4844-90dd-d29b7a7e113f
let
	plt.close()
	f, ax = plt.subplots(1, 1)
	patchdata = []
	patchdata_bnd = []

	u = -0
	ax.set_title("z = ui + $(u)")

	zu = maximum([maximum(p.yi) for p in patchMeta])
	for (i, patch) in enumerate(patchMeta)
		if maximum(patch.yi) >= zu
			li = patch.li
			ui = patch.ui
			xx, yy = MUST.meshgrid(patch.xi, patch.zi)
			data = read_patch(patchDataFiles[i], patchMeta[i], patch.offset, patch.n .+ 2 .*patch.ng)

			auxname = patchAuxFiles[i]
			aux = MUST.readaux(auxname)
			for a in aux
				asym = Symbol(a.name)
				if asym in variablesAndAux
					data[asym] = deepcopy(a.data)
				end
			end
			
			
			#d = data[:e][li[1]:ui[1],ui[2]+u,li[3]:ui[3]] ./ data[:qr][li[1]:ui[1],ui[2]+u,li[3]:ui[3]]
			d = data[:qr][li[1]:ui[1],ui[2]+u,li[3]:ui[3]]
			
			append!(patchdata, [(xx, yy, d)])

			xx, yy = MUST.meshgrid(patch.x, patch.z)
			#d = data[:e][:,ui[2]+u,:] ./ data[:qr][:,ui[2]+u,:]
			d = data[:qr][:,ui[2]+u,:]
			
			d[li[1]:ui[1],li[3]:ui[3]] .= NaN
			append!(patchdata_bnd, [(xx, yy, d)])
		end
	end
	
	vmin = minimum([minimum(d[3]) for d in patchdata])
	vmax = maximum([maximum(d[3]) for d in patchdata])
	for i in eachindex(patchdata)
		n =1
		xx, yy, d = patchdata[i]
		im = ax.pcolormesh(
			xx, 
			yy, 
			d./n, 
			vmin=vmin, vmax=vmax, rasterized=true, cmap="rainbow"
		)
		xx, yy, d = patchdata_bnd[i]
		im = ax.pcolormesh(
			xx, 
			yy, 
			d./n, 
			vmin=vmin, vmax=vmax,  rasterized=true, cmap="rainbow"
		)
	end
	
	f.colorbar(im, ax=ax)
	f
end

# ╔═╡ Cell order:
# ╠═8068b412-910a-11ef-2c08-e7a811075511
# ╠═82b0b107-7de4-4038-83d5-c4f603e59a86
# ╠═7751c1d4-15a0-4a6c-bb51-4a9e35e56f5c
# ╠═538710f8-9788-4a43-b838-ead75f2f3be8
# ╟─87d5bd1b-6a69-4b58-91d9-d3b022cd189d
# ╠═47a6d82b-8a91-4c86-bd3f-4c6d35d5badd
# ╠═247ef560-5e97-4f00-aa0f-7c5025c8caf3
# ╠═708834f3-e81e-4b98-bda8-3a30352b22a6
# ╠═7f856856-5f24-42cc-a271-0fcb51eba184
# ╠═4a621b8b-33ee-48e8-b353-6bbe0b91c32c
# ╠═b3b2b560-90b3-4844-90dd-d29b7a7e113f
