### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d69a0642-cb28-11ee-2eb3-d33c850d3799
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using PlutoUI
	using PythonPlot

	plt = matplotlib.pyplot
end;

# ╔═╡ 2189049a-20e6-43fb-96cb-20ae7994f96b
TableOfContents()

# ╔═╡ ac0a2cde-9228-47e8-b6cf-b866d5fd39d7
md"# Restarting DISPATCH
Currently in order to restart a DISPATCH simulation it is easiest to not restart at all. Currently your best option is to convert one of the later snapshots to the M3D format and start a fresh DISPATCH run from this snapshot. In this notebook you can perform all the needed steps automatically. The time of the snapshot you selected will be added to the name of the new run, to that you know where you started from."

# ╔═╡ a7983431-7ccc-4b06-bb00-92052b31065e
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ 9c3bba9f-bed4-41da-b06f-7e58c38ff657
datafolder = @in_dispatch "data"

# ╔═╡ 7871f99e-8495-4d27-a3ab-3d56f3846eef


# ╔═╡ 3c70fb6d-ff6d-4359-af7e-d7a91a7cef75
md"## Select Simulation"

# ╔═╡ 0a0c5350-2b15-46e6-b907-e580e6f2c3a1
availableRuns(folder) = begin
	allruns = MUST.glob("*/", folder)
	split.(allruns, "/", keepempty=false) .|> last
end

# ╔═╡ 44fa2200-f245-4e60-928c-69b7c7a1dd23
availableMonitor(folder) = begin
	allruns = MUST.glob("*/", folder)
	mask = isdir.(joinpath.(allruns, "monitoring"))
	split.(allruns[mask], "/", keepempty=false) .|> last
end

# ╔═╡ e595f80f-5ab4-421e-a47e-8b80d1612b29


# ╔═╡ d45457a4-9e33-4492-9a06-930f41f10fdc
md"""__Select one of the available runs__\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""

# ╔═╡ cbbe908e-01f2-485c-a6f1-a915c841838b
runfolder = joinpath(datafolder, selectedRun);

# ╔═╡ 263c6c0f-9310-4373-901e-ec0a74fc870c
md"""
Alternatively, you can select a list of run paths that should be restarted. Note that in this case, the namelist will be searched in the same folder as the data. This is convenient if you want to restart packed models. It will always restart from the second-to-last snapshot that is available.

	Enter path to the list you want to restart
$(@bind selectedRunList confirm(TextField(35)))
"""

# ╔═╡ 994d81c3-454f-4145-a32c-00d16dba6b05
if length(selectedRunList) > 0
	runNamesList = open(selectedRunList, "r") do f
		strip.(readlines(f))
	end
end;

# ╔═╡ 12653e83-2724-400f-b70e-ca9112e2e94a
findEoS(folder) = begin
	allfolders = MUST.glob("*/", folder)
	found = false
	fol = Ref("")
	for fo in allfolders
		allfiles = MUST.glob("*", fo)
		for f in allfiles
			if occursin("eos.hdf5", f)
				found=true
				fol[] = fo
				break
			end
		end

		found && break
	end

	if !found
		nothing
	else	
		abspath(fol[])
	end
end

# ╔═╡ f96bbc8a-ebae-4a50-b6b0-01e3d211fa7b
if length(selectedRunList) > 0
	eosNamesList = []
	folderNamesList = []
	runnameList = []
	for d in runNamesList
		eosfolder = findEoS(d)
		runname = last(split(eosfolder, '/', keepempty=false))
		
		if !isnothing(eosfolder)
			append!(eosNamesList, [String(eosfolder)])
			append!(folderNamesList, [String(d)])
			append!(runnameList, [String(runname)])
		end
	end

	runnameList
end

# ╔═╡ 1d51c950-69cb-4fa7-aa5e-8d4cce225516


# ╔═╡ b8a9a5ce-be8c-4ded-8914-96a9670b7f0e
md"## Select Snapshots"

# ╔═╡ 14aa40eb-cfa3-4d8e-b21e-ed9cb4e85b20
allsnaps = if length(selectedRunList) > 0
	sn = MUST.list_snapshots.(MUST.converted_snapshots.(folderNamesList))
	i = maximum(length.(sn))
	0:-1:-i+1 |> collect
else
	sort(MUST.list_of_snapshots(runfolder))
end;

# ╔═╡ 3f76e640-40a5-45a6-b46e-4b5a119fb180
if length(selectedRunList) > 0
	md"""
		Pick the snapshot (from the end of each run) from which you want to start
	$(@bind firstSnap Slider(allsnaps, show_value=true, default=allsnaps[1]))
	"""
else
	md"""
		Pick the snapshot from which you want to start
	$(@bind firstSnap Slider(allsnaps, show_value=true, default=allsnaps[end-1]))
	"""
end

# ╔═╡ a1c18af7-0d55-4c6c-b503-3950e13b6b17


# ╔═╡ 9a08c5ec-bdab-4323-895c-7fec358b5af8
md"## Prepare Restart
The restart is rather simple. A new namelist will be created as a copy of the old one. The new namelist will be modified such that the chosen snapshot will be used as initial condition in the M3D format."

# ╔═╡ f27c5eb0-6a6e-42bc-bf2e-f0aab0be2368
"""
	decay_newton_friction(time, nml) 

Decay the timescales to the appropriate time, so that one can continiue even during the newton or friction decay phases.
"""
decay_newton_friction!(time, nml) = begin
	time_scaling = nml.scaling_params["t_cgs"]
	newton_end_time = nml.newton_params["end_time"]
	friction_end_time = nml.friction_params["end_time"]
	current_time = time / time_scaling

	# the end time would be in the past
	new_newton_end_time = newton_end_time - current_time	
	new_friction_end_time = friction_end_time - current_time

	MUST.set!(
		nml,
		newton_params=(
			:end_time=>new_newton_end_time,
		),
		friction_params=(
			:end_time=>new_friction_end_time,
		)
	)
end

# ╔═╡ 4b782275-ce60-4ca0-ac17-be52d3573925
"""
	prepare_restart(name, snapshot; decay_timescales=false, datadir=@in_dispatch("data"), kwargs...)

Prepare the restart of a DISPATCH simulation by creating the new namelist and converting the chosen snapshot to the M3D format (if not already done). Kwargs can be used to modify the namelist. If `decay_timescales` is set, the newton and friction timescales are decaying even after the restart.
"""
function prepare_restart(name, snapshot; decay_timescales=false, datadir=joinpath(@in_dispatch("data"), name), eos_path=nothing, nml_name=@in_dispatch(name)*".nml", extend=0.0, kwargs...)
	# already converted snapshots
	snappath = joinpath(datadir, "m3dis_$(snapshot).mesh")
	snaptpath = joinpath(datadir, "m3dis_$(snapshot).atmos")
	csnap = isfile(snappath) &isfile(snaptpath)

	# first load the snapshot to get its time
	b = if !csnap
		@info "[$(name)] Converting snapshot $(snapshot)."
		b, _ = snapshotBox(
			snapshot; 
			folder=datadir, 
			optical_depth_scale=false, 
			save_snapshot=true, 
			#add_selection=false, 
			#is_box=true, 
			to_multi=true,
			legacy=false,
			eos_path=eos_path
		)

		if isnothing(b)
			Box("box_sn$(snapshot)", folder=datadir)
		else
			b
		end
	else
		@info "[$(name)] Loading snapshot $(snapshot)."
		Box("box_sn$(snapshot)", folder=datadir)
	end

	t = b.parameter.time

	# round the time to an integer
	t = floor(Int, t/60.)

	# the new namelist
	saveat = @in_dispatch("$(name)_t$(t).nml")

    # Init namelist
    nml = MUST.StellarNamelist(nml_name)

	# data folder without StAt
	datawithoutstat = if occursin("stellar_atmospheres/", datadir)
		last(split(datadir, "stellar_atmospheres/", keepempty=false))
	else
		datadir
	end

    # change the relevant fields
	MUST.set!(
		nml, 
		stellar_params=(
			:initial_path=>joinpath(datadir, "m3dis_$(snapshot)"),
			:initial_model=>"must"
		)
	)

	if !isnothing(eos_path)
		MUST.set!(
			nml,
			eos_params=(
			:table_loc=>eos_path,
			)
		)
	end

	# Adjust the cube size if one wants to extend the cube
	if extend>0
		l_old = nml.scaling_params["l_cgs"]
		size_old = nml.cartesian_params["size"]
		shift_old = nml.cartesian_params["position"]

		l_new = l_old*(1. +extend)
		
		# we make the cube in general bigger
		MUST.set!(nml, scaling_params=("l_cgs"=>l_new,))

		# now it is bigger, however we need to also shift it down by a bit 
		# more to retain the optical surface location
		# we only want to extend above the surface
		zbottom = (-size_old[3]/2.0 + shift_old[3])  * l_old / l_new
		shift_new = abs(size_old[3]/2.0 - abs(zbottom))

		# add a 1% tolerance
		shift_new = round(shift_new - 0.01*shift_new, sigdigits=3)
		
		MUST.set!(nml, cartesian_params=("position"=>[0.0,0.0,-shift_new],))		
	end

	if decay_timescales
		decay_newton_friction!(b.parameter.time, nml)
	else
		MUST.set!(
			nml,
			newton_params=(
				:on=>false,
				:delay_rt=>false
			),
			friction_params=(
				:on=>false,
			)
		)
	end

	MUST.set!(nml; kwargs...)
	MUST.write(nml, saveat)

	@info "[$(name)] Preparing restart complete."
	@info "[$(name)] Saved at $(saveat)."
	
	saveat
end

# ╔═╡ 27f1e1ac-77cf-459b-820e-b167ee11f391
function parametersFromName(path; carbon=false, vmic=false)
	datadir = if occursin('/', path)
		split(path, '/', keepempty=false) .|> last
	else
		path
	end
	
	parts = split(datadir, '_', keepempty=false)
	irelevantPart = findfirst(
		x->(occursin("t", x)&occursin("g", x)&occursin("m", x)),
		parts
	)
	relevantPart = parts[irelevantPart]
	t1 = findfirst(x->x=='t', relevantPart) + 1
	t2 = findfirst(x->x=='g', relevantPart) - 1
	t3 = findfirst(x->x=='m', relevantPart) - 1

	t = parse(Float64, relevantPart[t1:t2]) .* 100
	g = parse(Float64, relevantPart[t2+2:t3]) ./10
	m = parse(Float64, relevantPart[t3+2:end])

	d = Dict("teff"=>t, "logg"=>g, "feh"=>m)
	
	# find other stuff
	if carbon
		relevantPart = parts[2]
		t1 = findfirst(x->x=='c', relevantPart) + 1
		t2 = findfirst(x->x=='v', relevantPart) - 1
		c = parse(Float64, relevantPart[t1:t2])
		d["cfe"] = c
	end

	if vmic
		relevantPart = parts[2]
		t1 = findfirst(x->x=='v', relevantPart) + 1
		t2 = length(relevantPart)
		c = parse(Float64, relevantPart[t1:t2])

		d["vmic"] = c
	end

	d
end

# ╔═╡ cf03da3f-8a94-4697-a4e4-0e0f9c13b803
begin
	 #  using Teff + logg (because scale height ~ T/g) for htop scale
	x0 = 5777.0 / exp10(4.44)
	y0 = 1.1
	x1 = 0.29
	y1 = 1.5
	m = (y1 - y0) / (x1 - x0)
	htop_scale(teff, logg) = m * (teff/exp10(logg)) + (y1-m*x1)
end

# ╔═╡ b64ad2d7-78d4-49dd-9998-4d58166a649d


# ╔═╡ 0f2b46ac-9317-4d23-be32-e7d0265f6258
md"You can modify parameters of the new namelist by giving them specifically here.
For example: Set some new htop_scale at the upper boundary with some initial damping of the velocity field."

# ╔═╡ 3efea52b-cafa-4440-90fe-d6b226847a2d
new_namelist_params = if length(selectedRunList) > 0
	# get teff logg from name
	paras = parametersFromName.(runnameList)
	htop = [htop_scale(d["teff"], d["logg"]) for d in paras]
	
	 [Dict(
		:friction_params=>(
			:on=>true,
			:end_time=>50,
			:decay_scale=>20
		),
		:boundary_params=>(
			:htop_scale=>htop[i],
		),
	) for i in eachindex(runnameList)]
else
	new_namelist_params = Dict(
		#:friction_params=>(
		#	:on=>true,
		#	:end_time=>50,
		#	:decay_scale=>5
		#),
		#:boundary_params=>(
		#	:htop_scale=>7.0,
		#),
		#:cartesian_params=>(
		#	:dims=>[12,12,6],
		#),
		#:patch_params=>(
		#	:n=>[20,20,20],
		#),
		#:sc_rt_params=>(
		#	:rt_res=>[-1,-1,40],
		#)
	)
end

# ╔═╡ ac0f2184-e817-432a-860b-36cd2a03af3e


# ╔═╡ 28d04061-0040-45b5-9490-54c4cdb943f4
md"__Enter factor to extend atmosphere (%):__ $(@bind exAt TextField(default=\"0.0\"))"

# ╔═╡ 1e1e380f-9684-496c-8bf4-158f281a679b
md"__Click to continue decaying timescales:__ $(@bind decay CheckBox(default=false))"

# ╔═╡ 946b586e-5c04-44c8-9300-0f2b7a8babf8
md"__Click to start preparations:__ $(@bind convert_given CheckBox(default=false))"

# ╔═╡ 4c7c7f16-cfa4-4547-8355-72744211219a
begin
	restartModels = []
	if length(selectedRunList) > 0
		if convert_given
			for (i, run) in enumerate(runnameList)
				snaps = MUST.list_snapshots(
					MUST.converted_snapshots(folderNamesList[i])
				)
				isnap = if length(snaps) == 0
					@warn "There are no snaps for $(run)."
					continue
				elseif length(snaps) < abs(firstSnap)
					1
				end
				nmlName = joinpath(folderNamesList[i], run*".nml")
				s = prepare_restart(
					run, 
					snaps[end+firstSnap]; 
					datadir=abspath(folderNamesList[i]), 
					decay_timescales=decay,
					extend=parse(Float64, exAt)/100.0,
					nml_name=nmlName,
					eos_path=eosNamesList[i],
					new_namelist_params[i]...
				)

				append!(restartModels, [s])
			end
		end
	else
		if convert_given
			prepare_restart(
				selectedRun, 
				firstSnap; 
				datadir=joinpath(datafolder, selectedRun),
				decay_timescales=decay,
				extend=parse(Float64, exAt)/100.0,
				new_namelist_params...
			)
		end
	end
end

# ╔═╡ 18941dbb-1eed-4836-82a8-68350208ffeb


# ╔═╡ 184e5c6a-6629-499b-a1af-b7a0542b2120
if length(selectedRunList) > 0
	md"""
		Pick a name for the new models submission
	
	$(@bind newModelsName confirm(TextField(35, default="failed_runsX")))
	"""
end

# ╔═╡ 8b17e383-2bc1-4362-9f60-ca9cf55ebc8b
if length(selectedRunList) > 0
	if length(newModelsName) > 0
		!isdir(@in_dispatch(newModelsName)) && mkdir(@in_dispatch(newModelsName))
		for model in restartModels
			nml_name = split(model, "/", keepempty=false) |> last
			cp(model, @in_dispatch(joinpath(newModelsName, nml_name)), force=true)
		end
		@info "Model namelists collected in $(@in_dispatch(newModelsName))."
	end
end

# ╔═╡ Cell order:
# ╠═d69a0642-cb28-11ee-2eb3-d33c850d3799
# ╟─2189049a-20e6-43fb-96cb-20ae7994f96b
# ╟─ac0a2cde-9228-47e8-b6cf-b866d5fd39d7
# ╠═a7983431-7ccc-4b06-bb00-92052b31065e
# ╠═9c3bba9f-bed4-41da-b06f-7e58c38ff657
# ╟─7871f99e-8495-4d27-a3ab-3d56f3846eef
# ╟─3c70fb6d-ff6d-4359-af7e-d7a91a7cef75
# ╟─0a0c5350-2b15-46e6-b907-e580e6f2c3a1
# ╟─44fa2200-f245-4e60-928c-69b7c7a1dd23
# ╟─e595f80f-5ab4-421e-a47e-8b80d1612b29
# ╟─d45457a4-9e33-4492-9a06-930f41f10fdc
# ╟─cbbe908e-01f2-485c-a6f1-a915c841838b
# ╟─263c6c0f-9310-4373-901e-ec0a74fc870c
# ╟─994d81c3-454f-4145-a32c-00d16dba6b05
# ╟─12653e83-2724-400f-b70e-ca9112e2e94a
# ╟─f96bbc8a-ebae-4a50-b6b0-01e3d211fa7b
# ╟─1d51c950-69cb-4fa7-aa5e-8d4cce225516
# ╟─b8a9a5ce-be8c-4ded-8914-96a9670b7f0e
# ╟─14aa40eb-cfa3-4d8e-b21e-ed9cb4e85b20
# ╟─3f76e640-40a5-45a6-b46e-4b5a119fb180
# ╟─a1c18af7-0d55-4c6c-b503-3950e13b6b17
# ╟─9a08c5ec-bdab-4323-895c-7fec358b5af8
# ╟─f27c5eb0-6a6e-42bc-bf2e-f0aab0be2368
# ╟─4b782275-ce60-4ca0-ac17-be52d3573925
# ╟─27f1e1ac-77cf-459b-820e-b167ee11f391
# ╟─cf03da3f-8a94-4697-a4e4-0e0f9c13b803
# ╟─b64ad2d7-78d4-49dd-9998-4d58166a649d
# ╟─0f2b46ac-9317-4d23-be32-e7d0265f6258
# ╠═3efea52b-cafa-4440-90fe-d6b226847a2d
# ╟─ac0f2184-e817-432a-860b-36cd2a03af3e
# ╟─28d04061-0040-45b5-9490-54c4cdb943f4
# ╟─1e1e380f-9684-496c-8bf4-158f281a679b
# ╟─946b586e-5c04-44c8-9300-0f2b7a8babf8
# ╟─4c7c7f16-cfa4-4547-8355-72744211219a
# ╟─18941dbb-1eed-4836-82a8-68350208ffeb
# ╟─184e5c6a-6629-499b-a1af-b7a0542b2120
# ╟─8b17e383-2bc1-4362-9f60-ca9cf55ebc8b
