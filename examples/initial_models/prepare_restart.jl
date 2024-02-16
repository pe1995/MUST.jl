### A Pluto.jl notebook ###
# v0.19.36

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
end

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
md"""Select one of the available runs:\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""

# ╔═╡ cbbe908e-01f2-485c-a6f1-a915c841838b
runfolder = joinpath(datafolder, selectedRun);

# ╔═╡ b8a9a5ce-be8c-4ded-8914-96a9670b7f0e
md"## Select Snapshots"

# ╔═╡ 14aa40eb-cfa3-4d8e-b21e-ed9cb4e85b20
allsnaps = sort(MUST.list_of_snapshots(runfolder));

# ╔═╡ 288fb1ad-529e-40e1-ae8e-5994d438157d
md"Pick the snapshot from which you want to start: $(@bind firstSnap Slider(allsnaps, show_value=true, default=allsnaps[end-11]))
"

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
function prepare_restart(name, snapshot; decay_timescales=false, datadir=@in_dispatch("data"), kwargs...)
	# already converted snapshots
	snappath = joinpath(datadir, name, "m3dis_$(snapshot).mesh")
	snaptpath = joinpath(datadir, name, "m3dis_$(snapshot).atmos")
	csnap = isfile(snappath) &isfile(snaptpath)
	
	# first load the snapshot to get its time
	b = if !csnap
		@info "[$(name)] Converting snapshot $(snapshot)."
		b, _ = snapshotBox(
			snapshot; 
			folder=joinpath(datadir, name), 
			optical_depth_scale=false, 
			save_snapshot=true, 
			add_selection=false, 
			is_box=true, 
			to_multi=true
		)

		b
	else
		@info "[$(name)] Loading snapshot $(snapshot)."
		Box("box_sn$(snapshot)", folder=joinpath(datadir, name))
	end

	t = b.parameter.time

	# round the time to an integer
	t = floor(Int, t/60.)

	# the new namelist
	saveat = @in_dispatch("$(name)_t$(t).nml")

	# Name of the namelist of the current folder
    nml_name = @in_dispatch name

    # Init namelist
    nml = MUST.StellarNamelist(nml_name*".nml")

	# data without StAt
	datawithoutstat = last(split(datadir, "stellar_atmospheres/", keepempty=false))

    # change the relevant fields
	MUST.set!(
		nml, 
		stellar_params=(
			:initial_path=>joinpath(datawithoutstat, name, "m3dis_$(snapshot)"),
			:initial_model=>"must"
		)
	)

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

	MUST.set!(nml, kwargs...)
	MUST.write(nml, saveat)

	@info "[$(name)] Preparing restart complete."
	
	nothing
end

# ╔═╡ 27f1e1ac-77cf-459b-820e-b167ee11f391


# ╔═╡ 1e1e380f-9684-496c-8bf4-158f281a679b
md"Click to continue decaying timescales: $(@bind decay CheckBox(default=false))"

# ╔═╡ 946b586e-5c04-44c8-9300-0f2b7a8babf8
md"Click to start preparations: $(@bind convert_given CheckBox(default=false))"

# ╔═╡ 4c7c7f16-cfa4-4547-8355-72744211219a
begin
	if convert_given
		prepare_restart(
			selectedRun, 
			firstSnap, 
			datadir=datafolder, 
			decay_timescales=decay
		)
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
# ╟─b8a9a5ce-be8c-4ded-8914-96a9670b7f0e
# ╟─14aa40eb-cfa3-4d8e-b21e-ed9cb4e85b20
# ╟─288fb1ad-529e-40e1-ae8e-5994d438157d
# ╟─a1c18af7-0d55-4c6c-b503-3950e13b6b17
# ╟─9a08c5ec-bdab-4323-895c-7fec358b5af8
# ╟─f27c5eb0-6a6e-42bc-bf2e-f0aab0be2368
# ╟─4b782275-ce60-4ca0-ac17-be52d3573925
# ╟─27f1e1ac-77cf-459b-820e-b167ee11f391
# ╟─1e1e380f-9684-496c-8bf4-158f281a679b
# ╟─946b586e-5c04-44c8-9300-0f2b7a8babf8
# ╟─4c7c7f16-cfa4-4547-8355-72744211219a
