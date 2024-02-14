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
md"# Packing DISPATCH
Collect the last N dispatch snapshots + the namelist + the binned opacity table together in one folder that can easily be stored for the future. If monitoring is available, it will also be included in the packing. Since this involves conversions of possibly a large number of snapshots it is advisable to do this in a parallized script."

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
md"""Select one of the available monitored runs:\
$(@bind selectedRun confirm(Select(availableRuns(datafolder))))"""

# ╔═╡ cbbe908e-01f2-485c-a6f1-a915c841838b
runfolder = joinpath(datafolder, selectedRun);

# ╔═╡ b8a9a5ce-be8c-4ded-8914-96a9670b7f0e
md"## Select Snapshots"

# ╔═╡ 14aa40eb-cfa3-4d8e-b21e-ed9cb4e85b20
allsnaps = sort(MUST.list_of_snapshots(runfolder));

# ╔═╡ 288fb1ad-529e-40e1-ae8e-5994d438157d
md"Pick the first snapshot: $(@bind firstSnap Slider(allsnaps, show_value=true, default=allsnaps[end-11]))
"

# ╔═╡ 503f5ffb-d19f-4e10-b6d3-80a09d6f20ab
md"Pick the last snapshot: $(@bind lastSnap Slider(allsnaps, show_value=true, default=allsnaps[end-2]))
"

# ╔═╡ b8d80646-8f91-4f39-ae0e-3a22aa3a1e29
snapsSelected = allsnaps[findfirst(x->x==firstSnap, allsnaps):findfirst(x->x==lastSnap, allsnaps)];

# ╔═╡ 9a08c5ec-bdab-4323-895c-7fec358b5af8
md"## Convert Snapshots
Convert the snapshots (if not already done so). If you do this here it will take a bit longer because they will be converted sequentially. Please use the convert.jl script if you want to do this in parallel."

# ╔═╡ 4b782275-ce60-4ca0-ac17-be52d3573925
"""
	shipBox(name, snapshots; copymonitoring=true, copysnaps=true, saveat="pack_\$(name)", datadir=@in_dispatch("data"), kwargs...)

Convert (if needed) the given snapshots to `MUST.Box` objects, collect the EoS + possible monitoring in a common folder `saveat`.
"""
function shipBox(name, snapshots; copymonitoring=true, copysnaps=true, saveat="pack_$(name)", datadir=@in_dispatch("data"), kwargs...)
	# make sure output folder exists
	!isdir(saveat) && mkdir(saveat)

	# Name of the namelist of the current folder
    nml_name = @in_dispatch name

    # Init namelist
    nml = MUST.StellarNamelist(nml_name*".nml")

    # find the Squaregas EOS dir
    eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
	eos_name = last(split(eos_path, "/", keepempty=false))

	# copy over
	cp(nml_name*".nml", joinpath(saveat, name*".nml"), force=true)
	cp(MUST.@in_dispatch(eos_path), joinpath(saveat, eos_name), force=true)

	# copy monitoring
	if copymonitoring
		if isdir(joinpath(datadir, name, "monitoring"))
			cp(
				joinpath(datadir, name, "monitoring"), 
				joinpath(saveat, "monitoring"),
				force=true
			)
		end
	end

	# check for converted snapshots
	csnaps = MUST.list_snapshots(MUST.converted_snapshots(joinpath(datadir, name)))
	
	for (i, snap) in enumerate(snapshots)
		if !(snap in csnaps)
			@info "[$(name)] Converting snapshot $(snap) ($(i)/$(length(snapshots)))"
			snapshotBox(
				snap; 
				folder=joinpath(datadir, name), 
				optical_depth_scale=true, 
				save_snapshot=true, 
				add_selection=false, 
				is_box=true, 
				kwargs...
			)
		else
			@info "[$(name)] snapshot $(snap) already converted ($(i)/$(length(snapshots)))"
		end

		# after conversion, move (or copy) the snapshot
		snappath = joinpath(datadir, name, "box_sn$(snap).hdf5")
		snaptpath = joinpath(datadir, name, "box_tau_sn$(snap).hdf5")
		if isfile(snappath)
			if copysnaps
				cp(snappath, joinpath(saveat, "box_sn$(snap).hdf5"), force=true)
				cp(snaptpath, joinpath(saveat, "box_tau_sn$(snap).hdf5"), force=true)
			else
				mv(snappath, joinpath(saveat, "box_sn$(snap).hdf5"), force=true)
				mv(snaptpath, joinpath(saveat, "box_tau_sn$(snap).hdf5"), force=true)
			end
		end

		snappath = joinpath(datadir, name, "m3dis_$(snap).mesh")
		snaptpath = joinpath(datadir, name, "m3dis_$(snap).atmos")
		if isfile(snappath)
			if copysnaps
				cp(snappath, joinpath(saveat, "m3dis_$(snap).mesh"), force=true)
				cp(snaptpath, joinpath(saveat, "m3dis_$(snap).atmos"), force=true)
			else
				mv(snappath, joinpath(saveat, "m3dis_$(snap).mesh"), force=true)
				mv(snaptpath, joinpath(saveat, "m3dis_$(snap).atmos"), force=true)
			end
		end
	end

	@info "[$(name)] packing complete."
	nothing
end

# ╔═╡ 27f1e1ac-77cf-459b-820e-b167ee11f391


# ╔═╡ 169240b9-90c2-4a5c-9f37-0ad566993278
md"Tick box to also convert to Multi3D: $(@bind convert_m3d CheckBox(default=true))"

# ╔═╡ 5c02c516-3a8e-48a7-8ac5-ebd7380fde43
md"Tick box to move snapshots instead of copy: $(@bind moveSnaps CheckBox(default=false))"

# ╔═╡ 1e1e380f-9684-496c-8bf4-158f281a679b


# ╔═╡ 946b586e-5c04-44c8-9300-0f2b7a8babf8
md"Tick box to start packing: $(@bind convert_given CheckBox(default=false))"

# ╔═╡ 4c7c7f16-cfa4-4547-8355-72744211219a
begin
	if convert_given
		shipBox(
			selectedRun, 
			snapsSelected, 
			datadir=datafolder, 
			to_multi=convert_m3d,
			copysnaps=!moveSnaps
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
# ╟─503f5ffb-d19f-4e10-b6d3-80a09d6f20ab
# ╟─b8d80646-8f91-4f39-ae0e-3a22aa3a1e29
# ╟─9a08c5ec-bdab-4323-895c-7fec358b5af8
# ╟─4b782275-ce60-4ca0-ac17-be52d3573925
# ╟─27f1e1ac-77cf-459b-820e-b167ee11f391
# ╟─169240b9-90c2-4a5c-9f37-0ad566993278
# ╟─5c02c516-3a8e-48a7-8ac5-ebd7380fde43
# ╟─1e1e380f-9684-496c-8bf4-158f281a679b
# ╟─946b586e-5c04-44c8-9300-0f2b7a8babf8
# ╟─4c7c7f16-cfa4-4547-8355-72744211219a
