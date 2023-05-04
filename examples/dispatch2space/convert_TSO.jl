### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 059e7070-e033-11ed-109f-01304e9a85af
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using Glob
	using Interpolations
	using Distributed
end

# ╔═╡ 232daac9-8c9f-4da5-8efb-60497fdb4b94
md"# Converting Snapshots - Patches to Box"

# ╔═╡ 2dbf66fa-c6b5-4951-921c-d617ef4dd48f
md"## Setup"

# ╔═╡ 45e30f2d-ed66-47b3-9c3a-448b22e1b6a1
begin
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" EOS select
end

# ╔═╡ ddf76c9d-1a8c-40b6-9f70-b5407c0047d0
md"Inputs"

# ╔═╡ 0ed4ed19-2390-4d02-a3d6-aa609cd2e629
begin
	add_selection = false
	folder = MUST.@in_dispatch ARGS[1]
end

# ╔═╡ 536ecaa7-e6e0-4e13-b573-47f98755642c
md"Available snapshots in this output folder"

# ╔═╡ 7dbfec27-4b70-40ce-b9d0-30cb4f1cc149
begin
	content_of_folder = glob("*/", folder)
	available_snapshots = sort(MUST.list_of_snapshots(content_of_folder))
end

# ╔═╡ 86b3945f-97c4-4c76-a996-65e00a660163
md"Pick sublist of snapshots to convert. Always has to be a list."

# ╔═╡ 9d803429-0ae7-4896-9db2-3c3042b3a000
snapshots = [available_snapshots[end-2]]

# ╔═╡ 85f644bd-8c4a-4bb3-bfd8-9176bfaa247e
md"## Conversion
We read the associated input and see if there is an effective temperature file
"

# ╔═╡ 3bc36de5-0042-4e05-b574-270b738ed52a
nml = MUST.StellarNamelist(MUST.@in_dispatch(splitpath(folder)[end])*".nml")

# ╔═╡ 6e380259-b0e8-4be8-bc2e-ec1fbc39d3a2
md"The EoS"

# ╔═╡ 3d112dcf-4a30-4659-b6a7-682dcc5deb90
begin
	eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
	@info "Loading EoS from $(eos_path)"

	eos_sq = reload(SqEoS, MUST.@in_dispatch(
								joinpath(eos_path, "eos.hdf5")))
	opa_sq = reload(SqOpacity, MUST.@in_dispatch(
								joinpath(eos_path, "binned_opacities.hdf5")))
end

# ╔═╡ 2c5ce9b6-777a-482b-aa5b-3eb66be8ca7c
begin
	save_info = ispath(joinpath(folder, "teff.dat")) 
	if save_info
	    # Read and interpolate teff data
	    teff_func = MUST.teff_interpolated(joinpath(folder, "teff.dat"))
	end
end

# ╔═╡ c7d771a7-cbf1-448c-81b4-9b3be8919868
"""
	convert_snapshots(snapshots)
Convert a list of snapshots from the Distpatch format to a simple, corss-platform 
format (HDF5). Collects the info from patches and pieces together 3D arrays on
regular grids.
"""
function convert_snapshots(snapshots, eos=eos_sq, opa=opa_sq; 
										folder=folder, 
										add_selection=add_selection, 
										do_teff=save_info)
    for i_s in eachindex(snapshots)
        @info "Converting snapshot $(snapshots[i_s]) on worker $(myid())"
        if true
            # The dispatch snapshot object (Python)
            snap = dispatch.snapshot(snapshots[i_s], data=folder)

            # Units for conversion to CGS
            units = MUST.StaggerCGS(snap)

            # Convert its content to pure Julia
            s = if add_selection
				try
                	MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg)
				catch
					MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)
				end
            else
                MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)
            end

			@assert isfinite.(s.data[:d])

            # Apply the conversion
            if add_selection
				try
	                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
	                                    qr=:qr, pg=:p,
	                                    ux=:u, uy=:u, uz=:u,
	                                    x=:l, y=:l, z=:l)
				catch
					MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l)
				end
            else
                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l)
            end

			@assert isfinite.(s.data[:d])
    
            # Check if Teff should be added
            if do_teff
                @info "setting teff for sn $(snapshots[i_s])"
                MUST.set!(s.parameter, teff, ini_nml)
            end
            
            b_s = MUST.Box(s)
            s = nothing


			# we lookup aux quantities directly in the EoS and then add them.
			# We use the TSO lookup functions for this, which is faster in some cases
			@assert TSO.is_internal_energy(@axed(eos))

			lgd = log.(b_s[:d])
			lge = log.(b_s[:ee])

			T  = lookup(eos, :lnT, lgd, lge)
			kr = lookup(eos, opa, :κ_ross, lgd, lge)

			b_s.data[:T]  = exp.(T)
			b_s.data[:kr] = kr
			
            # Add the optical depth
            τ = MUST.optical_depth(b_s, opacity=:kr, density=:d)
            MUST.add!(b_s, τ, :τ_ross)

            # First save
            MUST.save(b_s; name="box_sn$(snapshots[i_s])", folder=folder)

			
            # NEW: Convert the height scale from cm to optical depth
            min_plane = MUST.plane_statistic(minimum, b_s, :τ_ross)
            b_s       = MUST.height_scale(b_s, :τ_ross)
            MUST.save(b_s; name="box_tau_sn$(snapshots[i_s])", folder=folder)

		else
            @warn "snapshot $(snapshots[i_s]) could not be loaded."
            continue
        end
    end
    nothing
end

# ╔═╡ 26c61249-bfac-4461-a91f-399a5f7c1000
md"Run the snapshot conversion"

# ╔═╡ 9512893f-759c-4a94-9087-4ccf2f075215
convert_snapshots(snapshots)

# ╔═╡ Cell order:
# ╟─232daac9-8c9f-4da5-8efb-60497fdb4b94
# ╟─2dbf66fa-c6b5-4951-921c-d617ef4dd48f
# ╠═059e7070-e033-11ed-109f-01304e9a85af
# ╠═45e30f2d-ed66-47b3-9c3a-448b22e1b6a1
# ╟─ddf76c9d-1a8c-40b6-9f70-b5407c0047d0
# ╠═0ed4ed19-2390-4d02-a3d6-aa609cd2e629
# ╟─536ecaa7-e6e0-4e13-b573-47f98755642c
# ╠═7dbfec27-4b70-40ce-b9d0-30cb4f1cc149
# ╟─86b3945f-97c4-4c76-a996-65e00a660163
# ╠═9d803429-0ae7-4896-9db2-3c3042b3a000
# ╟─85f644bd-8c4a-4bb3-bfd8-9176bfaa247e
# ╠═3bc36de5-0042-4e05-b574-270b738ed52a
# ╟─6e380259-b0e8-4be8-bc2e-ec1fbc39d3a2
# ╟─3d112dcf-4a30-4659-b6a7-682dcc5deb90
# ╟─2c5ce9b6-777a-482b-aa5b-3eb66be8ca7c
# ╠═c7d771a7-cbf1-448c-81b4-9b3be8919868
# ╟─26c61249-bfac-4461-a91f-399a5f7c1000
# ╠═9512893f-759c-4a94-9087-4ccf2f075215
