### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ d5a09b78-36ae-11ee-0605-77cc2bb6e03d
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using TSO
	using Glob
	using Plots
	PythonCall = MUST.PythonCall
end;

# ╔═╡ d8e88318-8824-47dc-9a8b-46577dc06df5
md"# Patches to ```MUST.jl``` Boxes"

# ╔═╡ 747455cb-8984-4777-8ce1-c823d80a0839
md"## Import DISPATCH"

# ╔═╡ 1d1e9fa3-b0af-474a-8f50-63fa54986885
begin
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
	MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" EOS select
	add_selection = false
end;

# ╔═╡ 7a185a6f-a3c5-48e6-9bfe-32199fcc5ec0
convert2must = ingredients("convert2must.jl") 

# ╔═╡ a43acc06-bb39-4fdf-9b0f-03555f500742
md"## Simulation Information"

# ╔═╡ 2ebbbcf5-e0b4-4d31-922a-1899d2f03c09
md"Select the available snapshots from a given simulation"

# ╔═╡ 49e3c9ac-f2c0-459b-8242-60b9659efada
simulation = @in_dispatch "data/pretty_good_sun_new_magg6_strong"

# ╔═╡ b374cea0-c80e-45fc-bc87-ee1ba43b5551
begin
	content_of_folder = glob("*/", simulation)
    snapshots = sort(MUST.list_of_snapshots(content_of_folder))
end

# ╔═╡ d77bf831-2e5d-4968-9569-d1110506e480
md"Name of the namelist of the current folder"

# ╔═╡ 14da665d-501c-4862-b9e0-763f974ad4dc
nml_name = @in_dispatch splitpath(simulation)[end]

# ╔═╡ af821310-7a57-46c6-b1b1-e9bf369dc99d
nml = MUST.StellarNamelist(nml_name*".nml")

# ╔═╡ a2912124-faad-406d-bb73-1deecc0d379f
begin
	eos_path = @in_dispatch replace(nml.eos_params["table_loc"], "'"=>"")
	eos_sq = MUST.SquareGasEOS(eos_path)
end

# ╔═╡ 8257f871-01ad-4a2b-b2ba-9e9903a8211d
md"We can also load it in the ```TSO.jl``` format to use this infrastructure"

# ╔═╡ e5d2c232-0b45-4b1b-8401-f507981eafe6
begin
	eos = reload(SqEoS, joinpath(eos_path, "eos.hdf5"))
	opa = reload(SqOpacity, joinpath(eos_path, "binned_opacities.hdf5"))
end

# ╔═╡ cc478462-6c3b-4116-ac90-527b93ad6c42
md"## Convert the snapshots"

# ╔═╡ 09197272-a346-4401-a640-cfa637dbcc3b
b, bτ = Box(simulation, snapshots |> first, eos=eos, opacity=opa)

# ╔═╡ 9940aa4c-7507-423f-8061-abeefc2f12da
begin
	plot(framestyle=:box, grid=false)
	plot!(profile(MUST.mean, bτ, :log10τ_ross, :T)...)
end

# ╔═╡ Cell order:
# ╟─d8e88318-8824-47dc-9a8b-46577dc06df5
# ╠═d5a09b78-36ae-11ee-0605-77cc2bb6e03d
# ╟─747455cb-8984-4777-8ce1-c823d80a0839
# ╟─1d1e9fa3-b0af-474a-8f50-63fa54986885
# ╠═7a185a6f-a3c5-48e6-9bfe-32199fcc5ec0
# ╟─a43acc06-bb39-4fdf-9b0f-03555f500742
# ╟─2ebbbcf5-e0b4-4d31-922a-1899d2f03c09
# ╠═49e3c9ac-f2c0-459b-8242-60b9659efada
# ╠═b374cea0-c80e-45fc-bc87-ee1ba43b5551
# ╟─d77bf831-2e5d-4968-9569-d1110506e480
# ╟─14da665d-501c-4862-b9e0-763f974ad4dc
# ╟─af821310-7a57-46c6-b1b1-e9bf369dc99d
# ╟─a2912124-faad-406d-bb73-1deecc0d379f
# ╟─8257f871-01ad-4a2b-b2ba-9e9903a8211d
# ╟─e5d2c232-0b45-4b1b-8401-f507981eafe6
# ╟─cc478462-6c3b-4116-ac90-527b93ad6c42
# ╠═09197272-a346-4401-a640-cfa637dbcc3b
# ╠═9940aa4c-7507-423f-8061-abeefc2f12da
