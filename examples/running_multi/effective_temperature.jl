### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 433c5690-2705-11ee-0b29-f9ac4b2948e6
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	import MUST
	using Plots	
end

# ╔═╡ caf8295a-4fc7-4b48-b2d4-b1140b42c44b
md"# MUST & M3D"

# ╔═╡ 8ee46151-9205-4e8b-ab33-636003bc713c
pc(x) = MUST.pyconvert(Any, x)

# ╔═╡ 814bb5f5-5e9c-4fda-b04b-48485600bff3
md"## Dispatch & M3D
We load the Python modules for handling Dispatch, M3D and loading the correct module locations."

# ╔═╡ 59848375-b439-46b1-8180-f3323cb966eb
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 87f69065-1a3a-443c-a132-d8435aaa975f
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 27a33e87-425d-40ae-ae1c-15202fbe6d61
md"## Model Atmospheres
Pick the model atmosphere. It should either be located in input_multi3d/MUST or the correct folder should be given to the namelist function."

# ╔═╡ 766cfb9d-6084-4cde-8856-d3f3c6509abe
modelatmos = "m3dis_sun_magg22_10x10x280_5"

# ╔═╡ 1a9a7946-fda5-4e01-b715-c3c36c1bf6ea
modelatmosfolder = "input_multi3d/sun_magg_120x240/"

# ╔═╡ 573b2070-18e5-4e53-8393-76549f2efce5
md"## Line List
One can either use no line list to compute Teff, use a line list in the Turbospectrum format, or use an absmet file. The absmet file contains sampled backgournd opacities, and might be good enough for the spectrum computation in terms of Teff. To activate the absmet, one has to set the linelist to nothing, in order for the namelist to be created correctly."

# ╔═╡ 1bd3e93c-0f57-4cf6-ace5-1f2f4fe9d57f
linelist = "./input_multi3d/vald_2490-25540.list"

# ╔═╡ fc1859d9-4a45-4372-95ad-403134ce7370
absmet = "./input_multi3d/absmet"

# ╔═╡ 908e3e3c-cd57-4221-967a-256e46246055
md"## Running M3D
Running M3D is then straight forward. One can either run this within an slurm allocation as a job step, directly run it as a new slurm job, or turn off slurm entirely. In the latter M3D will be run interactively. Please make sure you have enough resources available in this case. If you are using a line list this might easily not be the case. So be aware of that.\
If you have compiled ```M3D``` for a MPI system, it might be possible that you need to recompile without MPI for your interactive runs!"

# ╔═╡ ac2ba741-52f9-4192-837e-79d760548805
#=m3load = MUST.whole_spectrum(
	modelatmos, 
	namelist_kwargs=(
		:model_folder=>modelatmosfolder,
		:linelist=>nothing,
		:absmet=>nothing
	),
	slurm=false
)=#

# ╔═╡ 87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
md"## Computing Teff
The effective temperature is computed from integrating the Flux."

# ╔═╡ 3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
@info "Teff of $(modelatmos): $(MUST.Teff(m3load))"

# ╔═╡ Cell order:
# ╟─caf8295a-4fc7-4b48-b2d4-b1140b42c44b
# ╠═433c5690-2705-11ee-0b29-f9ac4b2948e6
# ╟─8ee46151-9205-4e8b-ab33-636003bc713c
# ╟─814bb5f5-5e9c-4fda-b04b-48485600bff3
# ╠═59848375-b439-46b1-8180-f3323cb966eb
# ╠═87f69065-1a3a-443c-a132-d8435aaa975f
# ╟─27a33e87-425d-40ae-ae1c-15202fbe6d61
# ╠═766cfb9d-6084-4cde-8856-d3f3c6509abe
# ╠═1a9a7946-fda5-4e01-b715-c3c36c1bf6ea
# ╟─573b2070-18e5-4e53-8393-76549f2efce5
# ╠═1bd3e93c-0f57-4cf6-ace5-1f2f4fe9d57f
# ╠═fc1859d9-4a45-4372-95ad-403134ce7370
# ╟─908e3e3c-cd57-4221-967a-256e46246055
# ╠═ac2ba741-52f9-4192-837e-79d760548805
# ╟─87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
# ╠═3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
