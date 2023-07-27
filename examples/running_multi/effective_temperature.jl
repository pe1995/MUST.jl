### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 433c5690-2705-11ee-0b29-f9ac4b2948e6
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using PythonPlot	
	using LaTeXStrings
end

# ╔═╡ caf8295a-4fc7-4b48-b2d4-b1140b42c44b
md"# MUST & M3D"

# ╔═╡ a73a987f-aae3-484d-9a96-b5c79754d22a
visual = MUST.ingredients("visual.jl")

# ╔═╡ c4c30863-3676-42f8-b13b-37692c3c4d40
plt = PythonPlot.matplotlib.pyplot

# ╔═╡ 0b3db435-1e02-4852-b881-b7b997b78e4a
PythonPlot.matplotlib.rcParams["font.size"] = 12

# ╔═╡ 8ee46151-9205-4e8b-ab33-636003bc713c
begin
	pc(x) = MUST.pyconvert(Any, x)
	pc(T, x) = MUST.pyconvert(T, x)
end

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
modelatmos = "m3dis_sun_magg22_20x20x280_3"

# ╔═╡ 1a9a7946-fda5-4e01-b715-c3c36c1bf6ea
modelatmosfolder = "input_multi3d/magg2022_150x300/"
#modelatmosfolder = "input_multi3d/stagger_sun/"

# ╔═╡ 76bfd78b-0d58-4c59-9587-4b7627a2d357
snapshots = [
	"m3dis_sun_magg22_20x20x280_1",
	"m3dis_sun_magg22_20x20x280_2",
	"m3dis_sun_magg22_20x20x280_3",
	"m3dis_sun_magg22_20x20x280_4",
	"m3dis_sun_magg22_20x20x280_5",
	"m3dis_sun_magg22_20x20x280_6",
	"m3dis_sun_magg22_20x20x280_7",
	"m3dis_sun_magg22_20x20x280_8",
	"m3dis_sun_magg22_20x20x280_9"
]

# ╔═╡ 8dfa58c2-54e4-4ca1-bf80-4badadfaf2b9
#snapshots = ["m3dis_sun_stagger_10x10x230_1"]

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

# ╔═╡ ecc8099e-ead7-4826-b17a-bb5b1787bfd2
compute = true

# ╔═╡ ac2ba741-52f9-4192-837e-79d760548805
begin
	m3druns = if compute
		m3load = MUST.whole_spectrum(
			snapshots, 
			namelist_kwargs=(
				:model_folder=>modelatmosfolder,
				:linelist=>nothing,
				:absmet=>absmet,
				:atom_params=>(:atom_file=>"", ),
				:spectrum_params=>(
					:daa=>1.0, :aa_blue=>1000, :aa_red=>100000
				)
			),
			slurm=true
		)
		m3load
	else
		MUST.M3DISRun.(joinpath.("data", snapshots))
	end
end

# ╔═╡ 87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
md"## Computing Teff
The effective temperature is computed from integrating the Flux."

# ╔═╡ bc423ddb-8027-4b97-a88f-a87188cbcb5f
stagger = MUST.M3DISRun(joinpath.("data", "m3dis_sun_stagger_10x10x230_1"))

# ╔═╡ 1bfef52c-7c67-4029-8dde-9149583b28a6
teff_nan(run) = MUST.Teff(run.lam[.!isnan.(run.flux)], run.flux[.!isnan.(run.flux)])

# ╔═╡ 3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
for (i, m) in enumerate(m3druns)
	@info "Teff of $(snapshots[i]): $(teff_nan(m))"
end

# ╔═╡ d74a883f-39f6-4675-bf1f-a14032fca71d
begin
	plt.close()
	fS, axS = plt.subplots(1, 1, figsize=(5, 5))
	visual.basic_plot!(axS)

	axS.plot(
		teff_nan.(m3druns), 
		color="k", 
		marker="s",
		markerfacecolor="white",
		markersize=8,
		lw=1
	)

	axS.axhline(teff_nan(stagger), color="red")

	axS.set_xlabel(L"\rm snapshot")
	axS.set_ylabel(L"\rm T_{eff}\ [K]")
	
	gcf()
end

# ╔═╡ 7fb7723f-0960-4735-8442-5e25d8a60783
begin
	plt.close()
	fF, axF = plt.subplots(1, 1, figsize=(5, 5))
	visual.basic_plot!(axF)

	l,f = pc.(Array, m3druns[1].crop(per_aa=true))
	axF.plot(
		l, f, 
		color="k", 
		marker="",
		lw=1
	)
	
	axF.set_xlabel(L"\rm \lambda\ [\AA]")
	axF.set_ylabel(L"\rm flux\ [normalized]")
	axF.set_xscale("log")
	axF.set_yscale("log")
	
	gcf()
end

# ╔═╡ Cell order:
# ╟─caf8295a-4fc7-4b48-b2d4-b1140b42c44b
# ╠═433c5690-2705-11ee-0b29-f9ac4b2948e6
# ╠═a73a987f-aae3-484d-9a96-b5c79754d22a
# ╠═c4c30863-3676-42f8-b13b-37692c3c4d40
# ╠═0b3db435-1e02-4852-b881-b7b997b78e4a
# ╟─8ee46151-9205-4e8b-ab33-636003bc713c
# ╟─814bb5f5-5e9c-4fda-b04b-48485600bff3
# ╠═59848375-b439-46b1-8180-f3323cb966eb
# ╠═87f69065-1a3a-443c-a132-d8435aaa975f
# ╟─27a33e87-425d-40ae-ae1c-15202fbe6d61
# ╠═766cfb9d-6084-4cde-8856-d3f3c6509abe
# ╠═1a9a7946-fda5-4e01-b715-c3c36c1bf6ea
# ╠═76bfd78b-0d58-4c59-9587-4b7627a2d357
# ╠═8dfa58c2-54e4-4ca1-bf80-4badadfaf2b9
# ╟─573b2070-18e5-4e53-8393-76549f2efce5
# ╠═1bd3e93c-0f57-4cf6-ace5-1f2f4fe9d57f
# ╠═fc1859d9-4a45-4372-95ad-403134ce7370
# ╟─908e3e3c-cd57-4221-967a-256e46246055
# ╠═ecc8099e-ead7-4826-b17a-bb5b1787bfd2
# ╠═ac2ba741-52f9-4192-837e-79d760548805
# ╟─87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
# ╠═bc423ddb-8027-4b97-a88f-a87188cbcb5f
# ╠═1bfef52c-7c67-4029-8dde-9149583b28a6
# ╠═3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
# ╠═d74a883f-39f6-4675-bf1f-a14032fca71d
# ╟─7fb7723f-0960-4735-8442-5e25d8a60783
