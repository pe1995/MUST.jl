### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 433c5690-2705-11ee-0b29-f9ac4b2948e6
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using TSO
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
modelatmos = "m3dis_sun_magg22_20x20x299_3"

# ╔═╡ 1a9a7946-fda5-4e01-b715-c3c36c1bf6ea
modelatmosfolder = "input_multi3d/magg2022_150x300/"
#modelatmosfolder = "input_multi3d/stagger_sun/"
#modelatmosfolder = "input_multi3d/atmos"
#modelatmosfolder = "input_multi3d/magg2022_195x390/"

# ╔═╡ 76bfd78b-0d58-4c59-9587-4b7627a2d357
#=snapshots = [
	"m3dis_sun_magg22_20x20x299_1",
	"m3dis_sun_magg22_20x20x299_2",
	"m3dis_sun_magg22_20x20x299_3",
	"m3dis_sun_magg22_20x20x299_4",
	"m3dis_sun_magg22_20x20x299_5",
	"m3dis_sun_magg22_20x20x299_6",
	"m3dis_sun_magg22_20x20x299_7",
	"m3dis_sun_magg22_20x20x299_8",
	"m3dis_sun_magg22_20x20x299_9"
]=#
#snapshots = ["m3dis_sun_stagger_20x20x230_1"]
#snapshots = ["t5777g44m0005_20.5x5x230"]
#snapshots = ["atmos.sun_MARCS"]
snapshots = [
	"m3dis_sun_magg22_20x20x299_1",
	"m3dis_sun_magg22_20x20x299_2",
	"m3dis_sun_magg22_20x20x299_3",
	"m3dis_sun_magg22_20x20x299_4",
	"m3dis_sun_magg22_20x20x299_5",
	"m3dis_sun_magg22_20x20x299_6",
	"m3dis_sun_magg22_20x20x299_7",
	"m3dis_sun_magg22_20x20x299_8",
	"m3dis_sun_magg22_20x20x299_9"
]

# ╔═╡ 4843f8e2-b1ce-4926-94bc-e3781098fb97
snapshots_m3dis = [
	"m3dis_sun_magg22_20x20x299_1",
	"m3dis_sun_magg22_20x20x299_2",
	"m3dis_sun_magg22_20x20x299_3",
	"m3dis_sun_magg22_20x20x299_4",
	"m3dis_sun_magg22_20x20x299_5",
	"m3dis_sun_magg22_20x20x299_6",
	"m3dis_sun_magg22_20x20x299_7",
	"m3dis_sun_magg22_20x20x299_8",
	"m3dis_sun_magg22_20x20x299_9"
]

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
	 if compute
		m3load = MUST.whole_spectrum(
			snapshots, 
			namelist_kwargs=(
				:model_folder=>modelatmosfolder,
				:linelist=>nothing,
				:absmet=>absmet,
				:atom_params=>(:atom_file=>"", ),
				:spectrum_params=>(
					:daa=>1.0, :aa_blue=>1000, :aa_red=>100000
				),
				:composition_params=>(
					:absdat_file=>"./input_multi3d/TS_absdat.dat", :abund_file=>"./input_multi3d/abund_magg"
				),
				:atmos_params=>(
					:atmos_format=>"MUST", 
					:use_density=>true, 
					:use_ne=>false
				),
			),
			slurm=true
		)
	end
end

# ╔═╡ 0320767c-8bb4-4aa1-bf19-d44e0aa52ec5
m3druns = MUST.M3DISRun.(joinpath.("data", snapshots_m3dis))

# ╔═╡ 87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
md"## Computing Teff
The effective temperature is computed from integrating the Flux."

# ╔═╡ bc423ddb-8027-4b97-a88f-a87188cbcb5f
stagger = MUST.M3DISRun(joinpath.("data", "m3dis_sun_stagger_20x20x230_1"))

# ╔═╡ 6df8b739-f894-43d1-a4c6-4907b54e96a4
marcs = MUST.M3DISRun(joinpath.("data", "atmos.sun_MARCS"))

# ╔═╡ 1bfef52c-7c67-4029-8dde-9149583b28a6
teff_nan(run) = MUST.Teff(run.lam[.!isnan.(run.flux)], run.flux[.!isnan.(run.flux)])

# ╔═╡ 3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
for (i, m) in enumerate(m3druns)
	@info "Teff of $(snapshots_m3dis[i]): $(teff_nan(m))"
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
		markersize=8, label=L"\rm M3DIS",
		lw=1
	)

	axS.axhline(teff_nan(stagger), color="red", label=L"\rm Stagger")
	axS.axhline(teff_nan(marcs), color="blue", label=L"\rm MARCS")
	

	axS.set_xlabel(L"\rm snapshot")
	axS.set_ylabel(L"\rm T_{eff}\ [K]")

	axS.legend(framealpha=0, labelspacing=0.001, handlelength=4, loc="lower right")
	
	gcf()
end

# ╔═╡ 9495f7be-08bb-4805-bb98-9a91248799da
u(λ, T) = (8π*TSO.HPlanck * TSO.CLight) / λ^5 * 1 / (exp(TSO.hc_k / λ /T) - 1)

# ╔═╡ 7fb7723f-0960-4735-8442-5e25d8a60783
begin
	plt.close()
	fF, axF = plt.subplots(1, 1, figsize=(5, 5))
	visual.basic_plot!(axF)

	l,f = pc.(Array, m3druns[1].crop(per_aa=true, norm=false))
	axF.plot(
		l, f ./ 1e-8 , 
		color="k", 
		marker="",
		lw=1
	)
	lp = range(800, 300000, length=2000)
	fp = TSO.Bλ.(lp, 5777.0) .* 4π
	axF.plot(lp, fp)
	
	axF.set_xlabel(L"\rm \lambda\ [\AA]")
	axF.set_ylabel(L"\rm flux")
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
# ╠═4843f8e2-b1ce-4926-94bc-e3781098fb97
# ╟─573b2070-18e5-4e53-8393-76549f2efce5
# ╠═1bd3e93c-0f57-4cf6-ace5-1f2f4fe9d57f
# ╠═fc1859d9-4a45-4372-95ad-403134ce7370
# ╟─908e3e3c-cd57-4221-967a-256e46246055
# ╠═ecc8099e-ead7-4826-b17a-bb5b1787bfd2
# ╠═ac2ba741-52f9-4192-837e-79d760548805
# ╠═0320767c-8bb4-4aa1-bf19-d44e0aa52ec5
# ╟─87b0d3a3-1b62-49d6-bc1c-36ae26ef314f
# ╠═bc423ddb-8027-4b97-a88f-a87188cbcb5f
# ╠═6df8b739-f894-43d1-a4c6-4907b54e96a4
# ╠═1bfef52c-7c67-4029-8dde-9149583b28a6
# ╠═3f6a9327-35fc-4e82-bf5d-f9c1e0d658e8
# ╠═d74a883f-39f6-4675-bf1f-a14032fca71d
# ╠═9495f7be-08bb-4805-bb98-9a91248799da
# ╠═7fb7723f-0960-4735-8442-5e25d8a60783
