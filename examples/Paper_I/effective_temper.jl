### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 29f57c4c-b063-11ee-1771-119c66d0a87b
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate("."); 
	using MUST
	using Glob
	using PythonPlot
	using TSO
	using LaTeXStrings
	using Printf
	using DelimitedFiles
	using DataFrames
	using CSV

	plt = matplotlib.pyplot
end;

# ╔═╡ 36dbc744-5d71-4275-9cfa-a3bda2c88e5a
begin
	np = MUST.pyimport("numpy")
	matplotlib.style.use(joinpath(dirname(pathof(MUST)), "Bergemann2023.mplstyle"))
	hl = 4
end;

# ╔═╡ 7c9e05b7-25f5-43e7-b43c-c2e445c45dce
md"# Dispatch and Multi3D"

# ╔═╡ 4de81775-5fe3-4bb4-ac45-fb612fef7566
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 4cb1ca95-d8ea-403f-9eec-af8e02f71b07
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 167acebd-c974-43bf-8345-14851b2d57ef
md"# Dispatch models"

# ╔═╡ b9bf719c-bb48-4e2d-a198-f070fb374965
names = [
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t55g45m00_v0.1",
	"DIS_MARCS_E_t60g45m00_v0.1",
	"DIS_MARCS_E_t65g45m00_v0.1",
	#"DIS_MARCS_E_t45g40m00_v0.1"
]

# ╔═╡ ae5dfae7-b170-4547-a4ad-c42a785dacdb
out_folder = [
	"M3DIS",
	"M3DIS-lres",
	"M3DIS-ires",
	"M3DIS-hres",
	"M3DIS-t55g45m00",
	"M3DIS-t60g45m00",
	"M3DIS-t65g45m00",
	#"M3DIS-t45g40m00"
]

# ╔═╡ 50d2f42d-e0aa-4537-8f1d-8b22cfd8c8fe
labels = [
	"model D",
	"model A",
	"model B",
	"model C",
	"t55g45m00",
	"t60g45m00",
	"t65g45m00"
]

# ╔═╡ 11e18101-36b6-449e-b9f1-613b52b9e4ca
ls = [
	"-",
	"--",
	"-.",
	":",
	"-",
	"--",
	"-.",
]

# ╔═╡ 7c94e58f-28a6-46b4-9d92-03708a928f08
models = Dict(
	"best"      => 1,
	"lres"      => 2,
	"ires"      => 3,
	"hres"      => 4,
	"t55g45m00" => 5,
	"t60g45m00" => 6,
	"t65g45m00" => 7,
	#"t45g40m00" => 8
)

# ╔═╡ 37d0a236-41fb-4a41-bc55-f06c6c606ed3
md"# Running Multi3D
You can either run it with slurm (if you run this notebook on a cluster) or directly on this computer."

# ╔═╡ 96743ca1-db5a-44b9-8fad-0fbca7fcf2d5
modelatmosfolders = abspath.(out_folder)

# ╔═╡ 240110da-7b96-4f31-8e4e-b64169939327
eos_name = joinpath.(modelatmosfolders, "eos_opa")

# ╔═╡ c0d37d0b-f082-44a6-acdc-7607c9aac8ad
getmodelnames(f) = begin
	allfiles = glob("m3dis*", f)
	allfiles = first.(split.(basename.(allfiles), ".", keepempty=false))

	mask = (occursin.("20x20", allfiles)) .& (.!occursin.("sun", allfiles))
	unique(convert.(String, allfiles[mask]))
end

# ╔═╡ dbcdba93-6c53-4997-af48-098b50826a48
snapshots_m3dis = [getmodelnames(f) for f in modelatmosfolders]

# ╔═╡ 107fd679-1470-4bed-a7a8-db0e48b0fb2d


# ╔═╡ 5f55e6e0-b865-48a3-be64-c066466aadc3
compute = false

# ╔═╡ d89c555b-0138-4d91-94ef-8ddf0e70e54d
begin
	 if compute
		 for (k, i) in models
			MUST.heating(                                 
				snapshots_m3dis[i], 
				eos_name[i],
				name="_$(k)",
				namelist_kwargs=(
					:model_folder=>modelatmosfolders[i],
					:linelist=>nothing,
					:absmet=>nothing,
					:atom_params=>(:atom_file=>"", ),
					:atmos_params=>(
						:atmos_format=>"MUST", 
						:use_density=>true, 
						:use_ne=>false
					),
				), 
				slurm=false
			)
		 end
	end
end

# ╔═╡ fc8722cd-aa7a-46df-b32d-11163f15c573


# ╔═╡ 6c9892a3-43fe-407b-9209-361a244b74a4
md"We reload the results from the corresponding folders."

# ╔═╡ 315ce664-ee3e-4972-8105-427385f06f8e
runs = Dict(
	k => [
		MUST.M3DISRun(
			MUST.@in_m3dis joinpath("data", n*"_$(k)_binned")
		)
		for n in snapshots_m3dis[i]
	]
	for (k, i) in models
)

# ╔═╡ 2323051f-5570-49c6-860a-9d98383b7ab7
md"# Effective Temperature
We can now compute the effective temperature from the binned flux. We do this for every snapshot and then compute the average."

# ╔═╡ 3c9439bb-93cb-447c-a81c-d647e8084ee2
Tₑ = Dict(
	k => [MUST.Teff(run.flux) for run in runs[k]] 
	for (k, i) in models
)

# ╔═╡ ce49c1f7-6996-44fe-99cf-d11dde3da369
list_models = collect(keys(models))

# ╔═╡ b247d428-a791-4229-8262-4b6d3744cc32
data = (
	"model" => [labels[models[k]] for k in list_models],
	"<Teff>" => [MUST.mean(Tₑ[k]) for k in list_models],
	"sigma Teff" => [MUST.std(Tₑ[k]) for k in list_models],
	"max Teff" => [maximum(Tₑ[k]) for k in list_models],
	"min Teff" => [minimum(Tₑ[k]) for k in list_models],
)

# ╔═╡ dad1097d-8525-417c-9eaa-5f5a03598c89
sorted_names = sortperm(last(data[1]))

# ╔═╡ defbbb5b-77b4-4b51-9e21-a4c6bfbdae25
tab = DataFrame((k=>v[sorted_names] for (k, v) in data)...);

# ╔═╡ 1b1773b7-70c5-472c-9fd1-ebde6847f0a7


# ╔═╡ d28a6516-d591-42cb-873c-125e3aa3f570
md"# Effective Temperature Statistics
$(tab)"

# ╔═╡ 4a79fd2e-e2b5-4703-8bce-f10d751d5c8a
CSV.write("effective_temperature_stats.csv", tab);

# ╔═╡ 5785c134-29a3-4936-a4d8-dcd68a084861


# ╔═╡ 7e12567a-3a23-4648-83c7-fd18c2a64816
cmap = plt.cm.get_cmap("rainbow")

# ╔═╡ b567f7ee-e943-4120-9f94-39eece98552a
let
	plt.close()
	
	f, ax = plt.subplots(1, 1, figsize=(5, 6))

	modelnames = list_models[sorted_names]
	colors = [cmap(i/length(modelnames)) for i in eachindex(modelnames)]
	colors = ["k", "k", "k", "k", "r", "r", "r"]

	for (i, k) in enumerate(modelnames)
		x = range(1, length(runs[k])) |> collect
		ax.plot(
			x, Tₑ[k], 
			label=labels[models[k]], 
			color=colors[i], 
			marker="", 
			markerfacecolor="None", 
			lw=2, markeredgewidth=1, 
			ls=ls[models[k]]
		)
	end

	ax.set_ylabel(@L_str("\\rm T_{eff}\\ [K]"))
	ax.set_xlabel(@L_str("\\rm snapshot"))

	ax.legend(
		framealpha=0, labelspacing=0.01, handlelength=hl, ncol=2,
		bbox_to_anchor=(0.5, 0.6), loc="lower center"
	)
	

	gcf()
end

# ╔═╡ Cell order:
# ╠═29f57c4c-b063-11ee-1771-119c66d0a87b
# ╠═36dbc744-5d71-4275-9cfa-a3bda2c88e5a
# ╟─7c9e05b7-25f5-43e7-b43c-c2e445c45dce
# ╠═4de81775-5fe3-4bb4-ac45-fb612fef7566
# ╠═4cb1ca95-d8ea-403f-9eec-af8e02f71b07
# ╟─167acebd-c974-43bf-8345-14851b2d57ef
# ╠═b9bf719c-bb48-4e2d-a198-f070fb374965
# ╠═ae5dfae7-b170-4547-a4ad-c42a785dacdb
# ╠═50d2f42d-e0aa-4537-8f1d-8b22cfd8c8fe
# ╠═11e18101-36b6-449e-b9f1-613b52b9e4ca
# ╠═7c94e58f-28a6-46b4-9d92-03708a928f08
# ╟─37d0a236-41fb-4a41-bc55-f06c6c606ed3
# ╟─96743ca1-db5a-44b9-8fad-0fbca7fcf2d5
# ╟─240110da-7b96-4f31-8e4e-b64169939327
# ╟─c0d37d0b-f082-44a6-acdc-7607c9aac8ad
# ╟─dbcdba93-6c53-4997-af48-098b50826a48
# ╟─107fd679-1470-4bed-a7a8-db0e48b0fb2d
# ╠═5f55e6e0-b865-48a3-be64-c066466aadc3
# ╟─d89c555b-0138-4d91-94ef-8ddf0e70e54d
# ╟─fc8722cd-aa7a-46df-b32d-11163f15c573
# ╟─6c9892a3-43fe-407b-9209-361a244b74a4
# ╠═315ce664-ee3e-4972-8105-427385f06f8e
# ╟─2323051f-5570-49c6-860a-9d98383b7ab7
# ╟─3c9439bb-93cb-447c-a81c-d647e8084ee2
# ╟─ce49c1f7-6996-44fe-99cf-d11dde3da369
# ╟─b247d428-a791-4229-8262-4b6d3744cc32
# ╟─dad1097d-8525-417c-9eaa-5f5a03598c89
# ╟─defbbb5b-77b4-4b51-9e21-a4c6bfbdae25
# ╟─1b1773b7-70c5-472c-9fd1-ebde6847f0a7
# ╟─d28a6516-d591-42cb-873c-125e3aa3f570
# ╟─4a79fd2e-e2b5-4703-8bce-f10d751d5c8a
# ╟─5785c134-29a3-4936-a4d8-dcd68a084861
# ╟─7e12567a-3a23-4648-83c7-fd18c2a64816
# ╟─b567f7ee-e943-4120-9f94-39eece98552a
