### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 49a4e636-166c-11ee-03ed-c129e8cdc7fc
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".");
	import MUST
	using Plots
end

# ╔═╡ bb45e6a0-2b8b-4731-bd53-190c44a5ea23
md"# Prepare & Run M3DIS"

# ╔═╡ a6d00a44-bb7d-4a40-9b45-c40a6d7e1f1b
pc(x) = MUST.pyconvert(Any, x)

# ╔═╡ 21951d47-3aef-4297-9edd-cf7afcb4544f
MUST.@import_m3dis "/u/peitner/DISPATCH/Multi3D"

# ╔═╡ 4bb87ab7-94fa-460d-9c7b-4b7fe06c8500
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ 90818cad-34c7-4ad0-8a8a-530391a3e176
md"## The Namelists
Similar to the automatic execution of the stellar atmospheres code, the Multi code can be executed within a slurm environment. The overall code is the same, so first a namelist has to be created. One can choose a native default namelist, or just construct it from scratch. All parameters that are needed are available as fields. Note that this strictly speaking is not needed, but can easier to fit into the rest of the existing code."

# ╔═╡ b7d1105d-3a57-4425-83dd-20fa03cd4194
md"A couple of defaults however are usefull when it comes to compute the effective temperature of a model."

# ╔═╡ b8f2e1cf-8d1c-402e-a4ac-347cf60612d6
nml = MUST.whole_spectrum_namelist("m3dis_sun_12x12x120")

# ╔═╡ 0ab5616f-3b2d-48a4-a7ac-eebb0632fc7c
md"## The Submission
The following convenient shortcut will execute M3D and return the output in an convenient format for us to look at. This actually runs if you are connected to a cluster using the slurm system and it will take some time."

# ╔═╡ ec0efed0-6565-4b95-a82a-84d9662b6f5c
#=m3load = MUST.whole_spectrum(
	"t5777g44m0005_20.5x5x230", 
	namelist_kwargs=(
		:model_folder=>"./input_multi3d/atmos",
		:atmos_params=>(:atmos_format=>"Multi",)
	)
)=#

# ╔═╡ c9dc6751-4287-4962-8e4b-73e644548346
md"## The Output
If the code is executed, it will return the result directly. However, it is also possible to load a precomputed setup just like in the m3dis python package."

# ╔═╡ 71e71f24-348a-42dc-a639-ecebf78326ea
m3load = MUST.M3DISRun("data/t5777g44m0005_20.5x5x230")

# ╔═╡ cbd69346-5ab9-4ce1-8c5e-da69e08f9d87
f = m3load.flux

# ╔═╡ d49fd89a-79cd-4901-b4d8-e738c778e945
l = m3load.lam

# ╔═╡ 6a39476c-0124-4290-82f1-ec51831c0cb7
length(l)

# ╔═╡ 73144c07-6580-476e-b9a6-30fbd829ed58
begin
	x, y = pc.(m3load.line[8].crop())
	plot(x, y)
end

# ╔═╡ 825ec223-06c2-431e-abad-037d4792fef6
md"And then integrate it over the entire spectrum"

# ╔═╡ c970adc5-2b61-45c4-abf2-28764eb567cb
@info "The effective temperature is: $(MUST.Teff(m3load))"

# ╔═╡ 32c8243d-d310-4b0a-b9b7-80b39e181dc7
begin
	plot(framestyle=:box, grid=false)

	plot!(log10.(l), f, color=:black, label=nothing)
	plot!(ylabel="flux", xlabel="log10 λ")
	#plot!(xlim=log10.([6562.8-10, 6562.8+10]))
end

# ╔═╡ 24ae360b-64b7-486d-8533-2ab0160ff4d3
md"## Experiment: replace windows"

# ╔═╡ 95286c76-a5b7-486b-b789-89df7cbc94d1
function replacewindows(run)
	l, F = deepcopy(run.lam), deepcopy(run.flux)

	for (i, line) in enumerate(run.line)
		ll, ff = pc.(line.crop(LTE=true, norm=false))

		if length(ll) <= 1
			continue
		end
		
		window = [minimum(ll), maximum(ll)]
		maskwindow = (l .> first(window)) .& (l .< last(window))
		if count(maskwindow) > length(l)
			continue
		end

		i_first = findfirst(maskwindow)
		i_last  = findlast(maskwindow)

		l2 = similar(l, 0)
		append!(l2, l[1:i_first-1])
		append!(l2, ll)
		append!(l2, l[i_last+1:end])

		F2 = similar(F, 0)
		append!(F2, F[1:i_first-1])
		append!(F2, ff)
		append!(F2, F[i_last+1:end])

		F = F2
		l = l2
	end

	l, F
end

# ╔═╡ 63fed658-da7c-4fd2-9994-ee044f3b0a0e
#lnew, Fnew = replacewindows(m3load)

# ╔═╡ b0873543-404b-4665-a9ef-ab56f689fe64
md"## Planck function"

# ╔═╡ 9b59b9af-2897-4b2f-8627-98e6cc53ea01
Bλ = MUST.B(l, 5477.0)

# ╔═╡ 04f8df8a-d760-47cd-af07-decc992d7c63
begin
	plot(framestyle=:box, grid=false)

	colors=palette(:hot, 10)

	plot!(log10.(l), f ./ maximum(f), color=:black, alpha=0.5)
	for (i,T) in enumerate(range(5000, 6000, length=10))
		Bi = MUST.B(l, T)
		plot!(log10.(l), Bi ./ maximum(Bi), label="$T", color=colors[i])
		vline!([log10(MUST.wien(T))], label=nothing, color=colors[i])
	end

	plot!()
end

# ╔═╡ Cell order:
# ╟─bb45e6a0-2b8b-4731-bd53-190c44a5ea23
# ╠═49a4e636-166c-11ee-03ed-c129e8cdc7fc
# ╠═a6d00a44-bb7d-4a40-9b45-c40a6d7e1f1b
# ╠═21951d47-3aef-4297-9edd-cf7afcb4544f
# ╠═4bb87ab7-94fa-460d-9c7b-4b7fe06c8500
# ╟─90818cad-34c7-4ad0-8a8a-530391a3e176
# ╟─b7d1105d-3a57-4425-83dd-20fa03cd4194
# ╠═b8f2e1cf-8d1c-402e-a4ac-347cf60612d6
# ╟─0ab5616f-3b2d-48a4-a7ac-eebb0632fc7c
# ╠═ec0efed0-6565-4b95-a82a-84d9662b6f5c
# ╟─c9dc6751-4287-4962-8e4b-73e644548346
# ╠═71e71f24-348a-42dc-a639-ecebf78326ea
# ╠═cbd69346-5ab9-4ce1-8c5e-da69e08f9d87
# ╠═d49fd89a-79cd-4901-b4d8-e738c778e945
# ╠═6a39476c-0124-4290-82f1-ec51831c0cb7
# ╟─73144c07-6580-476e-b9a6-30fbd829ed58
# ╟─825ec223-06c2-431e-abad-037d4792fef6
# ╟─c970adc5-2b61-45c4-abf2-28764eb567cb
# ╠═32c8243d-d310-4b0a-b9b7-80b39e181dc7
# ╟─24ae360b-64b7-486d-8533-2ab0160ff4d3
# ╟─95286c76-a5b7-486b-b789-89df7cbc94d1
# ╠═63fed658-da7c-4fd2-9994-ee044f3b0a0e
# ╟─b0873543-404b-4665-a9ef-ab56f689fe64
# ╠═9b59b9af-2897-4b2f-8627-98e6cc53ea01
# ╠═04f8df8a-d760-47cd-af07-decc992d7c63
