### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 31bff5f4-9a71-11ee-1f13-b516023dc02d
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using MUST
	using Plots
	using Printf
	using DataFrames
	using DelimitedFiles 
	using LazySets
	using Polyhedra
end

# ╔═╡ 0f42da5a-b9ad-4307-bfec-4843e0692269
modelgrids = MUST.ingredients("modelgrids.jl")

# ╔═╡ fd69fd51-dcc7-4ea6-a892-e258017357f6
staggergrid = MUST.StaggerGrid("stagger_grid_full_o.mgrid")

# ╔═╡ 9f26c2b9-09e6-4cb9-9f2d-3458d02d1a49
md"What I want: I want to find a initial model, that is adiabatic in the interior and lies in between the the other grid nodes. Apparently just interpolating the models is giving you something that looks like interpolated, however the interior is not represented by an adiabat anymore."

# ╔═╡ c4e4285d-3024-4093-b261-ef67027906f9
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	scatter!(staggergrid["teff"], staggergrid["logg"])
	
	xflip!(true)
	yflip!(true)
end

# ╔═╡ f07ea04b-7e2e-4c48-a1fc-8189a272febc
md"Now we pick models in Stagger grid and try to find a model in the center."

# ╔═╡ a4dae5fa-4c7c-4d6b-a103-b63e83ca6c47
m1 = findfirst(staggergrid.info[!, "name"] .== "t60g40m00")

# ╔═╡ 5cc4258a-1c4e-4473-b8bb-d04cb87403c9
m2 = findfirst(staggergrid.info[!, "name"] .== "t60g45m00")

# ╔═╡ 72d441a2-f425-4d6e-9049-b597aa81857b
#mask = .![m in ["t55g45m00", "t60g45m00"] for m in staggergrid.info[!, "name"]]

# ╔═╡ 529bcb7b-d212-438a-9ebb-3b12e66a7717
#deleteat!(staggergrid.info, mask)

# ╔═╡ 1cebf806-ec11-43d8-b9bc-4176465ae7a7
md"What are those average models"

# ╔═╡ 26a85f8a-0e5f-41f3-b460-5409e87a79ca
model1 = TSO.Average3D(staggergrid.info[m1, "av_path"], logg=staggergrid.info[m1, "logg"])

# ╔═╡ b7ae42af-31a6-4eb4-bdb6-8c55aa172f10
model2 = TSO.Average3D(staggergrid.info[m2, "av_path"], logg=staggergrid.info[m2, "logg"])

# ╔═╡ 0e557f66-1488-485b-b483-38810913e2b5
"""
    optical_surface(model)

Return the z coordinate of the optical surface
"""
function optical_surface(model::TSO.OpticalModel1D)
    mask = sortperm(log10.(model.τ))
    TSO.linear_interpolation(
        TSO.Interpolations.deduplicate_knots!(log10.(model.τ[mask]), move_knots=true), 
        model.z[mask]
    )(0.0)
end

# ╔═╡ 140ac579-f9a5-4934-ad9e-b45cafa0ba7f


# ╔═╡ 96d5e25a-af45-4820-bbe6-3d76af86e078
md"And what are their adiabats"

# ╔═╡ 1d644b12-49e9-438b-b95a-05b14e0f7e99
eosE = reload(SqEoS, joinpath("DIS_MARCS_E_t60g45m00_v0.1", "eos.hdf5"))

# ╔═╡ a97d6e84-5c6b-400c-a46f-a54f9e1646b3
eos = reload(SqEoS, joinpath("DIS_MARCS_t60g45m00_v0.1", "eos.hdf5"))

# ╔═╡ 961141b9-84fe-481b-ba15-7c0bf9b0835d
model1optical = TSO.flip(@optical(model1, eos), depth=true)

# ╔═╡ 522a82da-f5ad-495d-913e-e2ac8441cbd6
z_test = TSO.rosseland_depth(eos, model1optical)

# ╔═╡ acc77539-20c3-4bf7-b94f-029e670437a9
begin
	model_constructed = deepcopy(model1)
	model_constructed.z .= z_test

	model_constructed = @optical model_constructed eos
end

# ╔═╡ 98e7ae6d-f32c-4d3a-b2c2-ef245fd5bb4f
let
	plot(model_constructed.z .- optical_surface(model_constructed), model_constructed.lnT)
	plot!(model1.z, model1.lnT)
end

# ╔═╡ c57b68e9-8edc-458b-be4c-751ec52617cb
adbt(model) = begin
	m = TSO.flip(model)
	s = TSO.pick_point(m, 1)
	e = TSO.pick_point(m, length(m.z))
	TSO.adiabat(
		s, e, eosE
	)
end

# ╔═╡ 723d5e33-b6f4-4ee5-b466-5b4830b68fdb
adiabat1 = adbt(model1)

# ╔═╡ 4ed1ca06-353c-44ae-8033-d1f0a6d2a736
adiabat2 = adbt(model2)

# ╔═╡ 9655ffae-1069-40e2-87f9-727bb0d22870
md"all other adiabats"

# ╔═╡ 481b99c7-90e5-44b9-aaa0-229d34e0f287
#=adiabats = [
	adbt(
		TSO.Average3D(staggergrid.info[i, "av_path"], logg=staggergrid.info[i, "logg"])
	)
	for i in 1:nrow(staggergrid.info)
]=#

# ╔═╡ d96e80c6-cfac-40dc-be03-2efca8ad7cb9
#=begin
	names = []
	for i in 1:nrow(staggergrid.info)
		name = "$(staggergrid.info[i, "name"])_ad.dat"
		open(name, "w") do f
			MUST.writedlm(f, [-reverse(adiabats[i].z) reverse(exp.(adiabats[i].lnT)) reverse(adiabats[i].lnρ)])
		end
		append!(names, [name])
	end

	staggergrid.info[!, "ad_path"] = names
end=#

# ╔═╡ 1e93c9e6-2c07-44b5-aa8a-525e5a9bcce1
#MUST.save(staggergrid, "stagger_grid_avad.mgrid")

# ╔═╡ 2113efb7-4dd1-407c-9b3b-28c2c8de8ff5
md"what would the interpolated model look like"

# ╔═╡ c77d6b1d-f1bf-4e0c-a7d3-de7079af1c83
m_ip = modelgrids.interpolate_average(
	staggergrid, eos,
	teff=6000.0, 
	logg=4.25, 
	feh=0.0,
)

# ╔═╡ 9aa74969-fa08-4575-8467-e5c2be915b0e


# ╔═╡ a89ff26f-22c9-44d5-9041-38f2ae698476
adiabat_m_ip = adbt(m_ip)

# ╔═╡ 9dad2a30-950f-447c-8c02-7d3a7f5282e1
adz1 = TSO.interpolate_to(TSO.flip(adiabat1), z=model1.z)

# ╔═╡ 4684f7c6-a30e-4435-b279-5997f4e3cfcb
adz2 = TSO.interpolate_to(TSO.flip(adiabat2), z=model2.z)

# ╔═╡ 9c43c0cc-556c-4fd9-9ca1-92018d7db988
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	plot!(model1.z, model1.lnT, label="t60g40m00", color=:red)
	plot!(model2.z, model2.lnT, label="t60g45m00", color=:black)

	plot!(-adiabat1.z, adiabat1.lnT, label="adiabat t60g45m00",color=:red, ls=:dash)
	plot!(-adiabat2.z, adiabat2.lnT, label="adiabat t55g45m00",color=:black, ls=:dash)

	plot!(-m_ip.z, m_ip.lnT, label="interpolated", color=:green)
	plot!(-adiabat_m_ip.z, adiabat_m_ip.lnT, label="adiabat interpolated", color=:green, ls=:dash)
	
end

# ╔═╡ 81ae7c00-ea49-4a34-b47b-82daa5bfab32
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	plot!(model1.z, model1.lnρ, label="t60g40m00", color=:red)
	plot!(model2.z, model2.lnρ, label="t60g45m00", color=:black)

	plot!(-adiabat1.z, adiabat1.lnρ, label="adiabat t60g45m00",color=:red, ls=:dash)
	plot!(-adiabat2.z, adiabat2.lnρ, label="adiabat t55g45m00",color=:black, ls=:dash)

	plot!(-m_ip.z, m_ip.lnρ, label="interpolated", color=:green)
	plot!(-adiabat_m_ip.z, adiabat_m_ip.lnρ, label="adiabat interpolated", color=:green, ls=:dash)
	
end

# ╔═╡ e6c255a3-0ed2-4d82-b4f6-ac7baf059d8d
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	plot!(model1.z, adz1.lnT, label="t60g40m00", color=:red)
	plot!(model2.z, adz2.lnT, label="t60g45m00", color=:black)

	#plot!(m_ip.z, m_ip.lnT, label="interpolated", color=:green)
	#plot!(-adiabat_m_ip.z, adiabat_m_ip.lnT, label="adiabat interpolated", color=:green, ls=:dash)
	
end

# ╔═╡ fcb94492-6366-48b7-8d1d-f24188f57557
let
	a1 = TSO.flip(adiabat1)
	a2 = TSO.flip(adiabat2)
	m1 = TSO.flip(model1)
	m2 = TSO.flip(model2)
	
	a1 = TSO.interpolate_to(a1, z=m1.z)
	a2 = TSO.interpolate_to(a2, z=m2.z)

	plot(m1.z, m1.lnT ./ a1.lnT)
	plot!(m2.z, m2.lnT ./ a2.lnT)
end

# ╔═╡ 4a184e23-a70b-4773-9416-0ae18ea0917f


# ╔═╡ 5451fc0c-b4fb-410b-9b85-89408808db85
md"# All boundary conditions"

# ╔═╡ cbac74b6-703f-434e-8619-65ddf7339d7d
allModels = [
	TSO.Average3D(staggergrid.info[i, "av_path"], logg=staggergrid.info[i, "logg"])
	for i in 1:nrow(staggergrid.info)
]

# ╔═╡ b68363dc-11a6-4b29-83ea-be11724e899b
dbot = [
	TSO.pick_point(m, length(m.z)).lnρ[1]
	for m in allModels
]

# ╔═╡ 0965b4bd-0817-4ca4-a1c3-438b3c85ca95
zbot = [
	TSO.pick_point(m, length(m.z)).z[1]
	for m in allModels
]

# ╔═╡ c6406993-54eb-4a9b-bb9f-1c2b755707ef
tbot = [
	TSO.pick_point(m, length(m.z)).lnT[1]
	for m in allModels
]

# ╔═╡ 6eb36fb3-544e-4625-87e4-d821448a5fe1
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	mask = staggergrid["feh"] .== 0.0
	scatter!(
		staggergrid["teff"][mask], dbot[mask], marker_z=staggergrid["logg"][mask], label=nothing
	)

	scatter!([6000.0], [maximum(m_ip.lnρ)], marker_z=m_ip.logg, label="interpolated")
	
	plot!(ylabel="dbot")
end

# ╔═╡ a53e3003-e05d-4ef7-b7c6-892f0d033309
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	mask = staggergrid["feh"] .== 0.0
	scatter!(
		staggergrid["teff"][mask], tbot[mask], marker_z=staggergrid["logg"][mask], label=nothing
	)

	scatter!([6000.0], [maximum(m_ip.lnT)], marker_z=m_ip.logg, label="interpolated")
	
	plot!(ylabel="tbot")
end

# ╔═╡ 70e269a2-66da-4712-ae9a-f4a8e599040a
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	mask = staggergrid["feh"] .== 0.0
	scatter!(
		staggergrid["teff"][mask], log10.(zbot[mask]), marker_z=staggergrid["logg"][mask], label=nothing
	)

	scatter!([6000.0], [log10.(maximum(m_ip.z))], marker_z=m_ip.logg, label="interpolated")
	
	plot!(ylabel="zbot")
end

# ╔═╡ ed69bf83-55ed-49bc-a875-30ccd1b3fd6a
md"# Different adiabats"

# ╔═╡ b02d330d-6427-4aac-a09a-e0a9238e4170
m_flipped = TSO.flip(m_ip)

# ╔═╡ 76c85fef-c966-42fe-8d8b-0b7b3874d476
drange = [sum([
	TSO.pick_point(model1, length(model1.z)).lnρ[1], 
	TSO.pick_point(model2, length(model2.z)).lnρ[1]
]) /2] #=range(
	sort([
		TSO.pick_point(model1, length(model1.z)).lnρ[1], 
		TSO.pick_point(model2, length(model2.z)).lnρ[1]
	])...,
	length=2
)=#

# ╔═╡ 10c3c868-af11-4236-a597-14162472a9d6
trange = [sum([
	TSO.pick_point(model1, length(model1.z)).lnT[1], 
	TSO.pick_point(model2, length(model2.z)).lnT[1]
]) /2] #=range(
	sort([
		TSO.pick_point(model1, length(model1.z)).lnT[1], 
		TSO.pick_point(model2, length(model2.z)).lnT[1]
	])...,
	length=2
)=#

# ╔═╡ 702ae97a-71af-46df-89ec-c1dc3d6f4c42
zrange = -exp.([sum([
	log(TSO.pick_point(model1, length(model1.z)).z[1]), 
	log(TSO.pick_point(model2, length(model2.z)).z[1])
]) /2]) #=range(
	sort([
		-TSO.pick_point(model1, length(model1.z)).z[1], 
		-TSO.pick_point(model2, length(model2.z)).z[1]
	])...,
	length=10
)=#

# ╔═╡ ac86085b-f4b6-4d69-ac1f-c3c439e40dfb
s = TSO.pick_point(m_flipped, 1)

# ╔═╡ 8e0eb790-1076-4290-911e-ae8483629e35
e = TSO.pick_point(m_flipped, length(m_flipped.z))

# ╔═╡ 5eabc126-0780-4546-8cb0-976355052cfb


# ╔═╡ 7ca20d09-90e3-431c-8001-7b0b40118567
md"Generate adiabats in between, and see if one of them looks good"

# ╔═╡ 95bea873-1121-4ff2-99aa-601f5c4e1c5e
begin
	adiabatsTest = []
	for (d, t, z) in Base.Iterators.product(drange, trange, zrange)
		@show d,t,z
		s.lnρ .= d
		s.lnT .= t
		s.z   .= z
		
		append!(adiabatsTest, [TSO.adiabat(
			s, e, eosE
		)])
	end
end

# ╔═╡ d6186e60-a5d4-4092-aa19-e438e6ddd0bf
begin
	labels = []
	for (d, t, z) in Base.Iterators.product(drange, trange, zrange)
		append!(labels, ["$d-$t-$z"])
	end
end

# ╔═╡ 0243ba50-582e-4d8d-9c5a-d477180abd36


# ╔═╡ f76287e4-cd03-49bc-a193-4e509b5dd065
md"How do the adiabats look like?"

# ╔═╡ ba765f8d-9eb4-42eb-92e9-565b36f664b0
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	plot!(model1.z, model1.lnρ, label="t60g40m00", color=:red)
	plot!(model2.z, model2.lnρ, label="t60g45m00", color=:black)

	plot!(-adiabat1.z, adiabat1.lnρ, label="adiabat t60g45m00",color=:red, ls=:dash)
	plot!(-adiabat2.z, adiabat2.lnρ, label="adiabat t55g45m00",color=:black, ls=:dash)

	plot!(m_ip.z, m_ip.lnρ, label="interpolated", color=:green, lw=3)

	for (i, a) in enumerate(adiabatsTest)
		plot!(-a.z, a.lnρ, label=nothing)
	end

	plot!(ylim=(-20, -10), xlim=(-2e8, 1e9))
end

# ╔═╡ cd3749bb-dd2d-4e26-959d-8750556b3e10
let
	plot(framestyle=:box, grid=false, legendforegroundcolor=nothing)
	
	plot!(model1.z, model1.lnT, label="t60g40m00", color=:red)
	plot!(model2.z, model2.lnT, label="t60g45m00", color=:black)

	plot!(-adiabat1.z, adiabat1.lnT, label="adiabat t60g45m00",color=:red, ls=:dash)
	plot!(-adiabat2.z, adiabat2.lnT, label="adiabat t55g45m00",color=:black, ls=:dash)

	plot!(m_ip.z, m_ip.lnT, label="interpolated", color=:green, lw=3)

	for (i, a) in enumerate(adiabatsTest)
		d, t, z = parse.(Float64, split(labels[i], "-", keepempty=false))
		#(-z!=zrange[2]) && continue
		plot!(-a.z, a.lnT, label=labels[i])
	end

	plot!(ylim=(8,10.1), xlim=(-2e8, 1e9))
end

# ╔═╡ Cell order:
# ╠═31bff5f4-9a71-11ee-1f13-b516023dc02d
# ╠═0f42da5a-b9ad-4307-bfec-4843e0692269
# ╠═fd69fd51-dcc7-4ea6-a892-e258017357f6
# ╟─9f26c2b9-09e6-4cb9-9f2d-3458d02d1a49
# ╟─c4e4285d-3024-4093-b261-ef67027906f9
# ╟─f07ea04b-7e2e-4c48-a1fc-8189a272febc
# ╠═a4dae5fa-4c7c-4d6b-a103-b63e83ca6c47
# ╠═5cc4258a-1c4e-4473-b8bb-d04cb87403c9
# ╠═72d441a2-f425-4d6e-9049-b597aa81857b
# ╠═529bcb7b-d212-438a-9ebb-3b12e66a7717
# ╟─1cebf806-ec11-43d8-b9bc-4176465ae7a7
# ╠═26a85f8a-0e5f-41f3-b460-5409e87a79ca
# ╠═b7ae42af-31a6-4eb4-bdb6-8c55aa172f10
# ╠═961141b9-84fe-481b-ba15-7c0bf9b0835d
# ╠═522a82da-f5ad-495d-913e-e2ac8441cbd6
# ╠═acc77539-20c3-4bf7-b94f-029e670437a9
# ╠═0e557f66-1488-485b-b483-38810913e2b5
# ╠═98e7ae6d-f32c-4d3a-b2c2-ef245fd5bb4f
# ╟─140ac579-f9a5-4934-ad9e-b45cafa0ba7f
# ╟─96d5e25a-af45-4820-bbe6-3d76af86e078
# ╠═1d644b12-49e9-438b-b95a-05b14e0f7e99
# ╠═a97d6e84-5c6b-400c-a46f-a54f9e1646b3
# ╠═c57b68e9-8edc-458b-be4c-751ec52617cb
# ╠═723d5e33-b6f4-4ee5-b466-5b4830b68fdb
# ╠═4ed1ca06-353c-44ae-8033-d1f0a6d2a736
# ╟─9655ffae-1069-40e2-87f9-727bb0d22870
# ╠═481b99c7-90e5-44b9-aaa0-229d34e0f287
# ╠═d96e80c6-cfac-40dc-be03-2efca8ad7cb9
# ╠═1e93c9e6-2c07-44b5-aa8a-525e5a9bcce1
# ╟─2113efb7-4dd1-407c-9b3b-28c2c8de8ff5
# ╠═c77d6b1d-f1bf-4e0c-a7d3-de7079af1c83
# ╟─9aa74969-fa08-4575-8467-e5c2be915b0e
# ╠═a89ff26f-22c9-44d5-9041-38f2ae698476
# ╠═9dad2a30-950f-447c-8c02-7d3a7f5282e1
# ╠═4684f7c6-a30e-4435-b279-5997f4e3cfcb
# ╠═9c43c0cc-556c-4fd9-9ca1-92018d7db988
# ╠═81ae7c00-ea49-4a34-b47b-82daa5bfab32
# ╠═e6c255a3-0ed2-4d82-b4f6-ac7baf059d8d
# ╠═fcb94492-6366-48b7-8d1d-f24188f57557
# ╟─4a184e23-a70b-4773-9416-0ae18ea0917f
# ╟─5451fc0c-b4fb-410b-9b85-89408808db85
# ╠═cbac74b6-703f-434e-8619-65ddf7339d7d
# ╠═b68363dc-11a6-4b29-83ea-be11724e899b
# ╠═0965b4bd-0817-4ca4-a1c3-438b3c85ca95
# ╠═c6406993-54eb-4a9b-bb9f-1c2b755707ef
# ╟─6eb36fb3-544e-4625-87e4-d821448a5fe1
# ╟─a53e3003-e05d-4ef7-b7c6-892f0d033309
# ╟─70e269a2-66da-4712-ae9a-f4a8e599040a
# ╟─ed69bf83-55ed-49bc-a875-30ccd1b3fd6a
# ╠═b02d330d-6427-4aac-a09a-e0a9238e4170
# ╠═76c85fef-c966-42fe-8d8b-0b7b3874d476
# ╠═10c3c868-af11-4236-a597-14162472a9d6
# ╠═702ae97a-71af-46df-89ec-c1dc3d6f4c42
# ╠═ac86085b-f4b6-4d69-ac1f-c3c439e40dfb
# ╠═8e0eb790-1076-4290-911e-ae8483629e35
# ╟─5eabc126-0780-4546-8cb0-976355052cfb
# ╟─7ca20d09-90e3-431c-8001-7b0b40118567
# ╠═95bea873-1121-4ff2-99aa-601f5c4e1c5e
# ╠═d6186e60-a5d4-4092-aa19-e438e6ddd0bf
# ╟─0243ba50-582e-4d8d-9c5a-d477180abd36
# ╟─f76287e4-cd03-49bc-a193-4e509b5dd065
# ╟─ba765f8d-9eb4-42eb-92e9-565b36f664b0
# ╠═cd3749bb-dd2d-4e26-959d-8750556b3e10
