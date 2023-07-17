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

# ╔═╡ dbf0121b-b717-4540-acc7-67a1cb7b3f94
begin
	Plots.default(fontfamily = ("Courier"), titlefont=("Courier"))
	
	function copy_ticks(sp::Plots.Subplot; minorticks=10)
		ptx = twinx(sp)
		plot!(ptx,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
		
		pty = twiny(sp)
		plot!(pty,
				xlims=xlims(sp),
				ylims=ylims(sp),
				xformatter=_->"",
				yformatter=_->"", 
				minorticks=minorticks)
	end
	
	copy_ticks(plt::Plots.Plot = current(); minorticks=10) = copy_ticks(plt[1], minorticks=minorticks)
	
	
	function basic_plot!(plot::Plots.Plot=current(); minorticks=10, 
						tickfontsize=8, 
						legendfontsize=10, 
						guidefontsize=12, 
						size=(600,400), 
						lm=2, bm=2, tm=2, rm=2, copyticks=true)
		plot!(plot, formatter=:auto)
		
		if copyticks
			copy_ticks(plot, minorticks=15)
		end
		
		plot!(plot, framestyle=:box, 
				minorticks=minorticks, 
				tickfontsize=tickfontsize, 
				legendfontsize=legendfontsize, 
				titlefontsize=legendfontsize,
				grid=false, 
				format=false,
				foreground_color_legend=nothing, 
				size=size, 
				leftmargin=lm*Plots.mm, rightmargin=rm*Plots.mm, 
				bottommargin=bm*Plots.mm, topmargin=tm*Plots.mm, 
				guidefontsize=guidefontsize)

		plot!(
			legendforegroundcolor=nothing, 
			legendbackgroundcolor=nothing
		)
		plot!()
	end;
end

# ╔═╡ 90818cad-34c7-4ad0-8a8a-530391a3e176
md"## The Namelists
Similar to the automatic execution of the stellar atmospheres code, the Multi code can be executed within a slurm environment. The overall code is the same, so first a namelist has to be created. One can choose a native default namelist, or just construct it from scratch. All parameters that are needed are available as fields. Note that this strictly speaking is not needed, but can easier to fit into the rest of the existing code."

# ╔═╡ b7d1105d-3a57-4425-83dd-20fa03cd4194
md"A couple of defaults however are usefull when it comes to compute the effective temperature of a model."

# ╔═╡ b8f2e1cf-8d1c-402e-a4ac-347cf60612d6
nml = MUST.spectrum_namelist(
	"m3dis_sun_HD90-RT180", 
	model_folder="./input_multi3d/from_stagger"
)

# ╔═╡ 14f77786-6eb4-49fa-9f92-0df3121478c8
nml

# ╔═╡ 0ab5616f-3b2d-48a4-a7ac-eebb0632fc7c
md"## The Submission
The following convenient shortcut will execute M3D and return the output in an convenient format for us to look at. This actually runs if you are connected to a cluster using the slurm system and it will take some time."

# ╔═╡ ec0efed0-6565-4b95-a82a-84d9662b6f5c
#=hd120_rt240_old = MUST.spectrum(
	"m3dis_sun_HD120-RT240-old", 
	namelist_kwargs=(:model_folder=>"./input_multi3d/from_stagger",
					 :atmos_params=>(:atmos_format=>"MUST",))
)=#

# ╔═╡ c9dc6751-4287-4962-8e4b-73e644548346
md"## The Output
If the code is executed, it will return the result directly. However, it is also possible to load a precomputed setup just like in the m3dis python package."

# ╔═╡ 179dc012-ed76-419d-8025-40b7c5f07508
stagger = MUST.M3DISRun("data/t5777g44m0005_20.5x5x230")

# ╔═╡ 71e71f24-348a-42dc-a639-ecebf78326ea
rt180 = MUST.M3DISRun("data/m3dis_sun_HD90-RT180")

# ╔═╡ 40e2bc60-f1a8-4ded-96b1-a3da9b79e6c0
rt360 = MUST.M3DISRun("data/m3dis_sun_HD90-RT360")

# ╔═╡ 60c391f3-7b4d-4150-824e-6ed218bf52e0
hd120_rt240_old = MUST.M3DISRun("data/m3dis_sun_HD120-RT240-old")

# ╔═╡ 73144c07-6580-476e-b9a6-30fbd829ed58
begin	
	# which line should be shown
	il = 8
	norm = true
	
	plot(
		pc(stagger.line[il].crop(norm=norm, LTE=true))..., 
		label="Stagger", 
		color=:black
	)
	plot!(
		pc(rt180.line[il].crop(norm=norm, LTE=true))..., 
		label="HD90 - RT180", 
		color=:red
	)
	plot!(
		pc.(rt360.line[il].crop(norm=norm))..., 
		label="HD90 - RT360", 
		color=:blue
	)
	#=plot!(
		pc.(hd120_rt240_old.line[il].crop(norm=norm))..., 
		label="HD120 - RT240 (old)", 
		color=:grey
	)=#

	plot!(xlim=(6545, 6580))
	plot!(xlabel="λ [Å]", ylabel="flux")
	basic_plot!(size=(600,600))

	plot!(legendposition=:bottomleft, legendcolumn=1)
end

# ╔═╡ Cell order:
# ╟─bb45e6a0-2b8b-4731-bd53-190c44a5ea23
# ╠═49a4e636-166c-11ee-03ed-c129e8cdc7fc
# ╠═a6d00a44-bb7d-4a40-9b45-c40a6d7e1f1b
# ╠═21951d47-3aef-4297-9edd-cf7afcb4544f
# ╠═4bb87ab7-94fa-460d-9c7b-4b7fe06c8500
# ╟─dbf0121b-b717-4540-acc7-67a1cb7b3f94
# ╟─90818cad-34c7-4ad0-8a8a-530391a3e176
# ╟─b7d1105d-3a57-4425-83dd-20fa03cd4194
# ╠═b8f2e1cf-8d1c-402e-a4ac-347cf60612d6
# ╠═14f77786-6eb4-49fa-9f92-0df3121478c8
# ╟─0ab5616f-3b2d-48a4-a7ac-eebb0632fc7c
# ╠═ec0efed0-6565-4b95-a82a-84d9662b6f5c
# ╟─c9dc6751-4287-4962-8e4b-73e644548346
# ╠═179dc012-ed76-419d-8025-40b7c5f07508
# ╠═71e71f24-348a-42dc-a639-ecebf78326ea
# ╠═40e2bc60-f1a8-4ded-96b1-a3da9b79e6c0
# ╠═60c391f3-7b4d-4150-824e-6ed218bf52e0
# ╠═73144c07-6580-476e-b9a6-30fbd829ed58
