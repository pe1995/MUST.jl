### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ b95854ea-f871-11ed-3e17-4101e15fb286
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using Glob
	using PyPlot
	using TSO
	using LaTeXStrings
	using Printf
	using KernelDensity
	using DelimitedFiles
end

# ╔═╡ def9d891-401b-44ec-a38a-fbd06c240e83
include_helper(name) = include(joinpath(dirname(pathof(MUST)), name))

# ╔═╡ f6b74677-55e2-4828-9a17-9430ec234eb3
md"# Paper I: Validation & the Sun"

# ╔═╡ 6fda69d3-ffb1-45c3-8693-fb550ae8cbe6
md"## Code Setup"

# ╔═╡ 2098b0b7-386e-4cc4-a2d6-1d9c953ba5d1
md"### Plotting defaults"

# ╔═╡ 21e9bc39-68b6-4eb3-be42-68b49ef82ac4
begin
	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["font.size"] = 10
end

# ╔═╡ 888e3baa-a039-4ad3-a4a5-8478bd9acd3d
include_helper("visual.jl")

# ╔═╡ fe3ef038-a656-4478-b399-064a53ccee38
md"### Dispatch module"

# ╔═╡ e4fb03de-b118-4c33-bfb2-fa34f7277767
mean = MUST.mean

# ╔═╡ 095fd801-1365-431f-bc47-4aab3de3cddc
rms(x) = √(sum(x .^2) / length(x))

# ╔═╡ 4eec9985-a3a5-4173-b445-310c67e39f38
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"

# ╔═╡ e3597171-bd01-41f5-ab38-59bc0482b9a1
mesh(m::MUST.Box) = MUST.meshgrid(MUST.axis(m, :x) ./1e8, 
									MUST.axis(m, :y) ./1e8, 
									MUST.axis(m, :z) ./1e8)

# ╔═╡ 0841981a-a16c-4b04-8b92-3ed1f9ac1f93
is_log(x) = begin
	sx = String(x)
	
	xnew, f = if occursin("log10", sx)
		Symbol(sx[findfirst("log10", sx)[end]+1:end]), log10
	elseif occursin("log", sx)
		Symbol(sx[findfirst("log", sx)[end]+1:end]), log
	else
		x, identity
	end

	xnew, f
end

# ╔═╡ 844889f4-1480-40cd-9969-4aefb0ab6f98
profile(f, model, x=:z, y=:T) = begin
	xs, logx = is_log(x)
	ys, logy = is_log(y)
	
	if xs == :τ_ross
		logx.(MUST.axis(model, xs, 3)), logy.(MUST.plane_statistic(f, model, ys)) 
	else
		logx.(MUST.axis(model, xs)), logy.(MUST.plane_statistic(f, model, ys))
	end
end

# ╔═╡ 3d4d31d0-279c-4917-80bc-a9566a20431e
md"## General Models"

# ╔═╡ e1ef2d97-6730-4e4f-ad58-beaafeca63a4
md"Solar stagger model (converted to optical depth with current EoS."

# ╔═╡ 7594d84e-1393-4842-88d6-319b799608e1
md"Stagger"

# ╔═╡ 23d83ebb-cd4c-4f58-a8a1-7cf5b9958b51
begin
	folder_stagger = "/ptmp/peitner/model_grid/MUST.jl/examples/stagger2bifrost"
	stagger = MUST.Box("box_solar_stagger_MARCS_v1.4.31", 
		folder=folder_stagger)
	stagger_τ = MUST.Box("box_solar_stagger_MARCS_v1.4.31_t",
		folder=folder_stagger)	
end

# ╔═╡ f81e64d8-7931-4005-bc34-5a6022348015
md"Muram"

# ╔═╡ 92df7339-dc9d-415e-b281-6dbb84da40b1
begin
	folder_muram = "/u/peitner/DISPATCH/examples/Muram"
	muram = MUST.Box("box_MURaM_cube_small.221000_HDSun", folder=folder_muram)
	muram_τ = MUST.Box("box_MURaM_cube_small.221000_HDSun_t", folder=folder_muram)
end

# ╔═╡ 7e30c19d-2092-4042-a3a9-e09b01a6f75b
keys(muram.data)

# ╔═╡ ca5346fd-5e7f-4817-91ab-7509c9eae37d
#muram_τ = MUST.height_scale(muram, :τ_ross)

# ╔═╡ c6e0d3f7-163b-45ca-8f5c-052c3beb2b60
#MUST.save(muram_τ, folder=folder_muram, name="box_MURaM_cube_small.221000_HDSun_t")

# ╔═╡ 13935fec-b609-468e-a623-79fbdf1dc1a6
md"## A) Resolution Comparison"

# ╔═╡ 5f2b25d3-4738-40b2-b934-1814e58b6423
md"### Data"

# ╔═╡ 133a2fb5-80d9-4ee5-a61d-ab34027264ac
names_res = [
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1",
	"DIS_MARCS_E_t5777g44m00_v0.1"
]

# ╔═╡ c79686f9-155f-4561-99da-1ab6229e7cf8
out_folder_res = [
	"models/sun_5bins_lres",
	"models/sun_5bins_ires",
	"models/sun_5bins_hres",
	"models/sun_5bins_rtres"
]

# ╔═╡ 0912c3c4-2fa4-4f16-82aa-d6582028ba26
in_folder_res = [MUST.@in_dispatch "input_data/grd/$(name)" for name in names_res]

# ╔═╡ 057b35ed-3a84-44a6-9332-f5632312ad4d
eos_folder_res = [
	MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35"),
	MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35"),
	MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35"),
	MUST.@in_dispatch("input_data/DIS_MARCS_E_v1.4.35")
]

# ╔═╡ a6668016-fbb0-43de-be05-879bb4dd4b26
colors_res = ["limegreen", "tomato", "turquoise", "violet"]

# ╔═╡ 54450d4a-90ef-4089-a73e-9ad01a19b813
ls_res = ["-", "-", "-", "-", "-"]

# ╔═╡ 908ba257-018e-4cd3-afe2-08498d0be47e
begin
	snapshots_res = []
	snapshots_τ_res = []
	
	for i in eachindex(names_res)
		snapshot, snapshot_τ = pick_snapshot(out_folder_res[i], :recent)
		
		append!(snapshots_res, [snapshot])
		append!(snapshots_τ_res, [snapshot_τ])
	end
end

# ╔═╡ 3e7ea0ef-f7b3-4a2b-887b-29c2ba331505
begin
	eos_res = [reload(SqEoS, joinpath(eos_folder_res[i], "eos.hdf5")) 
				for i in eachindex(eos_folder_res)]
	opa_res = [reload(SqOpacity, joinpath(eos_folder_res[i], "binned_opacities.hdf5"))
				for i in eachindex(eos_folder_res)]
end

# ╔═╡ f30715f9-918a-4e05-b2ce-50bd7df9d223
for i in eachindex(eos_res)
	@info "Opacity table size: $(size(opa_res[i].κ))"
end

# ╔═╡ fb73d428-e189-41db-861f-c0ea17217758
initial_model = [
	@optical(Average3D(eos_res[i], joinpath(in_folder_res[i], "inim.dat")), 
						eos_res[i], opa_res[i])
		for i in eachindex(eos_res)
]

# ╔═╡ ffa7c65d-0df5-412a-a3e9-0436dee447c2
resolution(snap) = @sprintf "%.1f km" first(diff(MUST.axis(snap, :z) ./1e5))

# ╔═╡ f53c1ebb-dbed-412e-803a-d010223c934b
for (i,snap) in enumerate(snapshots_res)
	@info "Resolution $(out_folder_res[i]): $(resolution(snap))"
end

# ╔═╡ aa03989d-7e6e-452d-86f3-ddc9c25c6e41
labels_res = ["M3DIS - $(resolution(snap))" for snap in snapshots_res]

# ╔═╡ 9b59e8f3-df44-43a5-940a-528215b8fc94
md"### Figures: Compare different resolutions"

# ╔═╡ 44d2a5ea-89ee-453b-8028-f4d093b6f93a
begin
	close()
	
	f_res1, ax_res1 = plt.subplots(3, 1, figsize=(5,8), sharex=true)
	basic_plot!.(ax_res1)
	plt.subplots_adjust(wspace=0, hspace=0)

	rms5(x) = rms(x) ./ 1e5

	ls_res_2 = [":", "--", "-", "-.", ":"]

	# Stagger model
	#=ax_res1[1].plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
		color="k", 
		ls="--",
		lw=1.5,
		label="Stagger")
	ax_res1[2].plot(profile(mean, stagger_τ, :log10τ_ross, :log10d)..., 
		color="k", 
		ls="--",
		lw=1.5,
		label="Stagger")
	ax_res1[3].plot(profile(rms5, stagger_τ, :log10τ_ross, :uz)..., 
		color="k", 
		ls="--",
		lw=1.5,
		label="Stagger")=#

	
	# DISPATCH models
	for (i, snap) in enumerate(snapshots_τ_res)
		ax_res1[1].plot(profile(mean, snap, :log10τ_ross, :T)..., 
			color="k", 
			label=labels_res[i],
			ls=ls_res_2[i],
			lw=1.4)
		ax_res1[2].plot(profile(mean, snap, :log10τ_ross, :log10d)..., 
			color="k", 
			label=labels_res[i],
			ls=ls_res_2[i],
			lw=1.4)
		ax_res1[3].plot(profile(rms5, snap, :log10τ_ross, :uz)..., 
			color="k", 
			label=labels_res[i],
			ls=ls_res_2[i],
			lw=1.4)
	end


	# Figure props
	ax_res1[3].set_xlabel(L"\rm \log\ \tau_{ross}")
	ax_res1[1].set_ylabel(L"\rm T\ [K]")
	ax_res1[2].set_ylabel(L"\rm \log\ \rho\ [g \times cm^{-3}]")
	ax_res1[3].set_ylabel(L"\rm \sqrt{< U_z ^2 >}\ [km \times s^{-1}]")
	
	ax_res1[1].legend(framealpha=0)
	
	ax_res1[1].set_xlim(-4, 2)
	ax_res1[1].set_ylim(3900, 11000)
	ax_res1[2].set_ylim(-9, -6.25)

	f_res1.savefig("res_average_comp_zoom.pdf", bbox_inches="tight")

	gcf()
end

# ╔═╡ 4735928b-6a80-4915-b23f-8852251cf77f
begin
	uz_res = []
	for (i,snapshot_τ) in enumerate(snapshots_τ_res)
		isurf = MUST.closest(log10.(MUST.axis(snapshot_τ, :τ_ross, 3)), 0)
		uz_τ_surf = snapshot_τ[:uz][:, :, isurf] #./1e5
		append!(uz_res, [uz_τ_surf])
	end

	uz_res
end

# ╔═╡ 66862a2b-b0cd-485a-9dc4-d74b20a1fe26
begin
	close()
	
	f_res2, ax_res2 = plt.subplots(1, length(names_res), 
			figsize=(4*length(names_res),4), sharey=true, sharex=true)
	
	basic_plot!.(ax_res2)
	plt.subplots_adjust(wspace=0.05)

	extent(snap) = begin
		x = MUST.axis(snap, :x) ./1e8
		y = MUST.axis(snap, :x) ./1e8

		[minimum(x), maximum(x), minimum(y), maximum(y)]
	end

	vmin = minimum([minimum(u) for u in uz_res])
	vmax = maximum([maximum(u) for u in uz_res])
	
	
	# DISPATCH models
	im = []
	for (i, snap) in enumerate(snapshots_τ_res)
		imi = ax_res2[i].imshow(uz_res[i], 
			origin="lower", extent=extent(snap), 
			vmin=vmin, vmax=vmax,
			cmap="hot")

		ax_res2[i].set_title(labels_res[i])
		ax_res2[i].set_xlabel("x [Mm]")
		
		append!(im, [imi])
	end


	# Figure props
	c = f_res2.colorbar(im[1], ax=ax_res2)
	c.set_label(L"\rm U_z\ [km \times s^{-1}]")
	ax_res2[1].set_ylabel("y [Mm]")


	f_res2.savefig("res_optical_surface_uz.pdf", bbox_inches="tight")
	
	gcf()
end

# ╔═╡ bb3695ef-eaac-4953-b129-1e219d5c2606
md"## B) Full 3D structure
From the resolution tests, we pick the best model and have a look at the full 3D structure."

# ╔═╡ edbc29f7-7310-4cb9-ae36-7ed64125a3a6
i_best = 4

# ╔═╡ e16a3f22-b778-4467-860f-0022648672bb
begin
	close()
	
	#ft, axt = cube_with_velocities(snapshots_res[i_best], vmax_3d=15500)

	gcf()
end

# ╔═╡ e62ef28a-c05e-4842-b64b-47503ec49c43
md"## C) Comparison with other models"

# ╔═╡ 9cdb8768-798a-449b-93ea-e1cf662ccd51
function get_density(xin, yin)
	x_d = reshape(xin, :)
	y_d = reshape(yin, :)
	
	d = kde((x_d, y_d))
	x = range(minimum(x_d), maximum(x_d), length=250)
	y = range(minimum(y_d), maximum(y_d), length=250)
	
	xx, yy = MUST.meshgrid(x, y)
	xx, yy = reshape(xx, :), reshape(yy, :)
	ip = InterpKDE(d)
	dd = similar(xx)
	for i in eachindex(xx)
		dd[i] = pdf(ip, xx[i], yy[i])
	end

	#dd .+= 1.0
	#dd .= log.(dd ./ maximum(dd))
	dd[dd .< 1e-12] .= 0.0

	dd .= sqrt.(dd) 
	dd ./= maximum(dd)
	xx, yy, dd
end

# ╔═╡ 2f75525a-3049-424e-9bf7-842d640a77ee
# ╠═╡ disabled = true
#=╠═╡
xx_T, yy_T, dd_T = get_density(
	log10.(snapshots_res[i_best][:τ_ross]),
	snapshots_res[i_best][:T]
)
  ╠═╡ =#

# ╔═╡ 83d82b5f-e132-4625-817e-b44e0cb6d3b8
# ╠═╡ disabled = true
#=╠═╡
xx_d, yy_d, dd_d = get_density(
	log10.(snapshots_res[i_best][:τ_ross]),
	log10.(snapshots_res[i_best][:d])
)
  ╠═╡ =#

# ╔═╡ e5d2bb57-c85a-4f41-b868-67c1c80af16e
# ╠═╡ disabled = true
#=╠═╡
xx_u, yy_u, dd_u = get_density(
	log10.(snapshots_res[i_best][:τ_ross]),
	abs.(snapshots_res[i_best][:uz] ./1e5)
)
  ╠═╡ =#

# ╔═╡ 799c9bb0-7f9d-4b3c-9c56-94789108f194
begin
	close()
	
	f_others1, ax_others1 = plt.subplots(3, 2, figsize=(10,8), sharex=true)
	basic_plot!.(ax_others1)
	plt.subplots_adjust(wspace=0.15, hspace=0)

	


	# Stagger model
	ax_others1[1, 1].plot(profile(mean, stagger_τ, :log10τ_ross, :T)..., 
		color="cornflowerblue", 
		ls="-",
		lw=2.0,
		label="Stagger")
	ax_others1[2, 1].plot(profile(mean, stagger_τ, :log10τ_ross, :log10d)..., 
		color="cornflowerblue", 
		ls="-",
		lw=2.0,
		label="Stagger")
	ax_others1[3, 1].plot(profile(rms5, stagger_τ, :log10τ_ross, :uz)..., 
		color="cornflowerblue", 
		ls="-",
		lw=2.0,
		label="Stagger")


	# Muram model
	ax_others1[1, 1].plot(profile(mean, muram_τ, :log10τ_ross, :T)..., 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM")
	ax_others1[2, 1].plot(profile(mean, muram_τ, :log10τ_ross, :log10d)..., 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM")
	ax_others1[3, 1].plot(profile(rms5, muram_τ, :log10τ_ross, :uz)..., 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM")


	# DISPATCH models
	#ax_others1[1, 1].scatter(xx_T, yy_T, c=dd_T, cmap="Blues")
	ax_others1[1, 1].plot(profile(mean, 
						snapshots_τ_res[i_best], :log10τ_ross, :T)..., 
			color="k", 
			label="M3DIS",
			ls="--",
			lw=1.7)
	
	#ax_others1[2, 1].scatter(xx_d, yy_d, c=dd_d, cmap="Blues")
	ax_others1[2, 1].plot(profile(mean, 
						snapshots_τ_res[i_best], :log10τ_ross, :log10d)..., 
			color="k", 
			ls="--",
			label="M3DIS",
			lw=1.7)

	#ax_others1[3].scatter(xx_u, yy_u, c=dd_u, cmap="Blues")
	ax_others1[3, 1].plot(profile(rms5, 
						snapshots_τ_res[i_best], :log10τ_ross, :uz)..., 
			color="k", 
			ls="--",
			label="M3DIS",
			lw=1.7)
	

	# Differences
	common_tau = range(-4, 3, length=200) |> collect
	ip(x, y) = begin
		m = sortperm(x)
		TSO.linear_interpolation(x[m], y[m]).(common_tau)
	end


	# DISPATCH - Stagger
	ax_others1[1, 2].plot(common_tau,
		ip(profile(mean, snapshots_τ_res[i_best], :log10τ_ross, :T)...) .- 
		ip(profile(mean, stagger_τ, :log10τ_ross, :T)...), 
		color="k", 
		ls="--",
		label="M3DIS - Stagger",
		lw=1.7)
	ax_others1[2, 2].plot(common_tau,
		ip(profile(mean, snapshots_τ_res[i_best], :log10τ_ross, :log10d)...) .- 
		ip(profile(mean, stagger_τ, :log10τ_ross, :log10d)...), 
		color="k", 
		ls="--",
		label="M3DIS - Stagger",
		lw=1.7)
	ax_others1[3, 2].plot(common_tau,
		ip(profile(rms5, snapshots_τ_res[i_best], :log10τ_ross, :uz)...) .- 
		ip(profile(rms5, stagger_τ, :log10τ_ross, :uz)...), 
		color="k", 
		ls="--",
		label="M3DIS",
		lw=1.7)


	# DISPATCH - Muram
	ax_others1[1, 2].plot(common_tau,
		ip(profile(mean, muram_τ, :log10τ_ross, :T)...) .- 
		ip(profile(mean, stagger_τ, :log10τ_ross, :T)...), 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM - Stagger")
	ax_others1[2, 2].plot(common_tau,
		ip(profile(mean, muram_τ, :log10τ_ross, :log10d)...) .- 
		ip(profile(mean, stagger_τ, :log10τ_ross, :log10d)...), 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM - Stagger")
	ax_others1[3, 2].plot(common_tau,
		ip(profile(rms5, muram_τ, :log10τ_ross, :uz)...) .- 
		ip(profile(rms5, stagger_τ, :log10τ_ross, :uz)...), 
		color="k", 
		ls=":",
		lw=1.7,
		label="MURaM")

	ax_others1[1, 2].axhline(0.0, color="k", ls=":", alpha=0.5)
	ax_others1[2, 2].axhline(0.0, color="k", ls=":", alpha=0.5)
	ax_others1[3, 2].axhline(0.0, color="k", ls=":", alpha=0.5)


	

	# Figure props
	ax_others1[3, 1].set_xlabel(L"\rm \log\ \tau_{ross}")
	ax_others1[3, 2].set_xlabel(L"\rm \log\ \tau_{ross}")
	ax_others1[1, 1].set_ylabel(L"\rm T\ [K]")
	ax_others1[2, 1].set_ylabel(L"\rm \log\ \rho\ [g \times cm^{-3}]")
	ax_others1[3, 1].set_ylabel(L"\rm \sqrt{< U_z ^2 >}\ [km \times s^{-1}]")
	
	ax_others1[2, 1].legend(framealpha=0, loc="lower center", ncol=3)
	ax_others1[2, 2].legend(framealpha=0, loc="lower center", ncol=3)
	
	
	ax_others1[1, 1].set_xlim(-4, 3)
	ax_others1[1, 1].set_ylim(3200, 12000)
	ax_others1[2, 1].set_ylim(-9.1, -5.75)
	ax_others1[3, 1].set_ylim(-0.5, 4.15)

	ax_others1[1, 2].set_ylim(-220, 220)
	ax_others1[2, 2].set_ylim(-0.065, 0.065)
	ax_others1[3, 2].set_ylim(-0.65, 0.65)
	

	f_others1.savefig("others_average_comp.png", bbox_inches="tight", dpi=400)
	f_others1.savefig("others_average_comp.pdf", bbox_inches="tight")

	gcf()
end

# ╔═╡ 17969398-3ba2-4f41-b98c-f39a04c2cdec
md"## Extra) Convert to M3D format
### Models from A)"

# ╔═╡ cf368d28-1b44-4bb4-b34c-5029187b7b40
labels_new_res = [last(split(f, "/")) for f in deepcopy(out_folder_res)]

# ╔═╡ 3128bf5c-151d-4f54-8a38-544833c65a37
downsample_res = [10, 20, 50, 20]

# ╔═╡ 8597a8f6-1fde-4bd5-867b-10f1a6cf20f4
output_names_res = [
	"m3dis_$(label)" 
	for (name, label) in zip(names_res, labels_new_res)
]

# ╔═╡ 508bbc3f-0cbb-49f5-809c-6335fa5e1252
aos_res = [@axed(eos_res[i]) for i in eachindex(eos_res)]

# ╔═╡ 050e49df-d5e4-47f1-8311-714924488bce
ne_res = [lookup(aos_res[i], :lnNe, log.(model_box[:d]), log.(model_box[:ee])) 
			for (i, model_box) in enumerate(snapshots_res)]

# ╔═╡ ffa3bb64-2e9e-4e89-9bc6-862a26e833af
function save_tau(model, path)
	x, y = profile(mean, model, :log10τ_ross, :T)
	x2, y2 = profile(mean, model, :log10τ_ross, :d)

	@assert x==x2
	
	open(path; write=true) do f
    	write(f, "# tau_ross,T,d\n")
        writedlm(f, [x y y2], ',')
    end
end

# ╔═╡ 01a11645-2464-4e4c-8ac4-e9adc40c9a52
function save_z(model, path)
	x, y = profile(mean, model, :z, :T)
	x2, y2 = profile(mean, model, :z, :d)

	@assert x==x2
	
	
	open(path; write=true) do f
    	write(f, "# z,T,d\n")
        writedlm(f, [x y y2], ',')
    end
end

# ╔═╡ 41c5b6b7-f168-4137-ae48-2fdce6b3e3ca
for (i, model_box) in enumerate(snapshots_res)
	model_box.data[:ne] = exp.(ne_res[i])
	
	MUST.multiBox(model_box, joinpath(out_folder_res[i], output_names_res[i]), 
						downsample_xy=downsample_res[i])

	save_tau(snapshots_τ_res[i], 
		joinpath(out_folder_res[i], output_names_res[i]*"_av_tau.txt"))
	save_z(snapshots_res[i], 
		joinpath(out_folder_res[i], output_names_res[i]*"_av_z.txt"))
	
	@info "New size ($(output_names_res[i])): $(size(
		@view(ne_res[i][1:downsample_res[i]:end, 1:downsample_res[i]:end, :]))
	)"
end

# ╔═╡ Cell order:
# ╟─f6b74677-55e2-4828-9a17-9430ec234eb3
# ╟─6fda69d3-ffb1-45c3-8693-fb550ae8cbe6
# ╠═b95854ea-f871-11ed-3e17-4101e15fb286
# ╟─2098b0b7-386e-4cc4-a2d6-1d9c953ba5d1
# ╠═21e9bc39-68b6-4eb3-be42-68b49ef82ac4
# ╠═def9d891-401b-44ec-a38a-fbd06c240e83
# ╠═888e3baa-a039-4ad3-a4a5-8478bd9acd3d
# ╟─fe3ef038-a656-4478-b399-064a53ccee38
# ╟─e4fb03de-b118-4c33-bfb2-fa34f7277767
# ╠═095fd801-1365-431f-bc47-4aab3de3cddc
# ╠═4eec9985-a3a5-4173-b445-310c67e39f38
# ╟─e3597171-bd01-41f5-ab38-59bc0482b9a1
# ╟─0841981a-a16c-4b04-8b92-3ed1f9ac1f93
# ╟─844889f4-1480-40cd-9969-4aefb0ab6f98
# ╟─3d4d31d0-279c-4917-80bc-a9566a20431e
# ╟─e1ef2d97-6730-4e4f-ad58-beaafeca63a4
# ╟─7594d84e-1393-4842-88d6-319b799608e1
# ╠═23d83ebb-cd4c-4f58-a8a1-7cf5b9958b51
# ╟─f81e64d8-7931-4005-bc34-5a6022348015
# ╠═92df7339-dc9d-415e-b281-6dbb84da40b1
# ╠═7e30c19d-2092-4042-a3a9-e09b01a6f75b
# ╠═ca5346fd-5e7f-4817-91ab-7509c9eae37d
# ╠═c6e0d3f7-163b-45ca-8f5c-052c3beb2b60
# ╟─13935fec-b609-468e-a623-79fbdf1dc1a6
# ╟─5f2b25d3-4738-40b2-b934-1814e58b6423
# ╠═133a2fb5-80d9-4ee5-a61d-ab34027264ac
# ╠═c79686f9-155f-4561-99da-1ab6229e7cf8
# ╠═0912c3c4-2fa4-4f16-82aa-d6582028ba26
# ╠═057b35ed-3a84-44a6-9332-f5632312ad4d
# ╠═a6668016-fbb0-43de-be05-879bb4dd4b26
# ╠═54450d4a-90ef-4089-a73e-9ad01a19b813
# ╟─908ba257-018e-4cd3-afe2-08498d0be47e
# ╟─3e7ea0ef-f7b3-4a2b-887b-29c2ba331505
# ╟─f30715f9-918a-4e05-b2ce-50bd7df9d223
# ╟─fb73d428-e189-41db-861f-c0ea17217758
# ╠═ffa7c65d-0df5-412a-a3e9-0436dee447c2
# ╟─f53c1ebb-dbed-412e-803a-d010223c934b
# ╠═aa03989d-7e6e-452d-86f3-ddc9c25c6e41
# ╟─9b59e8f3-df44-43a5-940a-528215b8fc94
# ╟─44d2a5ea-89ee-453b-8028-f4d093b6f93a
# ╟─4735928b-6a80-4915-b23f-8852251cf77f
# ╟─66862a2b-b0cd-485a-9dc4-d74b20a1fe26
# ╟─bb3695ef-eaac-4953-b129-1e219d5c2606
# ╠═edbc29f7-7310-4cb9-ae36-7ed64125a3a6
# ╠═e16a3f22-b778-4467-860f-0022648672bb
# ╟─e62ef28a-c05e-4842-b64b-47503ec49c43
# ╟─9cdb8768-798a-449b-93ea-e1cf662ccd51
# ╠═2f75525a-3049-424e-9bf7-842d640a77ee
# ╠═83d82b5f-e132-4625-817e-b44e0cb6d3b8
# ╠═e5d2bb57-c85a-4f41-b868-67c1c80af16e
# ╠═799c9bb0-7f9d-4b3c-9c56-94789108f194
# ╟─17969398-3ba2-4f41-b98c-f39a04c2cdec
# ╠═cf368d28-1b44-4bb4-b34c-5029187b7b40
# ╠═3128bf5c-151d-4f54-8a38-544833c65a37
# ╠═8597a8f6-1fde-4bd5-867b-10f1a6cf20f4
# ╠═508bbc3f-0cbb-49f5-809c-6335fa5e1252
# ╠═050e49df-d5e4-47f1-8311-714924488bce
# ╠═ffa3bb64-2e9e-4e89-9bc6-862a26e833af
# ╠═01a11645-2464-4e4c-8ac4-e9adc40c9a52
# ╠═41c5b6b7-f168-4137-ae48-2fdce6b3e3ca
