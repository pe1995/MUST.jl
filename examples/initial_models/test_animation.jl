### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 0aa65414-d6e3-11ee-3afb-a9840965f7bf
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot

	plt = matplotlib.pyplot
end

# ╔═╡ 24c513c9-d697-493a-bd86-5d4bc6fb7728
np = MUST.numpy

# ╔═╡ d3c7bee9-fcc4-4aa9-8327-ed0ec7581a3f
visual = MUST.ingredients("visual.jl")

# ╔═╡ 91c92c7e-ab27-48d7-938a-3d4bf1d6ea69
MUST.@import_m3dis "../../../Multi3D"

# ╔═╡ dcba5288-3b15-40ba-bf7a-a79dfd9c0bc2


# ╔═╡ c28ad961-1f7a-4b13-8b25-6ab7d3541f0f
run = M3DISRun("/mnt/beegfs/gemini/groups/bergemann/users/klevas/model_analysis/d3t43g15mm10n01/0632173/data")

# ╔═╡ 6ecf9080-6752-40d9-8528-c50a2efc9626
x, y, z = MUST.meshgrid(
	MUST.pyconvert(Array, run.run.xx), 
	MUST.pyconvert(Array, run.run.yy), 
	MUST.pyconvert(Array, run.run.zz)
)

# ╔═╡ a9f9a525-15b0-459f-9b71-49d129045ec5
data = Dict(
	:T => MUST.pyconvert(Array, run.run.temp),
	:ux => MUST.pyconvert(Array, run.run.vx),
	:uy => MUST.pyconvert(Array, run.run.vy),
	:uz => MUST.pyconvert(Array, run.run.vz)
)

# ╔═╡ e3953f2b-2193-4c05-b6cd-0b351b93d68f
b = Box(x, y, z, data, MUST.AtmosphericParameters())

# ╔═╡ 3755b272-a9b5-4b00-bb0c-1c4e3a3be284


# ╔═╡ 9b8472fd-7091-4515-b129-73f7667573c5
b.x[:, 1, 1]

# ╔═╡ 9707c21f-c8b1-4505-be6b-489df2efe541
xx, yy, zz = MUST.pyconvert.(Array, np.meshgrid(
	MUST.pyconvert(Array, run.run.xx), 
	MUST.pyconvert(Array, run.run.yy), 
	MUST.pyconvert(Array, run.run.zz), indexing="ij"
))

# ╔═╡ 904cbb04-f51b-4fe4-9a65-1a1035de4ce5
all(xx .== x)

# ╔═╡ 6ce6f905-3d61-4c02-9145-cb31bb0a6693
begin
	velCube_var = :T
	velCube_vmin_3d = 2000
	velCube_vmax_3d = 15500
	velCube_s_3d = 12
	velCube_arrow_length_ratio = 0.2
	velCube_skipv = 1
	velCube_xoff = 1
	velCube_yoff = 1
	velCube_zoff = 1
	velCube_len_vec = 100
	velCube_cmap = "RdYlBu"
	velCube_show_time = false
	norm=1
	limx=5e3
	limy=4e3
	limz=1.0e3
end;

# ╔═╡ 29f6904f-9140-444f-a304-d02dc3079a7c
begin
	plt.close()
	f, ax = visual.cube_with_velocities(b, 
		velCube_var,
		vmin_3d=velCube_vmin_3d,
		vmax_3d=velCube_vmax_3d,
		s_3d=velCube_s_3d,
		arrow_length_ratio=velCube_arrow_length_ratio,
		skipv=velCube_skipv,
		xoff=velCube_xoff,
		yoff=velCube_yoff,
		zoff=velCube_zoff,
		len_vec=velCube_len_vec,
		cmap=velCube_cmap,
		show_time=velCube_show_time,
		limx=limx,
		limy=limy,
		limz=limz,
		norm=norm
	)

	gcf()
end

# ╔═╡ Cell order:
# ╠═0aa65414-d6e3-11ee-3afb-a9840965f7bf
# ╠═24c513c9-d697-493a-bd86-5d4bc6fb7728
# ╠═d3c7bee9-fcc4-4aa9-8327-ed0ec7581a3f
# ╠═91c92c7e-ab27-48d7-938a-3d4bf1d6ea69
# ╟─dcba5288-3b15-40ba-bf7a-a79dfd9c0bc2
# ╠═c28ad961-1f7a-4b13-8b25-6ab7d3541f0f
# ╠═6ecf9080-6752-40d9-8528-c50a2efc9626
# ╠═a9f9a525-15b0-459f-9b71-49d129045ec5
# ╠═e3953f2b-2193-4c05-b6cd-0b351b93d68f
# ╟─3755b272-a9b5-4b00-bb0c-1c4e3a3be284
# ╠═9b8472fd-7091-4515-b129-73f7667573c5
# ╠═9707c21f-c8b1-4505-be6b-489df2efe541
# ╠═904cbb04-f51b-4fe4-9a65-1a1035de4ce5
# ╠═6ce6f905-3d61-4c02-9145-cb31bb0a6693
# ╠═29f6904f-9140-444f-a304-d02dc3079a7c
