### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 0473f6fe-e832-11ee-114a-dfd1eabad789
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using TSO
	using PythonPlot
	plt = matplotlib.pyplot
	
end

# ╔═╡ 9707519a-1f1f-4808-8c0c-2e8b8b6332c2
eos_path = "DIS_MARCS_t5777g44m00_v0.5/"

# ╔═╡ 2b2d65ff-8810-452a-995e-63552064c033
aos = @axed reload(SqEoS, joinpath(eos_path, "eos.hdf5"))

# ╔═╡ d0938bd6-f49d-4e77-b240-cfdfe5ce2984


# ╔═╡ d4bc1289-b564-4572-8ac6-3cf622983eb2
av_path = "av_models/t5777g44m0005_00070_av.dat"

# ╔═╡ 3d4f8452-4366-4eb2-9067-55aa21083568
model_orig = TSO.flip(
	@optical(Average3D(aos, av_path, logg=4.4), aos), 
	depth=true
)

# ╔═╡ 5a355b88-4442-475e-850f-70a0f4900378


# ╔═╡ 0b0d4d92-f6be-4a19-ba3a-5a3f3907a5ea
z_new = TSO.rosseland_depth(aos, model_orig)

# ╔═╡ cd273657-8fe9-4461-b5f9-38662d855f48
model_new = deepcopy(model_orig)

# ╔═╡ d5dcb199-30fb-45f6-94ef-e1a8bcae727d
TSO.optical_surface!(model_new)

# ╔═╡ 0567d717-6773-4b49-9932-9572a11d2171
model_new.z .= z_new

# ╔═╡ b2286204-0084-4b95-bb7b-f0622d8aef5e
# We now can interpolate everthing to this new z scale
m = TSO.flip(
	TSO.interpolate_to(
		model_new, 
		z=collect(range(maximum(model_new.z), minimum(model_new.z), length=200))
	), 
	depth=false
)

# ╔═╡ 3c043ef1-de59-4ac1-8779-d23d34fd27bc
begin
	mask_orig = sortperm(model_orig.τ)
	
	ip_z_orig = TSO.linear_interpolation(
		TSO.Interpolations.deduplicate_knots!(
			log10.(model_orig.τ[mask_orig]), move_knots=true
		),
		model_orig.z[mask_orig], 
		extrapolation_bc=TSO.Flat()
	)
end

# ╔═╡ f7ef1b6a-0740-4a36-a2e3-880fa8dbde9f
begin
	mask = sortperm(model_new.τ)
	
	ip_z = TSO.linear_interpolation(
		TSO.Interpolations.deduplicate_knots!(
			log10.(model_new.τ[mask]), move_knots=true
		),
		model_new.z[mask], 
		extrapolation_bc=TSO.Flat()
	)
end

# ╔═╡ 39041c2a-b39a-4875-ac90-02d5057432fa


# ╔═╡ 44174f2f-0839-4d71-9fce-462b2a8c0a33
mload =  TSO.flip(
	@optical(Average3D(aos, "DIS_MARCS_E_t5777g44m00_v0.5/inim.dat", logg=4.4), aos), 
	depth=true
)

# ╔═╡ e0ace017-c69c-4487-ad2b-622e939721df


# ╔═╡ 545e8d9e-3324-48cb-9008-ff9c022145e8
let
	plt.close()

	plt.plot(model_orig.z, model_orig.lnT, label="original")
	plt.plot(-m.z, m.lnT, label="new")
	plt.plot(mload.z, mload.lnT)

	#plt.axvline(ip_z.(0.0))
	#plt.axvline(ip_z_orig.(0.0))

	norm = 9.158e7
	s = 3 
	lower = -(-s/2 -0.799) * norm
	upper = -(s/2 -0.799) * norm
	plt.axvline(lower)
	plt.axvline(upper)
	
	plt.legend()
	
	gcf()
end

# ╔═╡ c0568d5e-69d7-4be0-91e3-6bf2996ed27c
let
	plt.close()

	plt.plot(log10.(model_orig.τ), model_orig.z, label="original")
	plt.plot(log10.(m.τ), -m.z, label="new")

	plt.axvline(6.)
	plt.axhline(1.7e8)
	
	plt.legend()
	
	gcf()
end

# ╔═╡ 178282f2-d4fa-4e1b-8c7f-2a8c57613390


# ╔═╡ Cell order:
# ╠═0473f6fe-e832-11ee-114a-dfd1eabad789
# ╠═9707519a-1f1f-4808-8c0c-2e8b8b6332c2
# ╠═2b2d65ff-8810-452a-995e-63552064c033
# ╟─d0938bd6-f49d-4e77-b240-cfdfe5ce2984
# ╠═d4bc1289-b564-4572-8ac6-3cf622983eb2
# ╠═3d4f8452-4366-4eb2-9067-55aa21083568
# ╟─5a355b88-4442-475e-850f-70a0f4900378
# ╠═0b0d4d92-f6be-4a19-ba3a-5a3f3907a5ea
# ╠═cd273657-8fe9-4461-b5f9-38662d855f48
# ╠═d5dcb199-30fb-45f6-94ef-e1a8bcae727d
# ╠═0567d717-6773-4b49-9932-9572a11d2171
# ╠═b2286204-0084-4b95-bb7b-f0622d8aef5e
# ╠═3c043ef1-de59-4ac1-8779-d23d34fd27bc
# ╠═f7ef1b6a-0740-4a36-a2e3-880fa8dbde9f
# ╟─39041c2a-b39a-4875-ac90-02d5057432fa
# ╠═44174f2f-0839-4d71-9fce-462b2a8c0a33
# ╟─e0ace017-c69c-4487-ad2b-622e939721df
# ╠═545e8d9e-3324-48cb-9008-ff9c022145e8
# ╠═c0568d5e-69d7-4be0-91e3-6bf2996ed27c
# ╟─178282f2-d4fa-4e1b-8c7f-2a8c57613390
