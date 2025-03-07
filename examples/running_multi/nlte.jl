### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ ae78c684-f41b-11ef-0759-8bad589b0ec7
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ 35ccc947-55a1-4a2a-bb18-81739ad5d93a
@import_dispatch "../../../dispatch2"

# ╔═╡ 052f3135-8139-4446-a341-b2a5f5ee789b
@import_m3dis "../../../m3dis/experiments/Multi3D"

# ╔═╡ 3bb5ea4c-9ff6-440e-a1f1-ab2f6f4346b0
plt = matplotlib.pyplot

# ╔═╡ 94f7d445-c801-4ec9-bae5-af97764e2459


# ╔═╡ 1d0ce36f-82f6-4f0c-bb3c-587c30623302
md"""
# 3D Atmospheres with 1D NLTE Opacities
Running a hydrodynamic model atmosphere in full NLTE is currently computationally not possible. One workaround is to use a 1D NLTE approximation, that is precomputed before the simulation. For this we utilize the initial 1D model, which is assumed to ressemble the average structure reasonably well (question is to what extent) to compute NLTE depatures of level populations and the corresponding opacities. Those opacity departures can then be applied to the monochromatic opacity table in a approximative way.
"""

# ╔═╡ 033c8796-7325-4022-b44a-a124b68830a8


# ╔═╡ 7f13d02d-d495-45b1-8273-2e72bb9692b2
md"First we load the monochromatic opacity tables, that cover every model at a given composition."

# ╔═╡ ff951a16-6389-4811-a815-d209284d2666
eos = reload(SqEoS, "../../../opacity_tables/magg_m0_a0_vmic1_v1.2.3/combined_eos_magg_m0_a0_vmic1.hdf5")

# ╔═╡ 4c137681-266b-4782-8f76-e0a604544984
opa = reload(SqOpacity, "../../../opacity_tables/magg_m0_a0_vmic1_v1.2.3/combined_opacities_magg_m0_a0_vmic1.hdf5", mmap=true)

# ╔═╡ c0ff69ac-884e-4970-b9cf-1fb867212903


# ╔═╡ 439ff8ad-ca72-49bb-a78d-457d449b9ff3
md"Then we select a model for computing the NLTE opacities. This has to be done after the model already has been interpolated. It can be saved in the correct format before it is stored later. Here we do the conversion manully now."

# ╔═╡ 5f4d4f4b-2a5a-43d8-a92f-c6a44c6e5567
model = TSO.flip(Average3D(
	@in_dispatch("MainSequenceInterpolated/ST20_E_t57.77g44.40m0.000_v1.0/inim.dat")
), depth=true)

# ╔═╡ d61c1630-0620-4586-933b-354023242925
let
	plt.close()

	f, ax = plt.subplots()

	ax.plot(model.lnρ, model.lnT)
	
	ax.set_ylabel("ln T")
	ax.set_xlabel("ln rho")
	f
end

# ╔═╡ 39c1f40d-8c96-4217-9e4f-9df82d52dd49
begin
	f_new = @in_dispatch("MainSequenceInterpolated/ST20_E_t57.77g44.40m0.000_v1.0/inim_m3d.dat")

	open(f_new, "w") do f
		write(f, "ST20_E_t57.77g44.40m0.000_v1.0\n")
		write(f, "$(length(model.z))\n")
		
		for i in eachindex(model.z)
			line = MUST.@sprintf "%.6E    %.1f    %.4E    %.4E    %.1f\n" model.z[i] exp(model.lnT[i]) 0.0 exp(model.lnρ[i]) 0.0
			write(f, line)
		end
	end
end

# ╔═╡ de05e2d1-2a27-4f82-a60c-b581533e6c61


# ╔═╡ a998cdb8-729b-40f8-9579-57ad5b4767fb
md"We can now use the `spectrum.jl` script to perform a NLTE spectrum synthesis with whatever parameters are available there. This will provide the opacities we need to adjust the monochromatic table."

# ╔═╡ 2c2d1340-0c33-4bdd-b46a-17d25d5fbdb5
result = MUST.M3DISRun("/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/dispatch2/experiments/stellar_atmospheres/MainSequenceInterpolated/ST20_E_t57.77g44.40m0.000_v1.0/spectra_inim_m3d.dat_lam_7-12_NLTE")

# ╔═╡ b61bcf9c-f8b8-4840-945e-3ce354d5c612
result_LTE = MUST.M3DISRun("/mnt/beegfs/gemini/groups/bergemann/users/eitner/StAt/dispatch2/experiments/stellar_atmospheres/MainSequenceInterpolated/ST20_E_t57.77g44.40m0.000_v1.0/spectra_inim_m3d.dat_lam_7-12_LTE")

# ╔═╡ cfc931e5-1ff1-40b6-9556-37270e80616c
begin
	λ = MUST.pyconvert(Array, result.lam)
	χ_NLTE = MUST.pyconvert(Array, 
		result.run.read_patch_save("chi_01", concat=true, fdim=0, lazy=false)[0]
	)
	S_NLTE = MUST.pyconvert(Array, 
		result.run.read_patch_save("Snu_01", concat=true, fdim=0, lazy=false)[0]
	)
end;

# ╔═╡ e10ea5c3-e0ce-480c-9fa2-6aefa072b55d
begin
	χ_LTE2 = MUST.pyconvert(Array, 
		result_LTE.run.read_patch_save("chi_01", concat=true, fdim=0, lazy=false)[0]
	)
	S_LTE2 = MUST.pyconvert(Array, 
		result_LTE.run.read_patch_save("Snu_01", concat=true, fdim=0, lazy=false)[0]
	)
end;

# ╔═╡ 88057dc9-d90d-49f4-9aca-0923a0f544ea
begin
	lnrho = log.(MUST.pyconvert(Array, result.run.rho))
	lnT = log.(MUST.pyconvert(Array, result.run.temp))
end;

# ╔═╡ 4b07d299-d438-485b-b5c7-65d5547f5a3e


# ╔═╡ 4e709357-3156-4ea7-bd0b-9eccaa715587
md"In principle we could now substiute opacities and source function in the monochromatic opacity table at the correct positions, however we may also want to apply corrections close the the actual profile. So for this, we compute the opacity and source function departure."

# ╔═╡ b4f72e3b-953e-436b-aa33-fba7ce0d7265
χ_LTE = log.(lookup(eos, opa, :κ, lnrho, lnT))

# ╔═╡ bf310bc3-3246-465a-933f-4551ba4ec8ed


# ╔═╡ 8dba1049-0e80-4c6d-9612-d566199f2903
begin
	χ_fraction = similar(χ_LTE) 
	S_fraction = similar(χ_LTE) 
	
	ν = TSO.CLight ./ (λ .* 1e-8)
	frac = similar(ν)
	T = exp.(lnT)
	
	for i in eachindex(lnrho)
		# interpolate in wavelength
		frac .= χ_NLTE[:, 1, 1, i] ./ χ_LTE2[:, 1, 1, i]
		f_χ_NLTE = MUST.linear_interpolation(λ, frac, extrapolation_bc=NaN)
		#χ_fraction[i, :] .= exp.(f_χ_NLTE.(opa.λ)) ./ (exp.(χ_LTE[i, :]) .* exp(lnrho[i]))
		χ_fraction[i, :] .= f_χ_NLTE.(opa.λ)
		χ_fraction[i, isnan.(@view(χ_fraction[i, :]))] .= 1.0

		#frac .= 1 ./ (exp.(MUST.HPlanck .* ν ./ (MUST.KBoltzmann .* T[i])) .-1)
		#frac .= exp.(S_NLTE[:, 1, 1, i]) ./ frac 
		frac .= S_NLTE[:, 1, 1, i] ./ S_LTE2[:, 1, 1, i]
		f_S_NLTE = MUST.linear_interpolation(λ, frac, extrapolation_bc=NaN)
		S_fraction[i, :] .= f_S_NLTE.(opa.λ)
		S_fraction[i, isnan.(@view(S_fraction[i, :]))] .= 1.0
	end
end

# ╔═╡ 199620c9-c217-416b-b3ca-f1ece6e7a73f


# ╔═╡ f8eb5a48-cb3a-4ede-9329-832416aba12e
irho = 221

# ╔═╡ 499f3c76-709b-4ea7-9424-095988f46de4
begin
	@info "ρ = $(exp(lnrho[irho]))"
	@info "T = $(exp(lnT[irho]))"
end

# ╔═╡ 7689f092-0fac-4005-862c-797a53157bfa
let
	plt.close()

	f, ax = plt.subplots()

	ax.plot(opa.λ, χ_fraction[irho, :])
	ax.set_xscale("log")
	#ax.set_yscale("log")

	imax = argmax(χ_fraction[irho, :])
	@show opa.λ[imax] imax

	ax.set_ylabel(L"\rm χ_{NLTE}\ /\ χ_{LTE}")
	ax.set_xlabel(L"\rm \lambda")

	f
end

# ╔═╡ b88a8e5b-3a0a-411b-adb8-5de7082d92cf
let
	plt.close()

	f, ax = plt.subplots()

	ax.plot(opa.λ, S_fraction[irho, :])
	ax.set_xscale("log")
	#ax.set_yscale("log")

	ax.set_ylabel(L"\rm S\nu_{NLTE}\ /\ S\nu_{LTE} ")
	ax.set_xlabel(L"\rm \lambda")

	f
end

# ╔═╡ 7da2d30a-fa18-4566-9981-11926a2a95db


# ╔═╡ a7e94add-5962-4303-a9db-facacc2d3bf3
md"We can not apply those differences to the opacity table and investigate the differences. One can either apply them for all temperatures at a given density, or e.g. as a gaussian around the profile."

# ╔═╡ 9485ab43-229d-4145-9fec-84442755768e
begin
	correction_χ = similar(opa.κ)
	correction_S = similar(opa.src)
	dist = similar(eos.lnPg)
	idist = zeros(Int, size(eos.lnPg)...)
	idistT = zeros(Int, size(eos.lnPg)...)
	idistR = zeros(Int, size(eos.lnPg)...)
	temp = similar(opa.λ)

	for j in eachindex(eos.lnRho)
		for i in eachindex(eos.lnT)
			idist[i, j] = argmin(
				(eos.lnRho[j]/4 .- lnrho/4).^2 .+ (eos.lnT[i] .- lnT).^2 
			)
			idistT[i, j] = argmin(
				(eos.lnT[i] .- lnT).^2 
			)
			idistR[i, j] = argmin(
				(eos.lnRho[j] .- lnrho).^2 
			)
			dist[i, j] = (eos.lnRho[j]/4 .- lnrho[idist[i, j]]/4).^2 .+ (eos.lnT[i] .- lnT[idist[i, j] ]).^2 
		end
	end

	σ_rho = (abs(maximum(eos.lnRho)/4 - minimum(eos.lnRho)/4) / 16) ^2
	σ_T = (abs(maximum(eos.lnT) - minimum(eos.lnT)) / 16) ^2
	σ = σ_rho + σ_T

	nr = length(lnrho)
	function average_and_add!(lnRho, lnT, lam, idist, dist, χ_fraction, correction_χ, S_fraction, correction_S)
		for j in eachindex(lnRho)
			for i in eachindex(lnT)
				#temp .= 
				correction_χ[i, j, :] .= (χ_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
				#temp .= 
				correction_S[i, j, :] .= (S_fraction[idistR[i, j], :] .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
				#minoff = max(1, idist[i, j]-1)
				#maxoff = min(nr, idist[i, j]+1)
				#=temp  .= @view χ_fraction[idist[i, j], :] 
				temp .+= @view χ_fraction[minoff, :] 
				temp .+= @view χ_fraction[maxoff, :] 
				temp ./= 3		
				correction_χ[i, j, :] .= (temp .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1
	
				temp  .= @view S_fraction[idist[i, j], :] 
				temp .+= @view S_fraction[minoff, :] 
				temp .+= @view S_fraction[maxoff, :] 
				temp ./= 3	
				correction_S[i, j, :] .= (temp .- 1) .* exp(-dist[i, j]/(2*σ)) .+ 1=#
				#=w = exp(-dist[i, j]/(2*σ)) 
				@inbounds for k in eachindex(lam)
					correction_χ[i, j, k] = (χ_fraction[idist[i, j], k] +
											 χ_fraction[minoff,      k] +
											 χ_fraction[maxoff,      k]) / 3.0
					correction_S[i, j, k] = (S_fraction[idist[i, j], k] +
											 S_fraction[minoff,      k] +
											 S_fraction[maxoff,      k]) / 3.0
					
					correction_χ[i, j, k] = (correction_χ[i, j, k] - 1) * w +1
					correction_S[i, j, k] = (correction_S[i, j, k] - 1) * w +1
				end=#
			end
		end
	end
	average_and_add!(eos.lnRho, eos.lnT, opa.λ, idist, dist, χ_fraction, correction_χ, S_fraction, correction_S)
end

# ╔═╡ fea23d8f-c25e-4e06-9eb2-058011719e7c


# ╔═╡ cb998401-7836-4a6b-821e-dd4a037ce307
md"We can do the reverse test now:"

# ╔═╡ 371d4a90-16c0-4960-9574-081b71d8b9c5
opa_new = SqOpacity(
	opa.κ .* correction_χ,
	opa.κ_ross,
	opa.src .* correction_S,
	opa.λ,
	false
)

# ╔═╡ 998abbaa-f2f3-41c7-a586-1717595b9001
orig_test = lookup(eos, opa, :κ, lnrho[irho], lnT[irho])

# ╔═╡ a7ce8ea2-be22-4670-9aa8-af081661cf14
corr_test = lookup(eos, opa_new, :κ, lnrho[irho], lnT[irho])

# ╔═╡ 091827d1-cfb3-4bbb-9e46-55214a73e41a
let
	plt.close()

	f, ax = plt.subplots()

	ax.plot(opa_new.λ, corr_test ./ orig_test)
	ax.set_xscale("log")
	#ax.set_yscale("log")

	ax.set_ylabel(L"\rm χ_{NLTE}\ /\ χ_{LTE}\  after\ interpolation")
	ax.set_xlabel(L"\rm \lambda")

	f
end

# ╔═╡ ca63d02a-4edc-4b5c-bde1-d9d6439d38e1
let
	plt.close()

	f, ax = plt.subplots()

	TT, rr = meshgrid(@axed(eos))
	
	#im = ax.scatter(TT, rr, c=opa_new.κ[:,:,88749] ./ opa.κ[:,:,88749], rasterized=true, s=1)
	#im = ax.scatter(TT, rr, c=correction_χ[:,:,88749], rasterized=true, s=1)
	im = ax.scatter(TT, rr, c=exp.(-dist./(2*σ)), rasterized=true, s=1)

	ax.plot(lnT, lnrho, color="k")
	#im=ax.scatter(lnT, lnrho, c=χ_fraction[:, 88749])

	plt.colorbar(im)
	
	f
end

# ╔═╡ 1d64385c-c7a4-4fbd-9125-deeeed521ada
let
	plt.close()

	f, ax = plt.subplots()

	TT, rr = meshgrid(@axed(eos))
	
	#im = ax.scatter(TT, rr, c=opa_new.κ[:,:,88749] ./ opa.κ[:,:,88749], rasterized=true, s=1)
	im = ax.scatter(TT, rr, c=exp.(-dist./(2*σ)), rasterized=true, s=1)
	ax.plot(lnT, lnrho, color="k")
	im=ax.scatter(lnT, lnrho, c=χ_fraction[:, 88749])

	plt.colorbar(im)
	
	f
end

# ╔═╡ d5a7b6b1-6cb6-41d8-85ea-3049d4fb87e7
let
	plt.close()

	f, ax = plt.subplots()

	TT, rr = meshgrid(@axed(eos))
	
	#im = ax.scatter(TT, rr, c=opa_new.κ[:,:,88749] ./ opa.κ[:,:,88749], rasterized=true, s=1)
	im = ax.scatter(TT, rr, c=correction_χ[:,:,88749], rasterized=true, s=1)

	ax.plot(lnT, lnrho, color="k")
	#im=ax.scatter(lnT, lnrho, c=χ_fraction[:, 88749])

	plt.colorbar(im)
	
	f
end

# ╔═╡ Cell order:
# ╠═ae78c684-f41b-11ef-0759-8bad589b0ec7
# ╠═35ccc947-55a1-4a2a-bb18-81739ad5d93a
# ╠═052f3135-8139-4446-a341-b2a5f5ee789b
# ╠═3bb5ea4c-9ff6-440e-a1f1-ab2f6f4346b0
# ╟─94f7d445-c801-4ec9-bae5-af97764e2459
# ╟─1d0ce36f-82f6-4f0c-bb3c-587c30623302
# ╟─033c8796-7325-4022-b44a-a124b68830a8
# ╟─7f13d02d-d495-45b1-8273-2e72bb9692b2
# ╠═ff951a16-6389-4811-a815-d209284d2666
# ╠═4c137681-266b-4782-8f76-e0a604544984
# ╟─c0ff69ac-884e-4970-b9cf-1fb867212903
# ╟─439ff8ad-ca72-49bb-a78d-457d449b9ff3
# ╠═5f4d4f4b-2a5a-43d8-a92f-c6a44c6e5567
# ╟─d61c1630-0620-4586-933b-354023242925
# ╠═39c1f40d-8c96-4217-9e4f-9df82d52dd49
# ╟─de05e2d1-2a27-4f82-a60c-b581533e6c61
# ╟─a998cdb8-729b-40f8-9579-57ad5b4767fb
# ╠═2c2d1340-0c33-4bdd-b46a-17d25d5fbdb5
# ╠═b61bcf9c-f8b8-4840-945e-3ce354d5c612
# ╠═cfc931e5-1ff1-40b6-9556-37270e80616c
# ╠═e10ea5c3-e0ce-480c-9fa2-6aefa072b55d
# ╠═88057dc9-d90d-49f4-9aca-0923a0f544ea
# ╟─4b07d299-d438-485b-b5c7-65d5547f5a3e
# ╟─4e709357-3156-4ea7-bd0b-9eccaa715587
# ╠═b4f72e3b-953e-436b-aa33-fba7ce0d7265
# ╟─bf310bc3-3246-465a-933f-4551ba4ec8ed
# ╠═8dba1049-0e80-4c6d-9612-d566199f2903
# ╟─199620c9-c217-416b-b3ca-f1ece6e7a73f
# ╠═f8eb5a48-cb3a-4ede-9329-832416aba12e
# ╟─499f3c76-709b-4ea7-9424-095988f46de4
# ╠═7689f092-0fac-4005-862c-797a53157bfa
# ╟─b88a8e5b-3a0a-411b-adb8-5de7082d92cf
# ╟─7da2d30a-fa18-4566-9981-11926a2a95db
# ╟─a7e94add-5962-4303-a9db-facacc2d3bf3
# ╠═9485ab43-229d-4145-9fec-84442755768e
# ╟─fea23d8f-c25e-4e06-9eb2-058011719e7c
# ╟─cb998401-7836-4a6b-821e-dd4a037ce307
# ╠═371d4a90-16c0-4960-9574-081b71d8b9c5
# ╠═998abbaa-f2f3-41c7-a586-1717595b9001
# ╠═a7ce8ea2-be22-4670-9aa8-af081661cf14
# ╟─091827d1-cfb3-4bbb-9e46-55214a73e41a
# ╟─ca63d02a-4edc-4b5c-bde1-d9d6439d38e1
# ╟─1d64385c-c7a4-4fbd-9125-deeeed521ada
# ╠═d5a7b6b1-6cb6-41d8-85ea-3049d4fb87e7
