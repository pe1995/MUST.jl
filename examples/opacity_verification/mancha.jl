### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d64c7fa0-8212-11ef-2afe-7f449eb43a08
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using PythonPlot
end

# ╔═╡ 42039451-14e3-4f79-94fa-47a2b6b0684b
plt = matplotlib.pyplot

# ╔═╡ 1383ae68-993e-442a-8893-137470156af7
MUST.@import_dispatch "../../../dispatch2"

# ╔═╡ 8c85a22a-7aa7-420e-a5dc-dbe59ab5d7c7


# ╔═╡ 8617f4a8-fecc-4186-a6c0-2d22a92a0183
b, bτ = pick_snapshot(@in_dispatch("data/ST1_E_t57.77g44.40m0.000_v1.0_fast2"), 3841)

# ╔═╡ 9bec4de7-1a19-49f7-a872-cfabddd620df


# ╔═╡ b33169b4-ca20-4d10-b737-845df685fe5b
begin
	eosdir = "input_data/grd/MainSequenceInterpolated/ST6_E_t57.77g44.40m0.000_v1.0"
	eospath = @in_dispatch(joinpath(eosdir, "eos.hdf5"))
	opapath = @in_dispatch(joinpath(eosdir, "binned_opacities.hdf5"))
	
	eos = reload(SqEoS, eospath)
	opa = reload(SqOpacity, opapath)
end

# ╔═╡ 28372005-4de3-4960-aea3-f03cff582f78
begin
	eospathT = @in_dispatch(joinpath(eosdir, "eos_T.hdf5"))
	opapathT = @in_dispatch(joinpath(eosdir, "binned_opacities_T.hdf5"))
	
	eosT = reload(SqEoS, eospathT)
	opaT = reload(SqOpacity, opapathT)
end

# ╔═╡ 58bba96b-6516-46a4-9993-cf52e722b3f9
begin
	eosManchadir = "input_data/MANCHA_M3D_E_magg_m0_a0_vmic1_v5.1"
	eospathMancha = @in_dispatch(joinpath(eosManchadir, "eos.hdf5"))
	opapathMancha = @in_dispatch(joinpath(eosManchadir, "binned_opacities.hdf5"))
	
	eosMancha = reload(SqEoS, eospathMancha)
	opaMancha = reload(SqOpacity, opapathMancha)
end

# ╔═╡ 228fae7a-3089-4331-8066-566862ac7864
begin
	eospathTMancha = @in_dispatch(joinpath(eosManchadir, "eos_T.hdf5"))
	opapathTMancha = @in_dispatch(joinpath(eosManchadir, "binned_opacities_T.hdf5"))
	
	eosTMancha = reload(SqEoS, eospathTMancha)
	opaTMancha = reload(SqOpacity, opapathTMancha)
end

# ╔═╡ b1d4e175-bb3f-44e9-898c-ab41a068bcf8


# ╔═╡ d55e97a8-4779-4834-a8fd-58ecd4f409ff
begin	
	eosStaggerdir = "input_data/DISSTAG_E_v0.1"
	eospathStagger = @in_dispatch(joinpath(eosStaggerdir, "eos.hdf5"))
	opapathStagger = @in_dispatch(joinpath(eosStaggerdir, "binned_opacities.hdf5"))
	
	eosStagger = reload(SqEoS, eospathStagger)
	opaStagger = reload(SqOpacity, opapathStagger)

	_, ldStag = profile(MUST.mean, b, :z, :logd)
	_, lTStag = profile(MUST.mean, b, :z, :logT)
	leeStag = lookup(eosT, :lnEi, ldStag, lTStag)
	eeStag = lookup(eosStagger, :lnEi, ldStag, lTStag)
	corr = MUST.mean(eeStag .- leeStag)
	eosStagger.lnEi .= eosStagger.lnEi .- corr
end

# ╔═╡ f9a52b86-9a40-473c-b5b8-7dc9e0de82e0


# ╔═╡ fd9999e7-29ff-4538-b1dd-5359a5efa6e7
md"We can replace the EoS from Stagger by looking up the opacities on the new EoS grid."

# ╔═╡ db200f5d-7ce5-4031-a552-dbd86458591c
begin
	# lookup the opacity on the new EoS from the old EoS
	# Use the temperature for this
	TConv = deepcopy(eosStagger.lnT)
	rhoConv = similar(TConv)
	for j in axes(TConv, 2)
		rhoConv[:, j] .= eosStagger.lnRho[j]
	end
	
	lnEiConv = lookup(eos, :lnEi, rhoConv, TConv)
	κConv = lookup(eos, opa, :κ, rhoConv, lnEiConv)
	sConv = lookup(eos, opa, :src, rhoConv, lnEiConv)
	krConv = lookup(eos, opa, :κ_ross, rhoConv, lnEiConv)

	opa_StaggerMix = deepcopy(opaStagger)
	opa_StaggerMix.κ .= κConv
	opa_StaggerMix.κ_ross .= krConv
	opa_StaggerMix.src .= sConv

	eosStaggerMix = deepcopy(eosStagger)
	eosStaggerMix.lnRoss .= log.(krConv)

	newdir = MUST.@in_dispatch "input_data/DISTSO_E_v0.1"
	for_dispatch(eosStaggerMix, opa_StaggerMix, newdir)

	TSO.save(eosStaggerMix, joinpath(newdir, "eos.hdf5"))
	TSO.save(opa_StaggerMix, joinpath(newdir, "binned_opacities.hdf5"))
end;

# ╔═╡ 3e8bb05f-600c-4efa-9a5f-22f8d5c876ea


# ╔═╡ 857a0709-706f-4763-a6c0-1627be747ac2


# ╔═╡ 12e2ff3b-3b0a-4976-8e27-60d7ee040c1b
md"Test: Cut the opacity table beyond 1e-11?"

# ╔═╡ 21d36112-42c2-4b1a-943f-def2d65a0eaf
begin
	eos_mask = eosMancha.lnRho .<= 1e-11
end

# ╔═╡ 5b925a3c-ec86-49a3-85d7-736436a58f21


# ╔═╡ a49805e9-1cc9-4da4-a7ac-d8cd17a31b6d
let
	eos_ev = deepcopy(eosT)
	for j in axes(eos_ev.lnEi, 2)
		eos_ev.lnEi[:, j] .+= eos_ev.lnRho[j]
	end

	ev_dir = "test_4bins_ev"
	ev_dir_E = "test_4bins_ev_E"
	!isdir(ev_dir) && mkdir(ev_dir)
	model = @in_dispatch(joinpath(eosdir, "inim.dat"))
	TSO.save(eos_ev, joinpath(ev_dir, "eos.hdf5"))
	TSO.save(opaT, joinpath(ev_dir, "binned_opacities.hdf5"))

	TSO.convert_fromT_toE(ev_dir, ev_dir_E, model, lnEi_stretch=0.5)
end

# ╔═╡ 74ad6b6d-13a2-48a6-b739-60e05172514b
begin
	eos_ev = reload(SqEoS, joinpath("test_4bins_ev_E", "eos.hdf5"))
	opa_ev = reload(SqOpacity, joinpath("test_4bins_ev_E", "binned_opacities.hdf5"))
end

# ╔═╡ d3061189-067a-4127-b027-b9e85a3c94ea


# ╔═╡ d04b9cc1-99ed-4d6b-b9d6-5ff270795cd2
bin = 4

# ╔═╡ 2cd909d3-76c1-4015-aa79-c99352a4a0be
ee, rr = TSO.meshgrid(@axed(eos))

# ╔═╡ 2cd3ef34-1478-42cb-ae46-5368a9744251
eeMancha, rrMancha = TSO.meshgrid(@axed(eosMancha))

# ╔═╡ 130a2880-5166-44a6-956a-e58fe112bcfb
eeStagger, rrStagger = TSO.meshgrid(@axed(eosStagger))

# ╔═╡ 14787aa7-de6a-4b09-be21-fe6d516d5a3a
let
	f, ax = plt.subplots(1, 1)

	TT, dd = TSO.meshgrid(@axed(eosT))

	mask = eosT.lnEi .> 31.0
	lnei = deepcopy(eosT.lnEi)
	lnei[mask] .= 31.0
	
	im = ax.contourf(
		log10.(exp.(TT)), log10.(exp.(dd)), log10.(exp.(lnei)), 
		rasterized=true, levels=40
	)
	c = f.colorbar(im, ax=ax)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10T)
	ax.plot(ee, d, color="k", lw=2, ls="-", label="mean")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :log10T)
	ax.plot(ee, d, color="magenta", lw=2, ls="-", label="min(ρ)-max(T)")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10T)
	ax.plot(ee, d, color="cyan", lw=2, ls="-", label="min(ρ)-mean(T)")
	
	c.set_label("log Ei")

	ax.set_xlabel("log T")
	ax.set_ylabel("log ρ")

	ax.legend()
	
	f
end

# ╔═╡ 1f2d74ac-64e3-4eaa-bb25-8fd884e3c6ab
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	im1 = ax.contourf(
		log10.(exp.(ee)), log10.(exp.(rr)), log10.(exp.(eos.lnT)), 
		rasterized=true, levels=40
	)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10ee)
	ax.plot(ee, d, color="k", lw=2, ls="-", label="mean")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :log10ee)
	ax.plot(ee, d, color="magenta", lw=2, ls="-", label="min(ρ)-max(E)")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10ee)
	ax.plot(ee, d, color="cyan", lw=2, ls="-", label="min(ρ)-mean(E)")

	c1 = f.colorbar(im1, ax=ax)

	c1.set_label("log T")
	ax.set_xlabel("log Ei")
	ax.set_ylabel("log ρ")
	ax.legend()
	
	f
end

# ╔═╡ 0077c25e-5ec6-4e49-9378-4cd4f9867170
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	im1 = ax.contourf(
		log10.(exp.(eeStagger)), log10.(exp.(rrStagger)), 
		log10.(exp.(eosStagger.lnT)), 
		rasterized=true, levels=40
	)

	_, ld = profile(MUST.mean, b, :z, :logd)
	_, lT = profile(MUST.mean, b, :z, :logT)
	ee = lookup(eosStagger, :lnEi, ld, lT)
	ax.plot(
		log10.(exp.(ee)), log10.(exp.(ld)), 
		color="k", lw=2, ls="-", label="mean"
	)

	_, ld = profile(minimum, b, :z, :logd)
	_, lT = profile(maximum, b, :z, :logT)
	ee = lookup(eosStagger, :lnEi, ld, lT)
	ax.plot(
		log10.(exp.(ee)), log10.(exp.(ld)), 
		color="magenta", lw=2, ls="-", label="min(ρ)-max(E)"
	)

	_, ld = profile(minimum, b, :z, :logd)
	_, lT = profile(MUST.mean, b, :z, :logT)
	ee = lookup(eosStagger, :lnEi, ld, lT)
	ax.plot(
		log10.(exp.(ee)), log10.(exp.(ld)), 
		color="cyan", lw=2, ls="-", label="min(ρ)-mean(E)"
	)

	c1 = f.colorbar(im1, ax=ax)

	c1.set_label("log T")
	ax.set_xlabel("log Ei")
	ax.set_ylabel("log ρ")
	ax.set_title("Legacy Stagger")
	ax.legend()
	
	f
end

# ╔═╡ 813da454-74a0-492e-8092-dddeff054ad1


# ╔═╡ 654f2da6-e5d0-4f37-8a3b-563f5e8f0753
begin
	eos_test = deepcopy(eosT)
	for j in axes(eos_test.lnEi, 2)
		eos_test.lnEi[:, j] .+= eos_test.lnRho[j]
	end
end

# ╔═╡ d040380d-e975-4667-90a1-5de76d02e41a
eev, rrv = TSO.meshgrid(@axed(eos_ev))

# ╔═╡ f8ff0bb8-04d4-4d02-8380-c2eec52c170f
let
	f, ax = plt.subplots(1, 1)

	TT, dd = TSO.meshgrid(@axed(eos_test))

	mask = eos_test.lnEi .> 28.0
	lnei = deepcopy(eos_test.lnEi)
	lnei[mask] .= 28.0
	
	im = ax.contourf(
		log10.(exp.(TT)), log10.(exp.(dd)), lnei, 
		rasterized=true, levels=40
	)
	c = f.colorbar(im, ax=ax)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10T)
	ax.plot(ee, d, color="k", lw=2, ls="-", label="mean")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :log10T)
	ax.plot(ee, d, color="magenta", lw=2, ls="-", label="min(ρ)-max(T)")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10T)
	ax.plot(ee, d, color="cyan", lw=2, ls="-", label="min(ρ)-mean(T)")
	
	c.set_label("log Ei")

	ax.set_xlabel("log T")
	ax.set_ylabel("log ρ")

	ax.legend()
	
	f
end

# ╔═╡ 54a13656-f36d-47eb-bec8-40781968c473
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	im1 = ax.contourf(
		log10.(exp.(eev)), log10.(exp.(rrv)), log10.(exp.(eos_ev.lnT)), 
		rasterized=true, levels=40
	)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10e)
	ax.plot(ee, d, color="k", lw=2, ls="-", label="mean")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :log10e)
	ax.plot(ee, d, color="magenta", lw=2, ls="-", label="min(ρ)-max(E)")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10e)
	ax.plot(ee, d, color="cyan", lw=2, ls="-", label="min(ρ)-mean(E)")

	c1 = f.colorbar(im1, ax=ax)

	c1.set_label("log T")
	ax.set_xlabel("log Ei (V)")
	ax.set_ylabel("log ρ")
	ax.legend()
	
	f
end

# ╔═╡ 38bbc026-1597-4e96-8d48-5432a4a78cf5
begin
	z, d = profile(MUST.mean, b, :z, :logd)
	_, e = profile(MUST.mean, b, :z, :loge)
	_, Treal = profile(MUST.mean, b, :z, :logT)
	
	T = lookup(eos_ev, :lnT, d, e)
	plt.close()
	
	plt.plot(z, T)
	plt.plot(z, Treal)

	gcf()
end

# ╔═╡ 77e9ae8c-eee7-4495-9f36-1d7e9c8bab40


# ╔═╡ 8decebd5-2716-4568-bd2e-447f63128c9f
Cv(eos) = begin
	lne_step = diff(eos.lnEi) |> first

	ee = similar(eos.lnPg)
	d = similar(eos.lnPg)

	for j in axes(ee, 2)
		for i in axes(ee, 1)
			ee[i, j] = eos.lnEi[i]
			d[i, j] = eos.lnRho[j]
		end
	end
	
	dee1 = ee .- lne_step
	dee2 = ee .+ lne_step

	T1 = exp.(lookup(eos, :lnT, d, dee1))
	T2 = exp.(lookup(eos, :lnT, d, dee2))

	(exp.(dee2) .- exp.(dee1)) .* exp.(d) ./ (T2 .- T1)	

	#1.5 .* exp.(d) .*MUST.KBoltzmann
end

# ╔═╡ 445e3b6e-0d48-43ac-9518-1c28a0725e1d
Cv2(eos) = begin
	e = similar(eos.lnPg)
	e1 = similar(eos.lnPg)
	e2 = similar(eos.lnPg)
	d = similar(eos.lnPg)
	de = similar(eos.lnEi)
	
	for j in axes(ee, 2)
		e[:, j] .= exp.(eos.lnEi) .* exp.(eos.lnRho[j])
		d[:, j] .= exp.(eos.lnRho[j])

		de[2:end] .= e[2:end, j] .- e[1:end-1, j]
		de[1] = de[2]
		e1[:, j] .= e[:, j] .- de
		
		de[1:end-1] .= e[2:end, j] .- e[1:end-1, j] 
		de[end] = de[end-1]
		e2[:, j] .= e[:, j] .+ de
	end

	T1 = exp.(lookup(eos, :lnT, log.(d), log.(e1 ./ d)))
	T2 = exp.(lookup(eos, :lnT, log.(d), log.(e2 ./ d)))

	(e2 .- e1) ./ (T2 .- T1)	

	#1.5 .* d .*MUST.KBoltzmann
end

# ╔═╡ 1b21ac59-4977-43ac-a8ff-64424b0f0bc9
gamma(eos) = begin
	d = similar(eos.lnPg)
	for j in axes(d, 2)
		for i in axes(d, 1)
			d[i, j] = exp.(eos.lnRho[j])
		end
	end
	
	d .* MUST.KBoltzmann ./ Cv2(eos) .+ 1
end

# ╔═╡ dcf83225-72c5-4938-ada3-0784301e97aa
csound(eos) = begin
	e = similar(eos.lnPg)
	e1 = similar(eos.lnPg)
	e2 = similar(eos.lnPg)
	d = similar(eos.lnPg)
	de = similar(eos.lnRho)
	
	for i in axes(e, 1)
		e[i, :] .= exp(eos.lnEi[i])
		d[i, :] .= exp.(eos.lnRho)

		de[2:end] .= d[i, 2:end] .- d[i, 1:end-1]
		de[1] = de[2]
		e1[i, :] .= d[i, :] .- de
		
		de[1:end-1] .= d[i, 2:end] .- d[i, 1:end-1] 
		de[end] = de[end-1]
		e2[i, :] .= d[i, :] .+ de
	end

	T1 = exp.(lookup(eos, :lnPg, log.(e1), log.(e)))
	T2 = exp.(lookup(eos, :lnPg, log.(e2), log.(e)))

	
	abs.((T2 .- T1) ./ (e2 .- e1)) .^0.5 ./ 1e5

	#1.5 .* d .*MUST.KBoltzmann
end

# ╔═╡ 7024d5ab-4684-4fbe-8076-31021a06c7cf
gamma2(eos) = begin
	d = similar(eos.lnPg)	
	for i in axes(d, 1)
		d[i, :] .= exp.(eos.lnRho)
	end
	
	g = (csound(eos) .*1e5) .^2 ./ exp.(eos.lnPg) .* d

	g[g .> 1.7] .= 1.7
	g[g .< 1] .= 1

	g
end

# ╔═╡ 9aa170e3-21d0-4fc5-b76c-bdd5be28bbba


# ╔═╡ 227a1613-7281-4978-90f2-9a7c0df6528a
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	im1 = ax.contourf(
		log10.(exp.(ee)), log10.(exp.(rr)), log10.(exp.(eos.lnT)), 
		rasterized=true, levels=40
	)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10ee)
	ax.plot(ee, d, color="k", lw=2, ls="-", label="mean")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :log10ee)
	ax.plot(ee, d, color="magenta", lw=2, ls="-", label="min(ρ)-max(E)")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :log10ee)
	ax.plot(ee, d, color="cyan", lw=2, ls="-", label="min(ρ)-mean(E)")

	c1 = f.colorbar(im1, ax=ax)

	c1.set_label("log T")
	ax.set_xlabel("log Ei")
	ax.set_ylabel("log ρ")
	ax.legend()
	
	f
end

# ╔═╡ e8490b9b-bdcb-4e19-b51a-c174767cfcfb


# ╔═╡ 1ded5ad1-a76d-4aeb-9d0b-b2a8b9e0c8aa
let
	f, ax = plt.subplots(3, 1, sharex=true, sharey=true)

	im1 = ax[0].pcolormesh(
		ee, log10.(exp.(rr)), log10.(opa.κ[:,:,bin]), 
		rasterized=true
	)
	im2 = ax[1].pcolormesh(
		eeMancha, log10.(exp.(rrMancha)), log10.(opaMancha.κ[:,:,bin]), 
		rasterized=true
	)
	#=im3 = ax[2].pcolormesh(
		eeStagger, log10.(exp.(rrStagger)), log10.(opaStagger.κ[:,:,bin]), 
		rasterized=true
	)=#
	im3 = ax[2].pcolormesh(
		eeStagger, log10.(exp.(rrStagger)), log10.(opa_StaggerMix.κ[:,:,bin]), 
		rasterized=true
	)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :logee)
	ax[0].plot(ee, d, color="k", lw=2, ls="-")
	ax[1].plot(ee, d, color="k", lw=2, ls="-")
	ax[2].plot(ee, d, color="k", lw=2, ls="-")
	

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :logee)
	ax[0].plot(ee, d, color="magenta", lw=2, ls="-")
	ax[1].plot(ee, d, color="magenta", lw=2, ls="-")
	ax[2].plot(ee, d, color="magenta", lw=2, ls="-")

	c1 = f.colorbar(im1, ax=ax[0])
	c2 = f.colorbar(im2, ax=ax[1])
	c3 = f.colorbar(im3, ax=ax[2])
	

	c1.set_label("opacity bin $(bin)")
	c2.set_label("opacity bin $(bin)")
	c3.set_label("opacity bin $(bin)")
	

	ax[2].set_xlabel("internal energy")
	ax[1].set_ylabel("density")
	

	ax[0].set_title("TSO opacities")
	ax[1].set_title("MANCHA opacities")
	ax[2].set_title("Stagger (4 bin) opacities")
	
	
	
	
	f
end

# ╔═╡ 437706ae-8307-4e4b-9bc2-2cdc93eaa5d9
let
	f, ax = plt.subplots(2, 1, sharex=true, sharey=true)

	im1 = ax[0].pcolormesh(
		ee, log10.(exp.(rr)), log10.(opa.src[:,:,bin]), 
		rasterized=true
	)
	im2 = ax[1].pcolormesh(
		eeMancha, log10.(exp.(rrMancha)), log10.(opaMancha.src[:,:,bin]), 
		rasterized=true
	)

	_, d = profile(MUST.mean, b, :z, :log10d)
	_, ee = profile(MUST.mean, b, :z, :logee)
	ax[0].plot(ee, d, color="k", lw=2, ls="-")
	ax[1].plot(ee, d, color="k", lw=2, ls="-")

	_, d = profile(minimum, b, :z, :log10d)
	_, ee = profile(maximum, b, :z, :logee)
	ax[0].plot(ee, d, color="magenta", lw=2, ls="-")
	ax[1].plot(ee, d, color="magenta", lw=2, ls="-")

	c1 = f.colorbar(im1, ax=ax[0])
	c2 = f.colorbar(im2, ax=ax[1])

	c1.set_label("S bin $(bin)")
	c2.set_label("S bin $(bin)")

	ax[1].set_xlabel("internal energy")
	ax[0].set_ylabel("density")
	ax[1].set_ylabel("density")
	

	ax[0].set_title("TSO source function")
	ax[1].set_title("MANCHA source function")
	
	
	
	f
end

# ╔═╡ 61ae8511-874f-4fe4-854e-d209ef194337


# ╔═╡ efa1b8bc-20fb-4759-8187-cc5f84bca6d7
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	rho_target = log(1e-12)
	
	rho = similar(eos.lnEi)
	fill!(rho, rho_target)
	lnei = eos.lnEi
	kappa = lookup(eos, opa, :κ, rho, lnei)
	T = exp.(lookup(eos, :lnT, rho, lnei))
	ax.plot(T, log10.(kappa[:, bin]), label="TSO")

	rho = similar(eosMancha.lnEi)
	fill!(rho, rho_target)
	lnei = eosMancha.lnEi
	kappaMancha = lookup(eosMancha, opaMancha, :κ, rho, lnei)
	T = exp.(lookup(eosMancha, :lnT, rho, lnei))
	ax.plot(T, log10.(kappaMancha[:, bin]), label="Mancha")


	_, T = profile(MUST.maximum, b, :z, :T)
	ax.axvline(last(T), color="cyan")
	_, T = profile(MUST.mean, b, :z, :T)
	ax.axvline(last(T), color="k")
	_, T = profile(MUST.minimum, b, :z, :T)
	ax.axvline(last(T), color="magenta")

	ax.set_ylabel("opacity (bin $(bin))")
	ax.set_xlabel("Temperature")

	ax.set_xscale("log")
	ax.set_title("rho $(exp(rho_target))")
	
	ax.legend()
	
	f
end

# ╔═╡ 1ef3f475-0423-4584-b6bc-55b2d85fc176
let
	f, ax = plt.subplots(1, 1, sharex=true, sharey=true)

	T_target = 4000.0
	
	rho = eos.lnRho
	T = similar(rho)
	fill!(T, T_target)
	lnei = lookup(eos, :lnEi, rho, log.(T))
	kappa = lookup(eos, opa, :κ, rho, lnei)
	ax.plot(log10.(exp.(rho)), log10.(kappa[:, bin]), label="TSO")

	rho = eosMancha.lnRho
	T = similar(rho)
	fill!(T, T_target)
	lnei = lookup(eosMancha, :lnEi, rho, log.(T))
	kappa = lookup(eosMancha, opaMancha, :κ, rho, lnei)
	ax.plot(log10.(exp.(rho)), log10.(kappa[:, bin]), label="Mancha")

	_, d = profile(MUST.maximum, b, :z, :log10d)
	ax.axvline(last(d), color="cyan")
	_, d = profile(MUST.mean, b, :z, :log10d)
	ax.axvline(last(d), color="k")
	_, d = profile(MUST.minimum, b, :z, :log10d)
	ax.axvline(last(d), color="magenta")

	ax.set_ylabel("opacity (bin $(bin))")
	ax.set_xlabel("density")
	ax.set_title("temperature $(T_target)")
	
	ax.legend()
	
	f
end

# ╔═╡ 31de6658-bcc6-4e71-a035-d93d99250782


# ╔═╡ 8fb63802-497f-42c4-a450-356c74905397
let
	f, ax = plt.subplots(1, 1)
	
	_, ee = profile(MUST.mean, b, :z, :logee)
	_, d = profile(MUST.mean, b, :z, :logd)
	_, T = profile(MUST.mean, b, :z, :T)
	kappamean = log10.(lookup(eos, opa, :κ, d, ee))
	kappameanMancha = log10.(lookup(eosMancha, opaMancha, :κ, d, ee))
	ax.plot(
		log10.(exp.(d)), kappamean[:, bin], 
		color="k", ls="--", label="TSO mean"
	)
	ax.plot(
		log10.(exp.(d)), kappameanMancha[:, bin], 
		color="r", ls="--", label="Mancha mean"
	)
	

	_, ee = profile(MUST.maximum, b, :z, :logee)
	_, d = profile(MUST.minimum, b, :z, :logd)
	_, T = profile(MUST.maximum, b, :z, :T)
	kappaextream = log10.(lookup(eos, opa, :κ, d, ee))
	kappaextreamMancha = log10.(lookup(eosMancha, opaMancha, :κ, d, ee))
	ax.plot(
		log10.(exp.(d)), kappaextream[:, bin], 
		color="k", label="TSO extream"
	)
	ax.plot(
		log10.(exp.(d)), kappaextreamMancha[:, bin], 
		color="r", label="Mancha extream"
	)

	ax.set_ylabel("opacity (bin $(bin))")
	ax.set_xlabel("density")

	
	ax.legend()
	
	f
end

# ╔═╡ Cell order:
# ╠═d64c7fa0-8212-11ef-2afe-7f449eb43a08
# ╠═42039451-14e3-4f79-94fa-47a2b6b0684b
# ╠═1383ae68-993e-442a-8893-137470156af7
# ╟─8c85a22a-7aa7-420e-a5dc-dbe59ab5d7c7
# ╠═8617f4a8-fecc-4186-a6c0-2d22a92a0183
# ╟─9bec4de7-1a19-49f7-a872-cfabddd620df
# ╠═b33169b4-ca20-4d10-b737-845df685fe5b
# ╠═28372005-4de3-4960-aea3-f03cff582f78
# ╠═58bba96b-6516-46a4-9993-cf52e722b3f9
# ╠═228fae7a-3089-4331-8066-566862ac7864
# ╟─b1d4e175-bb3f-44e9-898c-ab41a068bcf8
# ╠═d55e97a8-4779-4834-a8fd-58ecd4f409ff
# ╟─f9a52b86-9a40-473c-b5b8-7dc9e0de82e0
# ╟─fd9999e7-29ff-4538-b1dd-5359a5efa6e7
# ╠═db200f5d-7ce5-4031-a552-dbd86458591c
# ╟─3e8bb05f-600c-4efa-9a5f-22f8d5c876ea
# ╟─857a0709-706f-4763-a6c0-1627be747ac2
# ╟─12e2ff3b-3b0a-4976-8e27-60d7ee040c1b
# ╠═21d36112-42c2-4b1a-943f-def2d65a0eaf
# ╟─5b925a3c-ec86-49a3-85d7-736436a58f21
# ╠═a49805e9-1cc9-4da4-a7ac-d8cd17a31b6d
# ╠═74ad6b6d-13a2-48a6-b739-60e05172514b
# ╟─d3061189-067a-4127-b027-b9e85a3c94ea
# ╠═d04b9cc1-99ed-4d6b-b9d6-5ff270795cd2
# ╠═2cd909d3-76c1-4015-aa79-c99352a4a0be
# ╠═2cd3ef34-1478-42cb-ae46-5368a9744251
# ╠═130a2880-5166-44a6-956a-e58fe112bcfb
# ╟─14787aa7-de6a-4b09-be21-fe6d516d5a3a
# ╟─1f2d74ac-64e3-4eaa-bb25-8fd884e3c6ab
# ╟─0077c25e-5ec6-4e49-9378-4cd4f9867170
# ╟─813da454-74a0-492e-8092-dddeff054ad1
# ╠═654f2da6-e5d0-4f37-8a3b-563f5e8f0753
# ╠═d040380d-e975-4667-90a1-5de76d02e41a
# ╟─f8ff0bb8-04d4-4d02-8380-c2eec52c170f
# ╟─54a13656-f36d-47eb-bec8-40781968c473
# ╟─38bbc026-1597-4e96-8d48-5432a4a78cf5
# ╟─77e9ae8c-eee7-4495-9f36-1d7e9c8bab40
# ╟─8decebd5-2716-4568-bd2e-447f63128c9f
# ╟─445e3b6e-0d48-43ac-9518-1c28a0725e1d
# ╟─1b21ac59-4977-43ac-a8ff-64424b0f0bc9
# ╟─dcf83225-72c5-4938-ada3-0784301e97aa
# ╟─7024d5ab-4684-4fbe-8076-31021a06c7cf
# ╟─9aa170e3-21d0-4fc5-b76c-bdd5be28bbba
# ╟─227a1613-7281-4978-90f2-9a7c0df6528a
# ╟─e8490b9b-bdcb-4e19-b51a-c174767cfcfb
# ╟─1ded5ad1-a76d-4aeb-9d0b-b2a8b9e0c8aa
# ╟─437706ae-8307-4e4b-9bc2-2cdc93eaa5d9
# ╟─61ae8511-874f-4fe4-854e-d209ef194337
# ╟─efa1b8bc-20fb-4759-8187-cc5f84bca6d7
# ╟─1ef3f475-0423-4584-b6bc-55b2d85fc176
# ╟─31de6658-bcc6-4e71-a035-d93d99250782
# ╟─8fb63802-497f-42c4-a450-356c74905397
