### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 9b235fc0-aec8-11ee-36f2-cf6c51e21b43
begin
	using Pkg; Pkg.activate(".")
	using MUST
	using TSO
	using Plots
end

# ╔═╡ ce08f64d-c891-4a30-ba0c-7a01a89a32b6
eos = reload(SqEoS, "DIS_MARCS_t50g40m00_v0.4/eos.hdf5")

# ╔═╡ 4ce6989e-a7c8-40ec-a8c2-4aca53ea3209
sip = MUST.scipy_interpolate

# ╔═╡ de2c2b8d-fdc3-49a4-8329-70a9503e13f7
np = MUST.numpy

# ╔═╡ b954762f-a8dd-410b-a336-d8a843360938
begin
	x = eos.lnT
	y = eos.lnEi[:, 20]
	y_macc = reverse(accumulate(MUST.min, reverse(y)))

	ip = MUST.Interpolations.extrapolate(
		MUST.Interpolations.interpolate(
			x, y_macc, MUST.Interpolations.SteffenMonotonicInterpolation()
		), 
		 MUST.Interpolations.Flat()
	)
	#ip = sip.UnivariateSpline(x, y_macc, k=1, s=0.05)

	x_new = range(minimum(x), maximum(x), length=1000)
	y_ip = ip.(x_new) #MUST.pyconvert.(Any, ip(x_new))
end

# ╔═╡ e873d3a7-13b5-43ea-aa94-70401c99a030
y_ip

# ╔═╡ b78bce1f-f3e3-4565-a3de-8e3a2ee9a157
begin
	plot(x, y)
	plot!(x_new, y_ip)
	#plot!(x, y_macc)
end

# ╔═╡ Cell order:
# ╠═9b235fc0-aec8-11ee-36f2-cf6c51e21b43
# ╠═ce08f64d-c891-4a30-ba0c-7a01a89a32b6
# ╠═4ce6989e-a7c8-40ec-a8c2-4aca53ea3209
# ╠═de2c2b8d-fdc3-49a4-8329-70a9503e13f7
# ╠═b954762f-a8dd-410b-a336-d8a843360938
# ╠═e873d3a7-13b5-43ea-aa94-70401c99a030
# ╠═b78bce1f-f3e3-4565-a3de-8e3a2ee9a157
