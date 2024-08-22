# Metal-poor sun, CH G-Band. Both with and without CEMP-like enhancement. M2 -> [C/Fe] = 0.4, M1 -> [C/Fe] = 3.0
##julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0" -n 390 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 60 
##julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --linelists="CHonly" --move -t 60 


# CEMP models diagnostics: 1) CEMP model + CEMP spectrum 2) MP model + MP spectrum 3) CEMP model + MP spectrum 4) MP model + CEMP spectrum
# This is usefull to test the effect of model differences on the final spectrum
##julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0_long" --extension="CEMP-CEMP" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 64 --nnu=64
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0"      --extension="MP-MP"     -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4                     --linelists="CHonly" --move -t 32 
#julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0_long" --extension="CEMP-MP"   -n 357 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4                     --linelists="CHonly" --move -t 32 
##julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0"      --extension="MP-CEMP"   -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 64 --nnu=64


# For abundance corrections we need a couple of spectra
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C3.1" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.1 --linelists="CHonly" --move -t 32 --nnu=64
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.9" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.9 --linelists="CHonly" --move -t 32 --nnu=64
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.8" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.8 --linelists="CHonly" --move -t 32 --nnu=64
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.7" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.7 --linelists="CHonly" --move -t 32 --nnu=64
#julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.6" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.6 --linelists="CHonly" --move -t 32 --nnu=64


# CEMP models with MARCS grid parameters
##julia spectrum.jl -r "M3_E_t57.50g45.00m-4.000_v1.0" --extension="CEMP-CEMP" -n 398 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M4_E_t57.50g45.00m-4.000_v1.0" --extension="MP-CEMP"   -n 398 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64


# Test with 1D models
##julia spectrum.jl -r "M3_E_t57.50g45.00m-4.000_v1.0"      --onedimensional="av_sn398CEMP-CEMP" -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M4_E_t57.50g45.00m-4.000_v1.0"      --onedimensional="av_sn398MP-CEMP"   -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0_long" --onedimensional="av_sn399CEMP-CEMP" -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0"      --onedimensional="av_sn399MP-CEMP"   -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0_long" --onedimensional="av_sn399CEMP-CEMP" --extension="C3.1" -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.1 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r "M1_E_t57.77g44.40m-4.000_v1.0_long" --onedimensional="av_sn399CEMP-CEMP" --extension="C2.9" -n 398 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.9 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional="av_sn398MP-CEMP" --extension=_C3.1 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.1 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional="av_sn398MP-CEMP" --extension=_C2.9 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.9 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional="av_sn398MP-CEMP" --extension=_C2.8 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.8 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional="av_sn398MP-CEMP" --extension=_C2.7 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.7 --linelists="CHonly" --move -t 32 --nnu=64
##julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional="av_sn398MP-CEMP" --extension=_C2.6 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.6 --linelists="CHonly" --move -t 32 --nnu=64


# For effective temperature we can use the whole range, lower resoltution should be fine though
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=teff -n 399 -x 10 -s 1000 -e 100000 --lambda_n=50000 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --move -t 32 --nnu=64
#julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0      --extension=teff -n 399 -x 10 -s 1000 -e 100000 --lambda_n=50000 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --move -t 32 --nnu=64


# MARCS spectra
#julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional=p5500_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional=p5750_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M2_E_t57.77g44.40m-4.000_v1.0 --onedimensional=p6000_g+4.5_m0.0_t02_st_z-5.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5600g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5700g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5800g4.50z-5.00a+0.40c+3.00t02-v2008-pturb.mod -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.0 --linelists=CHonly --move -t 32 --nnu=64


# abundance corrections
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C3.1 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.1 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C3.2 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.2 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C3.3 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.3 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C3.4 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.4 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C3.5 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.5 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C4.1 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_4.1 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C4.3 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_4.3 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C4.5 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_4.5 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --onedimensional=5772g4.44z-5.00a+0.40c+3.00t02-v2008-pturb.mod --extension=_C4.7 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_4.7 --linelists=CHonly --move -t 32 --nnu=64

#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=_C0.0 -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_0.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=_C2.0 -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.0 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=_C2.5 -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.5 --linelists=CHonly --move -t 32 --nnu=64
#julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=_C2.75 -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.75 --linelists=CHonly --move -t 32 --nnu=64
julia spectrum.jl -r M1_E_t57.77g44.40m-4.000_v1.0_long --extension=_C2.25 -n 399 -x 20 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.25 --linelists=CHonly --move -t 32 --nnu=64
