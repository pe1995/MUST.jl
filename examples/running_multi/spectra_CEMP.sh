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
julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C3.1" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_3.1 --linelists="CHonly" --move -t 32 --nnu=64
julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.9" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.9 --linelists="CHonly" --move -t 32 --nnu=64
julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.8" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.8 --linelists="CHonly" --move -t 32 --nnu=64
julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.7" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.7 --linelists="CHonly" --move -t 32 --nnu=64
julia spectrum.jl -r "M2_E_t57.77g44.40m-4.000_v1.0" --extension="_C2.6" -n 399 -x 30 -s 4297 -e 4303 -d 0.001 --feh=-5.0 --alpha=0.4 --composition=C_2.6 --linelists="CHonly" --move -t 32 --nnu=64
