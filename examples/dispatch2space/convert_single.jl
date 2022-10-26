using Pkg; Pkg.activate("."); #Pkg.update(); 
using MUST
using Glob
using PyCall
using Interpolations

MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2"
MUST.@import_dispatch "/u/peitner/DISPATCH/dispatch2" EOS select

add_selection = true

folder = MUST.@in_dispatch ARGS[1]
content_of_folder = glob("*/", folder)
if length(ARGS) > 1
    snapshots = sort(parse.(Int, ARGS[2:end]))
    #snapshots = sort(MUST.list_of_snapshots(content_of_folder))[sidx]
else
    content_of_folder = glob("*/", folder)
    snapshots         = sort(MUST.list_of_snapshots(content_of_folder))#[1:3:end]
end

# Name of the namelist of the current folder
nml_name = MUST.@in_dispatch splitpath(folder)[end]

# Init namelist
const nml = MUST.StellarNamelist(nml_name*".nml")

# Use the new Squaregas EOS 
eos_path = replace(nml.eos_params["table_loc"], "'"=>"")
const eos_sq = MUST.SquareGasEOS(MUST.@in_dispatch(eos_path))

# Effective temperature present?
save_info = ispath(joinpath(folder, "teff.dat")) 
if save_info
    # Read and interpolate teff data
    teff_func = MUST.teff_interpolated(joinpath(folder, "teff.dat"))
end

function convert_snapshots(snapshots)
    for i_s in eachindex(snapshots)
        @info "Converting snapshot $(snapshots[i_s])"
        try
            # The dispatch snapshot object (Python)
            snap = dispatch.snapshot(snapshots[i_s], data=folder)

            # Units for conversion to CGS
            units = MUST.StaggerCGS(snap)

            # Convert its content to pure Julia
            s = if add_selection
                MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e, :qr, :tt, :pg)
            else
                MUST.Space(snap, :d, :ee, :ux, :uy, :uz, :e)
            end

            # Apply the conversion
            if add_selection
                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    qr=:qr, pg=:p,
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l)
            else
                MUST.convert!(s, units; d=:d, ee=:ee, e=:e, 
                                    ux=:u, uy=:u, uz=:u,
                                    x=:l, y=:l, z=:l)
            end

            # Add additional columns already in CGS after converting
            MUST.add_from_EOS!(s, eos_sq, :T)
            MUST.add_from_EOS!(s, eos_sq, :kr)

            # Check if Teff should be added
            if save_info
                @info "Setting Teff for snapshot $(snapshots[i_s])"
                MUST.set!(s.parameter, teff_func, nml)
            end

            # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
            b_s = MUST.Box(s)
            s = nothing

            # Add the optical depth
            τ = MUST.optical_depth(b_s, opacity=:kr, density=:d)
            MUST.add!(b_s, τ, :τ_ross)

            # First save
            MUST.save(b_s; name="box_sn$(snapshots[i_s])",   folder=folder)

            # NEW: Convert the height scale from cm to optical depth
            min_plane = MUST.plane_statistic(minimum, b_s, :τ_ross)
            b_s       = MUST.height_scale(b_s, :τ_ross, Float32[maximum(min_plane), 10^(-6.0)])

            # Write to HDF5 file. Can easily be read as a Memory map later with the build-in functions
            #MUST.save(s;   name="space_sn$(snapshots[i_s])", folder=folder)
            MUST.save(b_s; name="box_tau_sn$(snapshots[i_s])",   folder=folder)

        catch
            @warn "snapshot $(snapshots[i_s]) could not be loaded."
            continue
        end
    end
    nothing
end

convert_snapshots(snapshots)


