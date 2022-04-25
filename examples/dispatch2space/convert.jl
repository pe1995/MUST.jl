using Pkg; Pkg.activate("."); #Pkg.update(); 
using MUST
using Glob
using PyCall

MUST.@import_dispatch "../../../dispatch2"
MUST.@import_dispatch "../../../dispatch2" EOS select

folder            = MUST.@in_dispatch ARGS[1]
content_of_folder = glob("*/", folder)
snapshots         = sort(MUST.list_of_snapshots(content_of_folder))

# Different EOS
eos        = dispatch.EOS.bifrost(top=MUST.dispatch_location, dir="data/eos/stagger/");
eos_legacy = dispatch.EOS.stagger(top=MUST.dispatch_location, file="table.dat");

for i_s in 1:length(snapshots)
    @info "Converting snapshot $(snapshots[i_s])/$(snapshots[end])"

    try
        # The dispatch snapshot object (Python)
        snap = dispatch.snapshot(snapshots[i_s], data=folder)

        # Convert its content to pure Julia
        s = MUST.Space(snap, :d, :ee, :uz, :e)

        # Add quantities from the EOS
        MUST.add_from_EOS!(s, eos_legacy, :K1, EOS_paras=(:d,:ee), EOS_log=(:true,:false))
        MUST.add_from_EOS!(s, eos_legacy, :K2, EOS_paras=(:d,:ee), EOS_log=(:true,:false))
        MUST.add_from_EOS!(s, eos_legacy, :K3, EOS_paras=(:d,:ee), EOS_log=(:true,:false))
        MUST.add_from_EOS!(s, eos_legacy, :K4, EOS_paras=(:d,:ee), EOS_log=(:true,:false))
        MUST.add_from_EOS!(s, eos_legacy, :T,  EOS_paras=(:d,:ee), EOS_log=(:true,:false), switch_name_to=:temp)

        # Units for conversion to CGS
        units = MUST.StaggerCGS(snap)

        # Apply the conversion
        MUST.convert!(s, units; d=:d,   ee=:ee, uz=:u,  e=:e, 
                                K1=:rk, K2=:rk, K3=:rk, K4=:rk,
                                x=:l,   y=:l,   z=:l)

        # Add additional columns already in CGS after converting
        MUST.add_from_EOS!(s, eos, :T)
        MUST.add_from_EOS!(s, eos, :rkap)

        # Also save the snapshot as Box (a regular gridded 3D-cube) to save time later
        b_s = MUST.Box(s)

        # Write to HDF5 file. Can easily be read as a Memory map later with the build in functions
        MUST.save(s;   name="space_sn$(snapshots[i_s])", folder=folder)
        MUST.save(b_s; name="box_sn$(snapshots[i_s])",   folder=folder)
    catch
        @warn "snapshot $(snapshots[i_s]) could not be loaded."
        continue
    end

end