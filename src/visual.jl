using PyPlot
using Glob
using MUST
using PyCall

const visual = true

function basic_subplot(nrows,ncols;
                        figsize=(6,6),
                        wspace=.1,hspace=.2,
                        label_fs="x-large",
                        frame_linewidth=1.2,
                        tick_length=5.,
                        label_color="black",
                        keep_axis_color=false,
                        add_outer_axes=true,
                        turn_off=[],
                        sharex=false,sharey=false,
                        fig_kwargs=Dict())
    # Create the frame of a basic plot
    f,ax = PyPlot.subplots(nrows,ncols,figsize=(figsize[1]*ncols,figsize[2]*nrows),
                        sharex=sharex,sharey=sharey,fig_kwargs...)

    for (i,axi) in enumerate(ax)
        if i in turn_off
            axi.axis("off")
            continue
        end
        # frame width
        try
            for loc in ["top","left","right","bottom"]
                axi.spines[loc].set_linewidth(frame_linewidth)
                axi.spines[loc].set_color(label_color)
            end
        catch 
            try
                axi.grid(linewidth=frame_linewidth,color=label_color)
            catch
                nothing
            end
        end

        # Edit ticks
        axi.minorticks_on()
        axi.tick_params(which="major",direction="in",top=true,right=true,length=tick_length,labelsize=label_fs,width=frame_linewidth,colors=label_color)
        axi.tick_params(which="minor",direction="in",top=true,right=true,length=0.5*tick_length,width=frame_linewidth,colors=label_color)
        for lab in axi.get_xticklabels()
            if keep_axis_color
                lab.set_color("black")
            end
        end
        for lab in axi.get_yticklabels()
            if keep_axis_color
                lab.set_color("black")
            end
        end
    end
    plt.subplots_adjust(wspace=wspace,hspace=hspace)
    if add_outer_axes
        # add a big axis, hide frame
        ax_xframe = f.add_subplot(111, frameon=False)
        # hide tick and tick label of the big axis
        ax_xframe.tick_params(labelcolor="none", top=false, bottom=false, left=false, right=false)
        out_ax = ax_xframe

        return f,ax,out_ax
    end

    return f,ax
end

function basic_plot(;
                    figsize=(6,6),
                    label_fs="large",
                    frame_linewidth=1.1,
                    tick_length=5.,
                    label_color="black",
                    keep_axis_color=false,
                    tick_direction="in")

    # Create the frame of a basic plot
    f = plt.figure(figsize=figsize)
    ax = f.add_subplot(111)

    # Edit ticks
    ax.minorticks_on()
    ax.tick_params(which="major",direction=tick_direction,top=true,right=true,length=tick_length,labelsize=label_fs,width=frame_linewidth,colors=label_color)
    ax.tick_params(which="minor",direction=tick_direction,top=true,right=true,length=0.5*tick_length,width=frame_linewidth,colors=label_color)

    # frame width
    try
        for loc in ["top","left","right","bottom"]
            ax.spines[loc].set_linewidth(frame_linewidth)
            ax.spines[loc].set_color(label_color)
        end
    catch
        try
            ax.grid(linewidth=frame_linewidth,color=label_color)
        catch
            nothing
        end
    end

    for lab in ax.get_xticklabels()
        if keep_axis_color
            lab.set_color("black")
        end
    end
    for lab in ax.get_yticklabels()
        if keep_axis_color
            lab.set_color("black")
        end
    end

    return f,ax
end

function basic_plot!(ax;
        label_fs="large",
        frame_linewidth=1.1,
        tick_length=5.,
        label_color="black",
        keep_axis_color=false,
        tick_direction="in")

    # Edit ticks
    ax.minorticks_on()
    ax.tick_params(which="major",direction=tick_direction,top=true,right=true,length=tick_length,labelsize=label_fs,width=frame_linewidth,colors=label_color)
    ax.tick_params(which="minor",direction=tick_direction,top=true,right=true,length=0.5*tick_length,width=frame_linewidth,colors=label_color)

    # frame width
    try
        for loc in ["top","left","right","bottom"]
        ax.spines[loc].set_linewidth(frame_linewidth)
        ax.spines[loc].set_color(label_color)
        end
    catch
        try
            ax.grid(linewidth=frame_linewidth,color=label_color)
        catch
            nothing
        end
    end

    for lab in ax.get_xticklabels()
        if keep_axis_color
        lab.set_color("black")
        end
    end
    for lab in ax.get_yticklabels()
        if keep_axis_color
            lab.set_color("black")
        end
    end
end

function gif_by_plane(stat::Function, folder::Vector{String}, labels::Vector{String}; 
                                    color=["k" for _ in folder], ls=["-" for _ in folder],
                                    ylim=(4000, 22000), duration=0.2,
                                    ylabel="", xlabel="z [cm]",
                                    names=["box" for _ in folder], variable=:temp, path_ext="box_stat")
    n_curves = length(folder)
    @assert length(folder) == length(labels) 
    @assert length(folder) == length(names) 

    iw    = zeros(Int64, n_curves)
    stats = []
    zs    = []

    n_snaps   = argmax([length(list_of_snapshots(f, n)) for (f,n) in zip(folder, names)])

    filenames = String[]

    for i in list_of_snapshots(folder[n_snaps], names[n_snaps])
        plt.title("snapshot $(i)")
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.ylim(ylim...)

        for j in 1:n_curves
            try
                b     = MUST.Box("$(names[j])_sn$(i)", folder=MUST.@in_dispatch(folder[j]))
                stats = MUST.plane_statistic(stat, b, variable)
                zs    = b.z[1,1,:]

                iw[j] = i
            catch
                if iw[j] == 0
                    continue
                else
                    b     = MUST.Box("$(names[j])_sn$(iw[j])", folder=MUST.@in_dispatch(folder[j]))
                    stats = MUST.plane_statistic(stat, b, variable)
                    zs    = b.z[1,1,:]
                end
            end
        
            if iw[j] == 0
                continue
            else
                plt.plot(zs, stats, color=color[j],label=labels[j], ls=ls[j])
            end
        end
        plt.legend()
        plt.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close()
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gif_from_pngs(filenames, "$(path_ext).gif", duration=duration)
end

function gif_by_value(stat::Function, folder::String, label::String; 
                                    cmap="Greys_r",
                                    duration=0.2,
                                    names="box", variable=:temp, path_ext="box_val", 
                                    clabel="", 
                                    kwargs...)
    iw    = 0
    n_snaps   = length(list_of_snapshots(folder, names))
    filenames = String[]

    for i in list_of_snapshots(folder, names)
        plt.title("$(label), snapshot $(i)")

        try
            b     = MUST.Box("$(names)_sn$(i)", folder=MUST.@in_dispatch(folder))
            bv = MUST.reduce_by_value(stat, b; kwargs...)

            im   = plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap);
            cb   = plt.colorbar(im)
            cb.set_label(clabel)

            iw = i
        catch
            if iw == 0
                continue
            else
                b     = MUST.Box("$(names)_sn$(iw)", folder=MUST.@in_dispatch(folder))
                bv = MUST.reduce_by_value(stat, b; kwargs...)

                im   = plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap);
                cb   = plt.colorbar(im)
                cb.set_label(clabel)
            end
        end
    
        plt.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close()
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gif_from_pngs(filenames, "$(path_ext).gif", duration=duration)
end

"""
Look for files starting with name in folder
"""
list_of_snapshots(folder::String, name::String) = begin
    f    = glob("$(name)*.hdf5", folder)
    sort([parse(Int64, fi[last(findlast("sn", fi))+1:end-5]) for fi in f])
end

gif_from_pngs(list_of_filenames, save_path; duration=0.2) = py"""
import imageio
import glob
import os
import sys

images = []

for filename in $(list_of_filenames):
    images.append(imageio.imread(filename))

imageio.mimsave($(save_path), images, duration=$(duration))

for f in $(list_of_filenames):
    os.remove(f)
"""