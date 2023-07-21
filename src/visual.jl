using Glob
using MUST
using PythonPlot
using Printf
using PythonCall

plt = matplotlib.pyplot

#using Plots

const visual = true
const cbar_fraction = 0.046
const cbar_pad = 0.04

gifs = MUST.@get_help_py gifs

#function copy_ticks(sp::Plots.Subplot; minorticks=10)
#    ptx = twinx(sp)
#    plot!(ptx,xlims=xlims(sp),ylims=ylims(sp),xformatter=_->"",yformatter=_->"", minorticks=10)
#    pty = twiny(sp)
#    plot!(pty,xlims=xlims(sp),ylims=ylims(sp),xformatter=_->"",yformatter=_->"", minorticks=10)
#end
#copy_ticks(plt::Plots.Plot = current(); minorticks=10) = copy_ticks(plt[1], minorticks=10)
#
#
#function basic_plot!(plot::Plots.Subplot=current(); minorticks=10, tickfontsize=12, legendfontsize=12, guidefontsize=12, size=(900,600), lm=5, rm=5)
#    copy_ticks(plot, minorticks=15)
#    plot!(plot, framestyle=:box, minorticks=minorticks, tickfontsize=tickfontsize, legendfontsize=legendfontsize, 
#            grid=false, foreground_color_legend=nothing, size=size, leftmargin=lm*Plots.mm, bottommargin=rm*Plots.mm, guidefontsize=guidefontsize)
#end



"""
    cube_with_velocities(m_3d, var=:T; vmin_3d=minimum(m_3d[var]), 
                        vmax_3d=maximum(m_3d[var]),
                        s_3d=1.1,
                        figsize=(10, 7),
                        limx=-0.4,
                        limy=0.4,
                        limz=-1.05,
                        cmap="rainbow",
                        len_vec=0.4,
                        lw_vec=0.6,
                        cvec="k",
                        arrow_length_ratio=0.35,
                        skipv=8,
                        xoff=10, 
                        yoff=15, 
                        zoff=10, 
                        show_time=false)

Create a python plot containing the 3D cube, with a cutout for velocities.
"""
function cube_with_velocities(m_3d, var=:T; vmin_3d=minimum(m_3d[var]), 
        vmax_3d=maximum(m_3d[var]),
        s_3d=1.1,
        figsize=(10, 7),
        limx=-0.4,
        limy=0.4,
        limz=-1.05,
        cmap="rainbow",
        clabel="T [K]",
        len_vec=0.4,
        lw_vec=0.6,
        cvec="k",
        arrow_length_ratio=0.35,
        skipv=8,
        xoff=10, 
        yoff=15, 
        zoff=10, 
        show_time=false)

    # Define dimensions
    Nx, Ny, Nz = size(m_3d)
    X, Y, Z    = MUST.mesh(m_3d)

    xmin, xmax = minimum(X), maximum(X)
    ymin, ymax = minimum(Y), maximum(Y)
    zmin, zmax = minimum(Z), maximum(Z)

    data = deepcopy(m_3d[var])

    # Create a figure with 3D ax
    fig_3d = plt.figure(figsize=figsize)
    ax_3d = fig_3d.add_subplot(111, projection="3d")


    mask_x1 = (X .>= limx) .& (Z .>= limz) .& (Y .<= limy)
    data[mask_x1] .= NaN


    # The cube faces
    im_3d = ax_3d.scatter(X[:,   1,  :],   
    Y[:,   1, :],   
    Z[:,   1, :],   
    c=data[:,   1, :],   
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap)
    im_3d = ax_3d.scatter(X[end, :,  :],   
    Y[end, :, :],   
    Z[end, :, :],   
    c=data[end, :, :],   
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap)
    im_3d = ax_3d.scatter(X[:,   :,  end], 
    Y[:,   :, end], 
    Z[:,   :, end], 
    c=data[:,   :, end], 
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap)




    # Draw the cut out part
    i_startx = findfirst(X[:, 1, 1] .> limx) 
    i_starty = findlast(Y[1, :, 1]  .<= limy)
    i_startz = findfirst(Z[1, 1, :] .> limz) 

    im_3d = ax_3d.scatter(  X[i_startx,   1:i_starty, i_startz:end],   
    Y[i_startx,   1:i_starty, i_startz:end],   
    Z[i_startx,   1:i_starty, i_startz:end],   
    c=m_3d[var][i_startx,   1:i_starty, i_startz:end],   
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap,zorder=0)

    im_3d = ax_3d.scatter(  X[i_startx:end,   1:i_starty, i_startz],   
    Y[i_startx:end,   1:i_starty, i_startz],   
    Z[i_startx:end,   1:i_starty, i_startz],   
    c=m_3d[var][i_startx:end,   1:i_starty, i_startz],   
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap,zorder=0)

    im_3d = ax_3d.scatter(  X[i_startx:end,   i_starty, i_startz:end],   
    Y[i_startx:end,   i_starty, i_startz:end],   
    Z[i_startx:end,   i_starty, i_startz:end],   
    c=m_3d[var][i_startx:end,   i_starty, i_startz:end],   
    s=s_3d, vmin=vmin_3d, vmax=vmax_3d, cmap=cmap,zorder=0)



    # Add the velocity field on top
    ux = deepcopy(m_3d[:ux] ./ maximum(abs.(m_3d[:ux])))
    uy = deepcopy(m_3d[:uy] ./ maximum(abs.(m_3d[:uy])))
    uz = deepcopy(m_3d[:uz] ./ maximum(abs.(m_3d[:uz])))	



    zp = fill!(similar(uz), 0.0)
    #i_startx = i_startx +10
    #i_starty = i_starty -20
    #i_startz = i_startz +10


    im_3d_vec = ax_3d.quiver(   X[i_startx+xoff,1:skipv:i_starty, i_startz:skipv:end],
    Y[i_startx+xoff,1:skipv:i_starty, i_startz:skipv:end],
    Z[i_startx+xoff,1:skipv:i_starty, i_startz:skipv:end],
    zp[i_startx, 1:skipv:i_starty, i_startz:skipv:end],
    uy[i_startx, 1:skipv:i_starty, i_startz:skipv:end],
    uz[i_startx, 1:skipv:i_starty, i_startz:skipv:end],
    length=len_vec, normalize=false, lw=lw_vec, 
    pivot="tail", arrow_length_ratio=arrow_length_ratio,
    color=cvec, zorder=100)

    im_3d_vec = ax_3d.quiver( X[i_startx:skipv:end,  1:skipv:i_starty, i_startz+zoff],
    Y[i_startx:skipv:end,  1:skipv:i_starty, i_startz+zoff],
    Z[i_startx:skipv:end,  1:skipv:i_starty, i_startz+zoff],
    ux[i_startx:skipv:end, 1:skipv:i_starty, i_startz],
    uy[i_startx:skipv:end, 1:skipv:i_starty, i_startz],
    zp[i_startx:skipv:end, 1:skipv:i_starty, i_startz],
    length=len_vec, normalize=false, lw=lw_vec, 
    pivot="tail", arrow_length_ratio=arrow_length_ratio,
    color=cvec, zorder=100)

    im_3d_vec = ax_3d.quiver(X[i_startx:skipv:end, i_starty-yoff, i_startz:skipv:end],
    Y[i_startx:skipv:end,  i_starty-yoff, i_startz:skipv:end],
    Z[i_startx:skipv:end,  i_starty-yoff, i_startz:skipv:end],
    ux[i_startx:skipv:end, i_starty, i_startz:skipv:end],
    zp[i_startx:skipv:end, i_starty, i_startz:skipv:end],
    uz[i_startx:skipv:end, i_starty, i_startz:skipv:end],
    length=len_vec, normalize=false, lw=lw_vec, color=cvec,
    pivot="tail", arrow_length_ratio=arrow_length_ratio, 
    zorder=100)


    # Set limits of the plot from coord limits
    ax_3d.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])


    # Plot edges (only until start of the second cube)
    xmax_lim = X[i_startx, 1, 1]
    ymax_lim = Y[1, i_starty, 1]
    zmax_lim = Z[1, 1, i_startz]


    ax_3d.plot([xmax, xmax], [ymax_lim, ymax], zmax, 
    color="k", linewidth=1, zorder=1e3)
    ax_3d.plot([xmax, xmax], [ymax_lim, ymax], 0.0, 
    color="k", linewidth=1, zorder=1e3, alpha=0.5, ls=":")
    ax_3d.plot([xmax, xmax], [ymax_lim, ymax], zmax_lim, 
    color="k", linewidth=1, zorder=1e3)


    ax_3d.plot([xmin, xmax_lim], [ymin, ymin], zmax, 
    color="k", linewidth=1, zorder=1e3)
    ax_3d.plot([xmin, xmax_lim], [ymin, ymin], 0.0, 
    color="k", linewidth=1, zorder=1e3, alpha=0.5, ls=":")
    ax_3d.plot([xmin, xmax_lim], [ymin, ymin], zmax_lim, 
    color="k", linewidth=1, zorder=1e3)


    ax_3d.plot([xmax_lim, xmax_lim], [ymin, ymax_lim], zmax, 
    color="k", linewidth=1, zorder=1e3)
    ax_3d.plot([xmax_lim, xmax_lim], [ymin, ymax_lim], 0.0, 
    color="k", linewidth=1, zorder=1e3, alpha=0.5, ls=":")
    #ax_3d.plot([xmax_lim, xmax_lim], [ymin, ymax_lim], zmax_lim, 
    #	color="k", linewidth=1, zorder=1e3)


    ax_3d.plot([xmax_lim, xmax], [ymax_lim, ymax_lim], zmax, 
    color="k", linewidth=1, zorder=1e3)
    ax_3d.plot([xmax_lim, xmax], [ymax_lim, ymax_lim], 0.0, 
    color="k", linewidth=1, zorder=1e3, alpha=0.5, ls=":")
    #ax_3d.plot([xmax_lim, xmax], [ymax_lim, ymax_lim], zmax_lim, 
    #	color="k", linewidth=1, zorder=1e3)


    ax_3d.plot([xmax_lim, xmax], [ymin, ymin], zmax_lim, 
    color="k", linewidth=1, zorder=1e3)

    ax_3d.plot([xmax, xmax], [ymin, ymax_lim], zmax_lim, 
    color="k", linewidth=1, zorder=1e3)



    # Set labels and zticks
    ax_3d.set(
    ylabel="Y [Mm]",
    xlabel="X [Mm]",
    zlabel="Z [Mm]")

    # Set zoom and angle view
    #ax.view_init(40, -30, 0)
    #ax.set_box_aspect(None, zoom=0.9)

    norm = matplotlib.colors.Normalize(vmin=vmin_3d, vmax=vmax_3d)
    sm   = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    #sm.set_array([]) 

    # Colorbar
    plt.colorbar(sm, ax=ax_3d, label=clabel)

    ax_3d.grid(false)

    # make the panes transparent
    ax_3d.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax_3d.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax_3d.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    if show_time
        ts = @sprintf("%i", m_3d.parameter.time)
        ax_3d.set_title("solar time: $(ts) min", fontsize="x-large")
    end

    fig_3d, ax_3d
end

function colors_map(cmap, n)
    c = plt.get_cmap(cmap, n)
    [c((i-1)/n) for i in 1:n]
end

function basic_subplot(nrows,ncols;
                        figsize=(6,6),
                        wspace=.1,hspace=.2,
                        label_fs="medium",
                        frame_linewidth=1.1,
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
        label_fs="medium",
        frame_linewidth=1.05,
        tick_length=4.4,
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

basic_colorbar!(cax) = begin
    basic_plot!(cax.ax)
end

add_cbar(im, ax, args...; fraction=cbar_fraction, pad=cbar_pad, kwargs...) = begin
    cax = plt.colorbar(im, ax=ax, args..., fraction=fraction, pad=pad, kwargs...)
    basic_plot!(cax.ax)

    cax
end

function gif_by_plane(stat::Function, folder::Vector{String}, labels::Vector{String}; 
                                    color=["k" for _ in folder], ls=["-" for _ in folder],
                                    ylim=(4000, 22000), duration=0.2,
                                    ylabel="", xlabel="z [cm]",ax_start=nothing, f_start=nothing,
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
        
        if isnothing(ax_start)
            f, ax = plt.subplots(1,1)
        else
            f, ax = deepcopy(f_start), deepcopy(ax_start)
        end
        
        ax.set_title("snapshot $(i)")
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.set_ylim(ylim...)

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
                ax.plot(zs, stats, color=color[j],label=labels[j], ls=ls[j])
            end
        end
        ax.legend()
        f.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close(f)
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gifs.gif_from_png(filenames, "$(path_ext).gif", duration=duration)
end

function gif_by_column(f, boxes, variable; 
                        vmin=nothing, vmax=nothing, cmap="Greys_r", clabel="$(variable)",
                        path_ext="box_val", duration=0.2, transformers...)

    filenames = String[]

    for i in eachindex(boxes)
        plt.title("snapshot $(i)")
        bv = boxes[i]
        surface = MUST.reduce_by_column(f, bv) 
        v = view(surface.data[variable], :, :, 1)

        v = length(keys(transformers)) > 0 ? transformers[variable].(v) : v
        if !isnothing(vmin)
            im = plt.imshow(v, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
            cb = plt.colorbar(im, boundaries=[range(vmin, vmax; length=100)...])
        else
            im = plt.imshow(v, origin="lower", cmap=cmap)
            cb = plt.colorbar(im)
        end

        cb.set_label(clabel)
        plt.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close()
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gifs.gif_from_png(filenames, "$(path_ext).gif", duration=duration)
end

function gif_by_yindex(yindex, boxes, variable; 
    vmin=nothing, vmax=nothing, cmap="Greys_r", clabel="$(variable)",
    path_ext="box_val", duration=0.2, zrange=nothing, transformers...)

    filenames = String[]

    for i in eachindex(boxes)
        plt.title("snapshot $(i)")
        surface = boxes[i]

        v = isnothing(zrange) ? view(surface.data[variable], :, yindex, :)' : view(surface.data[variable], :, yindex, zrange)'
        z = isnothing(zrange) ? view(surface.z, 1, yindex, :) : view(surface.z, 1, yindex, zrange)
        x = view(surface.x, :, yindex, 1) 
    
        minx, maxx = minimum(x)*1e-8, maximum(x)*1e-8
        minz, maxz = minimum(z)*1e-8, maximum(z)*1e-8
        extent = [minx, maxx, minz, maxz]

        v = length(keys(transformers)) > 0 ? transformers[variable].(v) : v
        if !isnothing(vmin)
            im = plt.imshow(v, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto", extent=extent)
            cb = plt.colorbar(im, boundaries=[range(vmin, vmax; length=100)...])
        else
            im = plt.imshow(v, origin="lower", cmap=cmap, aspect="auto", extent=extent)
            cb = plt.colorbar(im)
        end

        cb.set_label(clabel, fontsize="large")
        plt.xlabel("x [Mm]", fontsize="large")
        plt.ylabel("z [Mm]", fontsize="large")
        plt.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close()
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gifs.gif_from_png(filenames, "$(path_ext).gif", duration=duration)
end

function gif_by_value(stat::Function, folder::String, label::String; 
                                    cmap="Greys_r",
                                    duration=0.2, vmin=-99999999, vmax=-99999999,
                                    names="box", variable=:temp, path_ext="box_val", 
                                    clabel="", 
                                    kwargs...)
    iw    = 0
    n_snaps   = length(list_of_snapshots(folder, names))
    filenames = String[]

    for i in list_of_snapshots(folder, names)
        plt.title("$(label), snapshot $(i)")

        try
            b  = MUST.Box("$(names)_sn$(i)", folder=MUST.@in_dispatch(folder))
            bv = MUST.reduce_by_value(stat, b; kwargs...)

            im = ((vmin != -99999999) & (vmax != -99999999)) ? 
                        plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax) :
                        plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap)

            cb   = plt.colorbar(im)
            cb.set_label(clabel)

            iw = i
        catch
            if iw == 0
                continue
            else
                b  = MUST.Box("$(names)_sn$(iw)", folder=MUST.@in_dispatch(folder))
                bv = MUST.reduce_by_value(stat, b; kwargs...)

                im = ((vmin != -99999999) & (vmax != -99999999)) ? 
                        plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax) :
                        plt.imshow(bv.data[variable][:,:,1], origin="lower", cmap=cmap)
                        
                cb   = plt.colorbar(im)
                cb.set_label(clabel)
            end
        end
    
        plt.savefig("$(path_ext)_sn$(i).png", bbox_inches="tight")        
        plt.close()
        append!(filenames, ["$(path_ext)_sn$(i).png"])
    end

    gifs.gif_from_png(filenames, "$(path_ext).gif", duration=duration)
end

"""
Look for files starting with name in folder
"""
list_of_snapshots(folder::String, name::String) = begin
    f    = glob("$(name)*.hdf5", folder)
    sort([parse(Int64, fi[last(findlast("sn", fi))+1:end-5]) for fi in f])
end

;