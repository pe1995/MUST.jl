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