struct AuxPatch
    version
    id
    name
    rnk
    shp
    tp
    data
end


struct PatchMeta
    id
    mv
    rec
    nv
    ncell
    li
    ui
    gn
    n
    ng
    offset
    s
    pos
    llc
    urc
    ds
    x
    y 
    z
    xi
    yi 
    zi
end






#= Construction =#

"""
    PatchMeta(patch, snap)

Read the meta data from the given patch and snap namelists.
"""
function PatchMeta(patch, snap)
    guards = nmlValue(snap, "guard_zones")
    n = parse_multiple(nmlValue(snap, "n"))
    
    li, ui = if !guards
        [1, 1, 1], n
    else
        (
            parse_multiple(nmlValue(snap, "li")),
            parse_multiple(nmlValue(snap, "ui"))
        )
    end

    mv = nmlValue(snap, "mv")
    nv = nmlValue(snap, "nv")
    rec = nmlValue(patch, "record")
    ncell = parse_multiple(nmlValue(patch, "ncell"))
    gn = parse_multiple(nmlValue(snap, "gn"))
    ng = parse_multiple(nmlValue(snap, "ng"))
    bts = if guards
        prod(gn)
    else
        prod(n)
    end
    offset = [iv + (rec-1) * mv for iv in 0:mv-1] .* (4*bts)
    id = nmlValue(patch, "id")

    # for axes
    s = parse_multiple(nmlValue(patch, "size"))
    pos = parse_multiple(nmlValue(patch, "position"))
    llc= pos .- s./2.0
    urc= pos .+ s./2.0

    ds = parse_multiple(nmlValue(patch, "ds"))

    x, y, z = _add_axes(llc, ng, gn, ds)
    xi, yi, zi = x[li[1]:ui[1]], y[li[2]:ui[2]], z[li[3]:ui[3]]

    PatchMeta(
        id,
        mv,
        rec,
        nv,
        ncell,
        li,
        ui,
        gn,
        n,
        ng,
        offset,
        s,
        pos,
        llc,
        urc,
        ds,
        x,
        y,
        z,
        xi|>collect,
        yi|>collect,
        zi|>collect
    )
end

"""
    readaux(fname)

Read auxiliar data from the given filename.
"""
function readaux(fname)
    auxs = []
    if isfile(fname)
        f = FortranFiles.FortranFile(fname)
        version, id = read(f, (Int32, 2))
        while true
            try			
                name = read_string(f)
                rnk = read_numbers(f)
                shp = read_numbers(f)
                tp = read_string(f)
        
                data = if first(tp) == 'r'
                    a = zeros(Float32, shp...)
                    read(f, a)
                    a
                elseif first(tp) == 'd'
                    a = zeros(Float64, shp...)
                    read(f, a)
                    a
                else
                    a = zeros(Int32, shp...)
                    read(f, a)
                    a
                end

                append!(auxs, [AuxPatch(version, id, name, rnk, shp, tp, data)])
            catch
                break
            end
        end
        close(f)
    end
    auxs
end







#= Helper functions =#

function parse_multiple(s::AbstractString)
    sComponents = []
    
    # split at commas
    sComma = strip.(split(s, ",", keepempty=false))

    # split for multiplication signs
    sMulti = [strip.(split(si, "*", keepempty=false)) for si in sComma]

    # flatten the list
    for component in sMulti
        c = if length(component) > 1
            repeat([last(component)], parse(Int, first(component)))
        else
            component
        end

        append!(sComponents, c)
    end

    # try to convert to number
    sComponentsParsed = []
    for component in sComponents
        c_new = if occursin(".", component)
            try 
                parse(Float64, component)
            catch
                component
            end
        else
            try 
                parse(Int, component)
            catch
                component
            end
        end

        append!(sComponentsParsed, [c_new])
    end
        
    sComponentsParsed
end

parse_multiple(s::AbstractArray) = s


function read_stepwise(f, dtype)
    rec = FortranFiles.Record(f)
    res = []
    while true
        try
            a = read(rec, dtype)
            append!(res, [a])
        catch
            break
        end
    end

    close(rec)
    res
end

function read_string(f)
    str_arr = read_stepwise(f, FortranFiles.FString{1})
    str_arr = String[FortranFiles.trimstring(a) for a in str_arr]
    join(str_arr)
end

function read_numbers(f, dtype=Int32)
    a = dtype[read_stepwise(f, dtype)...]
    if length(a) == 1
        first(a)
    else
        a
    end
end

function _add_axes(llc, ng, gn, ds)
	first = llc .- ng .* ds
	n = gn
   
	x=first[1] .+ ds[1] .* range(0, n[1]-1)
	y=first[2] .+ ds[2] .* range(0, n[2]-1)
	z=first[3] .+ ds[3] .* range(0, n[3]-1)

	x, y, z
end




