using DelimitedFiles

path = ARGS[1]

function read_marcs(path)
    name = basename(abspath(path))
    dir = dirname(abspath(path))
    barename = split(name, ".out", keepempty=false) |> first
    new_path = joinpath(dir, barename*".mod")


    lines = readlines(path)
    # find key lines
    data = read_atmos(lines)
    thermo = read_thermo(lines)

    z = data["Depth"]
    T = data["T"]
    pe = data["Pe"]
    r = thermo["Density"]
    vmic = 2.0
    open(new_path, "w") do f
        write(f, name*"\n")
        write(f, "$(length(z))\n")
        for i in eachindex(z)
            write(f, "$(z[i]) $(T[i]) $(pe[i]) $(r[i]) $(vmic) \n")
        end
    end

    data, thermo
end

function read_atmos(lines)
    header = ["K","TauRoss","Tau500","Depth","T","Pe","Pg","Prad","Pturb","KappaRoss"]
    read_quant(lines, "M o d e l   A t m o s p h e r e", header)
end

function read_thermo(lines)
    header = ["K","TauRoss","Density","Mu","Cp","Cv","AdGrad","Q","Sound vel.","Conv. vel.","Fconv/F"]
    read_quant(lines, "T h e r m o d y n a m i c a l   q u a n t i t i e s  and   C o n v e c t i o n", header)
end

function read_quant(lines, q, header)
    imodelatmos = findfirst(l->occursin(q, l), lines)
    data = Dict(k=>[] for k in header)
    i = imodelatmos + 4
    while i < length(lines)
        try
            l = split(lines[i], " ", keepempty=false)
            for q in eachindex(header)
                append!(data[header[q]], [parse(Float64, l[q])])
            end
        catch
            break
        end
        i=i+1
    end

    data
end

read_marcs(path)