struct MARCSModel
    k
    info
    abundances
    structure
    pressures
end

MARCSModel(path::String) = begin
    info, abundances, structure, pressure = open(path, "r") do f
        title = readline(f)
        info = []
        line = ""
        while !occursin("chemical number abundances", line)
            line = readline(f)
            append!(info, [line])
        end
        info = info[1:end-1]

        abundances = []
        line = readline(f)
        while !occursin("Number of depth points", line)
            try
                append!(abundances, split(line, " "))
            catch
                nothing
            end

            line = readline(f)
        end
        abundances = [a for a in abundances if a!=""]

        structure = read_marcs_sub_table(
            f, 
            start_at = "Model structure", 
            end_at = "Assorted logarithmic partial pressures"
        )

        pressures = read_marcs_sub_table(f)

        info, abundances, structure, pressures
    end

    MARCSModel(length(structure["k"]), info, abundances, structure, pressure)
end

function read_marcs_sub_table(f; start_at=nothing, end_at=nothing)
    data = Dict()
    if !isnothing(start_at)
        header = [""]
        while !occursin(start_at, last(header))
            append!(header, [readline(f)])
        end
    end

    use_end = if isnothing(end_at)
        end_at = "dummy end at"
        use_end=false
    else
        use_end=true
    end

    line = readline(f)
    is_new_header = false
    current_header = []
    while !occursin(end_at, line)
        # read until you find a new header
        line = try
            is_new_header = false
            l = parse.(Float64, split(line, " ", keepempty=false))
            l
        catch
            is_new_header = true
            l_raw = line
            l = split(l_raw, " ", keepempty=false)

            if occursin("H I", l_raw)
                for (i, v) in enumerate(l)
                    if v=="I"
                        l[i-1] = "HI"
                        deleteat!(l, i)
                        break
                    end
                end
            end

            l
        end

        if is_new_header
            current_header = line
            # this is a new header, add or remove the existing columns in data
            for (i, quantity) in enumerate(current_header)
                data[quantity] = Float64[]
            end
        else
            # read otherwise
            for (i, q) in enumerate(current_header)
                append!(data[q], line[i])
            end
        end

        line = readline(f)

        if !use_end & (line=="")
            break
        end
    end

    data
end

function marcsBox(path::String)
    model = MARCSModel(path)
    xx = zeros(1, 1, model.k)
    yy = zeros(1, 1, model.k)
    zz = reshape(model.structure["Depth"], 1, 1, :)

    d = reshape(model.structure["Density"], 1, 1, :) 
    T = reshape(model.structure["T"], 1, 1, :)    
    t = exp10.(reshape(model.structure["lgTauR"], 1, 1, :))    
    t5 = exp10.(reshape(model.structure["lgTau5"], 1, 1, :))    
    pr = reshape(model.structure["Prad"], 1, 1, :)
    pg = reshape(model.structure["Pg"], 1, 1, :)    

    ne = reshape(model.structure["Pe"], 1, 1, :) ./ (KBoltzmann .* T)

    data = Dict{Symbol, Any}(
        :Ne=>ne, :T=>T, :d=>d, :τ_ross=>t, :τ500=>t5, :Pg=>pg, :Pr=>pr
    )

    for (k, v) in model.pressures
        data[Symbol(k)] = reshape(v, 1, 1, :)
    end

    Box(xx, yy, zz, data, AtmosphericParameters())
end

"""
	marcs_parameters_from_name(marcs_path)

Extract stellar parameters from the name of the model.
"""
marcs_parameters_from_name(marcs_path) = begin
	name = first(split(marcs_path, ".mod"))
	name = split(name, "_", keepempty=false)
	T = parse(Float64, name[1][2:end])
	g = parse(Float64, name[2][2:end])
	feh = parse(Float64, name[6][2:end])
	turb = parse(Float64, name[4][2:end])
	alpha = parse(Float64, name[7][2:end])
	sorp = name[1][1]
	
    ModelInformation(
        teff=T,
        logg=g,
        feh=feh,
        extension="$(sorp),alpha=$(alpha),turb=$(turb)",
        original_name=marcs_path,
        category="MARCS"
    )
end