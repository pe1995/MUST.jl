function reverse_parse(value)
    val_str = ""
    if typeof(value) <: AbstractArray
        for val in value
            val_str = val_str * "$(reverse_parse(val)),"
        end
        val_str = val_str[1:end-1]
    elseif typeof(value) <:Bool
        val_str = (value) ? ".true." : ".false."
    elseif typeof(value) <: AbstractString
        v = strip(value)
        if lowercase(v) in ["t", "f", ".true.", ".false."]
            val_str = v
        elseif tryparse(Float64, v) !== nothing
            val_str = v
        else
            if occursin("'", v)
                val_str = "\"$(v)\""
            elseif occursin("\"", v)
                val_str = "'$(v)'"
            else
                val_str = "$(v)"
            end
        end
    else
        val_str = "$(value)"
    end
    return val_str
end

println("word -> ", reverse_parse("word"))
println("1.0 -> ", reverse_parse("1.0"))
println("'word' -> ", reverse_parse("'word'"))
println("O'Neil -> ", reverse_parse("O'Neil"))
println("quoted 'word' -> ", reverse_parse("this is 'quoted'."))
println("double \"quoted\" -> ", reverse_parse("double \"quoted\""))
