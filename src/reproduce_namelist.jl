
function reverse_parse(value)
    val_str = ""
    if typeof(value) <: AbstractArray
        for val in value
            val_str = val_str * "$(reverse_parse(val)),"
        end
        val_str = val_str[1:end-1]
    elseif typeof(value) <:String
        if (lowercase(value) in ["t", "f", ".true.", ".false."]) | (occursin('*', value))
            v = strip(value)
            v = (v[1] == '\'') ? v[2:end] : v
            v = (v[end] == '\'') ? v[1:end-1] : v
            val_str = "$(strip(v))"
        elseif (occursin("'", value)) | (occursin('"', value))
            v = strip(value)
            v = (v[1] == '\'') ? v[2:end] : v
            v = (v[end] == '\'') ? v[1:end-1] : v
            val_str = "$(strip(v))"
        else
            can_be_float = try 
                parse(Float64, value)
                true
            catch
                false
            end

            if can_be_float
                val_str = "$(value)"
            else
                val_str = "'$(value)'"
            end
        end
    elseif typeof(value) <:Bool
        val_str = ".$(value)."
    else
        val_str = "$(value)"
    end
    return val_str
end

println("word -> ", reverse_parse("word"))
println("1.0 -> ", reverse_parse("1.0"))
println("'word' -> ", reverse_parse("'word'"))
println("O'Neil -> ", reverse_parse("O'Neil"))
println("quoted 'word' internal -> ", reverse_parse("this is 'quoted'."))
println("double quoted \"word\" -> ", reverse_parse("double \"quoted\""))
println("Bool true -> ", reverse_parse(true))
println("String .true. -> ", reverse_parse(".true."))
