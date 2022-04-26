#### Distributed execution functions #####

function sendsync(p::Int; args...)
    @sync for (nm, val) in args
        @async @spawnat(p, Core.eval(Main, Expr(:(=), nm, val)))
    end
end

function sendsync(p::Vector{Int}; args...)
    for ip in p
        sendsync(ip; args...)
    end
end

function wait(f::Vector{Distributed.Future})
    for fut in f 
        fetch(fut)
    end

    nothing
end