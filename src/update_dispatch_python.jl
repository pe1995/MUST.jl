file = joinpath(ARGS[1], "utilities/python/dispatch/EOS/stagger.py")
f = open(file, "r") 
lines = readlines(f; keep=true)
close(f)

if lines[1] != "from copy import deepcopy\n"
    insert!(lines, 1, "from copy import deepcopy\n")
    insert!(lines, 2, "from ._EOS import interp\n")
end

g = open(file, "w") 
for line in lines
    write(g, line)
end
close(g)