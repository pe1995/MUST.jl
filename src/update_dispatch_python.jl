f = open(ARGS[1], "r") 
lines = readlines(f; keep=true)
close(f)

if lines[1] != "from copy import deepcopy\n"
    insert!(lines, 1, "from copy import deepcopy\n")
    insert!(lines, 2, "from ._EOS import interp\n")
end

g = open(ARGS[1], "w") 
for line in lines
    write(g, line)
end
close(g)