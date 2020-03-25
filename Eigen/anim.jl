# this is a script to animate hydrogen orbitals

function main()
    file = ARGS[1]
    io = open(file, "r")
    data = read(io, String)
    while true
        println(data)
    end
end

main()

