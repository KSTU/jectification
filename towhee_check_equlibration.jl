
function towhee_equlibrate(equlibrationCheck, equlibrationId)
#function towhee_equlibrate()    #debug editionion
    fileId = open("equlibration","r")
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) >1
        s = split(s,"+++++ end of markov chain +++++")[2]
    else
        println("towhee equlibration error")
        exit()
    end
    ars = split(s,"\n")
    for curString in ars
        if (findfirst("Total Classical", curString) === nothing) == false
            energy = parse(Float64, split(curString)[4])
            println(energy, typeof(energy))
        end
        if (findfirst("Specific Density", curString) === nothing) == false
            density = parse(Float64, split(curString)[4])
            println(density, typeof(density))
        end
        if (findfirst("Virial Pressure", curString) === nothing) == false
            pressure = parse(Float64, split(curString)[4])
            println(pressure, typeof(pressure))
        end
    end
#    ars = split(s,"\n")
#    println(typeof(ars),ars)
    close(fileId)
end

a = 1
b = 1
towhee_equlibrate(a,b)
