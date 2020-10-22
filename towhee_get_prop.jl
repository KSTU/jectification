
function towhee_get_out_prop(filename, prop, curId, oneblock)
    fileId = open(filename,"r")
    
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) > 1
        s = split(s,"+++++ end of markov chain +++++")[2]
        #get molecule numbers
        ars = split(s,"\n")
        for curString in ars
            
            if (findfirst("Molecule Number", curString) === nothing) == false
                #println(curString)
                #println(split(curString))
                arrayId = parse(Int64, split(curString)[3])
                prop[curId].NType[arrayId] = round(Int64, parse(Float64, split(curString)[4]))
            end
        end
        if !oneblock
            s = split(s,"Block Averages")[2]
        end
    else
        println("towhee equlibration error")
        exit()
    end
    s = split(s, "Averages               Units")[2]
    ars = split(s,"\n")

    if oneblock #if one block without averages
        for curString in ars
            if length(curString) > 15
                if cmp(curString[1:16], " Total Classical") == 0
                    prop[curId].totEnergy = parse(Float64, split(curString)[4])
                end
                if cmp(curString[1:16], " Specific Densit") == 0
                    prop[curId].density = parse(Float64, split(curString)[4])
                end
                if cmp(curString[1:16], " Virial Pressure") == 0
                    prop[curId].pressure = parse(Float64, split(curString)[4])
                end
                if cmp(curString[1:16], " Volume         ") == 0
                    prop[curId].volume = parse(Float64, split(curString)[3])
                end
            end
        end
    else
        for curString in ars
            if length(curString) > 15
                if cmp(curString[1:16], " Total Classical") == 0
                    prop[curId].totEnergy = parse(Float64, split(curString)[5])
                end
                if cmp(curString[1:16], " Specific Densit") == 0
                    prop[curId].density = parse(Float64, split(curString)[5])
                end
                if cmp(curString[1:16], " Virial Pressure") == 0
                    prop[curId].pressure = parse(Float64, split(curString)[5])
                end
                if cmp(curString[1:16], " Volume         ") == 0
                    prop[curId].pressure = parse(Float64, split(curString)[4])
                end
            end
        end
    end
    
    close(fileId)
end
