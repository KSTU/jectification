
function towhee_equlibrate(equlibrationProp, equlibrationCheck, equlibrationId, curFlow, equlibrationNumber)
#function towhee_equlibrate()    #debug editionion
    cd("$(cfg.dir)/flow$(curFlow)")
    fileId = open("equlibration","r")
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) >1
        s = split(s,"+++++ end of markov chain +++++")[2]
    else
        println("towhee equlibration error")
        exit()
    end
    ars = split(s,"\n")
    #set UP propertyes
    if equlibrationId == equlibrationNumber
        for i = 1:equlibrationNumber-1
            equlibrationProp[i] = equlibrationProp[i+1]
        end
    end
    
    #parse new propertyes
    for curString in ars
        if (findfirst("Total Classical", curString) === nothing) == false
            equlibrationProp[equlibrationId].totEnergy = parse(Float64, split(curString)[4])
        end
        if (findfirst("Specific Density", curString) === nothing) == false
            equlibrationProp[equlibrationId].density = parse(Float64, split(curString)[4])
        end
        if (findfirst("Virial Pressure", curString) === nothing) == false
            equlibrationProp[equlibrationId].pressure = parse(Float64, split(curString)[4])
        end
    end
    
    #check equlibration
    if equlibrationId == equlibrationNumber
        meanEn = 0.0
        meanDens = 0.0
        meanPres = 0.0
        for i = 1:equlibrationNumber #get mean
            meanEn += equlibrationProp[i].totEnergy
            meanDens += equlibrationProp[i].density
            meanPres += equlibrationProp[i].pressure
        end
        meanEn /= equlibrationNumber
        meanDens /= equlibrationNumber
        meanPres /= equlibrationNumber
        
        errEn = 0.0
        errDens = 0.0
        errPres = 0.0
        for i = 1:equlibrationNumber
            errEn += abs(equlibrationProp[i].totEnergy - meanEn)
            errDens += abs(equlibrationProp[i].density - meanDens)
            errPres += abs(equlibrationProp[i].pressure - meanPres)
        end
        errEn = errEn / meanEn
        errDens = errDens / meanDens
        errPres = errPres / meanPres
        
        if errEn < 0.1 && errPres < 0.1 && errDens < 0.1
            equlibrationCheck[curFlow] = 1
        end
    end
    
    #plot equlibrationons
    println("current flow $(curFlow)")
    println("energy   \t    pressure \t density")
    for i = 1:equlibrationId
        println("$(equlibrationProp[i].totEnergy) \t $(equlibrationProp[i].pressure) \t $(equlibrationProp[i].density)")
    end
    if equlibrationId == equlibrationNumber
        println("_____________________________________")
        println("$(errEn) \t $(errPres) \t $(errDens)")
    end
    
    close(fileId)
end

function read_input_prop(prop, curFlow)
    cd("$(cfg.dir)/flow$(curFlow)")
    fileId = open("production","r")
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) >1
        s = split(s,"+++++ end of markov chain +++++")[2]
    else
        println("towhee production error")
        exit()
    end
    ars = split(s,"\n")
    #parse new propertyes
    for curString in ars
        if (findfirst("Total Classical", curString) === nothing) == false
            prop[curFlow].totEnergy = parse(Float64, split(curString)[4])
        end
        if (findfirst("Specific Density", curString) === nothing) == false
            prop[curFlow].density = parse(Float64, split(curString)[4])
        end
        if (findfirst("Virial Pressure", curString) === nothing) == false
            prop[curFlow].pressure = parse(Float64, split(curString)[4])
        end
    end
    
end
