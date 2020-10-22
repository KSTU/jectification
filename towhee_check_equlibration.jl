
function towhee_equlibrate(equlibrationProp, equlibrationCheck, equlibrationId, curFlow, equlibrationNumber)
#function towhee_equlibrate()    #debug editionion
    cd("$(cfg.dir)/flow$(curFlow)")
    curId = equlibrationId
    #set UP propertyes
    if equlibrationId > equlibrationNumber
        for i = 1:equlibrationNumber-1
            #println("old " , equlibrationProp[i], " new ", equlibrationProp[i+1])
            #equlibrationProp[i] = equlibrationProp[i+1]
            #copyto!(equlibrationProp[i+1], equlibrationProp[i])
        end
        curId = equlibrationId % equlibrationNumber
        if curId == 0
            curId = 4
        end
    end
    
    towhee_get_out_prop("equlibration", equlibrationProp, curId, false) #3 blocks
    
    #check equlibration
    if equlibrationId > equlibrationNumber
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
        errEn = errEn / abs(meanEn)
        errDens = errDens / abs(meanDens)
        errPres = errPres / abs(meanPres)
        
        if errEn < 0.05 && errPres < 0.3 && errDens < 0.05
            equlibrationCheck[curFlow] = 1
        end
    end
    
    #plot equlibrationons
    println("current flow $(curFlow)")
    println("energy   \t pressure \t density")
    
    if equlibrationId > equlibrationNumber
        for i = 1:equlibrationNumber
            println("$(equlibrationProp[i].totEnergy) \t $(equlibrationProp[i].pressure) \t $(equlibrationProp[i].density)")
        end
        println("_____________________________________")
        println("$(errEn) \t $(errPres) \t $(errDens)")
    else
        for i = 1:curId
            println("$(equlibrationProp[i].totEnergy) \t $(equlibrationProp[i].pressure) \t $(equlibrationProp[i].density)")
        end
    end
    
end

function read_input_prop(prop, curFlow)
    cd("$(cfg.dir)/flow$(curFlow)")
    fileId = open("production","r")
    s = read(fileId, String)
    if length(split(s,"+++++ end of markov chain +++++")) > 1
        s = split(s,"+++++ end of markov chain +++++")[2]
        ars = split(s,"\n")
        for curString in ars
            if (findfirst("Molecule Number", curString) === nothing) == false
            arrayId = parse(Int64, split(curString)[3])
            prop[curFlow].NType[arrayId] = ceil(parse(Float64, split(curString)[4]))
        end
        end
        s = split(s,"Block Averages")[2]
    else
        println("towhee production error")
        exit()
    end
    ars = split(s,"\n")
    #parse new propertyes
    for curString in ars
        if (findfirst("Total Classical", curString) === nothing) == false
            prop[curFlow].totEnergy = parse(Float64, split(curString)[5])
        end
        if (findfirst("Specific Density", curString) === nothing) == false
            prop[curFlow].density = parse(Float64, split(curString)[5])
        end
        if (findfirst("Virial Pressure", curString) === nothing) == false
            prop[curFlow].pressure = parse(Float64, split(curString)[5])
        end
        if (findfirst("Volume", curString) === nothing) == false
            prop[curFlow].pressure = parse(Float64, split(curString)[5])
        end
        if (findfirst("Number Density", curString) === nothing) == false
            #println("cur string ", curString)
            curId =  parse(Int64, split(curString)[4])
            prop[curFlow].NDens[curId] = parse(Float64, split(curString)[6])
            println("t ", curId, " val ", prop[curFlow].NDens[curId])
        end
        
        #read molecules number
    end
    
end

function towhee_equlibrate_prop(eqProp, curProp, state, nPlate, check)
    #set propertyes
    if check == "energy"
        eqProp[state.curEqBlock[nPlate]]. 
    end

    if state.curEqBlock[nPlate] < state.eqNumber
        state.curEqBlock[nPlate] += 1
    end
end
