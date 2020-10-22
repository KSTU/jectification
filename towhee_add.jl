
function towhee_add(cfg, prop, inputProp, cplate)
    for cflow = 1:cfg.inputNumber
        if cplate == cfg.inputPlates[cflow]
            #println("add debug ", cfg.inputPlates[cflow]," " , cfg.inputMolecules[cflow,:], " ", prop[cfg.inputPlates[cflow]].liqN, " -- ", sum(cfg.inputMolecules[cflow,:]))
            prop[cplate].liqN += sum(cfg.inputMolecules[cflow,:])   #add all molecules
            for cmol = 1:cfg.subNum #add molecules bytype
                prop[cplate].liqNType[cmol] += cfg.inputMolecules[cflow, cmol]
            end
            #println("energy ", inputProp[cflow].totEnergy," sum ", sum(inputProp[cflow].NType))
            prop[cplate].refEnergy += sum(cfg.inputInitMolecules[cflow,:]) * inputProp[cflow].totEnergy / sum(inputProp[cflow].NType)
            prop[cplate].inDens[1] = sum(inputProp[cflow].NDens)
            prop[cplate].inNumber[1] = sum(cfg.inputMolecules[cflow,:])
            prop[cplate].inTemp[1] = cfg.inputTemperature[cflow]
            prop[cplate].inPres[1] = inputProp[cflow].pressure
        else
            prop[cplate].inDens[1] = 0.0
            prop[cplate].inNumber[1] = 0
            prop[cplate].inTemp[1] = 0.0
            prop[cplate].inPres[1]= 0.0
        end
    end
end

function towhee_box_prepare(cfg, prop, inputProp, cplate)
    
    prop[cplate].liqVol = 0.0   #zeroes  volumeinvalid
    prop[cplate].temp = 0.0     #zeroes temperature
    prop[cplate].pressure = 0.0 #zeroes pressure
    for i = 1:3
        if prop[cplate].inNumber[i] > 0 && prop[cplate].inDens[i] > 0
            prop[cplate].liqVol += prop[cplate].inNumber[i] / prop[cplate].inDens[i] * 1000 # nm to angstrems
            prop[cplate].temp += prop[cplate].inTemp[i] * prop[cplate].inNumber[i]
            prop[cplate].pressure += prop[cplate].inPres[i] * prop[cplate].inNumber[i]
        end
    end #add input volumes  [1] - input [2] - up  [3] - down
    prop[cplate].temp /= sum(prop[cplate].inNumber)
    prop[cplate].pressure /= sum(prop[cplate].inNumber)
    prop[cplate].vapVol = cfg.plateVolume[cplate] - prop[cplate].liqVol
    println("vapor volume ", prop[cplate].liqVol)
    
end


function towhee_box_generate(cfg, prop, inputProp, cplate)
    
    #generate vapor
    if prop[cplate].vapN  > 0
    fileId = open("vapor.inp", "w")
    print(fileId, "tolerance 3.0\n")
    print(fileId, "output vapor.pdb\n")
    print(fileId, "filetype pdb\n")
    vBox = prop[cplate].vapVol ^ (1.0/3.0)
    for cmol = 1:cfg.subNum
        if prop[cplate].vapNType[cmol] > 0
            print(fileId, "structure ../$(cfg.subName[cmol])\n")
            print(fileId, "number $(prop[cplate].vapNType[cmol])")
            print(fileId, "inside cube 0. 0. 0. $(vBox)\n")
            print(fileId, "end structure\n\n")
        end
    end
    close(fileId)
    run(pipeline(`packmol`, stdin = "vapor.inp"))
    end
    #generate liquid
    if prop[cplate].liqN > 0
    fileId = open("liquid.inp", "w")
    print(fileId, "tolerance 3.0\n")
    print(fileId, "output liquid.pdb\n")
    print(fileId, "filetype pdb\n")
    lBox = prop[cplate].liqVol ^ (1.0/3.0)
    for cmol = 1:cfg.subNum
        if prop[cplate].liqNType[cmol] > 0
            print(fileId, "structure ../$(cfg.subName[cmol])\n")
            print(fileId, "number $(prop[cplate].liqNType[cmol])\n")
            print(fileId, "inside cube 0. 0. 0. $(lBox)\n")
            print(fileId, "end structure\n\n")
        end
    end
    close(fileId)
    run(pipeline(`packmol`, stdin = "liquid.inp"))
    end

end

function towhee_pdb_to_cords(inputFile)
    fileId2 = open("towhee_coords", "w")
    for curFile = 1:length(inputFile)
        fileId = open(inputFile[curFile])
        for i in eachline(fileId)
            if length(i) > 6
                if i[1:6] == "HETATM" || i[1:6] == "ATOM  "
                    println(fileId2, i[31:54])
                end
            end
        end
        close(fileId)
    end
    close(fileId2)
    
end
