


function psf_generate(cfg, box, onemol, fileName)
    fileId = open(fileName, "w")
    print(fileId, "PSF\n\n")
    print(fileId, "4 !NTITLE")
    print(fileId, " REMARKS\n REMARKS\n REMARKS\n REMARKS\n\n")
    #atoms
    print(fileId, adj_right(box.atomNum,8),"!NATOM\n")
    for i = 1:box.atomNum
        print(fileId, adj_left(box.curAtom[i],8)," ")
        print(fileId, adj_left(" M1", 3), )
        print(fileId, adj_left(box.curMol[i], 6))
        print(fileId, adj_left(box.molName[i], 6))
        print(fileId, adj_left(box.atomName[i], 6))
        print(fileId, adj_left(box.nsFF[i], 6))
        print(fileId, adj_left(box.nsQ[i], 8, 3))
        print(fileId, "   ")
        print(fileId, adj_left(box.nsMass[i], 8, 3))
#        print(fileId, adj_left(0, 6))   ##?????
        print(fileId, "     0") #?????
        print(fileId,"\n")
    end
    
    #calculate bonds, angles and torsion
    nbond = zeros(Int64, cfg.subNum)    #numbers of bonds
    bondAtomA = Array{AbstractArray}(undef, cfg.subNum)
    bondAtomB = Array{AbstractArray}(undef, cfg.subNum)
    nangles = zeros(Int64, cfg.subNum)  #numbers of angles
    angleAtomA = Array{AbstractArray}(undef, cfg.subNum)
    angleAtomB = Array{AbstractArray}(undef, cfg.subNum)
    angleAtomC = Array{AbstractArray}(undef, cfg.subNum)
    ntorsion = zeros(Int64, cfg.subNum) #numbers of torsion angle
    torsionAtomA = Array{AbstractArray}(undef, cfg.subNum)
    torsionAtomB = Array{AbstractArray}(undef, cfg.subNum)
    torsionAtomC = Array{AbstractArray}(undef, cfg.subNum)
    torsionAtomD = Array{AbstractArray}(undef, cfg.subNum)
    for i = 1:cfg.subNum
        #set zeroes
        connectarray = zeros(Int64, onemol[i].atomNum, onemol[i].atomNum)
        bondAtomA[i] = zeros(Int64, onemol[i].atomNum*onemol[i].atomNum) # set arrays to zero
        bondAtomB[i] = zeros(Int64, onemol[i].atomNum*onemol[i].atomNum)
        for atoma = 1:onemol[i].atomNum
            for atomb =1:size(onemol[i].connect[atoma],1)
                #println("connettion ", onemol[i].connect[atoma][atomb], typeof(onemol[i].connect[atoma][atomb]))
                connectarray[atoma,onemol[i].connect[atoma][atomb]] = 1
            end
        end
        #println("molecule", i, " connectarray ", connectarray)
        #get bond numbers
        id = 0
        for atoma = 1:onemol[i].atomNum
            if(onemol[i].atomNum > 1)
                for atomb = atoma+1:onemol[i].atomNum
                    if connectarray[atoma, atomb] == 1
                        id += 1
                        nbond[i] += 1
                        bondAtomA[i][id] = atoma
                        bondAtomB[i][id] = atomb
                    end
                end
            end
        end
        println("numbers of bond ", nbond[i], bondAtomA[i] ,bondAtomB[i] )
        #get angles number
        for atoma = 1:onemol[i].atomNum
            for atomb = 1:onemol[i].atomNum
                for atomc = 1:onemol[i].atomNum
                    if connectarray[atoma, atomb] == 1 && connectarray[atomb, atomc] == 1  && atomc > atoma
                        println("aa ", atoma, " ab ", atomb, " ac ", atomc)
                        nangles[i] += 1
                    end
                end
            end
        end
        println("numbers af angles ", nangles[i])
        angleAtomA[i] = zeros(Int64, nangles[i])
        angleAtomB[i] = zeros(Int64, nangles[i])
        angleAtomC[i] = zeros(Int64, nangles[i])
        id = 0
        for atoma = 1:onemol[i].atomNum
            for atomb = 1:onemol[i].atomNum
                for atomc = 1:onemol[i].atomNum
                    if connectarray[atoma,atomb] == 1 && connectarray[atomb,atomc] == 1  && atomc > atoma
                        id += 1
                        angleAtomA[i][id] = atoma
                        angleAtomB[i][id] = atomb
                        angleAtomC[i][id] = atomc
                    end
                end
            end
        end
        # calculate numbers of torsion angle
        for atoma = 1:onemol[i].atomNum
            for atomb = 1:onemol[i].atomNum
                for atomc = 1:onemol[i].atomNum
                    for atomd = 1:onemol[i].atomNum
                        if connectarray[atoma,atomb] ==1 && connectarray[atomb,atomc] == 1 && connectarray[atomc,atomd] == 1 && atomd > atoma && atoma != atomc && atoma != atomd && atomb != atomc && atomb != atomd
                            ntorsion[i] +=1
                        end
                    end
                end
            end
        end
        torsionAtomA[i] = zeros(Int64, ntorsion[i])
        torsionAtomB[i] = zeros(Int64, ntorsion[i])
        torsionAtomC[i] = zeros(Int64, ntorsion[i])
        torsionAtomD[i] = zeros(Int64, ntorsion[i])
        id = 0
        for atoma = 1:onemol[i].atomNum
            for atomb = 1:onemol[i].atomNum
                for atomc = 1:onemol[i].atomNum
                    for atomd = 1:onemol[i].atomNum
                        if connectarray[atoma,atomb] ==1 && connectarray[atomb,atomc] == 1 && connectarray[atomc,atomd] == 1 && atomd > atoma && atoma != atomc && atoma != atomd && atomb != atomc && atomb != atomd
                            id +=1
                            torsionAtomA[i][id] = atoma
                            torsionAtomB[i][id] = atomb
                            torsionAtomC[i][id] = atomc
                            torsionAtomD[i][id] = atomd
                        end
                    end
                end
            end
        end
    end
    println(nbond, box.molType)
    println( sum(nbond .* box.molType))
    #println("angle arrays ", angleAtomA, angleAtomB, angleAtomC)
    for i = 1:cfg.subNum
        for j = 1:nangles[i]
            println("atoma ", angleAtomA[i][j], " atomb ", angleAtomB[i][j], " atomc ", angleAtomC[i][j])
        end
    end
    for i = 1:cfg.subNum
        println("substance ", i, "torsion angles ", ntorsion[i])
        for j = 1:ntorsion[i]
            println("atoma ", torsionAtomA[i][j], " atomb ", torsionAtomB[i][j], " atomc ", torsionAtomC[i][j], "atomd", torsionAtomD[i][j])
        end
    end
    
    
    #print bonds 
    if maximum(nbond) > 0
        print(fileId, "\n")
        print(fileId, adj_right(sum(nbond .* box.molType), 10) ," !NBOND: bonds\n")
        id = 0
        curZero = 0 #zero for current molecule
        for i = 1:cfg.subNum
            for j = 1:box.molType[i]
                for k = 1:nbond[i]
                    id += 1
                    print(fileId, "   ", adj_right(bondAtomA[i][k]+curZero, 8 ), "   " , adj_right(bondAtomB[i][k] + curZero,8) )
                    if id % 3 == 0
                        print(fileId, "\n")
                    end
                end
                curZero += onemol[i].atomNum
            end
        end
    end
    
    #print angles
    if maximum(nangles) > 0
        print(fileId, "\n")
        print(fileId, adj_right(sum(nangles .* box.molType), 10), " !NTHETA: angles \n")
        id = 0
        curZero = 0
        for i = 1:cfg.subNum
            for j = 1:box.molType[i]
                for k = 1:nangles[i]
                    id += 1
                    print(fileId, "   ", adj_right(angleAtomA[i][k] + curZero, 8), "   ", adj_right(angleAtomB[i][k] + curZero, 8), "   ", adj_right(angleAtomC[i][k] + curZero, 10))
                    if id % 3 == 0
                        print(fileId, "\n")
                    end
                end
                curZero += onemol[i].atomNum
            end
        end
    end
    
    #dihedral
    if maximum(ntorsion) > 0
        print(fileId, "\n")
        print(fileId, adj_right(sum(ntorsion .* box.molType), 10), " !NPHI: dihedrals \n")
        id = 0
        curZero = 0
        for i = 1:cfg.subNum
            for j = 1:box.molType[i]
                for k = 1:ntorsion[i]
                    id += 1
                    print(fileId, "   ", adj_right(torsionAtomA[i][k] + curZero, 8), "   ", adj_right(torsionAtomB[i][k] + curZero, 8), "   ", adj_right(torsionAtomC[i][k] + curZero, 10), "   ", adj_right(torsionAtomD[i][k] + curZero, 8))
                    if id % 2 == 0
                        print(fileId, "\n")
                    end
                end
                curZero += onemol[i].atomNum
            end
        end
    end
    
    
    
    
    close(fileId)
end

