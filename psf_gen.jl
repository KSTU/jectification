


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
        
    end
    println(nbond, box.molType)
    println( sum(nbond .* box.molType))
    
    #print bonds 
    print(fileId, "\n")
    print(fileId, sum(nbond .* box.molType) ," !NBOND: bonds\n")
    id = 0
    curZero = 0 #zero for current molecule
    for i = 1:cfg.subNum
        for j = 1:box.molType[i]
            for k = 1:nbond[i]
                id += 1
                print(fileId, "    ", bondAtomA[i][k]+curZero, "     " , bondAtomB[i][k] + curZero )
                if id % 3 == 0
                    print(fileId, "\n")
                end
            end
            curZero += onemol[i].atomNum
        end
    end
    
    #angles
    
    #dihedral
    
    
    
    close(fileId)
end

