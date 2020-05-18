


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
        print(fileId,"\n")
    end
    
#    for i = 1:cfg.subNum
#        for j = 1:onemol[i].atomNum
#            print(fileId, adj_left(,6))
#        end
#    end
    #bonds
    
    #angles
    
    #dihedral
    
    
    
    close(fileId)
end

