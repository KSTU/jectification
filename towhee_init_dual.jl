
function set_initial_state(state)
    if state == "VAK"
        
    end
end

function towhee_coords_dual(onemol, cfg, plateProp, ndx)

    #liquid phase
    fileId = open("packmol.inp", "w")
    print(fileId, "#\n#\n#\n\n")
    println(fileId, "tolerance 3.0")
    println(fileId, "filetype pdb")
    println(fileId, "output towhee_coords.pdb")
    h = plateProp[ndx].liqVol ^ (1.0/3.0)
    for i = 1:cfg.subNum
        println(fileId, "structure ../$(cfg.subName[i])")
        println(fileId, " number $(plateProp[ndx].liqNType[i])")
        println(fileId, "inside box 0. 0. 0.  $(h)  $(h)  $(h)")
        println(fileId, "end structure\n")
    end
    
    close(fileId)
    
    run(pipeline(`packmol`, stdin="packmol.inp"))
    pdb_2_xyz("towhee_coords.pdb", "towhee_coords_liquid")
    #vapor phase
    fileId = open("packmol.inp", "w")
    print(fileId, "#\n#\n#\n\n")
    println(fileId, "tolerance 3.0")
    println(fileId, "filetype pdb")
    println(fileId, "output towhee_coords.pdb")
    h = plateProp[ndx].vapVol ^ (1.0/3.0)
    for i = 1:cfg.subNum
        println(fileId, "structure ../$(cfg.subName[i])")
        println(fileId, " number $(plateProp[ndx].vapNType[i])")
        println(fileId, "inside box 0. 0. 0.  $(h)  $(h)  $(h)")
        println(fileId, "end structure\n")
    end
     
    run(pipeline(`packmol`, stdin="packmol.inp"))
    pdb_2_xyz("towhee_coords.pdb", "towhee_coords_vapor")
    
    
    close(fileId)
    return 0
end

function towhee_input_dual()
    
    return 0
end
