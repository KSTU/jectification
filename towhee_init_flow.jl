

function generate_towhee_flow(onemol, cfg)
    mainFolder = pwd()
    for ndx = 1:size(cfg.inputPlates,1)
        #if NVT ensable
        cd("flow$(ndx)")
        towhee_initial(onemol, cfg, ndx, 1, true)
        towhee_coords_initial(onemol, cfg, cfg.inputMolecules, ndx)
        cd("..")
    end
end

function towhee_initial(onemol, cfg, ndx, nloop, initial)
    fileId = open("towhee_input", "w")
    println(fileId, "ensemble")
    if cfg.inputEnsamble[ndx] == "NVT"
        println(fileId, "'", "nvt", "'")
        println(fileId, "temperature")
        println(fileId, cfg.inputTemperature[ndx])
    elseif cfg.inputEnsamble[ndx] == "NPT"
        println(fileId, "'", "npt", "'")
        println(fileId, "temperature")
        println(fileId, cfg.inputTemperature[ndx])
        println(fileId, "pressure")
        println(fileId, cfg.inputPressure[ndx] * 100)    #to kPa
    end
    println(fileId, "nmolty")
    println(fileId, cfg.subNum)
    println(fileId, "nmolectyp")
    scale = ceil(Int64, cfg.inputInitMolecules[ndx] / sum(cfg.inputMolecules[ndx,:]))
    #println( sum(cfg.inputMolecules[ndx,:]))
    #println(scale)
    for i = 1:cfg.subNum
        print(fileId, scale * cfg.inputMolecules[ndx, i], "  ")
    end
    print(fileId, "\n")
    if cfg.inputEnsamble[ndx] == "NVT" || cfg.inputEnsamble[ndx] == "NPT"
        println(fileId, "numboxes")
        println(fileId, "1")
    end
    println(fileId, "stepstyle")
    println(fileId, "'moves'")
    println(fileId, "nstep")
    println(fileId, sum(cfg.inputMolecules[ndx,:]) * scale * 2 * nloop)
    println(fileId, "printfreq")
    println(fileId, 1000)
    println(fileId, "blocksize")
    println(fileId, sum(cfg.inputMolecules[ndx,:]) * scale * 2 * nloop)
    println(fileId, "moviefreq")
    println(fileId, 1000)
    println(fileId, "backupfreq")
    println(fileId, 5000)
    println(fileId, "runoutput")
    println(fileId, "full")
    println(fileId, "pdb_output_freq")
    println(fileId, 0)
    println(fileId, "pressure_virial_freq")
    println(fileId, 10)
    println(fileId, "trmaxdispfreq")
    println(fileId, 1000)
    println(fileId, "volmaxdispfreq")
    println(fileId, 1000)
#    println(fileId, "chempotperstep")
#    println(fileId, 0)
    println(fileId, "linit")
    if initial == true
        println(fileId, ".true.")
    else
        println(fileId, ".false.")
    end
    
    if cfg.inputEnsamble[ndx] == "NVT" || cfg.inputEnsamble[ndx] == "NPT"
        println(fileId, "initboxtype")
        println(fileId, "'dimensions'")
        println(fileId, "initstyle")
        for i = 1:cfg.subNum
            print(fileId, "'full cbmc'", "    ")
        end
        print(fileId, "\n")
        println(fileId, "initlattice")
        for i = 1:cfg.subNum
            print(fileId, "'simple cubic'", "    ")
        end
        print(fileId, "\n")
        println(fileId, "initmol")
        for i = 1:cfg.subNum
            print(fileId, scale * cfg.inputMolecules[ndx,i], " ")
        end
        print(fileId, "\n")
        println(fileId, "inix iniy iniz")
        #calculate dimentions
        dim = ceil(Int64, (scale * sum(cfg.inputMolecules[ndx,:]))^ (1/3))
        println(fileId, adj_left(dim, 6), adj_left(dim, 6), adj_left(dim, 6) )
        #calculate system size
        h = (scale * sum(cfg.inputMolecules[ndx,:]) / (cfg.inputDensity[ndx] * 6.02214076 /10^4))^(1/3)
        println(fileId, "hmatrix")
        println(fileId, h, "   0   ", "    0    ")
        println(fileId, "   0   ", h,  "    0    ")
        println(fileId, "   0   ", "    0    ", h)
        
    end
    if cfg.inputEnsamble[ndx] == "NPT"
        println(fileId, "pmvol     ")
        println(fileId, "0.01d0   ")
        println(fileId, "          pmvlpr")
        println(fileId, "          1.0d0")
        println(fileId, "          rmvol")
        println(fileId, "          0.1d0")
        println(fileId, "          tavol")
        println(fileId, "          0.5d0")
    end
    if cfg.inputEnsamble[ndx] == "NVT" || cfg.inputEnsamble[ndx] == "NPT"
        if cfg.subDim == "1D"
            println(fileId, "pmtracm   ")
            println(fileId, "1.0d0")
            println(fileId, "          pmtcmt")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          rmtrac")
            println(fileId, "          0.5d0 ")
            println(fileId, "          tatrac")
            println(fileId, "          0.5d0")
        end
        if cfg.subDim == "3D"
            #transition
            println(fileId, "pmtracm   ")
            println(fileId, "0.8d0")
            println(fileId, "          pmtcmt")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          rmtrac")
            println(fileId, "          0.5d0 ")
            println(fileId, "          tatrac")
            println(fileId, "          0.5d0")
            #rotation
            println(fileId, "pmrotate ")
            println(fileId, "1.0d0")
            println(fileId, "          pmromt  ")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          rmrot")
            println(fileId, "          0.05d0")
            println(fileId, "          tarot")
            println(fileId, "          0.5d0")
        end
        if cfg.subDim == "4D"
            println(fileId, "pmcb      ")
            println(fileId, "0.3d0    ")
            println(fileId, "          pmcbmt")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          pmall")
            println(fileId, "          0.0d0 ")
            #transition
            println(fileId, "pmtracm   ")
            println(fileId, "0.8d0")
            println(fileId, "          pmtcmt")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          rmtrac")
            println(fileId, "          0.5d0 ")
            println(fileId, "          tatrac")
            println(fileId, "          0.5d0")
            #rotation
            println(fileId, "pmrotate ")
            println(fileId, "1.0d0")
            println(fileId, "          pmromt  ")
            print(fileId, "           ")
            for i = 1:cfg.subNum
                print(fileId, i / cfg.subNum, "   ")
            end
            print(fileId, "\n")
            println(fileId, "          rmrot")
            println(fileId, "          0.05d0")
            println(fileId, "          tarot")
            println(fileId, "          0.5d0")
        end
    end
    println(fileId,"cbmc_formulation")
    println(fileId, "'Martin and Siepmann 1999 + Martin and Thompson 2004'")
    println(fileId, "cbmc_setting_style")
    println(fileId, "'default ideal'")
    println(fileId, "ffnumber")
    println(fileId, 1)
    println(fileId, "ff_filename")
    println(fileId, "../",cfg.ffname)
    println(fileId, "classical_potential")
    println(fileId, "'Lennard-Jones'")
    println(fileId, "classical_mixrule")
    println(fileId, "'Lorentz-Berthelot'")
    println(fileId, "lshift")
    println(fileId, ".false.")
    println(fileId, "ltailc")
    println(fileId, ".true.")
    println(fileId, "rmin  ")
    println(fileId, 0.1)
    println(fileId, "rcut  ")
    println(fileId, 15.0)
    println(fileId, "rcutin  ")
    println(fileId, 15.0)
    println(fileId, "electrostatic_form")
    println(fileId, "'none'")
    for i = 1:cfg.subNum
        println(fileId, "input_style")
        println(fileId, "'basic connectivity map'")
        println(fileId, "nunit")
        println(fileId, onemol[i].atomNum)
        println(fileId, "nmaxcbmc")
        println(fileId, onemol[i].atomNum)
        println(fileId, "lpdbnames")
        println(fileId, ".false.")
        println(fileId, "forcefield")
        println(fileId, "TraPPE-UA ")
        println(fileId, "charge_assignment")
        println(fileId, "'bond increment'")
        for j = 1:onemol[i].atomNum
            println(fileId, "unit ntype")
            println(fileId, adj_left(j, 6), "'", adj_left(onemol[i].nsFF[j], 6), "'")
            println(fileId, "vibration")
            println(fileId, size(onemol[i].connect[j],1))
            for k = 1:size(onemol[i].connect[j],1)
                print(fileId, onemol[i].connect[j][k], "  ")
            end
            print(fileId, "\n")
            println(fileId, "improper torsion")
            println(fileId, 0)
        end
    end
    close(fileId)
    println(onemol[2].connect)
end

function towhee_coords_initial(onemol, cfg, nmol, ndx)
    scale = ceil(Int64, cfg.inputInitMolecules[ndx] / sum(cfg.inputMolecules[ndx,:]))
    h = (scale * sum(cfg.inputMolecules[ndx,:]) / (cfg.inputDensity[ndx] * 6.02214076 /10^4))^(1/3)
    h = h - 1
    fileId = open("packmol.inp", "w")
    print(fileId, "#\n#\n#\n\n")
    println(fileId, "tolerance 3.0")
    println(fileId, "filetype pdb")
    println(fileId, "output towhee_coords.pdb")
    
    for i = 1:cfg.subNum
        println(fileId, "structure ../$(cfg.subName[i])")
        println(fileId, " number $(nmol[i] * scale)")
        println(fileId, "inside box 0. 0. 0.  $(h)  $(h)  $(h)")
        println(fileId, "end structure\n")
    end
    
    close(fileId)
    
    run(pipeline(`packmol`, stdin="packmol.inp"))
    pdb_2_xyz("towhee_coords.pdb", "towhee_coords")
    
end

function pdb_2_xyz(infile, outfile)
    fileIdIn = open(infile)
    fileIdOut = open(outfile, "w")
    for i in eachline(fileIdIn)
        if length(i) >6
            if i[1:6] == "HETATM" || i[1:6] == "ATOM  "
                print(fileIdOut, parse(Float64, i[31:38]) + 0.5, "  ", parse(Float64, i[39:46]) + 0.5, "  ", parse(Float64, i[47:54])+0.5,"\n")
            end
        end
    end
    close(fileIdOut)
    close(fileIdIn)
end

