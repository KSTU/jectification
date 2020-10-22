function towhee_input_gen(onemol, cfg, ndx, nloop, initial, genType, prop)
    fileId = open("towhee_input", "w")
    println(fileId, "ensemble")
    if genType == "energy"
        println(fileId, "'", "nvt", "'")
        println(fileId, "temperature")
        println(fileId, prop[ndx].temp)
    elseif genType == "pressure"
        println(fileId, "'npt'")
        println(fileId, "temperature")
        println(fileId, prop[ndx].temp)
        println(fileId, "pressure")
        println(fileId, prop[ndx].pressure)
    elseif genType == "NVT-GEMC"
        println(fileId, "'nvt'")
        println(fileId, "temperature")
        println(fileId, prop[ndx].temp)
    end
    println(fileId, "nmolty")
    println(fileId, cfg.subNum)
    println(fileId, "nmolectyp")
    if genType == "energy" || genType == "pressure"
        for i = 1:cfg.subNum
            print(fileId, prop[ndx].liqNType[i], "  ")
        end
        print(fileId, "\n")
    else    #GEMC
        for i = 1:cfg.subNum
            #println(i, "  ", prop[ndx].liqNType[i], "  ", prop[ndx].vapNType[i])
            print(fileId, prop[ndx].liqNType[i] + prop[ndx].vapNType[i], "  ")
        end
        print(fileId, "\n")
    end
    println(fileId, "numboxes")
    if genType == "energy" || genType == "pressure"
        println(fileId, "1")
    else
        println(fileId, "2")
    end
    println(fileId, "stepstyle")
    println(fileId, "'moves'")
    if genType == "energy" || genType == "pressure"
        println(fileId, "nstep")
        println(fileId, prop[ndx].liqN * 30 * nloop)
        println(fileId, "printfreq")
        println(fileId, 1000)
        println(fileId, "blocksize")
        println(fileId, prop[ndx].liqN * 30 * nloop)
        println(fileId, "moviefreq")
        println(fileId, 1000)
        println(fileId, "backupfreq")
        println(fileId, 5000)
        println(fileId, "runoutput")
        println(fileId, "full")
        println(fileId, "pdb_output_freq")
        println(fileId, prop[ndx].liqN * 30 * nloop)
        println(fileId, "pressure_virial_freq")
        println(fileId, 10)
        println(fileId, "trmaxdispfreq")
        println(fileId, 10000)
        println(fileId, "volmaxdispfreq")
        println(fileId, 10000)
    elseif genType == "NVT-GEMC"
        println(fileId, "nstep")
        println(fileId, (prop[ndx].liqN + prop[ndx].vapN) * 30 * nloop)
        println(fileId, "printfreq")
        println(fileId, 1000)
        println(fileId, "blocksize")
        println(fileId, (prop[ndx].liqN + prop[ndx].vapN) * 30 * nloop)
        println(fileId, "moviefreq")
        println(fileId, 1000)
        println(fileId, "backupfreq")
        println(fileId, 5000)
        println(fileId, "runoutput")
        println(fileId, "full")
        println(fileId, "pdb_output_freq")
        println(fileId, (prop[ndx].liqN + prop[ndx].vapN) * 30 * nloop)
        println(fileId, "pressure_virial_freq")
        println(fileId, 10)
        println(fileId, "trmaxdispfreq")
        println(fileId, 10000)
        println(fileId, "volmaxdispfreq")
        println(fileId, 10000)
    end
    if genType == "energy" || genType == "pressure"
        println(fileId, "linit")
        if initial == true
            println(fileId, ".true.")
        else
            println(fileId, ".false.")
        end
        println(fileId, "initboxtype")
        println(fileId, "'dimensions'")
        println(fileId, "initstyle")
        for i = 1:cfg.subNum
            print(fileId, "'coords'", "    ")
        end
        print(fileId, "\n")
        println(fileId, "initlattice")
        for i = 1:cfg.subNum
            print(fileId, "'none'", "    ")
        end
        print(fileId, "\n")
        println(fileId, "initmol")
        for i = 1:cfg.subNum
            print(fileId, prop[ndx].liqNType[i], " ")
        end
        print(fileId, "\n")
        println(fileId, "inix iniy iniz")
        #calculate dimentions
        dim = ceil(Int64, (prop[ndx].liqN)^ (1/3))+1
        println(fileId, adj_left(dim, 6), adj_left(dim, 6), adj_left(dim, 6) )
        #calculate system size
        h = (prop[ndx].liqVol) ^ (1.0/3.0)
        println(fileId, "hmatrix")
        println(fileId, h, "   0   ", "    0    ")
        println(fileId, "   0   ", h,  "    0    ")
        println(fileId, "   0   ", "    0    ", h)
    elseif genType == "NVT-GEMC"
        println(fileId, "linit")
        if initial == true
            println(fileId, ".true.")
        else
            println(fileId, ".false.")
        end
        println(fileId, "initboxtype")
        println(fileId, "'dimensions'")
        println(fileId, "initstyle")
        for j = 1:2
            for i = 1:cfg.subNum
                print(fileId, "'coords'", "    ")
            end
            print(fileId, "\n")
        end
        println(fileId, "initlattice")
        for j =1:2
            for i = 1:cfg.subNum
                print(fileId, "'none'", "    ")
            end
            print(fileId, "\n")
        end
        println(fileId, "initmol")
        for i = 1:cfg.subNum
            print(fileId, prop[ndx].liqNType[i], "  ")
        end
        print(fileId, "\n")
        for i = 1:cfg.subNum
            print(fileId, prop[ndx].vapNType[i], "  ")
        end
        print(fileId, "\n")
        println(fileId, "inix iniy iniz")
        #calculate dimentions
        dim = ceil(Int64, (prop[ndx].liqN)^ (1/3))+1
        println(fileId, adj_left(dim, 6), adj_left(dim, 6), adj_left(dim, 6) )
        dim = ceil(Int64, (prop[ndx].vapN)^ (1/3))+3
        println(fileId, adj_left(dim, 6), adj_left(dim, 6), adj_left(dim, 6) )
        #calculate system size
        h = (prop[ndx].liqVol) ^ (1.0/3.0)
        println(fileId, "hmatrix")
        println(fileId, h, "   0   ", "    0    ")
        println(fileId, "   0   ", h,  "    0    ")
        println(fileId, "   0   ", "    0    ", h)
        h = (prop[ndx].vapVol) ^ (1.0/3.0)
        println(fileId, h, "   0   ", "    0    ")
        println(fileId, "   0   ", h,  "    0    ")
        println(fileId, "   0   ", "    0    ", h)
    end
    if genType == "pressure"
        println(fileId, "pmvol   ")
        println(fileId, "0.05d0")
        println(fileId, "pmvlpr   ")
        println(fileId, "1.0d0")
        println(fileId, "rmvol   ")
        println(fileId, "0.01d0")
        println(fileId, "tavol   ")
        println(fileId, "0.5d0")
    elseif genType == "NVT-GEMC" && initial == true
        println(fileId, "pmvol   ")
        println(fileId, "0.001d0")
        println(fileId, "pmvlpr   ")
        println(fileId, "1.0d0")
        println(fileId, "rmvol   ")
        println(fileId, "0.01d0")
        println(fileId, "tavol   ")
        println(fileId, "0.5d0")
    elseif genType == "NVT-GEMC" && initial == false
        println(fileId, "pmvol   ")
        println(fileId, "0.01d0")
        println(fileId, "pmvlpr   ")
        println(fileId, "1.0d0")
        println(fileId, "rmvol   ")
        println(fileId, "0.01d0")
        println(fileId, "tavol   ")
        println(fileId, "0.5d0")
    end
    if genType == "NVT-GEMC" && initial == true
        println(fileId, "pm2boxcbswap")
        println(fileId, "0.002d0")
        println(fileId, "pm2cbswmt")
        for i = 1:cfg.subNum
            print(fileId, i / cfg.subNum, "   ")
        end
        print(fileId, "\n")
        println(fileId, "pm2cbswpr")
        println(fileId, "1.0d0")
    elseif genType == "NVT-GEMC" && initial == false
        println(fileId, "pm2boxcbswap")
        println(fileId, "0.2d0")
        println(fileId, "pm2cbswmt")
        for i = 1:cfg.subNum
            print(fileId, i / cfg.subNum, "   ")
        end
        print(fileId, "\n")
        println(fileId, "pm2cbswpr")
        println(fileId, "1.0d0")
    end
    if genType == "energy" || genType == "pressure" || genType == "NVT-GEMC"
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
    println(fileId, "../", cfg.ffname)
    println(fileId, "classical_potential")
    println(fileId, "'Lennard-Jones'")
    println(fileId, "classical_mixrule")
    println(fileId, "'Lorentz-Berthelot'")
    println(fileId, "lshift")
    println(fileId, ".false.")
    println(fileId, "ltailc")
    println(fileId, ".true.")
    println(fileId, "rmin  ")
    println(fileId, 0.5)
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
end
