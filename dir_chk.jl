##
function folder_check(dirname)
    for i = 1:size(cfg.inputPlates,1)
        if(isdir("$(dirname)$(i)"))
            rm("$(dirname)$(i)", recursive=true)
        end
        mkdir("$(dirname)$(i)")
    end
end


function generate_flow(onemol, cfg)
    mainFolder = pwd()
    for ndx = 1:size(cfg.inputPlates,1)
        #if NVT ensable
        cd("flow$(ndx)")
        println("|",cfg.inputEnsamble[ndx],"|")
        #set atom coordinates
        box = init_box_generate(cfg, onemol, ndx)
        #write pdb
        write_pdb(box[1],"box.pdb")
        #set topology
        psf_generate(cfg, box[1], onemol, "box.psf")
        #set configuration
        configuration_set(mainFolder, cfg, ndx)
        cd("..")
    end
end

"""
generate 

"""
function init_box_generate(cfg, onemol, ndx)
    box = set_config
    #recalculate molecule numbers   ##set recalculate in configure
    
    if(cfg.inputInitMolecules[ndx] > 4000)
        molsum = sum(cfg.inputMolecules[ndx,:])
        ratio = ceil(Int64, cfg.inputInitMolecules[ndx] / molsum)
        molnum = cfg.inputMolecules[ndx,:] .* ratio
    else
        molsum = sum(cfg.inputMolecules[ndx,:])
        ratio = ceil(Int64, 4000 / molsum)
        molnum = cfg.inputMolecules[ndx,:] .* ratio
    end

    #calculate box size
    if cfg.inputEnsamble[ndx] == "NVT"
        cfg.inputVolume[ndx] = cfg.inputInitMolecules[ndx] / (cfg.inputDensity[ndx] * 6.02214076 / 10^4)
    else
        sumatom = maximum(collect(onemol[i].atomNum for i = 1:cfg.subNum))
        cfg.inputVolume[ndx] = cfg.inputInitMolecules[ndx] * sumatom * 3.0^2 * pi
        println(sumatom)
    end
    #calculate coordinates
    box1 = set_molecule()
    #get total atom numbers
    box1.molNum = sum(molnum)
    box1.atomNum = sum([onemol[i].atomNum * molnum[i] for i = 1:cfg.subNum])
    println("inp mol ", onemol[2].atomNum, " molnum", molnum, " bmol ", box1.molNum," batom ",  box1.atomNum)
    
    box1.curAtom = Array{Int64}(undef, box1.atomNum)
    box1.curMol = Array{Int64}(undef, box1.atomNum)
    box1.atomName = Array{String}(undef, box1.atomNum)
    box1.molName = Array{String}(undef, box1.atomNum)
    box1.x = Array{Float64}(undef, box1.atomNum)
    box1.y = Array{Float64}(undef, box1.atomNum)
    box1.z = Array{Float64}(undef, box1.atomNum)
    box1.element = Array{String}(undef, box1.atomNum)
    
    ncube = ceil(Int64, box1.molNum^(1/3))
    tempx = Array{Float64}(undef, ncube^3)
    tempy = Array{Float64}(undef, ncube^3)
    tempz = Array{Float64}(undef, ncube^3)
    tempfree = zeros(Int64, ncube^3)
    id = 1
    for i = 1:ncube
        for j = 1:ncube
            for k = 1:ncube
                tempx[id] = (i-0.5)/ncube * cfg.inputVolume[ndx]^(1.0/3.0)
                tempy[id] = (j-0.5)/ncube * cfg.inputVolume[ndx]^(1.0/3.0)
                tempz[id] = (k-0.5)/ncube * cfg.inputVolume[ndx]^(1.0/3.0)
                id += 1
            end
        end
    end
    id = 1
    id2  = 1
    for i = 1:cfg.subNum
        insert = 0
        while insert < molnum[i]
            rmol = ceil(Int64, rand() * ncube^3)
            if tempfree[rmol] == 0
                for j = 1:onemol[i].atomNum
                    box1.curAtom[id] = id
                    box1.curMol[id] = id2
                    box1.atomName[id] = onemol[i].atomName[j]
                    box1.molName[id] = onemol[i].molName[j]
                    box1.element[id] = onemol[i].element[j]
                    box1.x[id] = tempx[rmol] + onemol[i].x[j]
                    box1.y[id] = tempy[rmol] + onemol[i].y[j]
                    box1.z[id] = tempz[rmol] + onemol[i].z[j]
                    #
                    if box1.x[id] < 0
                        box1.x[id] += cfg.inputVolume[ndx]^(1.0/3.0)
                    elseif box1.x[id] > cfg.inputVolume[ndx]^(1.0/3.0)
                        box1.x[id] -= cfg.inputVolume[ndx]^(1.0/3.0)
                    end
                    if box1.y[id] < 0
                        box1.y[id] += cfg.inputVolume[ndx]^(1.0/3.0)
                    elseif box1.y[id] > cfg.inputVolume[ndx]^(1.0/3.0)
                        box1.y[id] -= cfg.inputVolume[ndx]^(1.0/3.0)
                    end
                    if box1.z[id] < 0
                        box1.z[id] += cfg.inputVolume[ndx]^(1.0/3.0)
                    elseif box1.z[id] > cfg.inputVolume[ndx]^(1.0/3.0)
                        box1.z[id] -= cfg.inputVolume[ndx]^(1.0/3.0)
                    end
                    
                    
                    id += 1
                end
                tempfree[rmol] = 1
                id2 += 1
                insert += 1
            end
        end
    end
    
    
    return [box1]
end

function configuration_set(mainFolder, cfg, ndx)
    cp("$(mainFolder)/in.conf","in.conf")
    fileId = open("in.conf", "a")
    comln(fileId, "Parameters", "$(mainFolder)/ff.inp", "GENERATED")
    if cfg.inputEnsamble[ndx] == "GEMC_NPT" || cfg.inputEnsamble[ndx] == "GEMC_NVT"
        comln(fileId, "Coordinates 0", "BOX0.pdb")
        comln(fileId, "Coordinates 1", "BOX1.pdb")
        comln(fileId, "Structure 0", "BOX0.psf")
        comln(fileId, "Structure 1", "BOX1.psf")
    else
        comln(fileId, "Coordinates 0", "BOX0.pdb")
        comln(fileId, "Structure 0", "BOX0.psf")
    end
    comln(fileId, "Temperature", cfg.inputTemperature[ndx])
    if cfg.inputEnsamble[ndx] == "NVT"
        comln(fileId, "Rcut", 10)
    else
        comln(fileId, "Rcut", 10)
    end
    if cfg.inputEnsamble[ndx] == "GEMC_NPT" || cfg.inputEnsamble[ndx] == "NPT"
        comln(fileId, "Pressure", cfg.inputPressure[ndx])
    end
    
    #output file
    
    close(fileId)
end

function comln(fileId, var, val, comment = "1")
    if length(comment) > 1
        println(fileId, """#########################
        # $(comment)
        #########################""")
    end
    println(fileId, "$(var)\t$(val)")
end

function write_pdb(box, filename)
    fileId = open(filename, "w")
    println(fileId, "REMARK jectification generated")
    for i = 1:box.atomNum
        print(fileId, adj_left("ATOM",6))               #1-6
        print(fileId, adj_right(box.curAtom[i],5))      #7-11
        print(fileId, " ")                              #12whitespace
        print(fileId, adj_right(box.atomName[i],4))     #13-16
        print(fileId, " ")                              #17
        print(fileId, adj_right(box.molName[i],3))      #18-20
        print(fileId, " ")                              #21
        print(fileId, "A")                              #22
        print(fileId, adj_right(box.curMol[i],4))       #23-26
        print(fileId, " "^4)                            #27-30
        print(fileId, adj_right(box.x[i],8,3))          #31-38
        print(fileId, adj_right(box.y[i],8,3))          #39-46
        print(fileId, adj_right(box.z[i],8,3))          #47-54
        print(fileId, adj_right(1.0, 6, 2))             #55-60
        print(fileId, adj_right(0.0, 6, 2))             #61-66
        print(fileId," "^6)                             #67-72
        print(fileId, adj_left("M1", 4))                #73-76
        print(fileId, adj_right(box.element[i],2))      #77-78
        #print(fileId, adj_right(0.0,2))                 #79-80
        
        print(fileId, "\n")
    end
    print(fileId,"END")
    close(fileId)
end

function adj_left(var, s, dot=2)
    if typeof(var) == String
        var = lstrip(rstrip(var))
        if s > length(var)
            whadd = s - length(var)
            return var * " "^whadd
        else
            return var
        end
    elseif typeof(var) == Int
        des = floor(Int64, log10(var)) + 1
        if des < s
            return string(var) * " "^(s-des)
        else
            return string(var)
        end
    elseif typeof(var) == Float64
        return format(var, precision=s)
    end
end

function adj_right(var, s, dot = 2)
    if typeof(var) == String
        var = lstrip(rstrip(var))
        if s > length(var)
            whadd = s - length(var)
            return  " "^whadd * var
        else
            return var
        end
    elseif typeof(var) == Int
        des = floor(Int64, log10(var)) + 1
        if des < s
            return  " "^(s-des) * string(var)
        else
            return string(var)
        end
    elseif typeof(var) == Float64
        if var > 1 
            des = floor(Int64, log10(var)) + 1
        else
            des = 1
        end
        if des < s - dot-1
            return " "^(s-dot-des-1) * string(format(var, precision=dot))
        else
            return format(var, precision=dot)
        end
    end
end

