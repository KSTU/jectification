mutable struct molecule
    atomNum::Int64  #numbers of atoms
    molNum::Int64   #numbers of molecules
    
    curAtom::AbstractArray  #current atom number
    curMol::AbstractArray   #current molele
    atomName::AbstractArray #atom name
    molName::AbstractArray  #mol name
    x::AbstractArray        #atom ccordinates
    y::AbstractArray
    z::AbstractArray
    
    #connectivity
    connect::AbstractArray  #
    element::AbstractArray  #element name
    
    #non standart
    nsQ::AbstractArray #charge
    nsMass::AbstractArray #mass
    nsFF::AbstractArray #atomname in force field
    #molType::AbstractArray #numbers of each type
end

function set_molecule()
    return molecule(0,0,[0],[0],["0"],["0"],[0.0],[0.0],[0.0],[0,0],["a", "a"], [0.0, 0.0], [0.0, 0.0], ["A", "B"])
end

function read_pdb(fileName)
    mol = set_molecule()
    #get atom numbers
    fileId = open(fileName, "r")
    for str in eachline(fileId)
        if length(str) > 6
            if str[1:6] == "ATOM  " || str[1:6] == "HETATM"
                mol.atomNum += 1
            end
        end
    end
    seekstart(fileId)
    temp = collect(str for str in eachline(fileId) if length(str) > 6)  #need 2 step if length < 6
    atomString = collect(str for str in temp if str[1:6] == "ATOM  " || str[1:6] == "HETATM")    #get atomstring
    
    mol.atomName = collect(lstrip(str[13:16]) for str in atomString)
    mol.molName = collect(lstrip(str[18:20]) for str in atomString)
    mol.x = collect(parse(Float64,str[31:38]) for str in atomString)
    mol.y = collect(parse(Float64,str[39:46]) for str in atomString)
    mol.z = collect(parse(Float64,str[47:54]) for str in atomString)
    
    mol.curAtom = collect(parse(Int64,str[7:11]) for str in atomString)
    mol.curMol = collect(parse(Int64,str[23:26]) for str in atomString)
    mol.element = collect(lstrip(str[77:78]) for str in atomString)
    
    mol.molNum = maximum(collect(parse(Int64,str[23:26]) for str in atomString))
    
    #nonstandart properties
    mol.nsQ = Array{Float64}(undef, mol.atomNum)
    mol.nsMass = Array{Float64}(undef, mol.atomNum)
    mol.nsFF = Array{String}(undef, mol.atomNum)
    id = 0
    for str in atomString
        println("deb ", str)
        nsString = split(str,"!")[2]
        strArray = split(nsString)
        id += 1
        for i = 1:size(strArray,1)
            if strArray[i] == "charge"
                mol.nsQ[id] = parse(Float64, strArray[i+1])
            end
            if strArray[i] == "mass"
                mol.nsMass[id] = parse(Float64, strArray[i+1])
            end
            if strArray[i] == "ffaname"
                mol.nsFF[id] = strArray[i+1]
            end
        end
    end
    
    #mol.nsQ = collect( parse(Float64, ) for str in temp if str[1:6] == "ATOM" || str[1:6] == "HETATM" )
    
    
    println(mol.atomName, " check", mol.nsQ, mol.nsMass, mol.nsFF)
    
    #connectivity
    mol.connect = Array{AbstractArray}(undef, mol.atomNum)
    consize = zeros(Int64, mol.atomNum)
    for i in temp
        if i[1:6] == "CONECT"
            conatom = parse(Int64, split(i)[2])
            consize[conatom] += size(split(i),1) - 2
        end
    end
    println("consioze", consize)
    for i = 1:mol.atomNum
        mol.connect[i] = Array{Int64}(undef, consize[i])
    end
    
    #
    for cm = 1:mol.atomNum
        #mol.connect[cm]
        mol.connect[cm] = collect(parse(Int64,split(i)[j]) for i in temp for j = 3:size(split(i),1) if i[1:6] == "CONECT" && parse(Int64,split(i)[2]) == cm)
        #test = (parse(Int64,split(i)[j]) for i in temp for j = 3:size(split(i),1) if i[1:6] == "CONECT" && parse(Int64,split(i)[2]) == cm )
        
    end
    #println(test)
#    
    #println(mol.atomNum)
    #println("test ", mol.connect)
    close(fileId)
    return mol
end
