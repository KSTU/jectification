mutable struct config
    #substance
    subNum::Int64
    subName::AbstractArray
    #flows
    totalPlates::Int64
    inputPlates::AbstractArray
    inputNumber::Int64
    inputMolecules::AbstractArray
    inputInitMolecules::AbstractArray
    inputEnsamble::AbstractArray
    initialState::String
    inputTemperature::AbstractArray
    inputPressure::AbstractArray
    inputDensity::AbstractArray
    inputVolume::AbstractArray
    
    plateVolume::AbstractArray
end

function set_config()
    return config(0,["0","0"],0,[0,0],0,[0,0],[0,0],["NVT","NVT"],"NONE",[1.0, 1.0],[1.0, 1.0], [1.0, 1.0],[1.0,1.0], [1.0, 1.0])
end

function read_config(fileName)
cfg = set_config()
fileId = open(fileName, "r")
    for str in eachline(fileId)
        rem = lstrip(str)
        if(SubString(rem,1) != '#')
            #read substance number
            if rstrip(split(rem, "=")[1]) == "sub_name"
                write = split(split(rem, "=")[2])
                cfg.subName = write
                cfg.subNum = size(write,1)
                #println("substances ", cfg.subNum, cfg.subName)
                
                
            end
            # read total plates
            if rstrip(split(rem, "=")[1]) == "total_plates"
                cfg.totalPlates = parse(Int64,split(rem, "=")[2])
            end
            # read input plates
            if rstrip(split(rem, "=")[1]) == "input_plates"
                write = collect(parse(Int64,i) for i in split(split(rem, "=")[2]))
                cfg.inputPlates = write
                cfg.inputNumber = size(write,1)
                #println("inputs ", cfg.inputPlates, cfg.inputNumber)
                cfg.inputMolecules = Array{Int64}(undef, cfg.inputNumber, cfg.subNum)
            end

            for flow = 1:cfg.inputNumber
                if rstrip(split(rem, "=")[1]) == string("input_molecules",flow)
                    if size(split(split(rem, "=")[2]),1) != cfg.subNum
                        get_error("set flows to all molecules", @__FILE__, @__LINE__)
                    end
                    write = collect(parse(Int64,i) for i in split(split(rem, "=")[2]))
                    cfg.inputMolecules[flow,:] = write
                    println("inputs ", write, " ", cfg.inputMolecules, cfg.inputNumber,"",  cfg.subNum)
                end
            end
            
            
            # read initial state
            if rstrip(split(rem, "=")[1]) == "initial_state"
                cfg.initialState = rstrip(lstrip(split(rem, "=")[2]))
            end
            #read initial molecules
            if rstrip(split(rem, "=")[1]) == "input_init_molecules"
                write = collect(parse(Int64,i) for i in split(split(rem, "=")[2]))
                cfg.inputInitMolecules = write
            end
            #read initial state
            if rstrip(split(rem, "=")[1]) == "input_ensamble"
                write = collect( lstrip(rstrip(i)) for i in (split(split(rem, "=")[2])))
                cfg.inputEnsamble = write
            end
            #read initial temperature
            if rstrip(split(rem, "=")[1]) == "input_temperature"
                write = collect(parse(Float64,i) for i in (split(split(rem, "=")[2])))
                cfg.inputTemperature = write
            end
            #read initial pressure
            if rstrip(split(rem, "=")[1]) == "input_pressure"
                write = collect(parse(Float64,i) for i in (split(split(rem, "=")[2])))
                cfg.inputPressure = write
            end
            #read initial density
            if rstrip(split(rem, "=")[1]) == "input_density"
                write = collect(parse(Float64,i) for i in (split(split(rem, "=")[2])))
                cfg.inputDensity = write
            end
            #read plate volume
            if rstrip(split(rem, "=")[1]) == "input_density"
                write = collect(parse(Float64,i) for i in (split(split(rem, "=")[2])))
                if(size(write,1) == 1)
                    cfg.plateVolume = fill(write, cfg.inputPlates)
                else
                    cfg.plateVolume = write
                end
            end
            
            
            #set calculated pllates parameters
            cfg.inputVolume = Array{Float64}(cfg.inputPlates)
            
        end
    end
close(fileId)
return cfg
println("config read done")
end

function get_error(error, file, line)
    println(error, " file name: $(file) line number $(line) ")
    exit()
end



