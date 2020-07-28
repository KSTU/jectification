using Format
using Dates


include("input.jl")
include("read_pdb.jl")
include("dir_chk.jl")
include("psf_gen.jl")
include("gomc_init_flow.jl")
include("towhee_init_flow.jl")
include("towhee_check_equlibration.jl")
include("prop.jl")
include("towhee_init_dual.jl")


cfg = read_config("config.in")
onemol = collect(read_pdb(i) for i in cfg.subName)

#println(onemol)


#generate input flows
println("condition ", cfg.condition)
if cfg.condition == "start"
    folder_check(length(cfg.inputPlates), "flow")   #creat and free folders
    folder_check(cfg.totalPlates, "plate")
end

if cfg.back == "gomc"
    generate_gomc_flow(onemol, cfg)
elseif cfg.back == "towhee"
    generate_towhee_flow(onemol, cfg)
end

#run input flows in parralel
equlibrationNumber = 4
productionNumber = 5

equlibrationCheck = zeros(Int64, cfg.inputNumber)   #array of equlibration for each thread
equlibrationId = zeros(Int64, cfg.inputNumber)
equlibrationProp = Array{AbstractArray}(undef, cfg.inputNumber)  #Array of flow properties
for i = 1:cfg.inputNumber
    equlibrationProp[i] = Array{properties}(undef, equlibrationNumber)
    for j = 1:equlibrationNumber
        equlibrationProp[i][j] = set_properties()
    end
end
println("debug ", typeof(equlibrationProp), length(equlibrationProp))


if cfg.condition == "start"
    @sync begin
        for nThread = 1:cfg.inputNumber
            @async begin
                println("start $(nThread)  flow calculation $(Dates.Time(Dates.now()))")
                cd("$(cfg.dir)/flow$(nThread)")
                equlibrationId[nThread] = 1
                #### EQULIBRATION
                while equlibrationCheck[nThread] == 0
                    run(pipeline(`towhee`, stdout = "equlibration"))
                    println("stop $(nThread) $(pwd()) $(Dates.Time(Dates.now()))")
                    towhee_equlibrate(equlibrationProp[nThread], equlibrationCheck, equlibrationId[nThread], nThread, equlibrationNumber)
                    #rm("towhee_initial")
                    cp("towhee_final", "towhee_initial", force = true)
                    if equlibrationId[nThread] == 1     #change initial to false
                        #cd("$(cfg.dir)/flow$(nThread)")
                        towhee_initial(onemol, cfg, nThread, 1, false)
                        println("changed $(nThread)")
                    end
                    equlibrationId[nThread] += 1    #add to current ID
                    if equlibrationId[nThread] > equlibrationNumber
                        equlibrationId[nThread] = equlibrationNumber
                    end
                    
                end
                #### PRODUCTION
                println("start productation $(nThread) $(pwd()) $(Dates.Time(Dates.now()))")
                towhee_initial(onemol, cfg, nThread, productionNumber, false)
                run(pipeline(`towhee`, stdout = "production"))
                #cd("..")
            end
        end
    end
end

#get input propertyes
inputProp = Array{properties}(undef, cfg.inputNumber)
for i = 1:cfg.inputNumber
    inputProp[i] = set_properties()
end



@sync begin
    for nThread = 1:cfg.inputNumber
        @async begin
            read_input_prop(inputProp, nThread)
        end
    end
end



plateProp = Array{properties}(undef, cfg.totalPlates)
for i = 1:cfg.totalPlates
    plateProp[i] = set_properties()
end
#------ set initial state
if cfg.condition == "start"
    set_initial_state(cfg.initial_state)
else
    #read last state
end


endloop = 1

println("inputProp $(inputProp[1])")

while endloop == 0
    @sync begin
        for nThread = 1:cfg.totalPlates
            #add molecules
            
            #calculate
            
            #swich molecules
            
            #delete molecules
            
            
        end
    end
end



#@async run(pipeline(`towhee`, `tee equlibration`))



