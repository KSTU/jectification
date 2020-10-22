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
include("towhee_add.jl")
include("towhee_input_gen.jl")
include("current_state.jl")
include("towhee_get_prop.jl")

cfg = read_config("config.in")
onemol = collect(read_pdb(i) for i in cfg.subName)

#println(onemol)


#generate input flows
println("condition ", cfg.condition)
if cfg.condition == "start"
    folder_check(length(cfg.inputPlates), "flow")   #creat and free folders
    if cfg.back == "gomc"
        generate_gomc_flow(onemol, cfg)
    elseif cfg.back == "towhee"
        generate_towhee_flow(onemol, cfg)
    end
end


#run input flows in parralel
equlibrationNumber = 4
productionNumber = 5

equlibrationCheck = zeros(Int64, cfg.inputNumber)   #array of equlibration for each thread
equlibrationId = zeros(Int64, cfg.inputNumber)
equlibrationProp = Array{AbstractArray}(undef, cfg.inputNumber)  #Array of flow properties
for i = 1:cfg.inputNumber
    equlibrationProp[i] = Array{Properties}(undef, equlibrationNumber)
    for j = 1:equlibrationNumber
        equlibrationProp[i][j] = set_properties(cfg)
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
#                    if equlibrationId[nThread] > equlibrationNumber
#                        equlibrationId[nThread] = equlibrationNumber + 1
#                    end
                    
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
inputProp = Array{Properties}(undef, cfg.inputNumber)
for i = 1:cfg.inputNumber
    inputProp[i] = set_properties(cfg)
end

@sync begin
    for nThread = 1:cfg.inputNumber
        @async begin
            read_input_prop(inputProp, nThread)
        end
    end
end

plateProp = Array{Properties}(undef, cfg.totalPlates)
eqProp = Array{Properties}(undef, cfg.totalPlates)
state = Array{CState}(undef, cfg.totalPlates)
for i = 1:cfg.totalPlates
    plateProp[i] = set_properties(cfg)
    state[i] = set_state(cfg)
    eqProp[i] = set_properties(cfg)
end

#------ set initial state
cd(cfg.dir) #go to initial path
if cfg.condition == "start" || cfg.condition == "rectstart"
    folder_check(cfg.totalPlates, "plate")
    set_initial_state(cfg.initialState)
elseif cfg.condition == "restart"
    #read last state
    
end

loopEnd = 0
equlibrationCheck = zeros(Int64, cfg.totalPlates)
#println("inputProp $(inputProp[1])")
println("debug")
#while loopEnd == 0
for loopEnd = 1:1
    @sync begin
        for nThread = 1:cfg.totalPlates
            #zeroes current state
            set_initial_state(state[nThread])    #set initial state
            #add molecules
            towhee_add(cfg, plateProp, inputProp, nThread)
            towhee_box_prepare(cfg, plateProp, inputProp, nThread)
            #println("reference energy ", plateProp[nThread])
            #calculate if moplecules exist
            if plateProp[nThread].liqN + plateProp[nThread].vapN > 0
                cd("plate$(nThread)")
                towhee_box_generate(cfg, plateProp, inputProp, nThread)
                #box convert pdb to cooords
                towhee_pdb_to_cords(["liquid.pdb"])
                #constant energy
                towhee_input_gen(onemol, cfg, nThread, state[nThread].stepCoef, true, "energy", plateProp)
                println("start energy $(plateProp[nThread].liqN * 30 * state[nThread].stepCoef) MC steps")
                run(pipeline(`towhee`, stdout = "energy"))
                #constant pressure
                towhee_pdb_to_cords(["box_01_step_$(string(plateProp[nThread].liqN * 30 * state[nThread].stepCoef, pad=14)).pdb"])
                println("start pressure $(plateProp[nThread].liqN * 30 * state[nThread].stepCoef) MC steps")
                towhee_input_gen(onemol, cfg, nThread, state[nThread].stepCoef, true, "pressure", plateProp)
                run(pipeline(`towhee`, stdout = "pressure"))
                #dual box cycle
                towhee_pdb_to_cords(["box_01_step_$(string(plateProp[nThread].liqN * 30 * state[nThread].stepCoef, pad=14)).pdb"])
                #get new vapor volume
                towhee_get_out_prop("pressure", eqProp, nThread, true)
                plateProp[nThread].vapVol = cfg.plateVolume[nThread] - eqProp[nThread].liqVol
                while state[nThread].energyInit == true || (plateProp[nThread].liqEnergy + plateProp[nThread].vapEnergy - plateProp[nThread].refEnergy) / plateProp[nThread].refEnergy > 0.05
                    #equlibrate energy
                    while state[nThread].energyEqulibrate == false
                        towhee_input_gen(onemol, cfg, nThread, state[nThread].stepCoef, state[nThread].firstrun, "NVT-GEMC", plateProp)
                        println("start energy GEMC $(plateProp[nThread].vapN) $((plateProp[nThread].liqN + plateProp[nThread].vapN) * 30 * state[nThread].stepCoef) MC steps")
                        run(pipeline(`towhee`, stdout = "energy-gemc"))
                        towhee_get_out_prop("energy-gemc", eqProp, nThread, true)
                        state[nThread].firstrun = false
                        mv("towhee_final", "towhee_initial", force = true)
                        #state[nThread].energyEqulibrate == true
                        
                    end
                    #check energy
                    
                    #equlibrate pressure and concentration
                    
                    state[nThread].energyInit = false
                end
                #productation cycle
            
                #swich molecules
            
                #delete molecules
                cd("..")
            end
            #
        end
    end
    #loopEnd = 1
end



#@async run(pipeline(`towhee`, `tee equlibration`))



