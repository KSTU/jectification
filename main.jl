using Format
using Dates


include("input.jl")
include("read_pdb.jl")
include("dir_chk.jl")
include("psf_gen.jl")
include("gomc_init_flow.jl")
include("towhee_init_flow.jl")
include("towhee_check_equlibration.jl")


cfg = read_config("config.in")
onemol = collect(read_pdb(i) for i in cfg.subName)

#println(onemol)


#generate input flows
folder_check("flow")

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

@sync begin
    for nThread = 1:cfg.inputNumber
        @async begin
            println("start $(nThread) flow calculation $(Dates.Time(Dates.now()))")
            cd("$(cfg.dir)/flow$(nThread)")
            equlibrationId[nThread] = 1
            while equlibrationCheck[nThread] == 0
                run(pipeline(`towhee`, `tee equlibration`))
                println("stop $(nThread) $(pwd()) $(Dates.Time(Dates.now()))")
                towhee_equlibrate(equlibrationCheck, equlibrationId)
                equlibrationId[nThread] += 1    #add to current ID
                if equlibrationId[nThread] > equlibrationNumber
                    equlibrationId = equlibrationNumber
                end
            end
            #cd("..")
        end
    end
end

#@async run(pipeline(`towhee`, `tee equlibration`))



