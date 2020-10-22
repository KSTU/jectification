mutable struct CState
    energyInit::Bool
    energyEqulibrate::Bool
    firstrun::Bool
    stepCoef::Int64
    curEqBlock::AbstractArray   #current block number for all plates
    eqNumber::Int64
end

function set_state(cfg)
    return CState(true, 
        true,
        true,
        1,
        fill(1, cfg.plate),
        3
    )
end

function set_initial_state(state::CState)
    state.energyInit = true
    state.energyEqulibrate = false
    state.firstrun = true
    state.stepCoef = 1
end


