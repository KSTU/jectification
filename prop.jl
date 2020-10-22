mutable struct Properties
    totEnergy::Float64
    refEnergy::Float64
    density::Float64
    pressure::Float64
    volume::Float64
    temp::Float64
    NType::AbstractArray
    NDens::AbstractArray
    #liquid molecules
    liqN::Int64
    liqNType::AbstractArray
    liqEnergy::Float64
    liqDens::Float64
    liqPres::Float64
    liqVol::Float64
    #vapor molecules
    vapN::Int64
    vapNType::AbstractArray
    vapEnergy::Float64
    vapDens::Float64
    vapPres::Float64
    vapVol::Float64
    #input properties
    inDens::AbstractArray
    inPres::AbstractArray
    inNumber::AbstractArray
    inTemp::AbstractArray
end

function set_properties(cfg)
    zint = zeros(Int64, cfg.subNum)
    zfloat = zeros(Float64, cfg.subNum)
    return Properties(0.0,  #total energy
                     0.0,   #reference energy
                     0.0,   #density
                     0.0,   #pressure
                     0.0,   #volume
                     0.0,   #temperature
                     zeros(Int64, cfg.subNum),  #n type
                     zeros(Float64, cfg.subNum),    #number density of substances
                     0,     #n liq
                     zeros(Int64, cfg.subNum),    #n type liq
                     0.0,       #liq energy
                     0.0,       #liq dens
                     0.0,        #liq press
                     0.0,       #liq volume
                     0,         #n vap
                     zeros(Int64, cfg.subNum),    #n type liq
                     0.0,       #vap energy
                     0.0,       #vap dens
                     0.0,        #vap pressure
                     0.0,        #vap volume
                     zeros(Float64, 3), #inputDensity
                     zeros(Float64, 3),  #inputPressure
                     zeros(Int64, 3),    #input molecules
                     zeros(Float64, 3)  #input temperatures
                     )
end
