mutable struct properties
    totEnergy::Float64
    refEnergy::Float64
    density::Float64
    pressure::Float64
    #liquid molecules
    liqN::Int64
    liqNType::AbstractArray
    liqEnergy::Float64
    liqDens::Float64
    liqPres::Float64
    #vapor molecules
    vapN::Int64
    vapNType::AbstractArray
    vapEnergy::Float64
    vapDens::Float64
    vapPres::Float64
end

function set_properties()
    return properties(0.0,  #total energy
                     0.0,   #reference energy
                     0.0,   #density
                     0.0,   #pressure
                     0,     #n liq
                     [0, 0],    #n type liq
                     0.0,       #liq energy
                     0.0,       #liq dens
                     0.0,        #liq press
                     0,         #n vap
                     [0, 0],    #n type liq
                     0.0,       #vap energy
                     0.0,       #vap dens
                     0.0        #vap pressure
                     )
end
