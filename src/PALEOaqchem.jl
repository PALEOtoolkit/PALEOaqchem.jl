module PALEOaqchem



"""
    O2AlkUptakeRemin(Corg, (NO3, TNH3, Ngas), TPO4, Ccarb; rO2Corg=1) -> (O2, Alk)

Oxygen and alkalinity assimilated for production (or released by remineralisation) of particulate matter
with specified Corg, Ccarb from specified nitrate NO3, total ammonia TNH3, gaseous N, total phosphate TPO4

NB sign: for Corg +ve, O2 is -ve (ie quantity to subtract from tracer sms for production, or add for remineralisation)
"""
function O2AlkUptakeRemin(Corg, (NO3, TNH3, Ngas), TPO4, Ccarb; rO2Corg=1)

    O2 = -(rO2Corg*PB.get_total(Corg) +2.0*PB.get_total(NO3) + 3.0/4.0*PB.get_total(Ngas))  # -ve uptake/remin, ie oxygen released by production and consumed by remin
    Alk = - PB.get_total(NO3) + PB.get_total(TNH3) + 2.0*PB.get_total(Ccarb)                # production leads to increase in Alk with NO3 source, decrease with NH3 source

    # Biomass organic carbon stoichiometry is CH2Ox, where x follows from k_O_rO2Corg
    # CorgOx = 3 - 2*rO2Corg

    return (O2, Alk)
end


"""
    StoichPars() -> PB.ParametersTuple

Parameters defining particulate organic matter stoichiometry
"""
StoichPars() = PB.ParametersTuple(
    PB.ParDouble("rCorgPO4",  106.0, units="", description="Corg:P Redfield ratio"),
    PB.ParDouble("rNPO4",     16.0, units="",  description="N:P Redfield ratio"),
    PB.ParDouble("rCcarbCorg", 0.0, units="",  description="ratio of Ccarb to Corg produced"),
)

include("PALEOcarbchem/PALEOcarbchem.jl") # PALEOcarbchem module

include("CarbChem.jl")

include("Remin.jl")

include("SecondaryRedox.jl")

include("Particle.jl")

include("Boron.jl")


end # module PALEOaqchem
