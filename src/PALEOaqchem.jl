module PALEOaqchem

import PALEOboxes as PB

"""
    O2AlkUptakeRemin(Corg, (NO3, TNH3, Ngas), TPO4, Ccarb; rO2Corg=1) -> (O2, Alk)

Oxygen and alkalinity released by remineralisation
of particulate matter with specified Corg, Ccarb to specified nitrate NO3, total ammonia TNH3, gaseous N, total phosphate TPO4

NB sign: 
- for remineralization, Corg, NO3, TNH3, Ngas, TPO4, Ccarb are +ve, O2 is -ve 
  (ie provides O2eq (-ve) required and Alk to add to ocean solute tracers sms, for remineralisation of POC with C:N:P stoichiometry defined by Corg:sum(NO3, TNH3, Ngas):TPO4)
- for production by oxygenic photosynthesis, Corg, NO3, TNH3, Ngas, Ccarb are -ve, O2 is +ve 
 (ie provides O2 and Alk to add to ocean solute tracers sms for production of POC with C:N:P stoichiometry defined by Corg:sum(NO3, TNH3, Ngas):TPO4)

# Examples
```jldoctest; setup=:(import PALEOaqchem)
julia> PALEOaqchem.O2AlkUptakeRemin(106.0, (0.0, 0.0, 0.0), 1.0, 0.0) # Corg:P = 106:1, no N
(-106.0, -1.0)

julia> PALEOaqchem.O2AlkUptakeRemin(106.0, (16.0, 0.0, 0.0), 1.0, 0.0) # Corg:NO3:P = 106:16:1
(-138.0, -17.0)

julia> PALEOaqchem.O2AlkUptakeRemin(106.0, (0.0, 16.0, 0.0), 1.0, 0.0) # Corg:TNH3:P = 106:16:1
(-106.0, 15.0)
```
"""
function O2AlkUptakeRemin(Corg, (NO3, TNH3, Ngas), TPO4, Ccarb; rO2Corg=1)

    O2 = -(rO2Corg*PB.get_total(Corg) +2.0*PB.get_total(NO3) + 3.0/4.0*PB.get_total(Ngas))  # -ve uptake/remin, ie oxygen released by production and consumed by remin          
    Alk = - PB.get_total(NO3) + PB.get_total(TNH3) - PB.get_total(TPO4) + 2.0*PB.get_total(Ccarb)  # production leads to increase in Alk with NO3 source, decrease with NH3 source

    # Biomass organic carbon stoichiometry is CH2Ox, where x follows from k_O_rO2Corg
    # CorgOx = 3 - 2*rO2Corg

    return (O2, Alk)
end

"""
    parse_number_name(nname::AbstractString) -> (number, name)

Parse a string of form "-1*A"

"""
function parse_number_name(nname::AbstractString; sep=['*', ' '], number_first=true)

    # parse multiplier
    svmn = split(nname, sep, keepempty=false)
    if length(svmn) == 1
        mult, name = 1, svmn[1]
    elseif length(svmn) == 2
        if !number_first
            tmp = svmn[1]
            svmn[1] = svmn[2]
            svmn[2] = tmp
        end
        mult = tryparse(Int64, svmn[1])
        if isnothing(mult)
            mult = tryparse(Float64, svmn[1])
        end
        name = svmn[2]
    end

    !isnothing(mult) || 
        error("invalid field in nname, not of form number*name: ", nname)
  
    return (mult, name)
end

parse_name_to_power_number(nname::AbstractString) = parse_number_name(nname; sep=['^', ' '], number_first=false)

const _R_conc_attributes_base = (
    # :field_data=>rj.pars.field_data[],
    :totalnames=>String[],
    # :advect=>false,
    :vertical_movement=>0.0,
    :specific_light_extinction=>0.0,
    :vphase=>PB.VP_Undefined,
    :diffusivity_speciesname=>"",
    :gamma=>missing,
)

const R_conc_attributes_advectfalse = (_R_conc_attributes_base..., :advect=>false,)
const R_conc_attributes_advecttrue = (_R_conc_attributes_base..., :advect=>true,)


include("PALEOcarbchem/PALEOcarbchem.jl") # PALEOcarbchem module

include("Reservoirs.jl")

include("GenericReactions.jl")

include("CarbChem.jl")

include("Remin.jl")

include("SecondaryRedox.jl")

include("Particle.jl")

include("RCmultiG.jl")

include("Boron.jl")

include("CoPrecip.jl")

include("CFA.jl")

include("MolecularDiffusion.jl")

include("FeS_deprecated.jl")

end # module PALEOaqchem
