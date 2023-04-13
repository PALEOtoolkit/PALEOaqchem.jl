module Boron

import PALEOboxes as PB
using PALEOboxes.DocStrings

# import Infiltrator # Julia debugger

"""
    ReactionBoronIsotope

Calculate d11B for aqueous B(OH)4- and B(OH)3 species from mass balance, given  d11B for total B, 
total B concentration `B_conc` and B(OH)4- concentration `BOH4_conc`.

See eg [Zeebe2001](@cite), p220.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionBoronIsotope{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("alphaB", 1.0272, units="",
            description="isotopic fractionation factor B(OH)4m <-> B(OH)3"),

        PB.ParType(PB.AbstractData, "BIsotope", PB.IsotopeLinear,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable boron isotopes and specify isotope type"),
    )

end

function PB.register_methods!(rj::ReactionBoronIsotope)

    BIsotopeType = rj.pars.BIsotope[]
    PB.setfrozen!(rj.pars.BIsotope)

    vars = [
        PB.VarDep( "B_conc",      "mol m-3", "B concentration",
            attributes=(:field_data=>BIsotopeType,)),
        PB.VarDep( "BOH4_conc",   "mol m-3", "B(OH)4m concentration"),
        PB.VarDep( "B_delta",     "per mil", "d11B delta for total B"),

        PB.VarProp( "BOH4_delta",  "per mil", "d11B delta for B(OH)4- species"),
    ]

    PB.add_method_do!(rj, do_boron_isotope, (PB.VarList_namedtuple(vars),))

    return nothing
end

"Calculate d11BOH4m given xBO4 = [B(OH)4-]/[B_total], fractionation alphaB between B(OH)4- and B(OH)3, and delta d11Btot"
function d11BOH4m(d11Btot, xBO4, alphaB)
    d11BOH4m = (d11Btot - 1000.0*(1.0-xBO4)*(alphaB-1.0))/(alphaB-xBO4*(alphaB-1.0))
    return d11BOH4m
end

function do_boron_isotope(m::PB.ReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        xBO4 = vars.BOH4_conc[i]/PB.get_total(vars.B_conc[i])
        vars.BOH4_delta[i] = d11BOH4m(vars.B_delta[i], xBO4, pars.alphaB[])
    end

    return nothing
end


end
