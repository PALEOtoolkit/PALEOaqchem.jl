module CFA

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionCFAsimple

CFA formation (minimal version dependent only on phosphorus concentration)

    P -> CFA

at a rate `R_CFA * (P_conc - Peq_conc)` where `Peq_conc` is the saturation concentration
of total dissolved phosphate (Van Cappellen & Berner 1988, Dale etal 2016).

CFA formation assumed to follow idealized stoichiometry, Dale etal GCA (2016) eqn. 10:

    10*Ca++ + (6-x)*PO_4--- + x*CO_3-- + (2+x)*F- --> Ca_10 (PO_4)_(6-x) (CO_3)_x F_(2+x)
    
with x = 0.26 (NB: rates are normalized to mol P yr-1 here ie divide by (6-x)).

# References
Van Cappellen & Berner (1988) AJS <https://dx.doi.org/10.2475/ajs.288.4.289>
Dale etal GCA (2016) GCA <https://dx.doi.org/10.1016/j.gca.2016.05.046>

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionCFAsimple{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_CFA", 1.0 , units="yr-1",
            description="CFA formation rate constant"),
        PB.ParDouble("Peq_conc", 10e-6*1e3, units="mol m-3",
            description="saturation concentration of total dissolved phosphate below which CFA formation is inhibited"),
    )

    x::Float64 = 0.26 
    stoich_CFA = PB.RateStoich(
        PB.VarProp("rate_CFA", "mol P yr-1", "rate of CFA formation",
                        attributes=(:calc_total=>true,)),
        ((-10/(6-x), "Ca"), (-1.0, "P"), (-x/(6-x), "DIC"), (-(2+x)/(6-x), "F"), (1.0, "CFA"), (-2*10/(6-x), "TAlk")),  # TAlk from Ca consumption
        sms_prefix="", 
        sms_suffix="_sms",
        processname="coprecip",
    )

end

function PB.register_methods!(rj::ReactionCFAsimple)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        rj.stoich_CFA.ratevartemplate,
        PB.VarDep("P_conc", "mol m-3", "P concentration"),
    ]

    PB.add_method_do!(
        rj,
        do_CFAsimple, 
        (
            PB.VarList_namedtuple(vars),
        );
    )

    PB.add_method_do!(rj, rj.stoich_CFA)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end

function do_CFAsimple(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
   
    for i in cellrange.indices
        vars.rate_CFA[i] = max(0.0, vars.P_conc[i] - pars.Peq_conc[])
    end
   
    return nothing
end


end # module