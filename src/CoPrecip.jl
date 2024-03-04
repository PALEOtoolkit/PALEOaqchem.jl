module CoPrecip

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionPACoPrecip

Co-precipitation of P (eg iron-bound phosphorus)  with A (eg Fe oxide) formation

    P -> P=A

at a rate `gamma * A_formation_rate_<n>`, with limitation at low P concentration

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionPACoPrecip{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("A_rate_stoich_factors", [1.0],
            description="stoichiometry factor to multiply each A formation rate variable to convert to mol A"),
        PB.ParDouble("gamma", 0.15 , units="mol/mol",
            description="P:A molar ratio"),
        PB.ParDouble("P_limit", 1e-6, units="mol m-3",
            description="limiting P concentration below which co-precipitation is inhibited"),
    )

    stoich_PA_coprecip = PB.RateStoich(
        PB.VarProp("rate_PA_coprecip", "mol P yr-1", "rate of P co-precipitation",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "P"), (1.0, "PA"),), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="coprecip",
    )

end

function PB.register_methods!(rj::ReactionPACoPrecip)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        rj.stoich_PA_coprecip.ratevartemplate,
        PB.VarDep("P_conc", "mol m-3", "P concentration"),
    ]

    rate_vars = [
        PB.VarDep("A_formation_rate_$i", "mol m-3 yr-1", "substance A formation rate") for i in 1:length(rj.pars.A_rate_stoich_factors.v)
    ]

    PB.add_method_do!(
        rj,
        do_PA_coprecip, 
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_tuple(rate_vars),
        );
    )

    PB.add_method_do!(rj, rj.stoich_PA_coprecip)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end

function do_PA_coprecip(m::PB.ReactionMethod, pars, (vars, rate_vars), cellrange::PB.AbstractCellRange, deltat)
   
    for i in cellrange.indices
        P_conc_checked = max(0.0, vars.P_conc[i])
        P_lim_fac = P_conc_checked / (P_conc_checked + pars.P_limit[])

        vars.rate_PA_coprecip[i] = 0.0

        PB.IteratorUtils.foreach_longtuple(rate_vars, pars.A_rate_stoich_factors.v) do rv, stoich_factor
            vars.rate_PA_coprecip[i] += P_lim_fac*pars.gamma[]*rv[i]*stoich_factor
        end
    end
   
    return nothing
end


"""
    ReactionPARelease

Release of P (eg iron-bound phosphorus) with A (eg Fe oxide) destruction

    P=A -> P

at a rate `Prelease = theta * A_destruction_rate_<n>`, where `theta = PA_conc / A_conc`

`Prelease` is then partioned into two fluxes:

    P_sms = (1 - partion_frac) * Prelease
    P2_sms = partition_frac * Prelease

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionPARelease{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("A_rate_stoich_factors", [1.0],
            description="stoichiometry factor to multiply each A destruction rate variable to convert to mol A"),
        PB.ParDouble("partition_frac", 0.0,
            description="fraction of P to partition into P2_sms")
    )

end

function PB.register_methods!(rj::ReactionPARelease)

    @info "register_methods! $(PB.fullname(rj))"

    stoich_PA_release = PB.RateStoich(
        PB.VarProp("rate_PA_release", "mol P yr-1", "rate of coprecipitated P dissolution",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "PA"), ((1.0 - rj.pars.partition_frac[]), "P"), (rj.pars.partition_frac[], "P2")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="coprecip",
    )
    PB.setfrozen!(rj.pars.partition_frac)

    vars = [
        stoich_PA_release.ratevartemplate,
        PB.VarDep("PA_conc", "mol m-3", "adsorbed P concentration"),
        PB.VarDep("A_conc", "mol m-3", "adsorbant concentration"),
        PB.VarProp("PA_theta", "mol/mol", "P:A molar ratio of adsorbed of coprecipitated P"),
    ]

    rate_vars = [
        PB.VarDep("A_destruction_rate_$i", "mol m-3 yr-1", "substance A destruction rate") for i in 1:length(rj.pars.A_rate_stoich_factors.v)
    ]

    PB.add_method_do!(
        rj,
        do_PA_release, 
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_tuple(rate_vars),
        );
    )

    PB.add_method_do!(rj, stoich_PA_release)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end

function do_PA_release(m::PB.ReactionMethod, pars, (vars, rate_vars), cellrange::PB.AbstractCellRange, deltat)
   
    for i in cellrange.indices
        vars.PA_theta[i] = clamp(vars.PA_conc[i] / vars.A_conc[i], 0.0, 1.0) # enforce +ve and not unreasonably high for numerical stability

        vars.rate_PA_release[i] = 0.0

        PB.IteratorUtils.foreach_longtuple(rate_vars, pars.A_rate_stoich_factors.v) do rv, stoich_factor
            vars.rate_PA_release[i] += vars.PA_theta[i]*rv[i]*stoich_factor
        end
    end
   
    return nothing
end


end # module