module Remin

import PALEOboxes as PB
using PALEOboxes.DocStrings

"Stoichiometry for organic matter oxidation by O2"
const default_reminOrgOxO2 = PB.RateStoich(
    PB.VarProp("reminOrgOxO2", "mol O2 yr-1", "oxygen consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((1.0, "O2"),),
    sms_prefix="soluteflux_", 
    sms_suffix="",
    processname="remin",
)

"Stoichiometry and fractionation for organic matter oxidation by SO4.
 NB: normalized to O2eq"
const default_reminOrgOxSO4 = PB.RateStoich(
    PB.VarProp("reminOrgOxSO4", "mol O2eq yr-1", "2 * sulphate consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((0.5, "SO4::Isotope"), (-0.5, "H2S::Isotope"), (-1.0, "TAlk")),
    deltavarname_eta = ("SO4_delta", -35.0),   # constant fractionation 35 per mil
    sms_prefix="soluteflux_",
    sms_suffix="",
    processname="remin",
)

"""
Stoichiometry and fractionation for organic matter fermentation and methanogenesis, represented for stoichiometric purposes
as DIC -> CH4.
For programming convenience,  all Corg -> DIC already, so with a rate corresponding to 1/2 Corg -> 1/2 CH4 
the two reactions together represent Corg -> 1/2 DIC + 1/2 CH4.
Methane is produced at constant fractionation ~-35 per mil relative to the input organic matter Corg, so for the idealized stoichiometry
Corg -> 1/2 DIC + 1/2 CH4,  DIC is produced at +35 and CH4 at -35 relative to the input Corg.
"""
const default_reminOrgOxCH4 = PB.RateStoich(
    PB.VarProp("reminOrgOxCH4", "mol O2eq yr-1", "2 * DIC -> methane (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((0.5, "DIC::Isotope"), (-0.5, "CH4::Isotope")),
    deltavarname_eta = ("remin_Corg_delta", -35.0),   # constant fractionation 35 per mil relative to input Corg
    sms_prefix="soluteflux_",
    sms_suffix="",
    processname="remin",
)   

"""
    ReactionReminPonly

Organic particulate matter remineralization (no oxidant use)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminPonly{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionReminPonly)

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, rj.pars.CIsotope[])
   
    PB.add_method_do!(rj, do_remin_Ponly, (input_particulate_fluxes, output_solute_fluxes))

    # PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_Ponly(
    m::PB.ReactionMethod,
    (input_particulate_fluxes, output_solute_fluxes),
    cellrange::PB.AbstractCellRange,
    deltat
)

    @inbounds for i in cellrange.indices
        O2_eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=false)
    end

    return nothing
end


"""
    ReactionReminO2

Organic particulate matter remineralization (O2 oxidant only)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionReminO2)

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, rj.pars.CIsotope[])
   
    vars = [default_reminOrgOxO2.ratevartemplate]

    PB.add_method_do!(
        rj, 
        do_remin_O2, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vars)),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2(
    m::PB.ReactionMethod,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange, deltat
)

    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=false)
        vars.reminOrgOxO2[i] = O2eq
    end
 
    return nothing
end


"""
    ReactionReminO2_SO4

Organic particulate matter remineralization (O2, SO4 oxidants)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2_SO4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("oxreminlimit", 8e-3, units="mol m-3", 
            description="oxygen concentration below which use of O2 for remineralisation is inhibited"),            
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
        PB.ParType(PB.AbstractData, "SIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
          
    )
end


function PB.register_methods!(rj::ReactionReminO2_SO4)
    
    SIsotopeType    = rj.pars.SIsotope[]
    @info "register_methods! $(PB.fullname(rj)) SIsotopeType=$(SIsotopeType)"

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, rj.pars.CIsotope[])

    vars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        default_reminOrgOxO2.ratevartemplate,
        default_reminOrgOxSO4.ratevartemplate,
    ]

    PB.add_method_do!(
        rj, 
        do_remin_O2_SO4, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vars)),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)
    PB.add_method_do!(rj, default_reminOrgOxSO4, isotope_data=SIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2_SO4(
    m::PB.ReactionMethod,
    pars,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)

    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=false)
        freminOrgO2 = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))
        freminOrgSO4 = 1.0 - freminOrgO2
        vars.reminOrgOxO2[i] = O2eq*freminOrgO2
        vars.reminOrgOxSO4[i] = O2eq*freminOrgSO4 # NB: rate is already normalized to O2eq yr-1
    end
    
    return nothing
end


"""
    ReactionReminO2_SO4_CH4

Organic particulate matter remineralization (O2, SO4 oxidants, remaining Corg to CH4)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2_SO4_CH4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("oxreminlimit", 8e-3, units="mol m-3", 
            description="oxygen concentration below which use of O2 for remineralisation is inhibited"),

        PB.ParDouble("SO4reminlimit", 1000e-3, units="mol m-3", 
            description="sulphate concentration below which use of SO4 for remineralisation is inhibited"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
        PB.ParType(PB.AbstractData, "SIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionReminO2_SO4_CH4)
    
    CIsotopeType    = rj.pars.CIsotope[]
    SIsotopeType    = rj.pars.SIsotope[]
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType), SIsotopeType=$(SIsotopeType)"

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, CIsotopeType, add_input_Corg_delta=true)

    vars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        # ScalarData as we only want total, not isotopes (if any)
        PB.VarDep("SO4_conc","mol m-3", "SO4 concentration"; attributes=(:field_data=>PB.ScalarData,)),
        default_reminOrgOxO2.ratevartemplate,
        default_reminOrgOxSO4.ratevartemplate,
        default_reminOrgOxCH4.ratevartemplate,
    ]

    PB.add_method_do!(
        rj, 
        do_remin_O2_SO4_CH4, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vars)),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)
    PB.add_method_do!(rj, default_reminOrgOxSO4, isotope_data=SIsotopeType)
    PB.add_method_do!(rj, default_reminOrgOxCH4, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2_SO4_CH4(
    m::PB.ReactionMethod,
    pars,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=false)

        freminOrgO2     = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))

        freminOrgSO4    = (1.0 - freminOrgO2)*max(vars.SO4_conc[i], 0.0)/(pars.SO4reminlimit[] + max(vars.SO4_conc[i], 0.0))

        freminOrgCH4    = 1.0 - freminOrgO2 - freminOrgSO4

        vars.reminOrgOxO2[i] = O2eq*freminOrgO2 
        vars.reminOrgOxSO4[i] = O2eq*freminOrgSO4  # NB: rate already normalized to O2eq yr-1
        # CO2 -> CH4 + 2 O2eq, so each O2eq corresponds to 0.5 CH4 released
        # this maintains redox balance, but as no N cycle imagines nitrate and CH4 are produced ! (the least bad compromise)
        vars.reminOrgOxCH4[i] = O2eq*freminOrgCH4  
    end

    return nothing
end



"""
    vars_ReminParticulateFluxes(rj::PB.AbstractReaction, CIsotopeType::Type{<: PB.AbstractData}) -> (input_particulate_fluxes, output_solute_fluxes)

Return Variables needed for particulate input and solute output, without considering oxidants.

Optionally add Corg_delta to `input_particulate_fluxes`.
"""
function vars_ReminParticulateFluxes(
    rj::PB.AbstractReaction, CIsotopeType::Type{<: PB.AbstractData};
    add_input_Corg_delta=false
)
    
    # named tuple of particulate input fluxes
    input_particulate_fluxes = PB.Fluxes.FluxTarget(
        "remin_", ["P", "N", "Corg::$CIsotopeType", "Ccarb::$CIsotopeType"],
    )

    # named tuple of solute output fluxes
    output_solute_fluxes = PB.Fluxes.FluxContrib(
        "soluteflux_", ["P", "(TNH3)", "(DIC::$CIsotopeType)", "(TAlk)"],
        alloptional=false
    )
   
    if add_input_Corg_delta && (CIsotopeType <: PB.AbstractIsotopeScalar)
        input_particulate_fluxes = merge(
            input_particulate_fluxes,
            (Corg_delta = PB.VarProp("remin_Corg_delta", "per mil", "d13C of Corg"),),
        )
    end
       
    return (
        PB.VarList_namedtuple_fields(input_particulate_fluxes),
        PB.VarList_namedtuple_fields(output_solute_fluxes),
    )
end

# TODO Return remin stoichiometry in form needed by a Reaction to implement PALEOboxes.get_rate_stoichiometry(rj::AbstractReaction)"
# function get_rate_stoichiometry(vars_rpf)
# (will require additional variables to track P, C, N fluxes)
    

"""
    do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes); Ncycle::Bool) -> O2eq

Add solute_fluxes for particulate input in cell i, return `O2eq` oxygen demand
"""
function do_ReminParticulateFluxes(
    i, (input_particulate_fluxes, output_solute_fluxes);
    Ncycle::Bool
)
    P = input_particulate_fluxes.P[i]
    N = input_particulate_fluxes.N[i]
    Corg = input_particulate_fluxes.Corg[i]
    Ccarb = input_particulate_fluxes.Ccarb[i]

    output_solute_fluxes.P[i] += P

    if Ncycle
        Nremin = (0.0, N, 0.0)  # all as NH3
        output_solute_fluxes.TNH3[i] += N
    else
        Nremin = (N, 0.0, 0.0)  # all implicitly as NO3, no N cycle in this configuration
    end
        
    (O2eq, TAlk) = PALEOreactions.BioGeoChem.O2AlkUptakeRemin(Corg, Nremin, P, Ccarb)

    PB.add_if_available(output_solute_fluxes.DIC, i, Corg + Ccarb)
    PB.add_if_available(output_solute_fluxes.TAlk, i, TAlk)

    if hasfield(typeof(input_particulate_fluxes), :Corg_delta)
        input_particulate_fluxes.Corg_delta[i] = PB.get_delta(Corg, eps())  # guard against 0/0
    end

    return O2eq

end


#= 
"""
    ReactionRemin

Organic particulate matter remineralization

# Implementation:
This is implemented as two (almost) separate pieces:
- transfer of remineralized C, N, P to solute, with corresponding oxidant demand `O2eq`
- reduction of oxidants to provide `O2eq`

"""
Base.@kwdef mutable struct ReactionRemin <: PB.AbstractReaction
    base::PB.ReactionBase

    par_reminpathway::PB.ParString = PB.ParString("reminpathway", "Ponly", allowed_values=["Ponly", "O2", "O2_SO4", "O2_SO4_CH4"])

    par_Ncycle::PB.ParBool= PB.ParBool("Ncycle",false,   description="true to remineralize to NH3 for configuration with N cycle, false to implicitly remin to NO3")

    vars_soluteflux_P = PB.VarContrib("soluteflux_P", "mol yr-1", "P remineralized from organic particulate", PB.RP_React) 
    vars_soluteflux_TNH3 = PB.VarContrib("(soluteflux_TNH3)", "mol yr-1", "total NH3 from remineralization of organic particulate", PB.RP_React)
    vars_soluteflux_DIC = PB.VarContrib("(soluteflux_DIC)", "mol yr-1", "DIC from remineralization of organic particulate", PB.RP_React)
    vars_soluteflux_TAlk = PB.VarContrib("(soluteflux_TAlk)", "mol yr-1", "TAlk from remineralization of organic particulate", PB.RP_React)

    vars_input = nothing

    oxidants = nothing

end


function PB.initialize_parameters(rj::ReactionRemin)

    isotope_type = PB.VC_Single  # no isotopes yet
    
    rj.vars_input = PALEOreactions.BioGeoChem.fluxBioProdTarget(rj.par_fluxlist.v, prefix="remin_", isotope_type=isotope_type)
    PB.add_var(rj, rj.vars_input)

    if rj.par_reminpathway.v == "Ponly"
        rj.oxidants = Oxidants_Ponly
    elseif rj.par_reminpathway.v == "O2"
        rj.oxidants = Oxidants_O2
    else
        error("invalid reminpathway='$(rj.par.reminpathway.v)")
    end

    return nothing
end


function PB.initialize_data(rj::ReactionRemin, modeldata::PB.ModelData)
    
    return (rj.oxidants,
            PB.create_accessors_oxidants(rj.oxidants, modeldata),
            PB.create_accessors_namedtuple(rj.vars_remin, modeldata),
            PB.create_accessors_namedtuple(rj.variables, modeldata))
end

function PB.do_react(rj::ReactionRemin, (oxidants, vars_oxidants, vars_remin, vars), cellrange::PB.AbstractCellRange, deltat)

    for i in cellrange.indices
        (P, N, Corg, Ccarb) = (vars_remin.P[i], vars_remin.N[i], vars_remin.Corg[i], vars_remin.Ccarb[i])

        vars.soluteflux_P[i] += vars_remin.P[i]

        if rj.par_Ncycle.v
            Nremin = (0.0, N, 0.0)  # all as NH3
            vars.soluteflux_TNH3[i] += N
        else
            Nremin = (N, 0.0, 0.0)  # all implicitly as NO3, no N cycle in this configuration
        end
            
        (O2eq, TAlk) = PALEOreactions.BioGeoChem.O2AlkUptakeRemin(Corg, Nremin, P, Ccarb)

        PB.add_if_available(vars.soluteflux_DIC, i, Corg + Ccarb)
        PB.add_if_available(vars.soluteflux_TAlk, i, TAlk)

        do_oxidants(rj, oxidants, (vars_oxidants, vars), i, O2eq)
    end

    return nothing
end


Base.@kwdef mutable struct Oxidants_Ponly
end

function create_Oxidants_Ponly(rj::ReactionRemin)
    return Oxidants_Ponly()
end

function create_accessors_oxidants(oxidants::Oxidants_Ponly, modeldata::PB.ModelData)
    return nothing
end

function do_oxidants(rj::ReactionRemin, oxidants::Oxidants_Ponly, (vars_oxidants, vars), idx, O2eq)
end


Base.@kwdef mutable struct Oxidants_O2
    soluteflux_O2 = PB.VarContrib("soluteflux_O2", "mol yr-1", "O2 from remineralization of organic particulate", PB.RP_React) 
end

function create_Oxidants_O2(rj::ReactionRemin)
    oxidants = Oxidants_O2()
    PB.add_vars_fields(rj, oxidants)
end

function do_oxidants(rj::ReactionRemin, oxidants::Oxidants_O2, (vars_oxidants, vars), idx, O2eq)
    vars_oxidants.oxidantflux_O2[idx] += O2eq
end

function create_accessors_oxidants(oxidants::Oxidants_O2, modeldata::PB.ModelData)
    return PB.create_accessors_namedtuple_fields(oxidants, modeldata)
end

 =#


end # module
