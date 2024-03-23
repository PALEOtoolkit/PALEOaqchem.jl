module Remin

import PALEOboxes as PB
import PALEOaqchem
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

"Stoichiometry for organic matter oxidation by NO2.
 NB: normalized to O2eq
    4/3 NO2- + 4/3 H+ -> 2/3 N2 + 2/3 H2O [+O2]
"
const default_reminOrgOxNO2 = PB.RateStoich(
    PB.VarProp("reminOrgOxNO2", "mol O2eq yr-1", "3/4 NO2 consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((4/3, "NO2"), (-2/3, "N2"), (-4/3, "TAlk")), # NB: stoichiometry is *-1
    sms_prefix="soluteflux_",
    sms_suffix="",
    processname="remin",
)

"Stoichiometry for organic matter oxidation by NO3.
 NB: normalized to O2eq
    2 NO3- -> 2 NO2- [+O2]
"
const default_reminOrgOxNO3NO2 = PB.RateStoich(
    PB.VarProp("reminOrgOxNO3NO2", "mol O2eq yr-1", "1/2 NO3 consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((2.0, "NO3"), (-2.0, "NO2")), # NB: stoichiometry is *-1
    sms_prefix="soluteflux_",
    sms_suffix="",
    processname="remin",
)

"Stoichiometry for organic matter oxidation by NO3 only (no NO2)
 NB: normalized to O2eq
    4/5 NO3- + 4/5 H+ -> 2/5 N2 [+O2]
"
const default_reminOrgOxNO3only = PB.RateStoich(
    PB.VarProp("reminOrgOxNO3only", "mol O2eq yr-1", "5/4 NO3 consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((4/5, "NO3"), (-2/5, "N2"), (-4/5, "TAlk")), # NB: stoichiometry is *-1
    sms_prefix="soluteflux_",
    sms_suffix="",
    processname="remin",
)

"Stoichiometry for organic matter oxidation by MnIVOx (MnO2)
    4 H+ + 2 MnIVOx -> 2 Mn++ + 2 H2O [+ O2]
"
const default_reminOrgOxMnIVOx = PB.RateStoich(
    PB.VarProp("reminOrgOxMnIVOx", "mol O2eq yr-1", "oxygen consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((2.0, "MnIVOx"), (-2.0, "MnII"), (-4.0, "TAlk")), # NB: stoichiometry is *-1
    sms_prefix="soluteflux_", 
    sms_suffix="",
    processname="remin",
)

"Stoichiometry for organic matter oxidation by FeIIIOx (0.5 Fe2O3)
    8H+ + 4 FeIIIOx -> 4 Fe++ + 4 H2O [+ O2]
"
const default_reminOrgOxFeIIIOx = PB.RateStoich(
    PB.VarProp("reminOrgOxFeIIIOx", "mol O2eq yr-1", "oxygen consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((4.0, "FeIIIOx"), (-4.0, "FeII"), (-8.0, "TAlk")), # NB: stoichiometry is *-1
    sms_prefix="soluteflux_", 
    sms_suffix="",
    processname="remin",
)

"Stoichiometry and fractionation for organic matter oxidation by SO4.
 NB: normalized to O2eq
    0.5 SO4-- + H+ -> 0.5 H2S [+O2]
"
const default_reminOrgOxSO4 = PB.RateStoich(
    PB.VarProp("reminOrgOxSO4", "mol O2eq yr-1", "2 * sulphate consumption (-ve) by remineralization",
        attributes=(:calc_total=>true,)),
    ((0.5, "SO4::Isotope"), (-0.5, "H2S::Isotope"), (-1.0, "TAlk")), # NB: stoichiometry is *-1
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

        fO2  = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))
        fSO4 = max(vars.SO4_conc[i], 0.0)/(pars.SO4reminlimit[] + max(vars.SO4_conc[i], 0.0))

        freminOrgO2 = fO2
        freminOrgSO4 = (1.0 - fO2)*fSO4
        freminOrgCH4 = 1.0 - freminOrgO2 - freminOrgSO4

        vars.reminOrgOxO2[i] = O2eq*freminOrgO2 
        vars.reminOrgOxSO4[i] = O2eq*freminOrgSO4  # NB: rate already normalized to O2eq yr-1
        # CO2 -> CH4 + 2 O2eq, so each O2eq corresponds to 0.5 CH4 released
        # this maintains redox balance, but as no N cycle imagines nitrate and CH4 are produced ! (the least bad compromise)
        vars.reminOrgOxCH4[i] = O2eq*freminOrgCH4  
    end

    return nothing
end

"""
    ReactionReminO2_Fe_SO4_CH4

Organic particulate matter remineralization (O2, FeIIIOx, SO4 oxidants, remaining Corg to CH4)

FeIIIOxreminlimit is from [VanCappellen1996a](@cite),
100 umol Fe(OH)3 g-1, assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3) = 100e-6*2.5e6 = 250 mol m-3 solid phase

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2_Fe_SO4_CH4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("oxreminlimit", 8e-3, units="mol m-3", 
            description="oxygen concentration below which use of O2 for remineralisation is inhibited"),

        PB.ParDouble("FeIIIOxreminlimit", 100e-6*2.5e6, units="mol m-3", 
            description="FeIII oxide concentration below which use of FeIII for remineralisation is inhibited"),

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

function PB.register_methods!(rj::ReactionReminO2_Fe_SO4_CH4)
    
    CIsotopeType    = rj.pars.CIsotope[]
    SIsotopeType    = rj.pars.SIsotope[]
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType), SIsotopeType=$(SIsotopeType)"

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, CIsotopeType, add_input_Corg_delta=true)

    vars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("FeIIIOx_conc", "mol m-3", "FeIII oxide concentration"),
        # ScalarData as we only want total, not isotopes (if any)
        PB.VarDep("SO4_conc","mol m-3", "SO4 concentration"; attributes=(:field_data=>PB.ScalarData,)),
        default_reminOrgOxO2.ratevartemplate,
        default_reminOrgOxFeIIIOx.ratevartemplate,
        default_reminOrgOxSO4.ratevartemplate,
        default_reminOrgOxCH4.ratevartemplate,
    ]

    PB.add_method_do!(
        rj, 
        do_remin_O2_Fe_SO4_CH4, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vars)),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)
    PB.add_method_do!(rj, default_reminOrgOxFeIIIOx)
    PB.add_method_do!(rj, default_reminOrgOxSO4, isotope_data=SIsotopeType)
    PB.add_method_do!(rj, default_reminOrgOxCH4, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2_Fe_SO4_CH4(
    m::PB.ReactionMethod,
    pars,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=false)

        fO2             = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))
        fFeIIIOx        = max(vars.FeIIIOx_conc[i], 0.0)/(pars.FeIIIOxreminlimit[] + max(vars.FeIIIOx_conc[i], 0.0))
        fSO4            = max(vars.SO4_conc[i], 0.0)/(pars.SO4reminlimit[] + max(vars.SO4_conc[i], 0.0))

        freminOrgO2     = fO2
        freminOrgFeIIIOx= (1.0 - fO2)*fFeIIIOx
        freminOrgSO4    = (1.0 - fO2)*(1.0 - fFeIIIOx)*fSO4
        freminOrgCH4    = 1.0 - freminOrgO2 - freminOrgFeIIIOx - freminOrgSO4

        # NB: stoichiometries already normalized to O2eq yr-1
        vars.reminOrgOxO2[i] = O2eq*freminOrgO2
        vars.reminOrgOxFeIIIOx[i] = O2eq*freminOrgFeIIIOx
        vars.reminOrgOxSO4[i] = O2eq*freminOrgSO4    
        # CO2 -> CH4 + 2 O2eq, so each O2eq corresponds to 0.5 CH4 released
        # this maintains redox balance, but as no N cycle imagines nitrate and CH4 are produced ! (the least bad compromise)
        vars.reminOrgOxCH4[i] = O2eq*freminOrgCH4  
    end

    return nothing
end


"""
    ReactionReminO2_NN_Mn_Fe_SO4_CH4

Organic particulate matter remineralization (succesive use of O2, NO2, NO3, MnIVOx, FeIIIOx, SO4 oxidants, remaining Corg to CH4)

- FeIIIOxreminlimit is from [VanCappellen1996a](@cite),
  100 umol Fe(OH)3 g-1, assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3) = 100e-6*2.5e6 = 250 mol m-3 solid phase
- MnIVOxreminlimit is from [VanCappellen1996a](@cite),
  16 umol MnO2 g-1, assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3) = 16e-6*2.5e6 = 40 mol m-3 solid phase

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2_NN_Mn_Fe_SO4_CH4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("oxreminlimit", 1e-3, units="mol m-3", 
            description="oxygen concentration below which use of O2 for remineralisation is inhibited"),

        PB.ParDouble("NO2reminlimit", 10e-3, units="mol m-3", 
            description="NO2 concentration below which use of NO2 for remineralisation is inhibited"),

        PB.ParDouble("NO3reminlimit", 10e-3, units="mol m-3", 
            description="NO3 concentration below which use of NO3 for remineralisation is inhibited"),

        PB.ParDouble("MnIVOxreminlimit", 16e-6*2.5e6, units="mol m-3", 
            description="MnIV oxide concentration below which use of MnIV oxide for remineralisation is inhibited"),

        PB.ParDouble("FeIIIOxreminlimit", 100e-6*2.5e6, units="mol m-3", 
            description="FeIII oxide concentration below which use of FeIII oxide for remineralisation is inhibited"),

        PB.ParDouble("SO4reminlimit", 500e-3, units="mol m-3", 
            description="sulphate concentration below which use of SO4 for remineralisation is inhibited"),

        PB.ParDouble("aSO4", 1.0,
            description="attenuation factor for sulphate reduction"),

        PB.ParBool("oxidant_dependent_rate", false,
            description="true to scale decomposition rates by factor 'freminOrgTot'"),

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

function PB.register_methods!(rj::ReactionReminO2_NN_Mn_Fe_SO4_CH4)
    
    CIsotopeType    = rj.pars.CIsotope[]
    SIsotopeType    = rj.pars.SIsotope[]
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType), SIsotopeType=$(SIsotopeType)"

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, CIsotopeType, add_input_Corg_delta=true)

    ox_limit_vars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("NO2_conc", "mol m-3", "NO2 concentration"),
        PB.VarDep("NO3_conc", "mol m-3", "NO3 concentration"),
        PB.VarDep("MnIVOx_conc", "mol m-3", "MnIV oxide concentration"),
        PB.VarDep("FeIIIOx_conc", "mol m-3", "FeIII oxide concentration"),
        # ScalarData as we only want total, not isotopes (if any)
        PB.VarDep("SO4_conc","mol m-3", "SO4 concentration"; attributes=(:field_data=>PB.ScalarData,)),
    ]

    ox_frac_vars = [
        PB.VarProp("freminOrgO2", "", "fraction of remineralization by O2"),
        PB.VarProp("freminOrgNO2", "", "fraction of remineralization by NO2"),
        PB.VarProp("freminOrgNO3", "", "fraction of remineralization by NO3"),
        PB.VarProp("freminOrgMnIVOx", "", "fraction of remineralization by MnIVOx"),
        PB.VarProp("freminOrgFeIIIOx", "", "fraction of remineralization by FeIIIOx"),
        PB.VarProp("freminOrgSO4", "", "fraction of remineralization by SO4"),
        PB.VarProp("freminOrgCH4", "", "fraction of remineralization by CH4"),
    ]
    if rj.pars.oxidant_dependent_rate[]
        push!(ox_frac_vars, PB.VarProp("freminOrgTot", "", "overall oxidant-dependent scaling factor for remineralization rate"))
    end

    PB.add_method_do!(
        rj, 
        do_remin_O2_NN_Mn_Fe_SO4_CH4_ox_frac, 
        (PB.VarList_namedtuple(vcat(ox_limit_vars, ox_frac_vars)),),
    )
    
    rate_vars = [
        default_reminOrgOxO2.ratevartemplate,
        default_reminOrgOxNO2.ratevartemplate,
        default_reminOrgOxNO3NO2.ratevartemplate,
        default_reminOrgOxMnIVOx.ratevartemplate,
        default_reminOrgOxFeIIIOx.ratevartemplate,
        default_reminOrgOxSO4.ratevartemplate,
        default_reminOrgOxCH4.ratevartemplate,
    ]

    PB.add_method_do!(
        rj, 
        do_remin_O2_NN_Mn_Fe_SO4_CH4_rates, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vcat(PB.VarDep.(ox_frac_vars), rate_vars))),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)
    PB.add_method_do!(rj, default_reminOrgOxNO3NO2)
    PB.add_method_do!(rj, default_reminOrgOxNO2)
    PB.add_method_do!(rj, default_reminOrgOxMnIVOx)
    PB.add_method_do!(rj, default_reminOrgOxFeIIIOx)
    PB.add_method_do!(rj, default_reminOrgOxSO4, isotope_data=SIsotopeType)
    PB.add_method_do!(rj, default_reminOrgOxCH4, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2_NN_Mn_Fe_SO4_CH4_ox_frac(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        fO2             = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))
        fNO2            = max(vars.NO2_conc[i], 0.0)/(pars.NO2reminlimit[] + max(vars.NO2_conc[i], 0.0))
        fNO3            = max(vars.NO3_conc[i], 0.0)/(pars.NO3reminlimit[] + max(vars.NO3_conc[i], 0.0))
        fMnIVOx         = max(vars.MnIVOx_conc[i], 0.0)/(pars.MnIVOxreminlimit[] + max(vars.MnIVOx_conc[i], 0.0))
        fFeIIIOx        = max(vars.FeIIIOx_conc[i], 0.0)/(pars.FeIIIOxreminlimit[] + max(vars.FeIIIOx_conc[i], 0.0))
        fSO4            = max(vars.SO4_conc[i], 0.0)/(pars.SO4reminlimit[] + max(vars.SO4_conc[i], 0.0))

        freminOrgO2         = fO2
        freminOrgNO2        = (1.0 - fO2)*fNO2
        freminOrgNO3        = (1.0 - fO2)*(1.0 - fNO2)*fNO3
        freminOrgMnIVOx     = (1.0 - fO2)*(1.0 - fNO2)*(1.0 - fNO3)*fMnIVOx
        freminOrgFeIIIOx    = (1.0 - fO2)*(1.0 - fNO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*fFeIIIOx
        freminOrgSO4        = pars.aSO4[]*(1.0 - fO2)*(1.0 - fNO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*(1.0 - fFeIIIOx)*fSO4        
        if pars.oxidant_dependent_rate[]
            freminOrgCH4    = pars.aSO4[]*(1.0 - fO2)*(1.0 - fNO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*(1.0 - fFeIIIOx)*(1.0 - fSO4)
            freminOrgTot    = freminOrgO2 + freminOrgNO2 + freminOrgNO3 + freminOrgMnIVOx + freminOrgFeIIIOx + freminOrgSO4 + freminOrgCH4
            vars.freminOrgTot[i]        = freminOrgTot
        else
            freminOrgCH4    = 1.0 - freminOrgO2 - freminOrgNO2 - freminOrgNO3 - freminOrgMnIVOx - freminOrgFeIIIOx - freminOrgSO4
            freminOrgTot    = 1.0
        end

        vars.freminOrgO2[i]         = freminOrgO2 / freminOrgTot
        vars.freminOrgNO2[i]        = freminOrgNO2 / freminOrgTot
        vars.freminOrgNO3[i]        = freminOrgNO3 / freminOrgTot
        vars.freminOrgMnIVOx[i]     = freminOrgMnIVOx / freminOrgTot
        vars.freminOrgFeIIIOx[i]    = freminOrgFeIIIOx / freminOrgTot
        vars.freminOrgSO4[i]        = freminOrgSO4 / freminOrgTot
        vars.freminOrgCH4[i]        = freminOrgCH4 / freminOrgTot
        
    end

    return nothing
end

function do_remin_O2_NN_Mn_Fe_SO4_CH4_rates(
    m::PB.ReactionMethod,
    pars,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=true) # remineralize N to ammonia
     
        # NB: stoichiometries already normalized to O2eq yr-1
        vars.reminOrgOxO2[i]        = O2eq*vars.freminOrgO2[i]
        vars.reminOrgOxNO2[i]       = O2eq*vars.freminOrgNO2[i]
        vars.reminOrgOxNO3NO2[i]    = O2eq*vars.freminOrgNO3[i]
        vars.reminOrgOxMnIVOx[i]    = O2eq*vars.freminOrgMnIVOx[i]
        vars.reminOrgOxFeIIIOx[i]   = O2eq*vars.freminOrgFeIIIOx[i]
        vars.reminOrgOxSO4[i]       = O2eq*vars.freminOrgSO4[i]
        # CO2 -> CH4 + 2 O2eq, so each O2eq corresponds to 0.5 CH4 released
        vars.reminOrgOxCH4[i]       = O2eq*vars.freminOrgCH4[i]
    end

    return nothing
end


"""
    ReactionReminO2_NO3_Mn_Fe_SO4_CH4

Organic particulate matter remineralization (succesive use of O2, NO3, MnIVOx, FeIIIOx, SO4 oxidants, remaining Corg to CH4)

- FeIIIOxreminlimit is from [VanCappellen1996a](@cite),
  100 umol Fe(OH)3 g-1, assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3) = 100e-6*2.5e6 = 250 mol m-3 solid phase
- MnIVOxreminlimit is from [VanCappellen1996a](@cite),
  16 umol MnO2 g-1, assuming dry density is 2.5 g/cm^3 (= 2.5e6 g m^-3) = 16e-6*2.5e6 = 40 mol m-3 solid phase

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReminO2_NO3_Mn_Fe_SO4_CH4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("oxreminlimit", 1e-3, units="mol m-3", 
            description="oxygen concentration below which use of O2 for remineralisation is inhibited"),

        PB.ParDouble("NO3reminlimit", 10e-3, units="mol m-3", 
            description="NO3 concentration below which use of NO3 for remineralisation is inhibited"),

        PB.ParDouble("MnIVOxreminlimit", 16e-6*2.5e6, units="mol m-3", 
            description="MnIV oxide concentration below which use of MnIV oxide for remineralisation is inhibited"),

        PB.ParDouble("FeIIIOxreminlimit", 100e-6*2.5e6, units="mol m-3", 
            description="FeIII oxide concentration below which use of FeIII oxide for remineralisation is inhibited"),

        PB.ParDouble("SO4reminlimit", 500e-3, units="mol m-3", 
            description="sulphate concentration below which use of SO4 for remineralisation is inhibited"),

        PB.ParDouble("aSO4", 1.0,
            description="attenuation factor for sulphate reduction"),

        PB.ParBool("oxidant_dependent_rate", false,
            description="true to scale decomposition rates by factor 'freminOrgTot'"),

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

function PB.register_methods!(rj::ReactionReminO2_NO3_Mn_Fe_SO4_CH4)
    
    CIsotopeType    = rj.pars.CIsotope[]
    SIsotopeType    = rj.pars.SIsotope[]
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType), SIsotopeType=$(SIsotopeType)"

    (input_particulate_fluxes, output_solute_fluxes) = vars_ReminParticulateFluxes(rj, CIsotopeType, add_input_Corg_delta=true)

    ox_limit_vars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("NO3_conc", "mol m-3", "NO3 concentration"),
        PB.VarDep("MnIVOx_conc", "mol m-3", "MnIV oxide concentration"),
        PB.VarDep("FeIIIOx_conc", "mol m-3", "FeIII oxide concentration"),
        # ScalarData as we only want total, not isotopes (if any)
        PB.VarDep("SO4_conc","mol m-3", "SO4 concentration"; attributes=(:field_data=>PB.ScalarData,)),
    ]

    ox_frac_vars = [
        PB.VarProp("freminOrgO2", "", "fraction of remineralization by O2"),
        PB.VarProp("freminOrgNO3", "", "fraction of remineralization by NO3"),
        PB.VarProp("freminOrgMnIVOx", "", "fraction of remineralization by MnIVOx"),
        PB.VarProp("freminOrgFeIIIOx", "", "fraction of remineralization by FeIIIOx"),
        PB.VarProp("freminOrgSO4", "", "fraction of remineralization by SO4"),
        PB.VarProp("freminOrgCH4", "", "fraction of remineralization by CH4"),
    ]
    if rj.pars.oxidant_dependent_rate[]
        push!(ox_frac_vars, PB.VarProp("freminOrgTot", "", "overall oxidant-dependent scaling factor for remineralization rate"))
    end

    PB.add_method_do!(
        rj, 
        do_remin_O2_NO3_Mn_Fe_SO4_CH4_ox_frac, 
        (PB.VarList_namedtuple(vcat(ox_limit_vars, ox_frac_vars)),),
    )
    
    rate_vars = [
        default_reminOrgOxO2.ratevartemplate,
        default_reminOrgOxNO3only.ratevartemplate,
        default_reminOrgOxMnIVOx.ratevartemplate,
        default_reminOrgOxFeIIIOx.ratevartemplate,
        default_reminOrgOxSO4.ratevartemplate,
        default_reminOrgOxCH4.ratevartemplate,
    ]

    PB.add_method_do!(
        rj, 
        do_remin_O2_NO3_Mn_Fe_SO4_CH4_rates, 
        (input_particulate_fluxes, output_solute_fluxes, PB.VarList_namedtuple(vcat(PB.VarDep.(ox_frac_vars), rate_vars))),
    )

    PB.add_method_do!(rj, default_reminOrgOxO2)
    PB.add_method_do!(rj, default_reminOrgOxNO3only)
    PB.add_method_do!(rj, default_reminOrgOxMnIVOx)
    PB.add_method_do!(rj, default_reminOrgOxFeIIIOx)
    PB.add_method_do!(rj, default_reminOrgOxSO4, isotope_data=SIsotopeType)
    PB.add_method_do!(rj, default_reminOrgOxCH4, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)
    return nothing
end

function do_remin_O2_NO3_Mn_Fe_SO4_CH4_ox_frac(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        fO2             = max(vars.O2_conc[i], 0.0)/(pars.oxreminlimit[] + max(vars.O2_conc[i], 0.0))
        fNO3            = max(vars.NO3_conc[i], 0.0)/(pars.NO3reminlimit[] + max(vars.NO3_conc[i], 0.0))
        fMnIVOx         = max(vars.MnIVOx_conc[i], 0.0)/(pars.MnIVOxreminlimit[] + max(vars.MnIVOx_conc[i], 0.0))
        fFeIIIOx        = max(vars.FeIIIOx_conc[i], 0.0)/(pars.FeIIIOxreminlimit[] + max(vars.FeIIIOx_conc[i], 0.0))
        fSO4            = max(vars.SO4_conc[i], 0.0)/(pars.SO4reminlimit[] + max(vars.SO4_conc[i], 0.0))

        freminOrgO2         = fO2
        freminOrgNO3        = (1.0 - fO2)*fNO3
        freminOrgMnIVOx     = (1.0 - fO2)*(1.0 - fNO3)*fMnIVOx
        freminOrgFeIIIOx    = (1.0 - fO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*fFeIIIOx
        freminOrgSO4        = pars.aSO4[]*(1.0 - fO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*(1.0 - fFeIIIOx)*fSO4        
        if pars.oxidant_dependent_rate[]
            freminOrgCH4    = pars.aSO4[]*(1.0 - fO2)*(1.0 - fNO3)*(1.0 - fMnIVOx)*(1.0 - fFeIIIOx)*(1.0 - fSO4)
            freminOrgTot    = freminOrgO2 + freminOrgNO3 + freminOrgMnIVOx + freminOrgFeIIIOx + freminOrgSO4 + freminOrgCH4
            vars.freminOrgTot[i]        = freminOrgTot
        else
            freminOrgCH4    = 1.0 - freminOrgO2 - freminOrgNO3 - freminOrgMnIVOx - freminOrgFeIIIOx - freminOrgSO4
            freminOrgTot    = 1.0
        end

        vars.freminOrgO2[i]         = freminOrgO2 / freminOrgTot
        vars.freminOrgNO3[i]        = freminOrgNO3 / freminOrgTot
        vars.freminOrgMnIVOx[i]     = freminOrgMnIVOx / freminOrgTot
        vars.freminOrgFeIIIOx[i]    = freminOrgFeIIIOx / freminOrgTot
        vars.freminOrgSO4[i]        = freminOrgSO4 / freminOrgTot
        vars.freminOrgCH4[i]        = freminOrgCH4 / freminOrgTot
        
    end

    return nothing
end

function do_remin_O2_NO3_Mn_Fe_SO4_CH4_rates(
    m::PB.ReactionMethod,
    pars,
    (input_particulate_fluxes, output_solute_fluxes, vars),
    cellrange::PB.AbstractCellRange,
    deltat
)
                     
    @inbounds for i in cellrange.indices
        O2eq = do_ReminParticulateFluxes(i, (input_particulate_fluxes, output_solute_fluxes), Ncycle=true) # remineralize N to ammonia
     
        # NB: stoichiometries already normalized to O2eq yr-1
        vars.reminOrgOxO2[i]        = O2eq*vars.freminOrgO2[i]
        vars.reminOrgOxNO3only[i]   = O2eq*vars.freminOrgNO3[i]
        vars.reminOrgOxMnIVOx[i]    = O2eq*vars.freminOrgMnIVOx[i]
        vars.reminOrgOxFeIIIOx[i]   = O2eq*vars.freminOrgFeIIIOx[i]
        vars.reminOrgOxSO4[i]       = O2eq*vars.freminOrgSO4[i]
        # CO2 -> CH4 + 2 O2eq, so each O2eq corresponds to 0.5 CH4 released
        vars.reminOrgOxCH4[i]       = O2eq*vars.freminOrgCH4[i]
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
        
    (O2eq, TAlk) = PALEOaqchem.O2AlkUptakeRemin(Corg, Nremin, P, Ccarb)

    PB.add_if_available(output_solute_fluxes.DIC, i, Corg + Ccarb)
    PB.add_if_available(output_solute_fluxes.TAlk, i, TAlk)

    if hasfield(typeof(input_particulate_fluxes), :Corg_delta)
        input_particulate_fluxes.Corg_delta[i] = PB.get_delta(Corg, eps())  # guard against 0/0
    end

    return O2eq

end



end # module
