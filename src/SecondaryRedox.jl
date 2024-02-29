module SecondaryRedox

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionRedoxH2S_O2

Sulphide oxidation by oxygen.

|Rate `R_H2S_O2`|Units            | Ref                        | Notes |
|:--------------|:----------------|:---------------------------|:------|
|1.6e5          |(mol l-1)-1 yr-1 |[VanCappellen1996a](@cite)  |       |
|3.65e6         |(mol l-1)-1 yr-1 |[Ozaki2011](@cite)          |       |
|54e6           |(mol l-1)-1 yr-1 |[Romaniello2010a](@cite)    |       |

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxH2S_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        #                                          (mol l-1)-1 yr-1 / l/m^3 -> (mol m-3)-1 yr-1
        PB.ParDouble("R_H2S_O2", 3.65e6/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
        PB.ParType(PB.AbstractData, "SIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable sulphur isotopes and specify isotope type"),
    )

    stoich_redox_H2S_O2 = PB.RateStoich(
        PB.VarProp("redox_H2S_O2", "mol O2 yr-1", "oxygen consumption (+ve) by H2S oxidation",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "O2"), (-0.5, "H2S::Isotope"), (0.5, "SO4::Isotope"), (-1.0, "TAlk")), 
        deltavarname_eta = ("H2S_delta", 0.0),  # no fractionation
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxH2S_O2)

    SIsotopeType = rj.pars.SIsotope[]
    PB.setfrozen!(rj.pars.SIsotope)
    @info "register_methods! $(PB.fullname(rj)) SIsotopeType=$(SIsotopeType)"

    ratevar = rj.stoich_redox_H2S_O2.ratevartemplate
    concvars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("H2S_conc", "mol m-3", "H2S concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_H2S_O2),
        name="do_redox_H2S_O2"
    )

    PB.add_method_do!(rj, rj.stoich_redox_H2S_O2, isotope_data=SIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxCH4_O2

Methane oxidation by oxygen.

|Rate `R_CH4_O2`|Units            | Ref                        | Notes |
|:--------------|:----------------|:---------------------------|:------|
|1e10           |(mol l-1)-1 yr-1 |[VanCappellen1996a](@cite)  |       |
|10e6           |(mol l-1)-1 yr-1 |TODO                        |       |

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxCH4_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        #                                         (mol l-1)-1 yr-1 / l/m^3 -> (mol m-3)-1 yr-1
        PB.ParDouble("R_CH4_O2", 10e6/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    stoich_redox_CH4_O2 = PB.RateStoich(
        PB.VarProp("redox_CH4_O2", "mol O2 yr-1", "oxygen consumption (+ve) by CH4 oxidation",
                attributes=(:calc_total=>true,)),
        ((-1.0, "O2"), (-0.5, "CH4::Isotope"), (0.5, "DIC::Isotope")), 
        deltavarname_eta = ("CH4_delta", 0.0),  # no fractionation
        sms_prefix="",
        sms_suffix="_sms",
        processname="redox"
    )
end

function PB.register_methods!(rj::ReactionRedoxCH4_O2)

    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType)"

    ratevar = rj.stoich_redox_CH4_O2.ratevartemplate
    concvars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("CH4_conc", "mol m-3", "CH4 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")
    
    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_CH4_O2),
        name="do_redox_CH4_O2"
    )

    PB.add_method_do!(rj, rj.stoich_redox_CH4_O2, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxCH4_SO4

Methane oxidation by sulphate (anaerobic methane oxidation).

|Rate `R_CH4_SO4`|Units           | Ref                        | Notes |
|:--------------|:----------------|:---------------------------|:------|
|1e4            |(mol l-1)-1 yr-1 |[VanCappellen1996a](@cite)  |   k17 |

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxCH4_SO4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        #                                           (mol l-1)-1 yr-1 / l/m^3 -> (mol m-3)-1 yr-1
        PB.ParDouble("R_CH4_SO4", 1e4/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),

        PB.ParType(PB.AbstractData, "CIsotope", PB.ScalarData,
            external=true,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable carbon isotopes and specify isotope type"),
    )

    stoich_redox_CH4_SO4 = PB.RateStoich(
        PB.VarProp("redox_CH4_SO4", "mol CH4/SO4 yr-1", "sulphate consumption (+ve) by CH4 oxidation",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "SO4"), (-1.0, "CH4::Isotope"), (1.0, "DIC::Isotope"),(1.0, "H2S"), (2.0, "TAlk")), 
        deltavarname_eta = ("CH4_delta", 0.0),  # no fractionation
        sms_prefix="",
        sms_suffix="_sms",
        processname="redox"
    )
end

function PB.register_methods!(rj::ReactionRedoxCH4_SO4)

    CIsotopeType = rj.pars.CIsotope[]
    PB.setfrozen!(rj.pars.CIsotope)
    @info "register_methods! $(PB.fullname(rj)) CIsotopeType=$(CIsotopeType)"

    ratevar = rj.stoich_redox_CH4_SO4.ratevartemplate
    concvars = [
        PB.VarDep("SO4_conc", "mol m-3", "SO4 concentration"),
        PB.VarDep("CH4_conc", "mol m-3", "CH4 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_CH4_SO4), 
        name="do_redox_CH4_SO4"
    )

    PB.add_method_do!(rj, rj.stoich_redox_CH4_SO4, isotope_data=CIsotopeType)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxH2_O2

H2 oxidation by oxygen, producing water

    H2 + 0.5 O2 -> H2O

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxH2_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        #                                          (mol l-1)-1 yr-1 / l/m^3 -> (mol m-3)-1 yr-1
        PB.ParDouble("R_H2_O2", 1e6 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_H2_O2 = PB.RateStoich(
        PB.VarProp("redox_H2_O2", "mol H2 yr-1", "H2 consumption (+ve) by H2 oxidation by O2",
                        attributes=(:calc_total=>true,)),
        ((-0.5, "O2"), (-1.0, "H2")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxH2_O2)

    @info "register_methods! $(PB.fullname(rj))"

    ratevar = rj.stoich_redox_H2_O2.ratevartemplate
    concvars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("H2_conc", "mol m-3", "H2 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_H2_O2), 
        name="do_redox_H2_O2",
    )

    PB.add_method_do!(rj, rj.stoich_redox_H2_O2)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxH2_SO4

H2 oxidation by SO4, producing sulphide

     H2 + 1/4 SO4-- + 2/4 H+ -> 1/4 H2S + H2O 

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxH2_SO4{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        #                                          (mol l-1)-1 yr-1 / l/m^3 -> (mol m-3)-1 yr-1
        PB.ParDouble("R_H2_SO4", 1e6 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_H2_SO4 = PB.RateStoich(
        PB.VarProp("redox_H2_SO4", "mol H2 yr-1", "H2 consumption (+ve) by H2 oxidation by SO4",
                        attributes=(:calc_total=>true,)),
        ((-0.25, "SO4"), (-1.0, "H2"), (0.25, "H2S"), (+0.5, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxH2_SO4)

    @info "register_methods! $(PB.fullname(rj))"

    ratevar = rj.stoich_redox_H2_SO4.ratevartemplate
    concvars = [
        PB.VarDep("SO4_conc", "mol m-3", "SO4 concentration"),
        PB.VarDep("H2_conc", "mol m-3", "H2 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_H2_SO4), 
        name="do_redox_H2_SO4",
    )

    PB.add_method_do!(rj, rj.stoich_redox_H2_SO4)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


##################################################################################################################
# Iron 
###############################################################################################################

"""
    ReactionRedoxFeII_O2

FeII oxidation by oxygen, producing generic FeIIIOx (0.5 Fe2O3)

    4 Fe++ + O2 + 4 H2O -> 4 FeIIIOx + 8 H+

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeII_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # 
        PB.ParDouble("R_FeII_O2", 3.65e6/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeII_O2 = PB.RateStoich(
        PB.VarProp("redox_FeII_O2", "mol O2 yr-1", "oxygen consumption (+ve) by Fe++ oxidation",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "O2"), (-4.0, "FeII"), (4.0, "FeIIIOx"), (-8.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxFeII_O2)

    @info "register_methods! $(PB.fullname(rj))"

    ratevar = rj.stoich_redox_FeII_O2.ratevartemplate
    concvars = [
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
        PB.VarDep("FeII_conc", "mol m-3", "FeII concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_FeII_O2), 
        name="do_redox_FeII_O2",
    )

    PB.add_method_do!(rj, rj.stoich_redox_FeII_O2)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxFeIII_H2S

Sulphate-mediated iron reduction (of generic iron oxide FeIIIOx ~ 0.5 Fe2O3)

    H2S + 8 FeIIIOx + 14 H+ -> 8 Fe++ + SO4-- + 8 H2O

NB: rate constant k [units (mol m-3)^-0.5 yr-1] = 277.2 * k [units (mol l-1)^-0.5 hr-1]

where van de Velde (2021) give k_SMIs = 1.98e0 (mol l-1)^-0.5 hr-1 for solid-phase FeIII

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeIII_H2S{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
       
        PB.ParDouble("R_FeIII_H2S", 277.2*1.98e0 , units="(mol m-3)^-0.5 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeIII_H2S = PB.RateStoich(
        PB.VarProp("redox_FeIII_H2S", "mol H2S yr-1", "sulphide consumption (+ve) by Fe+++ reduction",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "H2S"), (-8.0, "FeIIIOx"), (8.0, "FeII"), (1.0, "SO4"), (+14.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxFeIII_H2S)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        rj.stoich_redox_FeIII_H2S.ratevartemplate,
        PB.VarDep("H2S_conc", "mol m-3", "H2S concentration"),
        PB.VarDep("FeIIIOx_conc", "mol m-3", "FeIIIOx concentration"),
        PB.VarDep("volume", "m3", "box fluid volume"),
    ]

    PB.add_method_do!(
        rj,
        do_FeIII_H2S_rate, 
        (
            PB.VarList_namedtuple(vars),
        );
    )

    PB.add_method_do!(rj, rj.stoich_redox_FeIII_H2S)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


function do_FeIII_H2S_rate(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
   
    rateparval = pars.R_FeIII_H2S[]

    for i in cellrange.indices
        # mol yr-1                                                        
        vars.redox_FeIII_H2S[i] = (  rateparval                             #  (mol m-3)^-0.5 yr-1 
                        * sqrt(max(PB.get_total(vars.H2S_conc[i]), 0.0))    #   (mol m-3)^0.5
                        * max(PB.get_total(vars.FeIIIOx_conc[i]), 0.0)       #   mol m-3  
                        * vars.volume[i])                                         #   m3
    end
   
    return nothing
end


"""
    ReactionRedoxFeS_O2

FeS oxidation

    FeS + 2 O2 -> Fe++ + SO4--

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeS_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_FeS_O2", 1e5*1e-3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeS_O2 = PB.RateStoich(
        PB.VarProp("redox_FeS_O2", "mol FeS yr-1", "FeS consumption (+ve) by oxidation by O2",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "FeS"), (-2.0, "O2"), (1.0, "FeII"), (1.0, "SO4"),), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxFeS_O2)

    @info "register_methods! $(PB.fullname(rj))"

    ratevar = rj.stoich_redox_FeS_O2.ratevartemplate
    concvars = [
        PB.VarDep("FeS_conc", "mol m-3", "FeS concentration"),
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_FeS_O2), 
        name="do_redox_FeS_O2",
    )

    PB.add_method_do!(rj, rj.stoich_redox_FeS_O2)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


"""
    ReactionRedoxFeS2pyr_O2

Pyrite oxidation

    FeS2 + 3.5 O2 + H2O -> Fe++ + 2 SO4-- + 2 H+

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeS2pyr_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_FeS2pyr_O2", 1e3*1e-3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeS2pyr_O2 = PB.RateStoich(
        PB.VarProp("redox_FeS2pyr_O2", "mol FeS2pyr yr-1", "pyrite consumption (+ve) by oxidation by O2",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "FeS2pyr"), (-3.5, "O2"), (1.0, "FeII"), (2.0, "SO4"), (-2.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

function PB.register_methods!(rj::ReactionRedoxFeS2pyr_O2)

    @info "register_methods! $(PB.fullname(rj))"

    ratevar = rj.stoich_redox_FeS2pyr_O2.ratevartemplate
    concvars = [
        PB.VarDep("FeS2pyr_conc", "mol m-3", "FeS2pyr (pyrite) concentration"),
        PB.VarDep("O2_conc", "mol m-3", "O2 concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        do_redox_rate, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(:R_FeS2pyr_O2), 
        name="do_redox_FeS2pyr_O2",
    )

    PB.add_method_do!(rj, rj.stoich_redox_FeS2pyr_O2)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


######################################################################################################################
# generic functions
######################################################################################################################

# calculate generic first-order rate for two concentration defined by concvars[1:2]
function do_redox_rate(m::PB.ReactionMethod, pars, (ratevar, concvars, volume), cellrange::PB.AbstractCellRange, deltat)
    rateparnameval = m.p
    # rateparval = getfield(pars, rateparname)[]

    getpar(pars, ::Val{PN}) where PN = getfield(pars, PN)

    rateparval = getpar(pars, rateparnameval)[]

    @inbounds for i in cellrange.indices
        # mol yr-1                                                        
        ratevar[i] = (  rateparval                                  #  (mol m-3)-1 yr-1 
                        * max(PB.get_total(concvars[1][i]), 0.0)    #   mol m-3 
                        * max(PB.get_total(concvars[2][i]), 0.0)    #   mol m-3  
                        *volume[i])                                 #   m3
    end
   
    return nothing
end


end
