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
