module SecondaryRedox

import PALEOboxes as PB
using PALEOboxes.DocStrings


#####################################################################
# Nitrogen
#####################################################################

"""
    ReactionRedoxNH4_NO2

Annamox

    NH4+ + NO2- -> N2

Rate

    R_NH4_NO2 * NH4_conc * NO2_conc (mol N2 m-3 yr-1)

[Dale2015](@cite) 1e8 M-1 yr-1 (-> 1e5 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxNH4_NO2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_NH4_NO2", 1e8*1e-3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_NH4_NO2 = PB.RateStoich(
        PB.VarProp("redox_NH4_NO2", "mol N2 yr-1", "N2 production (+ve) by anammox",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "NH4"), (-1.0, "NO2"), (1.0, "N2")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxNH4_NO2) = register_redox_methods(rj, rj.stoich_redox_NH4_NO2, "NH4", "NO2", :R_NH4_NO2)


"""
    ReactionRedoxNH4_O2

Aerobic oxidation of ammonia

    NH4+ + 1.5 O2 -> NO2- + H2O + 2H+

Rate

    R_NH4_O2 * NH4_conc * O2_conc (mol N m-3 yr-1)

[Dale2015](@cite) 1e7 M-1 yr-1 (-> 1e4 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxNH4_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_NH4_O2", 1e7*1e-3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_NH4_O2 = PB.RateStoich(
        PB.VarProp("redox_NH4_O2", "mol N yr-1", "NO2 production (+ve) by aerobic oxidation of ammonia",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "NH4"), (-1.5, "O2"), (1.0, "NO2"), (-2.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxNH4_O2) = register_redox_methods(rj, rj.stoich_redox_NH4_O2, "NH4", "O2", :R_NH4_O2)


"""
    ReactionRedoxNO2_O2

Aerobic oxidation of NO2

    NO2- + 0.5 O2 -> NO3-

Rate

    R_NO2_O2 * NO2_conc * NO3_conc (mol N m-3 yr-1)

[Dale2015](@cite) 1e7 M-1 yr-1 (-> 1e4 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxNO2_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_NO2_O2", 1e7*1e-3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_NO2_O2 = PB.RateStoich(
        PB.VarProp("redox_NO2_O2", "mol N yr-1", "NO3 production (+ve) by aerobic oxidation of NO2",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "NO2"), (-0.5, "O2"), (1.0, "NO3"),), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxNO2_O2) = register_redox_methods(rj, rj.stoich_redox_NO2_O2, "NO2", "O2", :R_NO2_O2)


##################################################################################################################
# Manganese
################################################################################################################

"""
    ReactionRedoxMnII_O2

MnII oxidation by oxygen, producing generic MnIVOx (MnO2)

    Mn++ + 1/2 O2 + H2O -> MnIVOx + 2 H+

Rate

    R_MnII_O2 * MnII_conc * O2_conc (mol Mn m-3 yr-1)

[Dale2015](@cite) 5e6 M-1 yr-1 (-> 5e3 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxMnII_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # 
        PB.ParDouble("R_MnII_O2", 5e6/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_MnII_O2 = PB.RateStoich(
        PB.VarProp("redox_MnII_O2", "mol Mn yr-1", "Mn consumption by Mn++ oxidation",
                        attributes=(:calc_total=>true,)),
        ((-1/2, "O2"), (-1.0, "MnII"), (1.0, "MnIVOx"), (-2.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxMnII_O2) = register_redox_methods(rj, rj.stoich_redox_MnII_O2, "O2", "MnII", :R_MnII_O2)


"""
    ReactionRedoxMnIV_FeII

Generic MnIVOx (MnO2) reduction by FeII, producing generic FeIIIOx (0.5 Fe2O3)

    MnIVOx + 2 Fe++ + H2O -> Mn++ + 2 FeIIIOx + 2 H+

Rate

    R_MnIV_FeII * MnIVOx_conc * FeII_conc  (mol Mn m-3 yr-1)

[Dale2015](@cite) 1e7, 1e5 M-1 yr-1 for MnHR, MnMR (-> 1e4, 1e2 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxMnIV_FeII{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # 
        PB.ParDouble("R_MnIV_FeII", 1e7/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_MnIV_FeII = PB.RateStoich(
        PB.VarProp("redox_MnIV_FeII", "mol Mn yr-1", "Mn consumption (+ve) by Mn oxide reduction by Fe2++",
                        attributes=(:calc_total=>true,)),
        ((-1.0, "MnIVOx"), (-2.0, "FeII"), (2.0, "FeIIIOx"), (1.0, "MnII"), (-2.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxMnIV_FeII) = register_redox_methods(rj, rj.stoich_redox_MnIV_FeII, "MnIVOx", "FeII", :R_MnIV_FeII)


"""
    ReactionRedoxMnIV_H2S

Sulphate-mediated manganese reduction (of generic Mn oxide FeMnIVOx ~ MnO2)

    1/4 H2S + MnIVOx + 3/2 H+ -> Mn++ + 1/4 SO4-- + H2O

Rate:

    R_MnIV_H2S * H2S_conc * MnIVOx_conc  (mol Mn m-3 yr-1)

[Dale2015](@cite) 1e5, 1e3 M-1 yr-1 for MnHR, MnMR (-> 1e2, 1e0 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxMnIV_H2S{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
       
        PB.ParDouble("R_MnIV_H2S", 1e2 , units="(mol m-3)^-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_MnIV_H2S = PB.RateStoich(
        PB.VarProp("redox_MnIV_H2S", "mol Mn yr-1", "Mn consumption (+ve) by MnIV reduction",
                        attributes=(:calc_total=>true,)),
        ((-1/4, "H2S"), (-1.0, "MnIVOx"), (1.0, "MnII"), (1/4, "SO4"), (+3/2, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxMnIV_H2S) = register_redox_methods(rj, rj.stoich_redox_MnIV_H2S, "MnIVOx", "H2S", :R_MnIV_H2S)


##################################################################################################################
# Iron and iron sulphides
###############################################################################################################

"""
    ReactionRedoxFeII_O2

FeII oxidation by oxygen, producing generic FeIIIOx (0.5 Fe2O3)

    Fe++ + 1/4 O2 + H2O -> FeIIIOx + 2 H+

Rate

    R_FeII_O2 * FeII_conc * O2_conc (mol Fe m-3 yr-1)

[Dale2015](@cite) 5e8 M-1 yr-1 (-> 5e5 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeII_O2{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # 
        PB.ParDouble("R_FeII_O2", 5e8/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeII_O2 = PB.RateStoich(
        PB.VarProp("redox_FeII_O2", "mol Fe yr-1", "Fe consumption (+ve) by Fe++ oxidation",
                        attributes=(:calc_total=>true,)),
        ((-1/4, "O2"), (-1.0, "FeII"), (1.0, "FeIIIOx"), (-2.0, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxFeII_O2) = register_redox_methods(rj, rj.stoich_redox_FeII_O2, "O2", "FeII", :R_FeII_O2)



"""
    ReactionRedoxFeII_NO3

FeII oxidation by nitrate, producing generic FeIIIOx (0.5 Fe2O3)

    Fe++ + 2/10 NO3- + 9/10 H2O -> FeIIIOx + 1/10 N2 + 18/10 H+

Rate

    R_FeII_NO3 * FeII_conc * NO3_conc (mol Fe m-3 yr-1)

[Dale2015](@cite) 1e5 M-1 yr-1 (-> 1e2 (mol m-3)-1 yr-1)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRedoxFeII_NO3{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        # 
        PB.ParDouble("R_FeII_NO3", 1e5/1e3 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_redox_FeII_NO3 = PB.RateStoich(
        PB.VarProp("redox_FeII_NO3", "mol Fe yr-1", "Fe consumption (+ve) by Fe++ oxidation by NO3",
                        attributes=(:calc_total=>true,)),
        ((-2/10, "NO3"), (-1.0, "FeII"), (1.0, "FeIIIOx"), (-18/10, "TAlk")), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )

end

PB.register_methods!(rj::ReactionRedoxFeII_NO3) = register_redox_methods(rj, rj.stoich_redox_FeII_NO3, "NO3", "FeII", :R_FeII_NO3)


"""
    ReactionRedoxFeIII_H2S

Sulphate-mediated iron reduction (of generic iron oxide FeIIIOx ~ 0.5 Fe2O3)

    1/8 H2S + FeIIIOx + 14/8 H+ -> Fe++ + 1/8 SO4-- + H2O

Rate:

    R_FeIII_H2S * H2S_conc^0.5 * FeIIIOx_conc  (mol Fe m-3 yr-1)

NB: rate constant k [units (mol m-3)^-0.5 yr-1] = 277.2 * k [units (mol l-1)^-0.5 hr-1]
where van de Velde (2021) give k_SMIs = 1.98e0 (mol l-1)^-0.5 hr-1 for solid-phase FeIII

[Dale2015](@cite) 100, 0.1, 3e-4 M-0.5 yr-1 -> / sqrt(1e3) = 3.16, 3.16e-3, 9.49e06 (mol m-3)^-0.5 yr-1

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
        ((-1/8, "H2S"), (-1.0, "FeIIIOx"), (1.0, "FeII"), (1/8, "SO4"), (+14/8, "TAlk")), 
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

Rate:

    R_FeS_O2 * FeS_conc * O2_conc  (mol Fe m-3 yr-1)

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

PB.register_methods!(rj::ReactionRedoxFeS_O2) = register_redox_methods(rj, rj.stoich_redox_FeS_O2, "FeS", "O2", :R_FeS_O2)


"""
    ReactionRedoxFeS2pyr_O2

Pyrite oxidation

    FeS2 + 3.5 O2 + H2O -> Fe++ + 2 SO4-- + 2 H+

Rate:

    R_FeS2pyr_O2 * FeS2pyr_conc * O2_conc (mol Fe m-3 yr-1)

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

PB.register_methods!(rj::ReactionRedoxFeS2pyr_O2) = register_redox_methods(rj, rj.stoich_redox_FeS2pyr_O2, "FeS2pyr", "O2", :R_FeS2pyr_O2)


##############################################
# Sulphur
###############################################

"""
    ReactionRedoxH2S_O2

Sulphide oxidation by oxygen.

    O2 + 0.5 H2S -> 0.5 SO4 + H+

Rate:

    R_H2S_O2 * O2_conc * H2S_conc (mol O2 m-3 yr-1)

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
    
    register_redox_methods(
        rj, rj.stoich_redox_H2S_O2, "O2", "H2S", :R_H2S_O2; 
        stoich_redox_args=(isotope_data=SIsotopeType,),
        append_message="SIsotopeType=$(SIsotopeType)"
    )

    return nothing
end


###################################################################
# Methane
###################################################################

"""
    ReactionRedoxCH4_O2

Methane oxidation by oxygen.

    O2 + 0.5 CH4 -> 0.5 DIC (CO2) + H2O

Rate:

    R_CH4_O2 * O2_conc * CH4_conc (mol O2 m-3 yr-1)


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
    
    register_redox_methods(
        rj, rj.stoich_redox_CH4_O2, "O2", "CH4", :R_CH4_O2;
        stoich_redox_args=(isotope_data=CIsotopeType,),
        append_message="CIsotopeType=$(CIsotopeType)"
    )

    return nothing
end


"""
    ReactionRedoxCH4_SO4

Methane oxidation by sulphate (anaerobic methane oxidation).

    SO4-- + CH4 + 2 H+ -> DIC + H2S + 2H2O 

Rate

    R_CH4_SO4 * SO4_conc * CH4_conc (mol C m-3 yr-1)

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
    
    register_redox_methods(
        rj, rj.stoich_redox_CH4_SO4, "SO4", "CH4", :R_CH4_SO4;
        stoich_redox_args=(isotope_data=CIsotopeType,),
        append_message="CIsotopeType=$(CIsotopeType)"
    )

    return nothing
end


####################################################################
# Hydrogen
####################################################################

"""
    ReactionRedoxH2_O2

H2 oxidation by oxygen, producing water

    H2 + 0.5 O2 -> H2O

Rate

    R_H2_O2 * H2_conc * O2_conc (mol H2 m-3 yr-1)

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

PB.register_methods!(rj::ReactionRedoxH2_O2) = register_redox_methods(rj, rj.stoich_redox_H2_O2, "O2", "H2", :R_H2_O2)


"""
    ReactionRedoxH2_SO4

H2 oxidation by SO4, producing sulphide

    H2 + 1/4 SO4-- + 2/4 H+ -> 1/4 H2S + H2O

Rate

    R_H2_SO4 * H2_conc * SO4_conc (mol H2 m-3 yr-1)

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

PB.register_methods!(rj::ReactionRedoxH2_SO4) = register_redox_methods(rj, rj.stoich_redox_H2_SO4, "SO4", "H2", :R_H2_SO4)


######################################################################################################################
# generic functions
######################################################################################################################

function register_redox_methods(
    @nospecialize(rj::PB.AbstractReaction), stoich_redox::PB.RateStoich, R1name::String, R2name::String, Rname::Symbol;
    stoich_redox_args=NamedTuple(), append_message="", rate_function=do_redox_rate,
)
    @info "register_methods! $(PB.fullname(rj)) $append_message"

    ratevar = stoich_redox.ratevartemplate
    concvars = [
        PB.VarDep(R1name*"_conc", "mol m-3", R1name*" concentration"),
        PB.VarDep(R2name*"_conc", "mol m-3", R2name*" concentration"),
    ]
    volume = PB.VarDep("volume", "m3", "box fluid volume")

    PB.add_method_do!(
        rj,
        rate_function, 
        (
            PB.VarList_single(ratevar),
            PB.VarList_tuple(concvars),
            PB.VarList_single(volume), 
        );
        p=Val(Rname), 
        name="do_redox_"*R1name*"_"*R2name,
    )

    PB.add_method_do!(rj, stoich_redox; stoich_redox_args...)

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
