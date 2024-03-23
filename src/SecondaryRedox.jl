module SecondaryRedox

import PALEOboxes as PB
using PALEOboxes.DocStrings



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
