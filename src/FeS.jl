module FeS

import PALEOboxes as PB
using PALEOboxes.DocStrings

"""
    ReactionFeSaq

Equilibrium chemistry of FeII, H2S, FeS system

Represents two equilibrium reactions as algebraic constraints to
define primary species [Fe2++] and [HS-] in terms of:
- total [S-II] = [HS-] + [H2S] + [FeSaq]
- total [FeII] = [Fe++] + [FeSaq] 
given a provided fixed value of [H+]

    H2S        <--> H+ + HS-         eqb. const. K_H2S
    FeSaq + H+ <--> Fe++ + HS-       eqb. const K2_FeSaq


Species concentration       |  PALEO variable name    | PALEO variable type
----------------------------|-------------------------|---------------------
[H+]                        | Hp_conc                 | Dependency 
total [S-II]                | SmIIaqtot_conc          | Dependency 
[HS-]                       | Hsm_conc                | State
constraint on total [S-II]  | SmIIaqtot_constraint    | Constraint
[H2S]                       | H2Ssp_conc              | Property
[HS-] + [H2S]               | H2Stot_conc             | Property
total [FeII]                | FeIIaqtot_conc          | Dependency
[Fe2++]                     | FeII_conc               | State
constraint on total [FeII]  | FeIIaqtot_constraint    | Constraint


(usually SmIIaqtot and FeIIaqtot are ODE variables and should be defined by ReactionReservoir)

See:
- Rickard (2006) GCA https://10.1016/j.gca.2006.02.029
- Rickard and Luther (2007) Chemical Reviews https://dx.doi.org/10.1021/cr0503658
- Lenton and Daines (2017) Ann. Rev. Mar. Sci. https://dx.doi.org/10.1146/annurev-marine-010816-060521
- van de Velde etal (2021) GMD https://10.5194/gmd-14-2713-2021

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionFeSaq{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParDouble("K_H2S", 2.15e-7 , units="mol kg-1",
            description="K_H2S = [H+][HS-]/[H2S]  units mol kg-1 eqb const for H2S <--> HS- + H+ (Hofmann (2010) Aquatic Geochemistry https://10.1007/s10498-009-9084-1)"),
        PB.ParDouble("K2_FeSaq", 10.0^2.2 , units="",
            description="K2_FeSaq = {Fe++}[HS-]/([FeSaq][H+]) eqb constant for FeSaq + H+ <--> Fe++ + HS- (Rickard (2006) eq 8)"),
        PB.ParDouble("activity_Fe2p", 0.23, units="",
            description="Fe++ activity (eg van de Velde 2012 Table 3 from Davies eqn)")

    )


end


# K_H2S from Hofmann (2010) Aquatic Geochemistry https://10.1007/s10498-009-9084-1
# Sal = 35
# TempK = 273.15 + 15
# lnKH2S= (225.838 + 0.3449*sqrt(Sal) - 0.0274*Sal) - 13275.3/TempK - 34.6435*log(TempK)
# -15.35415072896447
# exp(lnKH2S)
# 2.146728213653719e-7 # mol kg-1

function PB.register_methods!(rj::ReactionFeSaq)

    @info "register_methods! $(PB.fullname(rj))"

    setup_vars = [
        PB.VarState("HSm_conc", "mol m-3", "aqueous HS- concentration"),
        PB.VarState("FeII_conc", "mol m-3", "aqueous Fe++ concentration"),        
        PB.VarConstraint("SmIIaqtot_constraint", "mol m-3", "algebraic constraint on SmIIaqtot_conc (= 0)"),
        PB.VarConstraint("FeIIaqtot_constraint", "mol m-3", "algebraic constraint on FeIIaqtot_conc (= 0)"),
    ]

    vars_FeSaq_eqb = vcat(
        setup_vars, 
        [
            PB.VarDep("Hp_conc", "mol m-3", "aqueous H+ concentration"),
            PB.VarDep("rho_ref", "kg m-3", "reference density"),
            PB.VarDep("SmIIaqtot_conc", "mol m-3", "total aqueous S-II concentration"),
            PB.VarDep("FeIIaqtot_conc", "mol m-3", "total aqueous FeII concentration"),
            PB.VarProp("H2Ssp_conc", "mol m-3", "aqueous H2S species concentration"),
            PB.VarProp("H2Stot_conc", "mol m-3", "aqueous H2S + HS- concentration"),
            PB.VarProp("FeSaq_conc", "mol m-3", "aqueous FeS concentration"),            
        ]
    )

    PB.add_method_setup!(
        rj, 
        setup_FeSaq_eqb,
        (PB.VarList_namedtuple(setup_vars), ),
    )

    PB.add_method_do!(
        rj,
        do_FeSaq_eqb, 
        (
            PB.VarList_namedtuple(vars_FeSaq_eqb),
        );
    )

    return nothing
end

function setup_FeSaq_eqb(
    m::PB.ReactionMethod,
    pars,
    (setup_vars, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name in (:initial_value, :norm_value) || return
    
    norm_value = 1e-3 # for solver scaling only

    # constraint variables just need norm_value setting
    if attribute_name == :norm_value           
        setup_vars.SmIIaqtot_constraint[cellrange.indices] .= norm_value
        @info "ReactionFeSaq initialize_state:$(rpad(attribute_name,20)) $(rpad("SmIIaqtot_constraint", 30)) = $norm_value"
        setup_vars.FeIIaqtot_constraint[cellrange.indices] .= norm_value
        @info "ReactionFeSaq initialize_state:$(rpad(attribute_name,20)) $(rpad("FeIIaqtot_constraint", 30)) = $norm_value"
    end

    # state variables need both norm_value and initial_value setting
    initialnorm_value = Dict(:initial_value=>1e-9, :norm_value=>norm_value)[attribute_name]
    setup_vars.HSm_conc[cellrange.indices] .= initialnorm_value
    @info "ReactionFeSaq initialize_state:$(rpad(attribute_name,20)) $(rpad("HSm_conc", 30)) = $initialnorm_value"
    setup_vars.FeII_conc[cellrange.indices] .= initialnorm_value
    @info "ReactionFeSaq initialize_state:$(rpad(attribute_name,20)) $(rpad("FeII_conc", 30)) = $initialnorm_value"

    return nothing
end

function do_FeSaq_eqb(m::PB.ReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)
   
    K_H2S = pars.K_H2S[] # K_H2S = [H+][HS-]/[H2S]  units mol kg-1
    K2_FeSaq = pars.K2_FeSaq[]
    activity_Fe2p = pars.activity_Fe2p[]

    for i in cellrange.indices
        # H2S <--> HS- + H+ eqb
        #   mol m-3        =  mol m-3         *  mol m-3         / (mol kg-1 * kg m-3)
        vars.H2Ssp_conc[i] = vars.HSm_conc[i] * vars.Hp_conc[i] / (K_H2S * vars.rho_ref[i])
        vars.H2Stot_conc[i] = vars.HSm_conc[i] + vars.H2Ssp_conc[i]

        # FeSaq + H+ <--> Fe++ + HS- eqb              
        vars.FeSaq_conc[i] = (activity_Fe2p * vars.FeII_conc[i] * vars.HSm_conc[i] / vars.Hp_conc[i]) / K2_FeSaq

        # H2S budget constraint for DAE solver
        vars.SmIIaqtot_constraint[i] = vars.SmIIaqtot_conc[i] - (vars.H2Stot_conc[i] + vars.FeSaq_conc[i])

        # FeII budget constraint for DAE solver
        vars.FeIIaqtot_constraint[i] = vars.FeIIaqtot_conc[i] - (vars.FeII_conc[i] + vars.FeSaq_conc[i])       
    end
   
    return nothing
end


"""
    ReactionFeSm

FeS precipitation/dissolution

Represents 

    [FeSaq] <--> [FeSm]

as a fast precipitation-dissolution reaction, using the formulation from Van Cappellen & Wang (1996)

See:
- Van Cappellen & Wang (1996) AJS https://dx.doi.org/10.2475/ajs.296.3.197
- Rickard (2006) GCA https://dx.doi.org/10.1016/j.gca.2006.02.029
- Rickard and Luther (2007) Chemical Reviews https://dx.doi.org/10.1021/cr0503658

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionFeSm{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("K0_FeSm", 10.0^-5.7 , units="mol l-1",
            description="solubility constant for FeSm"),
        PB.ParDouble("rate_FeS_prec", 40.0 , units="mol m-3 yr-1",
            description="rate of FeS precipitation"),
        PB.ParDouble("rate_FeS_diss", 1e3 , units="yr-1",
            description="rate of FeS dissolution"),
    )

    stoich_FeS_precdiss = PB.RateStoich(
        PB.VarProp("rate_FeSm", "mol yr-1", "FeSm precipitation - dissolution rate"; 
            attributes=(:calc_total=>true,)),
        ((-1.0, "SmIIaqtot"), (-1.0, "FeIIaqtot"), (1.0, "FeSm"), ), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="mineral",
    )
end


function PB.register_methods!(rj::ReactionFeSm)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        rj.stoich_FeS_precdiss.ratevartemplate,
        PB.VarDep("volume_solute", "m3", "solute fluid volume"),
        PB.VarDep("volume_solid", "m3", "solid-phase volume"),
        PB.VarDep("FeSaq_conc", "mol m-3", "aqueous FeS concentration"),
        PB.VarProp("Omega_FeSaq", "", "aqueous FeS saturation"),
        PB.VarDep("FeSm_conc", "mol m-3", "FeS solid phase concentration"),
    ]

    PB.add_method_do!(
        rj,
        do_FeSm, 
        (
            PB.VarList_namedtuple(vars),
        );
    )

    PB.add_method_do!(rj, rj.stoich_FeS_precdiss)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_FeSm(m::PB.ReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)
    # FeSm saturation and precipitation / dissolution

    K0_FeSm = pars.K0_FeSm[]

    for i in cellrange.indices    
        # FeS saturation
        # mol m-3    = mol l-1 * l m-3
        FeS_sat_conc = K0_FeSm * 1000
        vars.Omega_FeSaq[i] = vars.FeSaq_conc[i] / FeS_sat_conc

        # FeSm precipitation and dissolution
        # NB: bodge so AD sparsity detection sees a dependency on both FeSaq_conc and FeSm_conc
        # set rate to zero first with a fake dependency on both vars, then add/subtract to set value
        vars.rate_FeSm[i] = 0*vars.FeSm_conc[i] + 0*vars.FeSaq_conc[i]
        
        if vars.Omega_FeSaq[i] >= 1.0 # FeSm precipitation
            # Van Cappellen & Wang (1996) formulation
            # mol yr-1      =            mol m-3 yr-1                           * m3        
            vars.rate_FeSm[i] += pars.rate_FeS_prec[]*(vars.Omega_FeSaq[i] - 1) * vars.volume_solute[i]          
        else # FeSm dissolution
            # mol yr-1      =  yr-1          *  mol m-3                          * m3   
            # Van Cappellen & Wang (1996) formulation
            vars.rate_FeSm[i] -= pars.rate_FeS_diss[]*(1 - vars.Omega_FeSaq[i])*max(vars.FeSm_conc[i] - 1e-12, 0.0)*vars.volume_solid[i] # limit to minimum conc of 1e-12            
        end

    end
   
    return nothing
end


"""
    ReactionPyrH2S

Pyrite formation by the "H2S" (Berzelius) mechanism

    FeSm + H2S -> FeS2pyr + H2

with rate

    R_Pyr_H2S * [FeSm] * [H2S]

Where R_Pyr_H2S is highly uncertain:

[Dale2015](@cite) use R_Pyr_H2S = 1e5 M-1 yr-1 (= 1e5 *1e-3 = 100 (mol m-3-1) yr-1)

[van_de_velde_iron_2021](@cite) use R_Pyr_H2S = 0.3708 M-1 h-1 (= 0.3708 * 8.766 (mol m-3)-1 yr-1),
where [H2S] is taken to be total H2S (HS- + H2S)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionPyrH2S{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("R_Pyr_H2S", 0.3708e0 * 8.766 , units="(mol m-3)-1 yr-1",
            description="rate constant"),
    )

    stoich_FeS2pyr = PB.RateStoich(
        PB.VarProp("rate_FeS2pyr_H2S", "mol yr-1", "pyrite formation rate via the H2S mechanism"; 
            attributes=(:calc_total=>true,)),
        ((-1.0, "FeSm"), (-1.0, "SmIIaqtot"), (1.0, "FeS2pyr"), (1.0, "H2"), ), 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="mineral",
    )
end

function PB.register_methods!(rj::ReactionPyrH2S)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        rj.stoich_FeS2pyr.ratevartemplate,
        PB.VarDep("FeSm_conc", "mol m-3", "FeSm concentration"),
        PB.VarDep("H2Ssp_conc", "mol m-3", "H2S (species) concentration"),
        PB.VarDep("volume", "m3", "box fluid volume"),
    ]

    PB.add_method_do!(rj, do_PyrH2S, (PB.VarList_namedtuple(vars), ))

    PB.add_method_do!(rj, rj.stoich_FeS2pyr)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end


function do_PyrH2S(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
   
    for i in cellrange.indices
        # mol yr-1                                                        
        vars.rate_FeS2pyr_H2S[i] = ( pars.R_Pyr_H2S[]           #  (mol m-3)-1 yr-1 
                * max(PB.get_total(vars.FeSm_conc[i]), 0.0)     #   mol m-3 
                * max(PB.get_total(vars.H2Ssp_conc[i]), 0.0)    #   mol m-3  
                * vars.volume[i])                               #   m3
    end
   
    return nothing
end


end
