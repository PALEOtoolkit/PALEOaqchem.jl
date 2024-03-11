module GenericReactions

import PALEOboxes as PB
using PALEOboxes.DocStrings

import ..PALEOaqchem

import Infiltrator

"""
    ReactionAqEqb

Define a new equilibrium species or omega N  

    N + a A <--> b B + c C

    [N] = (1/K_eqb) ([B]^b [C]^c) / (A^a)


# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAqEqb{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParStringVec("Reactants", ["H2S_conc"];
            description="concentrations or activities of new species followed by other reactants, write powers as X^2 etc"),
        PB.ParStringVec("Products", [ "HSm_conc", "Hp_conc"];
            description="concentrations or activities of products, write powers as X^2 etc"),
        PB.ParDouble("K_eqb", 1.814e-7*1027, units="mol m-3 or mol kg-1",
            description="equilibrium constant"),
        PB.ParDouble("K_density_power", 0.0;
            description="multiple K_eqb * rho_ref^K_density_power to convert units: 0.0 for K_eqb in mol m-3, 1.0 for K_eqb in mol kg-1, etc"),
        
        PB.ParStringVec("N_components", ["SmIIaqtot_calc_conc"] ,
            description="contribution of new species to element or component total concentrations, empty vector to just define an Omega"),
    )

    component_stoichs::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionAqEqb)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        PB.VarDep("volume", "m3", "cell solute volume"),
    ]
    K_density_power = rj.pars.K_density_power[]
    if !iszero(K_density_power)
        push!(vars, PB.VarDep("rho_ref", "kg m-3", "reference density"))
    end
    # convert to integer if possible, as integer powers are much faster
    K_density_power = isinteger(K_density_power) ? Int64(K_density_power) : K_density_power

    N_conc_name = first(rj.pars.Reactants.v)
    N_conc_var = PB.VarProp(N_conc_name, "mol m-3", "aqueous concentration or activity"; attributes=PALEOaqchem.R_conc_attributes)

    other_reactant_conc_vars = PB.VariableReaction[]
    other_reactant_powers = []
    for or_name_power in rj.pars.Reactants.v[2:end]
        orp, orn = PALEOaqchem.parse_name_to_power_number(or_name_power)
        push!(other_reactant_conc_vars, PB.VarDep(orn, "mol m-3", "aqueous concentration or activity"))
        push!(other_reactant_powers, orp)
    end

    product_conc_vars = PB.VariableReaction[]
    product_powers = []
    for pr_name_power in rj.pars.Products
        pp, pn = PALEOaqchem.parse_name_to_power_number(pr_name_power)
        push!(product_conc_vars, PB.VarDep(pn, "mol m-3", "aqueous concentration or activity"))
        push!(product_powers, pp)
    end

    empty!(rj.component_stoichs)
    component_vars = PB.VariableReaction[]
    for cv in rj.pars.N_components
        stoich, name = PALEOaqchem.parse_number_name(cv)
        push!(rj.component_stoichs, stoich)
        push!(component_vars, PB.VarContrib(name, "mol", "total moles"))
    end

    PB.add_method_do!(
        rj,
        do_aqeqb, 
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_single(N_conc_var),
            PB.VarList_tuple(other_reactant_conc_vars),
            PB.VarList_tuple(product_conc_vars),
            PB.VarList_tuple(component_vars),
        );
        p = (Tuple(other_reactant_powers), Tuple(product_powers), K_density_power)
    )

    return nothing
end

function do_aqeqb(m::PB.ReactionMethod, pars, (vars, N_conc, other_reactant_concs, product_concs, components), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    other_reactant_powers, product_powers, K_density_power = m.p

    for i in cellrange.indices
        sc_denom = mapreduce(x->max(PB.get_total(x[1][i]), 0.0)^x[2], *, zip(other_reactant_concs, other_reactant_powers); init=1.0)
        sc_num = mapreduce(x->max(PB.get_total(x[1][i]), 0.0)^x[2], *, zip(product_concs, product_powers); init=1.0)

        Keqb = pars.K_eqb[]
        if !iszero(K_density_power)
            Keqb *= vars.rho_ref[i]^K_density_power
        end

        N_conc[i] = sc_num / (sc_denom * Keqb)

        PB.IteratorUtils.foreach_longtuple(components, rj.component_stoichs) do c, s
            c[i] += s * N_conc[i] * vars.volume[i]
        end
    end
   
    return nothing
end


"""
    ReactionAqKinetic

Define a kinetic reaction with rate dependent on concentrations

    a A + b B --> c C + d D

Rate (default) is:

    R = K * [A] * [B]

where this can be modified to different power-law dependency.

Parameters `Reactants` and `Products` should be the vectors of stoichiometry * name of (total) species to accumulate fluxes into `_sms` variables.

`Reactant_conc` should be an empty vector to take concentration names from `Reactants`, or a Vector of name[^power] to calculate rate with explicit
species names (eg where `Reactants` refers to totals which are partioned into multiple species by equilibrium reactions) or non-first-order rate dependence.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAqKinetic{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParStringVec("Reactants", ["A", "2*B"];
            description="reactant species"),
        PB.ParStringVec("Products", [ "2*C", "D"];
            description="product species"),
        PB.ParStringVec("Reactant_concs", String[],
            description="names of concentration variables to calculate rate eg '[`\"A_conc^0.5\"]' etc, empty vector to used defaults from 'Reactants' eg 'A_conc', 'B_conc' ..."),
        PB.ParDouble("K", NaN;
            description="rate constant"),
    )
end


function PB.register_methods!(rj::ReactionAqKinetic)

    @info "register_methods! $(PB.fullname(rj))"

    reactant_names = String[]
    all_stoich = Float64[]
    for rns in rj.pars.Reactants
        stoich, name = PALEOaqchem.parse_number_name(rns)
        push!(all_stoich, -1*stoich) # -1 for reactant 
        push!(reactant_names, String(name))
    end
    product_names = String[]
    for pns in rj.pars.Products
        stoich, name = PALEOaqchem.parse_number_name(pns)
        push!(all_stoich, stoich)
        push!(product_names, String(name))
    end

    all_names = vcat(reactant_names, product_names)
    all_stoich_names = Tuple((s, n) for (s,n) in zip(all_stoich, all_names))
    ratevarname = replace("redox_"*join(all_names, "_"), "_TAlk"=>"")
    stoich_redox = PB.RateStoich(
        PB.VarProp(ratevarname, "mol yr-1", "rate variable",
                        attributes=(:calc_total=>true,)),
        all_stoich_names, 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="redox",
    )
    
    reactant_concs = isempty(rj.pars.Reactant_concs.v) ? reactant_names .* "_conc" : rj.pars.Reactant_concs.v 
    reactant_conc_vars = PB.VariableReaction[]
    reactant_powers = []
    for rnp in reactant_concs
        rp, rn = PALEOaqchem.parse_name_to_power_number(rnp)
        push!(reactant_conc_vars, PB.VarDep(rn, "mol m-3", "aqueous concentration or activity"))
        push!(reactant_powers, rp)
    end

    volume = PB.VarDep("volume", "m3", "cell solute volume")
    PB.add_method_do!(
        rj,
        do_aqkinetic, 
        (
            PB.VarList_single(stoich_redox.ratevartemplate),
            PB.VarList_tuple(reactant_conc_vars),
            PB.VarList_single(volume),
        );
        p=Tuple(reactant_powers)
    )

    PB.add_method_do!(rj, stoich_redox)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    

    return nothing
end

function do_aqkinetic(m::PB.ReactionMethod, pars, (ratevar, reactant_concs, volume), cellrange::PB.AbstractCellRange, deltat)
    reactant_powers = m.p

    rateparval = pars.K[]

    for i in cellrange.indices
        ratevar[i] = mapreduce(x->max(PB.get_total(x[2][i]), 0.0)^x[1], *, zip(reactant_powers, reactant_concs); init=rateparval*volume[i])
        # ratevar[i] = max(PB.get_total(reactant_concs[1][i]), 0.0) * max(PB.get_total(reactant_concs[2][i]), 0.0)*rateparval*volume[i]
    end
   
    return nothing
end


"""
    ReactionAqPrecipDissol

Define a precipitation and dissolution reaction for solid S

    a A + b B <--> s S + d D

Rate for the precipitation and dissolution reactions are:

    R_precip = K_precip * (Ω - 1)               (Ω > 1)
    R_dissol = K_dissol * S_conc * (1 - Ω)      (Ω < 1)

Parameters `Reactants` and `Products` should be the vectors of stoichiometry * name of (total) species to accumulate fluxes into `_sms` variables.

`Solid_conc` should be the name of the concentration variable for S, or an empty string to use default 'S_conc'.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAqPrecipDissol{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParStringVec("Reactants", ["A", "2*B"];
            description="reactant species"),
        PB.ParStringVec("Products", [ "S", "0.5*D"];
            description="product species (solid S first)"),
        PB.ParString("Solid_conc", "",
            description="name of solid S concentration variable (empty string to default to 'S_conc')"),
        PB.ParDouble("K_precip", 0.0; units="mol m-3 yr-1",
            description="rate constant for precipitation reaction"),
        PB.ParDouble("K_dissol", 0.0; units="yr-1",
            description="rate constant for dissolution reaction"),
        PB.ParDouble("dissol_rolloff_conc", 0.0; units="mol m-3",
            description="limiting concentration below which dissolution rate is rolled off to zero as Solid_conc^2")
    )
end


function PB.register_methods!(rj::ReactionAqPrecipDissol)

    @info "register_methods! $(PB.fullname(rj))"

    reactant_names = String[]
    all_stoich = Float64[]
    for rns in rj.pars.Reactants
        stoich, name = PALEOaqchem.parse_number_name(rns)
        push!(all_stoich, -1*stoich) # -1 for reactant 
        push!(reactant_names, String(name))
    end
    product_names = String[]
    for pns in rj.pars.Products
        stoich, name = PALEOaqchem.parse_number_name(pns)
        push!(all_stoich, stoich)
        push!(product_names, String(name))
    end

    all_names = vcat(reactant_names, product_names)
    all_stoich_names = Tuple((s, n) for (s,n) in zip(all_stoich, all_names))
    ratevarname = replace("precipdissol_"*join(all_names, "_"), "_TAlk"=>"")
    stoich_precip = PB.RateStoich(
        PB.VarProp(ratevarname, "mol yr-1", "rate variable",
                        attributes=(:calc_total=>true,)),
        all_stoich_names, 
        sms_prefix="", 
        sms_suffix="_sms",
        processname="precipdissol",
    )
    
    omega_var = PB.VarDep("Omega", "", "saturation state")
    solid_conc_name = isempty(rj.pars.Solid_conc[]) ? first(product_names) .* "_conc" : rj.pars.Solid_conc[]
    solid_conc_var = PB.VarDep(solid_conc_name, "mol m-3", "solid concentration or activity")
    volume = PB.VarDep("volume", "m3", "cell solid phase volume")

    PB.add_method_do!(
        rj,
        do_aqprecipdissol, 
        (
            PB.VarList_single(stoich_precip.ratevartemplate),
            PB.VarList_single(solid_conc_var),
            PB.VarList_single(omega_var),
            PB.VarList_single(volume),
        );
    )

    PB.add_method_do!(rj, stoich_precip)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    

    return nothing
end

function do_aqprecipdissol(m::PB.ReactionMethod, pars, (ratevar, solid_conc, omega, volume), cellrange::PB.AbstractCellRange, deltat)
  
    K_precip = pars.K_precip[]
    K_dissol = pars.K_dissol[]

    dissol_rolloff_conc = pars.dissol_rolloff_conc[]    
    f_rolloff(x, xmin) = x/(xmin + x)

    for i in cellrange.indices
        # set ratevar to 0, with fake dependencies on variables for AD sparsity detection 
        ratevar[i] = 0.0*omega[i] + 0.0*solid_conc[i]
        if omega[i] > 1
            # mol yr-1 =  mol m-3 yr-1          * m^3      * bodge
            ratevar[i] = K_precip*(omega[i] - 1)*volume[i] # * f_rolloff(omega[i]-1, 1e-1)
        elseif omega[i] < 1
            sc_checked = max(solid_conc[i], 0.0)
            # mol yr-1 = yr-1     mol m-3                     * m^3     * additional rolloff at low concentrations (to help numerics)
            ratevar[i] = -K_dissol*sc_checked*(1 - omega[i])*volume[i] *f_rolloff(sc_checked, dissol_rolloff_conc) # * f_rolloff(1 - omega[i], 1e-1)
        end
    end
   
    return nothing
end




end
