module GenericReactions

import PALEOboxes as PB
using PALEOboxes.DocStrings

import ..PALEOaqchem

import Infiltrator

"""
    ReactionAqEqb

Define a new equilibrium species `N`  

    N + a A <--> b B + c C

    [N] = K_eqb'^K_power ([B]^b [C]^c) / ([A]^a)

where to convert density units for `K_eqb`:

    K_eqb' = K_eqb * rho_ref^K_density_power

The first name in the `Reactants` list is the new species concentration `N`: 
other species concentrations in `Reactants` and `Products` lists must already be
defined elsewhere in the model configuration.

The contribution of the new species to element totals or components is defined by the `N_components` vector, which may be empty eg to
just calculate an Omega or a gas partial pressure etc.

# Examples

## Gas partial pressure from concentration

    solubility_H2:
        class: ReactionAqEqb
        parameters:
            Reactants:  ["pH2"]
            Products:  ["H2_conc"]
            K_eqb:      7.8e-1   # mol m-3 atm-1 at 298.15 K (Henry's law coefficent)
            K_power:   -1.0  #  pH2 = H2_conc * K_eqb^-1
        variable_attributes:
            pH2%units:   atm

Compilation of Henry's law coefficients: https://www.henrys-law.org/  which is
R. Sander: Compilation of Henry's law constants (version 5.0.0) for water as solvent, Atmos. Chem. Phys., 23, 10901-12440 (2023), doi:10.5194/acp-23-10901-2023

Unit conversions:  1 mol m^-3 Pa-1 = 1.01325e5 mol m-3 atm-1

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAqEqb{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(

        PB.ParStringVec("Reactants", ["N_conc", "A_conc^2"];
            description="concentrations or activities of new species followed by other reactants, write powers as X^2 etc"),
        PB.ParStringVec("Products", [ "B_conc", "C_conc"];
            description="concentrations or activities of products, write powers as X^2 etc"),
        PB.ParDouble("K_eqb", 1.814e-7*1027, units="",
            description="equilibrium constant"),
        PB.ParDouble("K_density_power", 0.0;
            description="multiple K_eqb * rho_ref^K_density_power to convert units: 0.0 for K_eqb in mol m-3, 1.0 for K_eqb in mol kg-1, etc"),
        PB.ParDouble("K_power", -1.0, 
            description="exponent of K_eqb"),
        PB.ParStringVec("N_components", String[],
            description="contribution of new species to element or component total eg '[\"2*TN_calc\"]' to add `2*N_conc*volume` to `TN_calc` (or empty vector to just define an Omega)"),
    )

    component_stoichs::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionAqEqb)

    @info "register_methods! $(PB.fullname(rj))"

    vars = PB.VariableReaction[]

    K_density_power = rj.pars.K_density_power[]
    if !iszero(K_density_power)
        push!(vars, PB.VarDep("rho_ref", "kg m-3", "reference density"))
    end
    # convert to integer if possible, as integer powers are much faster
    K_density_power = isinteger(K_density_power) ? Int(K_density_power) : K_density_power

    K_power = rj.pars.K_power[]
    K_power = isinteger(K_power) ? Int(K_power) : K_power

    N_conc_name = first(rj.pars.Reactants.v)
    N_conc_var = PB.VarProp(N_conc_name, "mol m-3", "aqueous concentration or activity"; attributes=PALEOaqchem.R_conc_attributes_advectfalse)

    other_reactant_conc_vars = PB.VariableReaction[]
    other_reactant_powers = []
    for or_name_power in rj.pars.Reactants.v[2:end]
        orp, orn = PB.parse_name_to_power_number(or_name_power)
        push!(other_reactant_conc_vars, PB.VarDep(orn, "mol m-3", "aqueous concentration or activity"))
        push!(other_reactant_powers, orp)
    end

    product_conc_vars = PB.VariableReaction[]
    product_powers = []
    for pr_name_power in rj.pars.Products
        pp, pn = PB.parse_name_to_power_number(pr_name_power)
        push!(product_conc_vars, PB.VarDep(pn, "mol m-3", "aqueous concentration or activity"))
        push!(product_powers, pp)
    end

    empty!(rj.component_stoichs)
    component_vars = PB.VariableReaction[]
    for cv in rj.pars.N_components
        stoich, name = PB.parse_number_name(cv)
        push!(rj.component_stoichs, stoich)
        push!(component_vars, PB.VarContrib(name, "mol", "total moles"))
    end
    if !isempty(rj.component_stoichs)
        push!(vars, PB.VarDep("volume", "m3", "cell solute volume"))
    end

    PB.setfrozen!(rj.pars.Reactants, rj.pars.Products, rj.pars.N_components)

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
        p = (Tuple(other_reactant_powers), Tuple(product_powers), K_density_power, K_power)
    )

    return nothing
end

function do_aqeqb(m::PB.ReactionMethod, pars, (vars, N_conc, other_reactant_concs, product_concs, components), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    other_reactant_powers, product_powers, K_density_power, K_power = m.p

    for i in cellrange.indices
        sc_denom = mapreduce(x->max(PB.get_total(x[1][i]), 0.0)^x[2], *, zip(other_reactant_concs, other_reactant_powers); init=1.0)
        sc_num = mapreduce(x->max(PB.get_total(x[1][i]), 0.0)^x[2], *, zip(product_concs, product_powers); init=1.0)

        Keqb = pars.K_eqb[]
        if !iszero(K_density_power)
            Keqb *= vars.rho_ref[i]^K_density_power
        end
        Keqb = Keqb^K_power

        N_conc[i] = Keqb * sc_num / sc_denom

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

where this can be modified to different functional form by defining a vector of `Rate_functions` to apply to each concentration variable.

Parameters `Reactants` and `Products` should be the vectors of stoichiometry * <name> of (total) species to accumulate fluxes into `<name>_sms` variables.

Parameter `Reactant_concs` should be an empty vector to take default concentration variable names from `Reactants`, or a Vector of names to specify 
concentration species names explicitly (required when eg where `Reactants` refers to totals which are partioned into multiple species).

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
            description="names of concentration variables to calculate rate eg '[`\"A_conc\"]' etc, empty vector to used defaults from 'Reactants' eg 'A_conc', 'B_conc' ..."),
        PB.ParStringVec("Rate_functions", String[];
            allowed_values=["linear", "sqrt", "monod"],
            description="functional form for rate function of each concentration (empty vector for default 'linear')"),
        PB.ParDouble("K", NaN;
            description="rate constant"),
        PB.ParDouble("K_lim", NaN, units="mol m-3";
            description="limiting concentration for 'monod' rate function"),
    )
end


function PB.register_methods!(rj::ReactionAqKinetic)

    @info "register_methods! $(PB.fullname(rj))"

    reactant_names = String[]
    all_stoich = Float64[]
    for rns in rj.pars.Reactants
        stoich, name = PB.parse_number_name(rns)
        push!(all_stoich, -1*stoich) # -1 for reactant 
        push!(reactant_names, String(name))
    end
    product_names = String[]
    for pns in rj.pars.Products
        stoich, name = PB.parse_number_name(pns)
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
    reactant_conc_vars = [PB.VarDep(rn, "mol m-3", "aqueous concentration or activity") for rn in reactant_concs]
    
    rate_functions = isempty(rj.pars.Rate_functions.v) ? fill("linear", length(reactant_concs)) : rj.pars.Rate_functions.v
    length(rate_functions) == length(reactant_concs) ||
        error("number of 'Rate_functions' != number of 'Reactant_concs")
    @inline f_linear_rate(c, klim) = max(PB.get_total(c), 0.0)
    @inline f_sqrt_rate(c, klim) = sqrt(max(PB.get_total(c), 0.0))
    @inline f_monod_rate(c, klim) = max(PB.get_total(c), 0.0)/(max(PB.get_total(c), 0.0) + klim)
    f_dict = Dict(
        "linear"=>f_linear_rate,
        "sqrt"=>f_sqrt_rate,
        "monod"=>f_monod_rate,
    )
    rate_funcs = [f_dict[fn] for fn in rate_functions]

    volume = PB.VarDep("volume", "m3", "cell solute volume")
    PB.add_method_do!(
        rj,
        do_aqkinetic, 
        (
            PB.VarList_single(stoich_redox.ratevartemplate),
            PB.VarList_tuple(reactant_conc_vars),
            PB.VarList_single(volume),
        );
        p=Tuple(rate_funcs) # supply a tuple so that code is fully typed
    )

    PB.setfrozen!(rj.pars.Reactants, rj.pars.Products, rj.pars.Reactant_concs, rj.pars.Rate_functions)

    PB.add_method_do!(rj, stoich_redox)

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    

    return nothing
end

function do_aqkinetic(m::PB.ReactionMethod, pars, (ratevar, reactant_concs, volume), cellrange::PB.AbstractCellRange, deltat)
    rate_functions = m.p

 
    # 0 reactant_concs
    function _do_aqkinetic(K, K_lim, ratevar, _, volume, cellrange)
        for i in cellrange.indices
            ratevar[i] = K*volume[i]    
        end       
        return nothing
    end

    # 1 reactant_concs
    function _do_aqkinetic(K, K_lim, ratevar, reactant_concs, volume, cellrange, f1)
        for i in cellrange.indices
            ratevar[i] = K*f1(reactant_concs[1][i], K_lim)*volume[i]
        end       
        return nothing
    end

    # 2 reactant_concs
    function _do_aqkinetic(K, K_lim, ratevar, reactant_concs, volume, cellrange, f1, f2)
        for i in cellrange.indices
            ratevar[i] = K*f1(reactant_concs[1][i], K_lim)*f2(reactant_concs[2][i], K_lim)*volume[i]
        end       
        return nothing
    end

    _do_aqkinetic(pars.K[], pars.K_lim[], ratevar, reactant_concs, volume, cellrange, rate_functions...)
   
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
        stoich, name = PB.parse_number_name(rns)
        push!(all_stoich, -1*stoich) # -1 for reactant 
        push!(reactant_names, String(name))
    end
    product_names = String[]
    for pns in rj.pars.Products
        stoich, name = PB.parse_number_name(pns)
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

    PB.setfrozen!(rj.pars.Reactants, rj.pars.Products, rj.pars.Solid_conc)

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
    f_rolloff(x, xmin) = iszero(xmin) ? x : x/(xmin + x)

    for i in cellrange.indices
        # set ratevar to 0, with fake dependencies on variables for AD sparsity detection 
        ratevar[i] = 0.0*omega[i] + 0.0*solid_conc[i]
        if omega[i] > 1
            # mol yr-1 =  mol m-3 yr-1          * m^3      * bodge
            ratevar[i] += K_precip*(omega[i] - 1)*volume[i] # * f_rolloff(omega[i]-1, 1e-1)
        elseif omega[i] < 1
            sc_checked = max(solid_conc[i], 0.0)
            # mol yr-1 = yr-1     mol m-3                     * m^3     * additional rolloff at low concentrations (to help numerics)
            ratevar[i] -= K_dissol*sc_checked*(1 - omega[i])*volume[i] *f_rolloff(sc_checked, dissol_rolloff_conc) # * f_rolloff(1 - omega[i], 1e-1)
        end
    end
   
    return nothing
end




end
