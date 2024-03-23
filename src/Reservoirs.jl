module ReservoirsAq

import PALEOboxes as PB
using PALEOboxes.DocStrings

import ..PALEOaqchem

import Infiltrator


"""
    ReactionConstraintReservoir

A primary species and (algebraic) constraint on a corresponding total or component.

The primary species concentration or amount is defined as a PALEO `State` Variable, which depending on the `primary_variable` parameter, may be:
- `Primary_conc`: (mol m-3)
- `Primary`: (mol)
- `Primary_pconc`: -log 10 (concentration (mol kg-1))

The corresponding `R_constraint_conc` or `R_constraint` (mol) defining the algebraic constraint on the corresponding total
(for use by the numerical solver) is defined as a PALEO `Constraint` Variable.

This `ReactionConstraintReservoir` would usually be used in combination with a `ReactionReservoir` that provides the required total component concentration or amount
as an ODE variable (where as usual reaction source and sink fluxes are applied to the 
corresponding `_sms` variable). Depending on the `constraint_variable` parameter, 
the total component may be supplied as either a per-cell concentration or amount:
- `R_conc`: (mol m-3)
- `R`: (mol)

Equilibrium reactions defining secondary species should add their contributions to the total to `R_calc` (mol).
A primary species contribution `R_calc += primary_total_stoich * Primary_conc * primary_volume` is added to `R_calc`
(where for the usual case parameter `primary_total_stoich` should be set to `1.0`). Primary species contributions to other
totals can be included by setting the `primary_other_components` parameter.

The numerical solver then solves for the primary species (and hence the secondary species concentrations)
that (depending on the `constraint_variable` parameter) satisfy one of:

    0 = R_constraint_conc = R_conc - R_calc/volume
    0 = R_constraint = R - R_calc

# Volume conversions
The total species concentration `R_conc` and primary species concentration `Primary_conc` use (potentially 
different) volume conversions provided in `volume` and `primary_volume` respectively.

This allows for cases eg equilibrium partitioning between solute and solid phases by surface complexation,
where `R_conc` refers to a cell total volume, and `Primary_conc` to a solute concentration.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionConstraintReservoir{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("primary_total_stoich", 1.0,
            description="stoichiometric factor R_calc_conc += primary_total_stoich * Primary_conc"),
        PB.ParStringVec("primary_other_components", String[] ,
            description="contribution of primary species to other element or component total concentrations"),
        PB.ParString("primary_variable", "concentration";
            allowed_values=["concentration", "amount", "p_concentration"],
            description="units for primary variable"),
        PB.ParString("constraint_variable", "concentration";
            allowed_values=["concentration", "amount"],
            description="units for constraint variable"),

    )

    component_stoichs::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionConstraintReservoir)

    R_calc = PB.VarTarget("R_calc", "m-3", "contributions to total R_calc_conc (NB: a total, not concentration, to generalize to multiphase eqb)")

    default_setup_vars = PB.VariableReaction[] # any state Variables that need "standard" setup (no volume conversions etc)

    # define the variables we need for Primary_conc etc
    primary_volume = PB.VarDep("primary_volume"=>"volume", "m3", "cell volume (as used by Primary_conc)")
    primary_vars = PB.VariableReaction[PB.VarContrib(R_calc), primary_volume]
    if rj.pars.primary_variable[] == "concentration"
        Primary_conc = PB.VarState("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse)
        push!(default_setup_vars, Primary_conc)
        append!(primary_vars, [
            Primary_conc,
        ])
    elseif rj.pars.primary_variable[] == "amount"        
        Primary = PB.VarState("Primary", "mol", "concentration of primary species")       
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [Primary],
            convertvars = [primary_volume],
            convertfn = ((primary_volume,), i) -> primary_volume[i],
            convertinfo = " * primary_volume",
        )
        append!(primary_vars, [
            Primary,
            PB.VarProp("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse),
        ])
    elseif rj.pars.primary_variable[] == "p_concentration"
        Primary_pconc = PB.VarState("Primary_pconc", "", "-log10(concentration of primary species)")
        push!(default_setup_vars, Primary_pconc)
        append!(primary_vars, [
            Primary_pconc,
            PB.VarProp("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse),
            PB.VarDep("rho_ref", "kg m-3", "reference density to convert mol m-3 to mol kg-1"),
        ])
    else
        error("invalid parameter 'primary_variable = $(rj.pars.primary_variable[])")
    end

    empty!(rj.component_stoichs)
    component_vars = PB.VariableReaction[]
    for cv in rj.pars.primary_other_components
        stoich, name = PALEOaqchem.parse_number_name(cv)
        push!(rj.component_stoichs, stoich)
        push!(component_vars, PB.VarContrib(name, "mol", "total moles"))
    end

    PB.add_method_do!(rj, do_constraintreservoir_primary, (PB.VarList_namedtuple(primary_vars), PB.VarList_tuple(component_vars),))


    # define the variables we need to calculate the constraint
    
    volume = PB.VarDep("volume", "m3", "cell volume (as used by total variable)")
    constraint_vars = PB.VariableReaction[R_calc]
    if rj.pars.constraint_variable[] == "concentration"
        R_constraint_conc = PB.VarConstraint("R_constraint_conc", "mol m-3", "algebraic constraint on R_conc (= 0)")
        push!(default_setup_vars, R_constraint_conc)
        append!(constraint_vars, [
            R_constraint_conc,
            PB.VarDep("R_conc", "mol m-3", "total R_conc"),
            volume,
        ])
        
    elseif rj.pars.constraint_variable[] == "amount"
        R_constraint = PB.VarConstraint("R_constraint", "mol", "algebraic constraint on R (= 0)")
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [R_constraint],
            convertvars = [volume],
            convertfn = ((volume,), i) -> volume[i],
            convertinfo = " * volume",
        )
        append!(constraint_vars, [
            R_constraint,
            PB.VarDep("R", "mol", "total R"),
        ])

    else
        error("invalid parameter 'constraint_variable' = $(pars.constraint_variable[])")
    end
    
    PB.add_method_do!(rj, do_constraintreservoir_constraint, (PB.VarList_namedtuple(constraint_vars),))
    
    PB.add_method_setup_initialvalue_vars_default!(rj, default_setup_vars)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_constraintreservoir_primary(m::PB.AbstractReactionMethod, pars, (vars, components), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    do_amount = (pars.primary_variable[] == "amount")
    do_p_concentration = (pars.primary_variable[] == "p_concentration")
   
    for i in cellrange.indices
        if do_amount
            vars.Primary_conc[i] = vars.Primary[i]/vars.primary_volume[i]
        end
        if do_p_concentration
            # mol m-3                    mol kg-1              kg m-3
            vars.Primary_conc[i] = 10^(-vars.Primary_pconc[i])*vars.rho_ref[i]
        end
        
        vars.R_calc[i] += pars.primary_total_stoich[]*vars.Primary_conc[i]*vars.primary_volume[i]

        PB.IteratorUtils.foreach_longtuple(components, rj.component_stoichs) do c, s
            c[i] += s * vars.Primary_conc[i]*vars.primary_volume[i]
        end
    end

    return nothing
end

function do_constraintreservoir_constraint(m::PB.AbstractReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)

    if pars.constraint_variable[] == "concentration"
        for i in cellrange.indices
            vars.R_constraint_conc[i]  = vars.R_conc[i] - vars.R_calc[i]/vars.volume[i]
        end
    elseif pars.constraint_variable[] == "amount"
        for i in cellrange.indices
            vars.R_constraint[i]  = vars.R[i] - vars.R_calc[i]
        end
    else
        error("do_constraintreservoir_constraint: unsupported 'constraint_variable' = $(pars.constraint_variable[])")
    end

    return nothing
end



"""
    ReactionImplicitReservoir

A primary species and corresponding total or component as an 'implicit' ODE variable.

This provides an implementation of the 'Direct Substitution Approach' to chemical speciation, where the total or component is a function
of the primary species concentration.

The primary species concentration or amount is defined as a PALEO `StateTotal` Variable, which depending on the `primary_variable` parameter, may be:
- `Primary_conc`: (mol m-3)
- `Primary`: (mol)
- `Primary_pconc`: -log 10 (concentration (mol kg-1))

The corresponding total component `R_conc` or `R` is defined as a PALEO `Total` Variable, which depending on the `constraint_variable` parameter, 
may be provided to the solver either as a per-cell concentration or amount:
- `R_conc = R_calc/volume`: (mol m-3)
- `R = R_calc`: (mol)

Equilibrium reactions defining secondary species should add their contributions to the total to `R_calc` (mol).
A primary species contribution `R_calc += primary_total_stoich * Primary_conc * primary_volume` is added to `R_calc`
(where for the usual case parameter `primary_total_stoich` should be set to `1.0`). Primary species contributions to other
totals can be included by setting the `primary_other_components` parameter.

Source - sink fluxes eg kinetic reactions should be added to `R_sms` (mol yr-1) defined as a PALEO `Deriv` Variable.

# Volume conversions
The total species concentration `R_conc` and primary species concentration `Primary_conc` use (potentially 
different) volume conversions provided in `volume` and `primary_volume` respectively.

This allows for cases eg equilibrium partitioning between solute and solid phases by surface complexation,
where `R_conc` refers to a cell total volume, and `Primary_conc` to a solute concentration.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionImplicitReservoir{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("primary_total_stoich", 1.0,
            description="stoichiometric factor R_calc_conc += primary_total_stoich * Primary_conc"),
        PB.ParStringVec("primary_other_components", String[] ,
            description="contribution of primary species to other element or component total concentrations"),
        PB.ParString("primary_variable", "concentration";
            allowed_values=["concentration", "amount", "p_concentration"],
            description="units for primary variable (specifies Primary_conc, Primary, Primary_pconc as StateTotal variable)"),
        PB.ParString("total_variable", "concentration";
            allowed_values=["concentration", "amount"],
            description="units for total variable (specifies R_conc, R as Total variable)"),
        PB.ParBool("total", false,
            description="true to calculate R_total = sum(R)"),
    )

    component_stoichs::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionImplicitReservoir)

    R_calc = PB.VarTarget("R_calc", "m-3", "contributions to total R_calc_conc (NB: a total, not concentration, to generalize to multiphase eqb)")

    default_setup_vars = PB.VariableReaction[] # any state Variables that need "standard" setup (no volume conversions etc)
    
    # define the variables we need for Primary_conc etc
    primary_volume = PB.VarDep("primary_volume"=>"volume", "m3", "cell volume (as used by Primary_conc)")
    primary_vars = PB.VariableReaction[PB.VarContrib(R_calc), primary_volume]
    if rj.pars.primary_variable[] == "concentration"
        Primary_conc = PB.VarStateTotal("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse)
        push!(default_setup_vars, Primary_conc)
        append!(primary_vars, [
            Primary_conc,
        ])
    elseif rj.pars.primary_variable[] == "amount"        
        Primary = PB.VarStateTotal("Primary", "mol", "concentration of primary species")       
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [Primary],
            convertvars = [primary_volume],
            convertfn = ((primary_volume,), i) -> primary_volume[i],
            convertinfo = " * primary_volume",
        )
        append!(primary_vars, [
            Primary,
            PB.VarProp("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse),
        ])
    elseif rj.pars.primary_variable[] == "p_concentration"
        Primary_pconc = PB.VarStateTotal("Primary_pconc", "", "-log10(concentration of primary species)")
        push!(default_setup_vars, Primary_pconc)
        append!(primary_vars, [
            Primary_pconc,
            PB.VarProp("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes_advectfalse),
            PB.VarDep("rho_ref", "kg m-3", "reference density to convert mol m-3 to mol kg-1"),
        ])
    else
        error("invalid parameter 'primary_variable = $(rj.pars.primary_variable[])")
    end

    empty!(rj.component_stoichs)
    component_vars = PB.VariableReaction[]
    for cv in rj.pars.primary_other_components
        stoich, name = PALEOaqchem.parse_number_name(cv)
        push!(rj.component_stoichs, stoich)
        push!(component_vars, PB.VarContrib(name, "mol", "total moles"))
    end

    PB.add_method_do!(rj, do_implicitreservoir_primary, (PB.VarList_namedtuple(primary_vars), PB.VarList_tuple(component_vars),))

    # define the variables we need to calculate the total and total_sms    
    volume = PB.VarDep("volume", "m3", "cell volume (as used by total variable)")
    total_vars = PB.VariableReaction[R_calc, volume]
    if rj.pars.total_variable[] == "concentration"
        R = PB.VarProp("R", "mol", "total or component R")
        R_conc = PB.VarTotal("R_conc", "mol m-3", "total or component R_conc"; attributes=PALEOaqchem.R_conc_attributes_advecttrue)
        push!(default_setup_vars, R_conc)

        append!(total_vars, [R, R_conc,])
        
        # provide R_sms to accumulate source - sink fluxes, then convert to R_conc_sms for solver
        total_sms_vars = [
            PB.VarTarget("R_sms", "mol yr-1", "total or component R source - sink"),
            PB.VarDeriv("R_conc_sms", "mol m-3 yr-1", "total or component R_conc source - sink"),
            volume,
        ]

        PB.add_method_do!(rj, do_implicitreservoir_sms, (PB.VarList_namedtuple(total_sms_vars),))

    elseif rj.pars.total_variable[] == "amount"
        R = PB.VarTotal("R", "mol", "total or component R")
        R_conc = PB.VarProp("R_conc", "mol m-3", "total or component R_conc"; attributes=PALEOaqchem.R_conc_attributes_advecttrue)
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [R],
            convertvars = [volume],
            convertfn = ((volume,), i) -> volume[i],
            convertinfo = " * volume",
        )
        
        append!(total_vars, [R, R_conc,])

        R_sms = PB.VarDeriv("R_sms", "mol yr-1", "total or component R source - sink")        
        PB.add_method_do_nothing!(rj, [R_sms]) # we don't need to access this, but must be registered

    else
        error("invalid parameter 'total_variable' = $(pars.total_variable[])")
    end
    
    PB.add_method_do!(rj, do_implicitreservoir_total, (PB.VarList_namedtuple(total_vars),))
    
    PB.add_method_setup_initialvalue_vars_default!(rj, default_setup_vars)

    if rj.pars.total[]
        # workaround: can't calculate total for a VarTotal R so create a VarDep that will link to R
        R_dep = PB.VarDep("R", "mol", "total or component R"; attributes=(:calc_total=>true, :field_data=>PB.get_attribute(R, :field_data),))
        PB.add_method_do_totals_default!(rj, [R_dep])
    end
    PB.setfrozen!(rj.pars.total)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_implicitreservoir_primary(m::PB.AbstractReactionMethod, pars, (vars, components), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    do_amount = (pars.primary_variable[] == "amount")
    do_p_concentration = (pars.primary_variable[] == "p_concentration")
   
    for i in cellrange.indices
        if do_amount
            vars.Primary_conc[i] = vars.Primary[i]/vars.primary_volume[i]
        end
        if do_p_concentration
            # mol m-3                    mol kg-1              kg m-3
            vars.Primary_conc[i] = 10^(-vars.Primary_pconc[i])*vars.rho_ref[i]
        end
        
        vars.R_calc[i] += pars.primary_total_stoich[]*vars.Primary_conc[i]*vars.primary_volume[i]

        PB.IteratorUtils.foreach_longtuple(components, rj.component_stoichs) do c, s
            c[i] += s * vars.Primary_conc[i]*vars.primary_volume[i]
        end
    end

    return nothing
end

function do_implicitreservoir_total(m::PB.AbstractReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)

    for i in cellrange.indices
        vars.R[i]  = vars.R_calc[i]
        vars.R_conc[i] = vars.R[i]/vars.volume[i]
    end
    
    return nothing
end

# only used for pars.total_variable[] == "concentration"
function do_implicitreservoir_sms(m::PB.AbstractReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)

    for i in cellrange.indices
        vars.R_conc_sms[i] += vars.R_sms[i]/vars.volume[i]
    end

    return nothing
end


"""
    ReactionAqConcSum

A sum of concentration variables (eg to get an element total)


# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAqConcSum{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("vars_to_add", ["2*myvar", "myothervar", "-1*mythirdvar"],
            description="vector of variable names to add, eg [2*myvar, myothervar, -1*mythirdvar]"),
        PB.ParBool("add_to_sum_volume", false;
            description="true to also add to a 'sum' Variable  += 'sum_conc * volume"),
        PB.ParBool("define_sum_volume", false;
            description="only if 'add_to_sum_volume == true': true to also define the 'sum' Variable"),
    )

    var_multipliers::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionAqConcSum)

    io = IOBuffer()
    println(io, "register_methods: $(PB.fullname(rj)) $(PB.typename(rj))")

    vars_to_add = []
    empty!(rj.var_multipliers)

    for varmultname in rj.pars.vars_to_add
        # parse multiplier
        mult, varname = PALEOaqchem.parse_number_name(varmultname)
        println(io, "    add $mult * $varname")
        push!(rj.var_multipliers, mult)

        v =  PB.VarDep(varname, "", "")
        push!(vars_to_add, v)
    end

    vars = PB.VariableReaction[
        PB.VarProp("sum_conc", "mol m-3", "sum of specified variables"; attributes=PALEOaqchem.R_conc_attributes_advectfalse)
    ]
    if rj.pars.add_to_sum_volume[]
        push!(vars, PB.VarDep("volume", "m^3", "cell or phase volume"))
        if rj.pars.define_sum_volume[]
            sum = PB.VarTarget("sum", "mol", "sum of specified variables * volume")
        else
            sum = PB.VarContrib("sum", "mol", "sum of specified variables * volume")
        end
        push!(vars, sum)
    end    

    PB.add_method_do!(
        rj,
        do_aqconcsum,
        (
            PB.VarList_namedtuple(vars),
            PB.VarList_tuple(vars_to_add),
        ),
    )

    PB.add_method_initialize_zero_vars_default!(rj) # only actually needed if sum = PB.VarTarget

    @info String(take!(io))

    return nothing
end

function do_aqconcsum(
    m::PB.ReactionMethod,
    pars,
    (vars, vars_to_add),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    for i in cellrange.indices
        vars.sum_conc[i] = 0.0
        PB.IteratorUtils.foreach_longtuple(vars_to_add, rj.var_multipliers) do d, m
            vars.sum_conc[i] += d[i]*m
        end

        if hasfield(typeof(vars), :sum)
            vars.sum[i] += vars.sum_conc[i]*vars.volume[i]
        end
    end

    return nothing
end



end # module