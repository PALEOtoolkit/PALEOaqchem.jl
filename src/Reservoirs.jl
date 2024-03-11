module ReservoirsAq

import PALEOboxes as PB
using PALEOboxes.DocStrings

import ..PALEOaqchem

import Infiltrator


"""
    ReactionTReservoirConc



# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionTReservoirConc{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("primary_total_stoich", 1.0,
            description="stoichiometric factor R_calc_conc += primary_total_stoich * Primary_conc"),
        PB.ParStringVec("primary_other_components", String[] ,
            description="contribution of primary species to other element or component total concentrations"),
    )

    component_stoichs::Vector{Float64} = Float64[]
end


function PB.register_methods!(rj::ReactionTReservoirConc)

    state_vars = [
        PB.VarState("Primary_conc", "mol m-3", "concentration of primary species"; attributes=PALEOaqchem.R_conc_attributes),   
        PB.VarConstraint("R_constraint_conc", "mol m-3", "algebraic constraint on R_conc (= 0)"),
    ]
    vars = vcat(
        state_vars, 
        PB.VarDep("R_conc", "mol m-3", "total R_conc"),            
        PB.VarTarget("R_calc", "m-3", "contributions to total R_calc_conc (NB: a total, not concentration, to generalize to multiphase eqb)"),
        PB.VarDep("volume", "m3", "cell volume (as used by R_conc)"),
        PB.VarDep("primary_volume", "m3", "cell volume (as used by Primary_conc)"),
    )

    empty!(rj.component_stoichs)
    component_vars = PB.VariableReaction[]
    for cv in rj.pars.primary_other_components
        stoich, name = PALEOaqchem.parse_number_name(cv)
        push!(rj.component_stoichs, stoich)
        push!(component_vars, PB.VarContrib(name, "mol", "total moles"))
    end

    PB.add_method_do!(rj, do_treservoirconc, (PB.VarList_namedtuple(vars), PB.VarList_tuple(component_vars),))

    PB.add_method_setup_initialvalue_vars_default!(rj, state_vars)
    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_treservoirconc(m::PB.AbstractReactionMethod, pars, (vars, components), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # @Infiltrator.infiltrate
    for i in cellrange.indices
        vars.R_calc[i] += pars.primary_total_stoich[]*vars.Primary_conc[i]*vars.primary_volume[i]
        vars.R_constraint_conc[i]  = vars.R_conc[i] - vars.R_calc[i]/vars.volume[i]

        PB.IteratorUtils.foreach_longtuple(components, rj.component_stoichs) do c, s
            c[i] += s * vars.Primary_conc[i]*vars.primary_volume[i]
        end
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
        PB.VarProp("sum_conc", "mol m-3", "sum of specified variables"; attributes=PALEOaqchem.R_conc_attributes)
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