module RCmultiG

import PALEOboxes as PB
using PALEOboxes.DocStrings
import SpecialFunctions
import Printf
# import Infiltrator

"""
    ReactionRCmultiG

Discrete representation of organic matter reaction-continuum model from Dale (2015)

The fraction of POC having reactivity `k (yr-1)` in interval `dk` at time zero is given by the PDF:

``f(k, 0) dk = \\frac{a^{\\nu} k^{\\nu-1} e^{-a k}}{\\Gamma(\\nu)} dk``

This is represented by a discrete distribution in n-2 bins spaced logarithmically in k, with additional
bins at the start and end to represent POC with reactivity below and above the first and last bin.

This Reaction: 
- Defines n state variables `Corg_<n>` with corresponding concentrations `Corg_<n>_conc` and source-sink `Corg_<n>_sms`,
  and also calculates total concentration (summed over POC fractions) `Corg_conc`.
- Applies POC degradation to each discrete fraction and accumulates to variable `POC_decay`.
- If `do_rc_distribution = true`, distributes input POC defined by `Corg_sms` into the `Corg_<n>` state variables,
  according to the discretised RC continuum distribution, optionally modified by multipliers specified by parameter `k_dist_modifier`.

# References
- Dale (2015) GBC <https://dx.doi.org/10.1002/2014GB005017>

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRCmultiG{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("k_bin_edges", [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2], units="yr-1",
            description="n-1 rate-bin edges. Organic matter reactivity is represented an initial bin for k < k_bin_edges[1], n-2 bins centered in log(k), and a final bin for k > bin_edges[end]"),

        PB.ParBool("do_RC_distribution", true,
            description="true to transfer Corg_sms to Corg_<n>_sms according to RC initial reactivity distribution, false to just define fractions and decay"),
        PB.ParDoubleVec("k_dist_modifier", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], units="",
            description="n multipliers for distribution of Corg input to reactivity bins (1.0 to use unmodified RC distribution)"),
        PB.ParDouble("a", 3e-4 , units="yr-1",
            description="reactive-continuum model average lifetime"),
        PB.ParDouble("v", 0.125, units="",
            description="reactive-continuum model shape parameter"),

        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
    )

    num_bins::Int64 = -1
    k_bin::Vector{Float64} = Float64[]
    initial_frac_bin::Vector{Float64} = Float64[]
end

function calc_rc_bin_fractions(
    a,
    v,
    k_bin_edges,
)
    # cumulative organic matter vs k is given by incomplete gamma function

    # tuple of (γ(v, a*k)/Γ(v), Γ(v, a*k)/Γ(v)) where γ is the lower incomplete gamma function (integral 0 - x) and Γ the upper incomplete gamma function (integral x - Inf)
    #  lower incomplete gamma function 
    gammas = [SpecialFunctions.gamma_inc(v, a*k) for k in k_bin_edges]

    k_bin, frac_bin = Float64[], Float64[]

    # first bin - everything up to k_bin_edges[1]
    push!(k_bin, k_bin_edges[begin])
    push!(frac_bin, first(gammas[begin]))

    # middle bins - everything between bin edges
    for i in 1:(length(k_bin_edges)-1)
        push!(k_bin, sqrt(k_bin_edges[i]*k_bin_edges[i+1])) # mid point of logarithmicaly spaced bins
        push!(frac_bin, first(gammas[i+1]) - first(gammas[i]))
    end

    # last bin - everything above k_bin_edges[end]
    push!(k_bin, k_bin_edges[end])
    push!(frac_bin, last(gammas[end]))

    return k_bin, frac_bin

end


function PB.register_methods!(rj::ReactionRCmultiG)

    @info "register_methods! $(PB.fullname(rj))"

    rj.num_bins = length(rj.pars.k_bin_edges.v) + 1

    volume  = PB.VarDep("volume",   "m3",       "cell volume (or cell phase volume eg for a sediment with solid and liquid phases)")

    state_vars = [
        PB.VarStateExplicit("Corg_$i",        "mol",      "organic carbon reservoir discrete RC fraction",
            attributes=(:field_data=>rj.pars.field_data[], :calc_total=>true,)
        )
        for i in 1:rj.num_bins
    ]
    PB.add_method_setup_initialvalue_vars_default!(
        rj, state_vars,
        convertvars = [volume],
        convertfn = ((volume,), i) -> volume[i],
        convertinfo = " * volume",
    )
     
    sms_vars = [
            PB.VarDeriv("Corg_$(i)_sms",    "mol yr-1", "organic carbon reservoir discrete RC fraction source - sink",
                attributes=(:field_data=>rj.pars.field_data[], ),
        )
        for i in 1:rj.num_bins
    ]

    conc_vars = [
        PB.VarProp("Corg_$(i)_conc",   "mol m-3",  "concentration",
            attributes=(
                :field_data=>rj.pars.field_data[], 
                :advect=>true,
                :vertical_movement=>0.0,
                :specific_light_extinction=>0.0,
                :vphase=>PB.VP_Undefined,
                :diffusivity_speciesname=>""
            )
        )
        for i in 1:rj.num_bins
    ]
    conc_sum_var = PB.VarProp("Corg_conc",   "mol m-3",  "total Corg concentration from multi G reservoirs",
                        attributes=(:field_data=>rj.pars.field_data[], ))

    PB.add_method_do!(
        rj,
        do_RCmultiG_conc, 
        (
            PB.VarList_single(volume),
            PB.VarList_single(conc_sum_var),
            PB.VarList_vector(state_vars),
            PB.VarList_vector(conc_vars),
        );
    )

    PB.add_method_setup!(
        rj,
        setup_RCmultiG, 
        (
            PB.VarList_vector(state_vars),
        );
    )
    if rj.pars.do_RC_distribution[]
        sms_input_var = PB.VarTarget("Corg_sms",    "mol yr-1", "total organic carbon reservoir source - sink (to distribute to RC fractions)",
                            attributes=(:field_data=>rj.pars.field_data[], ))
        PB.add_method_do!(
            rj,
            do_RCmultiG_distribute_input_flux, 
            (
                PB.VarList_single(sms_input_var),
                PB.VarList_vector(sms_vars),
            );
        )
    end
    PB.setfrozen!(rj.pars.do_RC_distribution)

    PB.add_method_do!(
        rj,
        do_RCmultiG_POC_decay, 
        (
            PB.VarList_single(
                PB.VarContrib("POC_decay", "mol C yr-1", "decay flux";
                    attributes=(:field_data=>rj.pars.field_data[], )),
            ),
            PB.VarList_vector(state_vars),
            PB.VarList_vector(sms_vars),
        );
    )

    PB.add_method_do_totals_default!(rj)
    PB.add_method_initialize_zero_vars_default!(rj)    
    return nothing
end

function setup_RCmultiG(m::PB.ReactionMethod, pars, (state_vars,), cellrange::PB.AbstractCellRange, attribute_name)
    rj = m.reaction

    if attribute_name == :setup 
        io = IOBuffer()
        println(io, "setup_RCmultiG $(PB.fullname(rj)) ReactionRCmultiG")

        # Vectors of VariableReactions corresponding to the (volume, Corg_conc, state_vars) data Arrays
        (v_state_vars, ) = PB.get_variables_tuple(m)

        # calculate bin reactivity and initial fractions for discrete multi-G reactive continuum model
        rj.num_bins == length(pars.k_bin_edges.v) + 1 ||
            error("$(PB.fullname(rj)) num_bins $(rj.num_bins) != length(k_bin_edges) + 1 $(length(pars.k_bin_edges.v) + 1) - length of k_bin_edges has been changed ?")
        
        k_bin, initial_frac_bin = calc_rc_bin_fractions(pars.a[], pars.v[], pars.k_bin_edges.v)  

        if pars.do_RC_distribution[]
            # bin reactivity and distribution          

            rj.num_bins == length(pars.k_dist_modifier.v) ||
                error("$(PB.fullname(rj)) num_bins $(rj.num_bins) != length(k_dist_modifier) $(length(pars.k_dist_modifier.v))")

            modified_initial_frac_bin = [ifb*kdm for (ifb, kdm) in zip(initial_frac_bin, pars.k_dist_modifier.v)]
            modified_initial_frac_bin /= sum(modified_initial_frac_bin) # renormalize to 1.0

            rj.k_bin, rj.initial_frac_bin = k_bin, modified_initial_frac_bin

            Printf.@printf(io, "    %20s%20s%20s%20s%30s\n", "Variable", "k_bin", "initial_frac_bin", "k_dist_modifier", "modified initial_frac_bin")
            for (v_sv, kb, ifb, kdm, mifb) in PB.IteratorUtils.zipstrict(v_state_vars, rj.k_bin, initial_frac_bin, pars.k_dist_modifier.v, rj.initial_frac_bin)
                Printf.@printf(io, "    %20s%20e%20g%20g%30g\n", PB.fullname(v_sv.linkvar), kb,  ifb,  kdm, mifb)
            end
        else
            # bin reactivity only
            rj.k_bin = k_bin
            Printf.@printf(io, "    %20s%20s\n", "Variable", "k_bin")
            for (v_sv, kb) in PB.IteratorUtils.zipstrict(v_state_vars, rj.k_bin)
                Printf.@printf(io, "    %20s%20e\n", PB.fullname(v_sv.linkvar), kb)
            end
        end
        @info String(take!(io))

    # bespoke initialisation from Corg_conc variable - not needed ! (just use standard PB.add_method_setup_initialvalue_vars_default!)
    # elseif attribute_name in (:initial_value, :norm_value)
    #     # initialise each of the Corg_<n> state variables using :initial_value, :initial_delta, :norm_value read from the Corg_conc (Domain) variable
    #     domv_Corg_conc = only(v_Corg_conc).linkvar  # linked Domain variable so we are guaranteed to see updated attribute values set in yaml file
    #     convertfn = (volume, i) -> volume[i]
    #     convertinfo = " * volume"
    #     for (v_sv, sv) in PB.IteratorUtils.zipstrict(v_state_vars, state_vars)
    #         PB.init_values!(sv, pars.field_data[], (), PB.CellSpace, attribute_name, domv_Corg_conc, convertfn, volume, cellrange, (PB.fullname(v_sv.linkvar), convertinfo, ""))
    #     end

    end
   
    return nothing
end

function do_RCmultiG_conc(m::PB.ReactionMethod, pars, (volume, Corg_conc, state_vars, conc_vars), cellrange::PB.AbstractCellRange, deltat)
   
    for i in cellrange.indices
        Corg_conc[i]  = zero(Corg_conc[i])
    end

    for (sv, cv) in PB.IteratorUtils.zipstrict(state_vars, conc_vars)
        for i in cellrange.indices
            cv[i]  = sv[i]/volume[i]
            Corg_conc[i] += cv[i]
        end
    end
   
    return nothing
end

function do_RCmultiG_distribute_input_flux(m::PB.ReactionMethod, pars, (sms_input_var, sms_vars), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # distribute input_flux to RC fractions
    for (sms_v, i_f_b, k_d_m) in PB.IteratorUtils.zipstrict(sms_vars, rj.initial_frac_bin, pars.k_dist_modifier.v)
        for (icol, colindices) in cellrange.columns
            for i in colindices
                sms_v[i]  += sms_input_var[i]*i_f_b*k_d_m
            end
        end
    end
   
    return nothing
end

function do_RCmultiG_POC_decay(m::PB.ReactionMethod, pars, (decay_flux, state_vars, sms_vars), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # first-order POC decay
    for (state_v, sms_v, k_bin) in PB.IteratorUtils.zipstrict(state_vars, sms_vars, rj.k_bin)
        for i in cellrange.indices
            # for numerical stability, set rate to zero for -ve values
            # (multiply by 0 or 1 instead of using a conditional test so that AD sparsity detection always sees a dependency on state_v)
            above_threshold = (PB.get_total(state_v[i]) > 0.0)
            # mol yr-1 =   mol        yr-1
            df = above_threshold*state_v[i] * k_bin
            sms_v[i] -= df
            decay_flux[i] += df
        end
    end
   
    return nothing
end

end # module