module Particle

import PALEOboxes as PB
using PALEOboxes.DocStrings

# import Infiltrator # Julia debugger

"""
    ReactionParticleDecay

Decay (eg organic matter to remineralization) at rate `1.0/decay_timescale` of eg an organic matter dissolved/particulate phase in Reservoir `Particle`,
to `decayflux`. `Particle` may have `isotope_type`.
The Reservoir `Particle` may have the `:vsink` attribute set to represent a marine sinking particulate phase.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionParticleDecay{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("decay_timescale", 0.5, units="yr",
            description="particle decay timescale"),

        PB.ParDouble("decay_threshold", -Inf, units="mol m-3",
            description="particle decay concentration below which decay stops"),

        PB.ParBool("show_decayrate", false,
            description="true to create a 'decayrate' variable to record decay rate"),

        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
 
    )



end

function PB.register_methods!(rj::ReactionParticleDecay)

    vars = [
        PB.VarDep(    "Particle",      "mol", "Particle amount", attributes=(:field_data=>rj.pars.field_data[],)),
        PB.VarContrib("Particle_sms",  "mol yr-1", "Particle source-sink", attributes=(:field_data=>rj.pars.field_data[],)),

        PB.VarContrib("decayflux",     "mol yr-1", "Particle decay flux", attributes=(:field_data=>rj.pars.field_data[],)),
    ]
    if rj.pars.decay_threshold[] > -Inf
        push!(vars, PB.VarDep("Particle_conc", "mol m-3", "Particle concentration"),)
    end
    PB.setfrozen!(rj.pars.decay_threshold)
    if rj.pars.show_decayrate[]
        push!(vars, PB.VarProp("decayrate", "mol yr-1", "Particle decay flux"),)
    end
    PB.setfrozen!(rj.pars.show_decayrate)

    PB.add_method_do!(rj, do_particle_decay, (PB.VarList_namedtuple(vars),))

    return nothing
end

function do_particle_decay(m::PB.ReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        above_threshold = true
        if pars.decay_threshold[] > -Inf
            above_threshold = (vars.Particle_conc[i] > pars.decay_threshold[])
        end
        # mol yr-1
        decayflux = above_threshold*vars.Particle[i]/pars.decay_timescale[]
        vars.Particle_sms[i] -= decayflux
        vars.decayflux[i] += decayflux
        if pars.show_decayrate[]
            vars.decayrate[i] = decayflux
        end
    end

    return nothing
end


"""
    ReactionFluxToComponents

Distribute a single `input_flux` (eg an organic matter flux)
to `output_fluxnames` components with stoichiometry `output_fluxstoich`. `input_flux` may have
an isotope type (set by `field_data`) in which case the isotopic composition will be sent to (usually one) `output_fluxname` with `::Isotope` suffix.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionFluxToComponents{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("outputflux_prefix", "",
            description="Prefix for output flux names"),
        PB.ParStringVec("outputflux_names", ["Corg", "N", "P"],
            description="Suffixes for output flux names.  Use ::Isotope suffix to identify a flux with 'isotope_type'"),
        PB.ParDoubleVec("outputflux_stoich", [106.0, 16.0, 1.0],
            description="Stoichiometry for output fluxes relative to input flux"),

        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

    )

end

function PB.register_methods!(rj::ReactionFluxToComponents)

    var_inputflux = PB.VarTarget(
        "inputflux",  "mol m-3", "input flux", attributes=(:field_data=>rj.pars.field_data[],)
    )

    vars_outputflux = PB.Fluxes.FluxContrib(
        rj.pars.outputflux_prefix[], rj.pars.outputflux_names,
        isotope_data=Dict("Isotope"=>rj.pars.field_data[]),
    )

    PB.add_method_initialize_zero_vars_default!(rj, [var_inputflux])

    PB.add_method_do!(
        rj,
        do_flux_to_components,
        (PB.VarList_single(var_inputflux), PB.VarList_tuple(vars_outputflux)),
    )

    return nothing
end


function do_flux_to_components(
    m::PB.ReactionMethod,
    pars,
    (var_inputflux, vars_outflux),
    cellrange::PB.AbstractCellRange,
    deltat
)
    # inputflux and outflux are the same type (both isotopes, or both scalars)
    function do_outflux(outflux::AbstractVector{T}, stoich, inputflux::AbstractVector{T})  where {T}
        @inbounds for i in cellrange.indices
            outflux[i] += stoich*var_inputflux[i]
        end
        return nothing
    end
    # outflux is a scalar and inputflux is an isotope
    function do_outflux(outflux, stoich, _)
        @inbounds for i in cellrange.indices
            outflux[i] += stoich*PB.get_total(var_inputflux[i])
        end
        return nothing
    end

    PB.IteratorUtils.foreach_longtuple_p(
        do_outflux, vars_outflux, pars.outputflux_stoich, var_inputflux
    )

    return nothing
end


end
