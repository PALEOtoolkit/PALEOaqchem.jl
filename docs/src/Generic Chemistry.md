# Generic Chemistry

```@meta
CurrentModule = PALEOaqchem
```

Generic chemical equilibrium and kinetic reactions following the
standard approach used by reaction-transport codes such as PHREEQ
and CrunchFlow, see eg [steefel_reactive_2015](@cite). This exploits timescale
separation between "fast" (assumed instantaneous) chemical equilibrium reactions, and "slow" kinetic reactions or transport.

- The chemical system is represented by a small number of `totals`
(or `components`) and an equal number of `primary` species concentrations, with `secondary` species
concentrations calculated from the primary species via a set of equilibrium reactions. Primary species concentrations
are determined by solving the set of algebraic equations given by the constraints on total concentrations.
- Kinetic reactions (with any species, `primary` or `secondary` as reactants) are then written with `totals` as products.
- Bulk transport (eg ocean advection or eddy diffusivity) transports `totals`.  Molecular diffusivity (eg in a sediment) transports `primary` or `secondary` species and accumulates fluxes into `totals`.

## Reservoirs

```@docs
ReservoirsAq.ReactionConstraintReservoir
ReservoirsAq.ReactionImplicitReservoir
ReservoirsAq.ReactionAqConcSum
```

## Equilibrium reactions

```@docs
GenericReactions.ReactionAqEqb
```

## Kinetic reactions
```@docs
GenericReactions.ReactionAqKinetic
```

## Precipitation-dissolution reactions

```@docs
GenericReactions.ReactionAqPrecipDissol
```

## Particulate fluxes
```@docs
Particle.ReactionParticleDecay
Particle.ReactionFluxToComponents
```


## Co-precipitation
```@docs
CoPrecip.ReactionPACoPrecip
CoPrecip.ReactionPARelease
```
