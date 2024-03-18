# PALEOaqchem Generic Chemical Reactions

```@meta
CurrentModule = PALEOaqchem
```

Generic chemical equilibrium and kinetic reactions following the
standard approach used by reaction-transport codes such as PHREEQ
and CrunchFlow, see eg [steefel_reactive_2015](@cite).

The chemical system is represented by a small number of totals
or components and associated primary species, with secondary species
in chemical equilibrium. Kinetic reactions are then written in terms
of primary species alone.

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
