# Secondary redox

```@meta
CurrentModule = PALEOaqchem
```

NB: this is just a minimal set of predefined secondary redox reactions for H2S and CH4.
Use [`GenericReactions.ReactionAqKinetic`](@ref) to define additional reactions.

## Sulphur
```@docs
SecondaryRedox.ReactionRedoxH2S_O2
```

## Methane
```@docs
SecondaryRedox.ReactionRedoxCH4_O2
SecondaryRedox.ReactionRedoxCH4_SO4
```
