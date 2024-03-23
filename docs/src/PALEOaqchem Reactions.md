# PALEOaqchem Reactions

```@meta
CurrentModule = PALEOaqchem
```

## Particulate organic carbon fluxes
```@docs
RCmultiG.ReactionRCmultiG
```

## Remin
```@docs
Remin.ReactionReminPonly
Remin.ReactionReminO2
Remin.ReactionReminO2_SO4
Remin.ReactionReminO2_SO4_CH4
Remin.ReactionReminO2_Fe_SO4_CH4
Remin.ReactionReminO2_NN_Mn_Fe_SO4_CH4
Remin.ReactionReminO2_NO3_Mn_Fe_SO4_CH4
```

## Secondary redox

NB: this is just a minimal set of predefined secondary redox reactions for H2S and CH4.
Use [`GenericReactions.ReactionAqKinetic`](@ref) to define additional reactions.

### Sulphur
```@docs
SecondaryRedox.ReactionRedoxH2S_O2
```

### Methane
```@docs
SecondaryRedox.ReactionRedoxCH4_O2
SecondaryRedox.ReactionRedoxCH4_SO4
```

## Carbonate chemistry
```@docs
CarbChem.ReactionCO2SYS
```

## Isotope systems
### Boron
```@docs
Boron.ReactionBoronIsotope
```
