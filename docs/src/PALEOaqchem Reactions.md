# PALEOaqchem Reactions

```@meta
CurrentModule = PALEOaqchem
```

## Particulate fluxes
```@docs
Particle.ReactionParticleDecay
Particle.ReactionFluxToComponents
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
```

## Secondary redox

### Nitrogen
```@docs
SecondaryRedox.ReactionRedoxNH4_NO2
SecondaryRedox.ReactionRedoxNH4_O2
SecondaryRedox.ReactionRedoxNO2_O2
```

### Manganese
```@docs
SecondaryRedox.ReactionRedoxMnII_O2
SecondaryRedox.ReactionRedoxMnIV_FeII
SecondaryRedox.ReactionRedoxMnIV_H2S
```

### Iron and iron sulphides
```@docs
SecondaryRedox.ReactionRedoxFeII_O2
SecondaryRedox.ReactionRedoxFeIII_H2S
SecondaryRedox.ReactionRedoxFeS_O2
SecondaryRedox.ReactionRedoxFeS2pyr_O2
```

### Sulphur
```@docs
SecondaryRedox.ReactionRedoxH2S_O2
```

### Methane
```@docs
SecondaryRedox.ReactionRedoxCH4_O2
SecondaryRedox.ReactionRedoxCH4_SO4
```

### Hydrogen
```@docs
SecondaryRedox.ReactionRedoxH2_O2
SecondaryRedox.ReactionRedoxH2_SO4
```

## Iron-sulphur system
```@docs
FeS.ReactionFeSaq
FeS.ReactionFeSm
FeS.ReactionPyrH2S
```


## Co-precipitation
```@docs
CoPrecip.ReactionPACoPrecip
CoPrecip.ReactionPARelease
```

## Phosphorus
```@docs
CFA.ReactionCFAsimple
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
