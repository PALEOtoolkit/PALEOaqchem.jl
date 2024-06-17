# Carbonate chemistry
```@meta
CurrentModule = PALEOaqchem
```
Marine carbonate chemistry is implemented by the PALEO [`CarbChem.ReactionCO2SYS`](@ref)
which uses an implementation provided by the [`PALEOcarbchem`](@ref) module.

```@docs
CarbChem.ReactionCO2SYS
```

# PALEOcarbchem
```@docs
PALEOcarbchem
```

```@meta
CurrentModule = PALEOaqchem.PALEOcarbchem
```

## Inputs and equilibrium constants

```@docs
ComponentsAll
get_components_inputs
ConstNames
CNameIdx
getConstDescription
calc_constants!
calc_limits
calc_modern_default_concs
```

## Solvers and outputs

```@docs
ResultNames
RNameIdx
getResultDescription
calculateTAfromTCpHfree!
calculatepHfromTATC!
calculateOmega
mappHscale
```

