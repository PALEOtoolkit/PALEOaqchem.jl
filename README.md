# PALEOaqchem.jl

[![CI](https://github.com/PALEOtoolkit/PALEOaqchem.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOaqchem.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOaqchem.jl/dev)

Aquatic biogeochemistry components for the [PALEOtoolkit](https://github.com/PALEOtoolkit) biogeochemical model framework.

The PALEOaqchem package implements biogeochemistry components that provide:
- representation of organic matter phases (as discrete reactivity fractions or as a reactive-continuum model) and remineralization by oxygen, nitrogen, manganese, iron, sulphur and disproportionation to methane.
- generic equilibrium, kinetic, and precipitation-dissolution reactions, which can be configured into reaction networks using the PALEO YAML configuration file.
- pre-packaged reactions to simplify implementation of marine carbonate chemistry and a subset of secondary redox reactions, including support for carbon and sulphur isotope systems.

It can be used in two main ways:
- to construct arbitrary aqueous biogeochemical models, including equilibrium, kinetic and precipitation-dissolution reactions for a configurable set of solute and solid-phase species.
- to construct more-or-less standard ocean or sediment biogeochemical models (with commonly used carbonate chemistry, primary and secondary redox reactions) out of predefined components.

## Installation

The PALEOaqchem Reactions are available to the [PALEOtoolkit](https://github.com/PALEOtoolkit) framework when the registered PALEOaqchem package is installed and loaded:

    julia> Pkg.add("PALEOaqchem")  # install PALEOaqchem in currently active Julia environment
    julia> import PALEOaqchem 

## Documentation

Documentation is available online at <https://paleotoolkit.github.io/PALEOaqchem.jl/>

## Examples of usage

The [PALEOtutorials](https://github.com/PALEOtoolkit/PALEOtutorials.jl) repository includes a [minimal model of the marine carbonate system ](https://paleotoolkit.github.io/PALEOtutorials.jl/stable/collated_examples/configurable_chemistry/README/) implementated using PALEOaqchem generic chemistry.

The [PALEOocean](https://github.com/PALEOtoolkit/PALEOocean.jl) and [PALEOsediment](https://github.com/PALEOtoolkit/PALEOsediment.jl) include examples of water column and sediment biogeochemistry.

## Credits

The PALEOcarbchem carbonate chemistry equilibrium code originally based on Matlab CO2SYS v1.1 ([van Heuven et al 2011](https://doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1), [Lewis and Wallace 1998](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)) and refactored for speed and extensibility.

Please see <https://github.com/jamesorr/CO2SYS-MATLAB> for the current version of CO2SYS-MATLAB, history, and full citation information.

The license for CO2SYS-MATLAB is included as file LICENSE_CO2SYS.md.