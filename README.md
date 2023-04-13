# PALEOaqchem.jl

[![CI](https://github.com/PALEOtoolkit/PALEOaqchem.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOaqchem.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOaqchem.jl/dev)

Aquatic biogeochemistry components for the PALEO biogeochemical model framework.


**NB: work-in-progress - this repo contains an initial minimal example only to test infrastructure.**

## Using PALEOaqchem Reactions from other models

The PALEOaqchem Reactions are available to the PALEO framework when the registered PALEOaqchem package is loaded (without downloading the repository), ie

    julia> Pkg.add("PALEOaqchem")
    julia> import PALEOaqchem