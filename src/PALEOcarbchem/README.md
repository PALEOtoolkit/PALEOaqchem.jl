# PALEOcarbchem.jl

Carbonate chemistry equilibrium originally based on Matlab CO2SYS v1.1 (van Heuven et al 2011, Lewis and Wallace 1998) and refactored for speed and extensibility.

Please see <https://github.com/jamesorr/CO2SYS-MATLAB> for the current version of CO2SYS-MATLAB, history, and full citation information.

The license for CO2SYS-MATLAB is included as file LICENSE_CO2SYS.md.

## Differences from CO2SYS-MATLAB
- Refactored for speed, to allow use in numerical models.
- Based on an older (2011) version 1.1 of CO2SYS Matlab
- Uses CO2SYS Matlab equilibrium constants, with a reimplemented equilibrium calculation. NB: Although constants can be calculated for different pH scales, TAlk or pH calculations require constants to be on free pH scale.
- H2S and NH3 added from Hofman etal (2010)

## REFERENCES

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2
System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf.  Anal. Cent.,
Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp.,
https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf

van Heuven, S., D. Pierrot, J.W.B. Rae, E. Lewis, and D.W.R. Wallace (2011)
MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b.  Carbon
Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S.
Department of Energy, Oak Ridge, Tennessee. https://doi.org/10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1

Hofmann AF etal  (2010)
AquaEnv : An Aquatic Acid–Base Modelling Environment in R.
Aquatic Geochemistry, https://doi.org/10.1007/s10498-009-9084-1