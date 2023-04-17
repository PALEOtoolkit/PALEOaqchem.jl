"""
    PALEOcarbchem

Carbonate chemistry equilibrium translated from Matlab CO2SYS v1.1 [VanHeuven2011](@cite), [Lewis1998](@cite) and refactored for speed.

Implementation based on CO2SYS constants, with reimplemented equilibrium calculation.

NB: 
- Although effective equilibrium constants can be calculated for different pH scales, TAlk or pH calculations using [`PALEOcarbchem.calculateTAfromTCpHfree!`](@ref)
require constants to be on the free pH scale.
- Effective equilibrium constants are empirical measurements hence are only accurate for modern seawater composition.

# Example usage with default choice of constants, all components enabled:
```jldoctest; setup = :(import ..PALEOcarbchem; using Printf)
julia> PALEOcarbchem.ComponentsAllStrings  # all available components
("Ci", "S", "F", "B", "P", "Si", "H2S", "NH3", "Omega")

julia> comps, concinputs = PALEOcarbchem.get_components_inputs(["Ci", "S", "F", "B", "P", "Si", "H2S", "NH3", "Omega"])
((Ci = Val{true}(), S = Val{true}(), F = Val{true}(), B = Val{true}(), P = Val{true}(), Si = Val{true}(), H2S = Val{true}(), NH3 = Val{true}(), Omega = Val{true}()), ["TCi", "TS", "TF", "TB", "TP", "TSi", "TH2S", "TNH3", "Ca"])

julia> println(concinputs) # concentrations required (in addition to TAlk)
["TCi", "TS", "TF", "TB", "TP", "TSi", "TH2S", "TNH3", "Ca"]  

julia> options = (; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=comps)
(WhichKs = Val{10}(), WhoseKSO4 = Val{1}(), pHScale = Val{3}(), Components = (Ci = Val{true}(), S = Val{true}(), F = Val{true}(), B = Val{true}(), P = Val{true}(), Si = Val{true}(), H2S = Val{true}(), NH3 = Val{true}(), Omega = Val{true}()))

julia> C = zeros(length(PALEOcarbchem.ConstNames));

julia> PALEOcarbchem.calc_constants!(C, 25.0, 1000.0, 35.0, Options=options) # Temp(C), P(dbar), Sal(psu)

julia> C_NT = NamedTuple{PALEOcarbchem.ConstNames}(C);  # convert to NamedTuple for convenience

julia> map(x -> @sprintf("%.14e", x), C_NT)
(TF = "6.83258396883673e-05", TS = "2.82354341328601e-02", fH = "7.13404318750000e-01", VPFac = "9.69344700036820e-01", KW = "5.12594224560177e-14", KF = "2.46750587115740e-03", KS = "1.07228495518292e-01", K0 = "2.83918818040155e-02", K1 = "1.23204547404511e-06", K2 = "9.14637929020938e-10", FugFac = "9.96810458692103e-01", KB = "2.23064975910959e-09", KP1 = "2.01400858730558e-02", KP2 = "9.31686820250764e-07", KP3 = "1.40307711726248e-09", KSi = "3.62200612635021e-10", KH2S = "2.60000514374855e-07", KNH3 = "4.74226503827862e-10", KCa = "4.92763187414538e-07", KAr = "7.39194729649679e-07")

julia> modern_concs = PALEOcarbchem.calc_modern_default_concs(35.0, Options=options);

julia> map(x -> @sprintf("%.14e", x), modern_concs)
(TF = "6.83258396883673e-05", TS = "2.82354341328601e-02", TB = "4.15700000000000e-04", Ca = "1.02845697008497e-02")

julia> input_concs = (TCi=[2000e-6], TS=[modern_concs.TS], TF=[modern_concs.TF], TSi=[1e-3], TP=[1e-6],  TB=[modern_concs.TB], TH2S=[1e-6], TNH3=[1e-6], Ca=[modern_concs.Ca]);  # all in mol kg-1

julia> res = zeros(length(PALEOcarbchem.ResultNames));

julia> (pHfree, steps) = PALEOcarbchem.calculatepHfromTATC!(res, C, options, 2300e-6, input_concs);

julia> @printf("%.14f", pHfree)
8.04695972423890

julia> steps
5

julia> res_NT = NamedTuple{PALEOcarbchem.ResultNames}(res); # convert to NamedTuple for convenience

julia> map(x -> @sprintf("%.14e", x), res_NT)
(pHfree = "8.04695972423890e+00", H = "8.97512024441008e-09", OH = "5.71127974446284e-06", TA = "2.30000000000000e-03", dTAdpH = "6.71952166583155e-04", TS = "2.82354341328601e-02", HSO4 = "2.36333069917159e-09", TF = "6.83258396883673e-05", HF = "2.48522365702671e-10", TCi = "2.00000000000000e-03", CO2 = "1.31351930129921e-05", HCO3 = "1.80311290118315e-03", CO3 = "1.83751905803853e-04", CAlk = "2.17061671279086e-03", fCO2 = "4.62639042514413e-04", pCO2 = "4.64119370418157e-04", xCO2dryinp = "4.78797037215479e-04", TB = "4.15700000000000e-04", BAlk = "8.27503245712347e-05", TP = "1.00000000000000e-06", H3PO4 = "3.68183719676091e-15", H2PO4 = "8.26200822875433e-09", HPO4 = "8.57660283752523e-07", PO4 = "1.34077704336885e-07", PengCorrection = "0.00000000000000e+00", PAlk = "1.12581568874446e-06", TSi = "1.00000000000000e-03", SiAlk = "3.87906357916080e-05", TH2S = "1.00000000000000e-06", HSAlk = "9.66632218352735e-07", TNH3 = "1.00000000000000e-06", NH3Alk = "5.01861680455629e-08", Ca = "1.02845697008497e-02", OmegaCA = "3.83512675291204e+00", OmegaAR = "2.55657840498853e+00")

julia> pHtot = PALEOcarbchem.mappHscale(C, pHfree, Val(3), Val(1), res_NT.TS, res_NT.TF);

julia> @printf("%.14f", pHtot)
7.94544626702045
```
"""
module PALEOcarbchem

include("Constants.jl")

include("Solvers.jl")

end