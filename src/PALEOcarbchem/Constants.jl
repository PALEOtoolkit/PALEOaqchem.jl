
"""
    ComponentsAll::NamedTuple
    ComponentsAllStrings::Vector{String}

All available components

```jldoctest; setup = :(import ..PALEOcarbchem; using Printf)
julia> PALEOcarbchem.ComponentsAll
(Ci = Val{true}(), S = Val{true}(), F = Val{true}(), B = Val{true}(), P = Val{true}(), Si = Val{true}(), H2S = Val{true}(), NH3 = Val{true}(), Omega = Val{true}())

julia> PALEOcarbchem.ComponentsAllStrings  # all available components
("Ci", "S", "F", "B", "P", "Si", "H2S", "NH3", "Omega")
```
"""
const ComponentsAll = (
    Ci=Val(true),
    S=Val(true),
    F=Val(true),
    B=Val(true),
    P=Val(true),
    Si=Val(true),
    H2S=Val(true),
    NH3=Val(true),
    Omega=Val(true)
)

const ComponentsAllStrings = String.(keys(ComponentsAll))

"""
    get_components_inputs(compstoenable) -> (components::NamedTuple, concinputs::Vector{String})

Returns `components` and `concinputs` required for specified `compstoenable`

`compstoenable` is a list (Vector or Tuple) with component names as Strings, eg ["Ci", "B", "Si"] 
    
Returns:
- a NamedTuple `components` with those components present in `compstoenable` as `Val{true}` and others as `Val{false}`
- a Vector `concinputs` of input concentrations required (P requires TP etc).

# Examples
All components:
```jldoctest; setup = :(import ..PALEOcarbchem; using Printf)
julia> compsall, concinputsall = PALEOcarbchem.get_components_inputs(["Ci", "S", "F", "B", "P", "Si", "H2S", "NH3", "Omega"])
((Ci = Val{true}(), S = Val{true}(), F = Val{true}(), B = Val{true}(), P = Val{true}(), Si = Val{true}(), H2S = Val{true}(), NH3 = Val{true}(), Omega = Val{true}()), ["TCi", "TS", "TF", "TB", "TP", "TSi", "TH2S", "TNH3", "Ca"])

julia> compsall
(Ci = Val{true}(), S = Val{true}(), F = Val{true}(), B = Val{true}(), P = Val{true}(), Si = Val{true}(), H2S = Val{true}(), NH3 = Val{true}(), Omega = Val{true}())

julia> println(concinputsall) # concentrations required (in addition to TAlk)
["TCi", "TS", "TF", "TB", "TP", "TSi", "TH2S", "TNH3", "Ca"]
```

All optional components disabled (just H2O)
```jldoctest; setup = :(import ..PALEOcarbchem; using Printf)
julia> compsminimal, concinputsminimal = PALEOcarbchem.get_components_inputs([])
((Ci = Val{false}(), S = Val{false}(), F = Val{false}(), B = Val{false}(), P = Val{false}(), Si = Val{false}(), H2S = Val{false}(), NH3 = Val{false}(), Omega = Val{false}()), String[])

julia> println(concinputsminimal)
String[]
```
"""
function get_components_inputs(compstoenable)
    allcomponents = keys(ComponentsAll) # Symbols ie :B, :P etc
    enabled = Set()
    inputs = String[]
    for comp in compstoenable
        symcomp = Symbol(comp)
        if symcomp in allcomponents
            push!(enabled, symcomp)
        else
            error("invalid component $comp")
        end
        if comp in ("Ci", "S", "F", "B", "P", "Si", "H2S", "NH3")
            push!(inputs, "T"*comp)
        elseif comp == "Omega"
            push!(inputs, "Ca")
        end
    end

    return (
        NamedTuple{allcomponents}(
            [
                comp in enabled ? Val(true) : Val(false)
                for comp in allcomponents
            ]
        ),
        inputs,
    )
end


"""
    ConstNames::Tuple

All available CO2SYS effective equilibrium constant names
"""
const ConstNames = (
    :TF,
    :TS,
    :fH,
    :VPFac,
    :KW,
    :KF,
    :KS,
    :K0,
    :K1,
    :K2,
    :FugFac,
    # :TB,
    :KB,
    :KP1,
    :KP2,
    :KP3,
    :KSi,
    :KH2S,
    :KNH3, 
    :KCa,
    :KAr,
    # :Ca
)

""" 
    CNameIdx::NamedTuple

Indices of constants in a Vector of length([`ConstNames`](@ref))
"""
const CNameIdx = NamedTuple{ConstNames}((1:length(ConstNames)))

"""
    getConstDescription(constName::AbstractString) -> (constName, units, description)

Get human-readable description of `constName`
"""
function getConstDescription(constName::AbstractString)
    ConstDescriptions = (
        ("TF",          "mol kg-1",  "Total Fluorine modern value from Sal (NB: used only to convert tabulated constants to SW pH scale)"),
        ("TS",          "mol kg-1",  "Total Sulphate modern value from Sal (NB: used only to convert tabulated constants to SW pH scale)"),
        ("fH",          "",          "Hydrogen ion fugacity modern value from Sal and Temp (NB: used only to convert tabulated constants on NBS scale to SW pH scale)"),
        ("VPFac",       "",         "1 - (saturated vapour pressure of H2O)" ),
        # KW              = KW, # water        
        # KF              = KF,   # hydrogen fluoride
        # KS              = KS,   # bisulfate        
        #    FugFac          = get_FugFac(E, WhichKs),
        
        #    K0              = K0,   # carbonic acid
        #    K1              = K1,
        #    K2              = K2,
        #    KW              = KW, # water        
        #    KF              = KF,   # hydrogen fluoride
        #    KS              = KS,   # bisulfate  

        #
    )
    for constdesc in ConstDescriptions
        if constName == constdesc[1]
            return constdesc
        end
    end
    return nothing
end
        


"""
    calc_constants!(
        Cout, TempC_input, Pdbar, Sal_input; 
        Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=ComponentsAll))

Calculate carbonate chemistry constants, results are returned in `Cout`.

# Arguments:
- `Cout::Vector`:  (output) Vector of length length([`ConstNames`](@ref)) with calculated equilibrium constants
- `TempC_input`:  temperature deg C. Will be limited to valid range for constants in use (from [`calc_limits`](@ref))
- `Pdbar`:  pressure, dbar. 
- `Sal_input`:      salinity.   Will be limited to valid range for constants in use (from [`calc_limits`](@ref))
- `Options`: NamedTuple with fields:
    - `WhichKs`: Val(Int), default Val(10) 
    - `WhoseKSO4`:  Val(Int), default Val(1) 
    - `pHScale`:  Val(Int), default Val(3) (free - required for [`calculateTAfromTCpHfree!`](@ref))
    - `Components`: components to enable as generated by [`get_components_inputs`](@ref)      


## WhichKs: 
K1 K2 dissociation constants that are to be used.
Val(Int) where Int is:
-   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
-   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
-   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
-   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
-   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0.
-   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
-  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
-  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
-  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
-  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
-  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.

### `WhichKs=Val(10)` (Lueker etal (2000) parameters.)
Recommended as 'best practice' by Dickson (2007) 
(as cited in Orr, J. C., Epitalon, J-M., & Gattuso, J-P. (2014). 
Comparison of seven packages that compute ocean carbonate chemistry. 
Biogeosciences Discussions, 11(4), 5327�5397. doi:10.5194/bgd-11-5327-2014 )

## WhoseKSO4:
KSO4 dissociation constants that are to be used,
in combination with the formulation of the borate-to-salinity ratio to be used.
Having both these choices in a single argument is somewhat awkward,
but it maintains syntax compatibility with the previous version.
Val(Int) where Int is:
-  1 = KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED)
-  2 = KSO4 of Khoo    & TB of Uppstrom 1979
-  3 = KSO4 of Dickson & TB of Lee 2010
-  4 = KSO4 of Khoo    & TB of Lee 2010

## pHScale: Set pH scale.
Val(Int) where Int is:
-  1 = Total scale
-  2 = Seawater scale
-  3 = Free scale   [required for PALEOcarbchem solvers]
-  4 = NBS scale

# Implementation
Modified CO2SYS Constants, split into functions for maintainability

## Comments from the original code:

    % SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
    % Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
    % Outputs: K0, K(), T(), fH, FugFac, VPFac
    % This finds the Constants of the CO2 system in seawater or freshwater,
    % corrects them for pressure, and reports them on the chosen pH scale.
    % The process is as follows: the Constants (except KS, KF which stay on the
    % free scale - these are only corrected for pressure) are
    %       1) evaluated as they are given in the literature
    %       2) converted to the SWS scale in mol/kg-SW or to the NBS scale
    %       3) corrected for pressure
    %       4) converted to the SWS pH scale in mol/kg-SW
    %       5) converted to the chosen pH scale
    %
    %       PROGRAMMER'S NOTE: all logs are log base e
    %       PROGRAMMER'S NOTE: all Constants are converted to the pH scale
    %               pHScale% (the chosen one) in units of mol/kg-SW
    %               except KS and KF are on the free scale
    %               and KW is in units of (mol/kg-SW)^2
    
## Julia-specific details:
`WhichKs`, `WhoseKSO4`, `pHScale` are passed as Types (using Val(Int)), not integer values,
so that Julia can work out which constants and functions to call at compile time.
Similarly, `Comps` is a Type which encodes the selection of components, so the appropriate
code is generated (components included or excluded) at compile time.
"""
function calc_constants!(
    Cout, TempC_input, Pdbar, Sal_input; 
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=ComponentsAll),
)
    TempC = clamp_TempC(TempC_input, Options.WhichKs)
    Sal = clamp_Sal(Sal_input, Options.WhichKs)
 
    RGasConstant        = 83.1451  # ml bar-1 K-1 mol-1, DOEv2
    # RGasConstant      = 83.14472; # ml bar-1 K-1 mol-1, DOEv3
    # C.RGasConstant      = RGasConstant;
    TempK               = TempC + 273.15

    # Collect functions of environmental conditions
    E = (
        TempC               = TempC,
        Sal                 = Sal,
        sqrSal              = sqrt(Sal),
        Pbar                = Pdbar/10.0,         
        RT                  = RGasConstant*TempK,
        TempK               = TempK,
        logTempK            = log(TempK),
        volkg               = (1 - 0.001005*Sal),  # l / kg  volume of 1 kg sw at 1 atm
        # CalculateIonS:
        # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
        IonS                = 19.924 * Sal / (1000 - 1.005   * Sal),
    )

    
    # Use modern default concentrations to convert constants to SW pH scale 
    TF = default_TF(E.Sal)
    TS = default_TS(E.Sal)  
    
    data_KS = getK_KS(E, Options.WhoseKSO4)
    data_KF = getK_KF(E)

    # pH scale conversion factors to convert constant data to SW pHScale,
    # These are not pressure-corrected and assume modern TF, TS to Sal correlations.
    # free = [H+] 
    # tot  = [H+] + [HSO4-]
    # pHtot < pHfree
    # see CO2SYS FindpHOnAllScales  
    pHScaleFacsConstsToSWS = (
        SWStoTOT  = (1 + TS/data_KS.K)/(1 + TS/data_KS.K + TF/data_KF.K),
        FREEtoTOT =  1 + TS/data_KS.K,
        fH = get_fH(E, Options.WhichKs),
    )

    # pressure correct KS, KF
    KS = pressure_correct(data_KS, E.Pbar, E.RT)
    KF = pressure_correct(data_KF, E.Pbar, E.RT)

    # pressure-corrected pH scale conversion factors to place output on chosen pH scale
    pHScaleFacsPcorrected = (
        SWStoTOT  = (1 + TS/KS)/(1 + TS/KS + TF/KF),
        FREEtoTOT =  1 + TS/KS,    
    )
    
    # conversion factor to place output on chosen pH scale
    if Options.pHScale == Val(1)
        #Total
        pHfactor = pHScaleFacsPcorrected.SWStoTOT;            
    elseif Options.pHScale == Val(2)
        #SWS, they are all on this now
        pHfactor = 1.0;
    elseif Options.pHScale == Val(3)
        #pHfree
        pHfactor = pHScaleFacsPcorrected.SWStoTOT/pHScaleFacsPcorrected.FREEtoTOT;
    elseif Options.pHScale == Val(4) #pHNBS
        pHfactor = pHfacs.fH;
    else
        error("invalid pHScale ", Options.pHScale)
    end

    # KS, KF are already pressure-corrected, and stay on the free pH scale
    
    # Water is always included
    data_KW = getK_KW(E, Options.WhichKs)
    KW = pressure_correct(data_KW.KW, E.Pbar, E.RT)*pHfactor

    # Inorganic carbon
    if Options.Components.Ci == Val(true)      
        data_K1K2 = getK_K1K2(E, Options.WhichKs, pHScaleFacsConstsToSWS)

        Cout[CNameIdx.K0]       = getK_K0(E).K0.K   # no pressure correction or pH scale conversion
        Cout[CNameIdx.K1]       = pressure_correct(data_K1K2.K1, E.Pbar, E.RT)*pHfactor
        Cout[CNameIdx.K2]       = pressure_correct(data_K1K2.K2, E.Pbar, E.RT)*pHfactor
        Cout[CNameIdx.FugFac]   = get_FugFac(E, Options.WhichKs)
    end

    # Boron
    if Options.Components.B == Val(true)
        data_KB = getK_KB(E, Options.WhichKs, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KB]       = pressure_correct(data_KB.KB, E.Pbar, E.RT)*pHfactor
    end

    # Phosphate    
    if Options.Components.P == Val(true)   # NB: Options.Components is a Type, so this is a compile-time (not run time) enable/disable.
        data_KP = getK_KP(E, Options.WhichKs, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KP1]      = pressure_correct(data_KP.KP1, E.Pbar, E.RT)*pHfactor
        Cout[CNameIdx.KP2]      = pressure_correct(data_KP.KP2, E.Pbar, E.RT)*pHfactor
        Cout[CNameIdx.KP3]      = pressure_correct(data_KP.KP3, E.Pbar, E.RT)*pHfactor
    end

    if Options.Components.Si == Val(true)
        data_KSi = getK_KSi(E, Options.WhichKs, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KSi]      = pressure_correct(data_KSi.KSi, E.Pbar, E.RT)*pHfactor
    end    

    if Options.Components.H2S == Val(true)
        data_KH2S = getK_KH2S(E, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KH2S]     = pressure_correct(data_KH2S.KH2S, E.Pbar, E.RT)*pHfactor
    end 
    if Options.Components.NH3 == Val(true)
        data_KNH3 = getK_KNH3(E, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KNH3]     = pressure_correct(data_KNH3.KNH3, E.Pbar, E.RT)*pHfactor
    end

    if Options.Components.Omega == Val(true)
        data_KCaAr = getK_KCaAr(E, Options.WhichKs, pHScaleFacsConstsToSWS)
        Cout[CNameIdx.KCa]      = pressure_correct(data_KCaAr.KCa, E.Pbar, E.RT)  # no pH correction
        Cout[CNameIdx.KAr]      = pressure_correct(data_KCaAr.KAr, E.Pbar, E.RT)
    end

           
    Cout[CNameIdx.TF]       = TF     # contemporary value from Sal
    Cout[CNameIdx.TS]       = TS  # contemporary value
    Cout[CNameIdx.fH]       = pHScaleFacsConstsToSWS.fH
    Cout[CNameIdx.VPFac]    = get_VPFac(E)
    Cout[CNameIdx.KW]       = KW # water        
    Cout[CNameIdx.KF]       = KF   # hydrogen fluoride
    Cout[CNameIdx.KS]       = KS   # bisulfate        
      
    return nothing

end

"""
    calc_modern_default_concs(Sal_input, Options=(; WhichKs=Val(10), WhoseKSO4=Val(1)) -> (; TF, TS, TB, Ca)

Calculate modern seawater default concentrations as a function of salinity `Sal_input`

See [`calc_constants!`](@ref) for documentation for `Options`.

# Returns
    - TF: (mol kg-1) Total Fluorine modern value from salinity
    - TS: (mol kg-1) Total Sulphate modern value from salinity
    - TB: (mol kg-1) Total Boron modern value from salinity
    - Ca: (mol kg-1) Ca modern value from salinity

# Example
```jldoctest; setup = :(import ..PALEOcarbchem; using Printf)
julia> options = (; WhichKs=Val(10), WhoseKSO4=Val(1));

julia> modern_concs = PALEOcarbchem.calc_modern_default_concs(35.0, Options=options);

julia> map(x -> @sprintf("%.14e", x), modern_concs)
(TF = "6.83258396883673e-05", TS = "2.82354341328601e-02", TB = "4.15700000000000e-04", Ca = "1.02845697008497e-02")
```
"""
function calc_modern_default_concs(
    Sal_input; 
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1)),
)
    Sal = clamp_Sal(Sal_input, Options.WhichKs)

    return (
        TF = default_TF(Sal),
        TS = default_TS(Sal),  
        TB = default_TB(Sal, Options.WhichKs, Options.WhoseKSO4),
        Ca = default_Ca(Sal, Options.WhichKs),
    )
end
        
"""
    calc_limits( Val{WhichKs::Int}) -> (TminC, TmaxC, salmin, salmax)

Return range limits for  constant set `WhichKs`:
-   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
-   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
-   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
-   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
-   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
-   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0.
-   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
-  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
-  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
-  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
-  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
-  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
"""
function calc_limits( ::Val{WhichKs}) where WhichKs
    # `where WhichKs` has found WhichKs::Int from the Val{WhichKs} Type

    WhichKs >= 1 && WhichKs <=14 || error("invalid WhichKs=", WhichKs)
    #        1      2       3       4       5       6       7       8       9       10      11      12      13      14
    cTminC = (0,    -1,     2,      2,      2,      2,      2,      0,      2,      2,      0,      -1.6,   0,      0)
    cTmaxC = (45,   40,     35,     35,     35,     35,     35,     50,     35,     35,     45,     35,     50,     50)
    csalmin= (5,    10,     20,     20,     20,     19,     19,     0,      0,      19,     5,      34,     1,      1)
    csalmax= (45,   50,     40,     40,     40,     43,     43,     0,      49,     43,     42,     37,     50,     50)
    
    TminC = cTminC[WhichKs]
    TmaxC = cTmaxC[WhichKs]
    salmin = csalmin[WhichKs]
    salmax = csalmax[WhichKs]

    return (TminC, TmaxC, salmin, salmax)
end

"""
    clamp_Sal(Sal_input, WhichKs) -> Sal

Limit salinity to valid range for constants in use
"""
function clamp_Sal(Sal_input, WhichKs)
    (_, _, salmin, salmax) = calc_limits(WhichKs)
    Sal = min(max(Sal_input, salmin), salmax)

    return Sal
end

"""
    clamp_TempC(TempC_input, WhichKs) -> TempC

Limit temperature (Celsius) to valid range for constants in use
"""
function clamp_TempC(TempC_input, WhichKs)
    
    (TminC, TmaxC, _, _) = calc_limits(WhichKs)
    TempC = min(max(TempC_input, TminC), TmaxC)

    return TempC
end

"modern default Total Borate"
function default_TB(Sal, WhichKs, WhoseKSO4)
    # CalculateTB - Total Borate:
    if WhichKs == Val(8)
        # Pure water case.
        TB = 0.0
    elseif WhichKs == Val(6) || WhichKs == Val(7)
        TB = 0.0004106*Sal/35; # in mol/kg-SW
        # this is .00001173*Sali
        # this is about 1% lower than Uppstrom's value
        # Culkin, F., in Chemical Oceanography,
        # ed. Riley and Skirrow, 1965:
        # GEOSECS references this, but this value is not explicitly
        # given here
    else
        if WhoseKSO4 == Val(1) || WhoseKSO4 == Val(2) # If user opted for Uppstrom's values:
            # Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
            # this is .000416*Sali/35. = .0000119*Sali
            # TB(FF) = (0.000232/10.811)*(Sal(FF)/1.80655); % in mol/kg-SW
            TB =  0.0004157*Sal/35; # in mol/kg-SW
        elseif WhoseKSO4 == Val(3) || WhoseKSO4 == Val(4) # If user opted for the new Lee values:
                
            # Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.
            # Geochimica Et Cosmochimica Acta 74 (6): 1801�1811.
            TB =  0.0004326*Sal/35; # in mol/kg-SW
        else
            error("unrecognized WhoseKSO4 ",WhoseKSO4)
        end
    end

    return TB
end

"modern default Total Fluoride"
function default_TF(Sal)
   # CalculateTF;
    # Riley, J. P., Deep-Sea Research 12:219-220, 1965:
    # this is .000068*Sali/35. = .00000195*Sali
    TF = (0.000067/18.998)*(Sal/1.80655); # in mol/kg-SW

    return TF
end

"modern default Total Sulphate"
function default_TS(Sal)
    # CalculateTS ;
    # Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    # this is .02824*Sali/35. = .0008067*Sali
    TS = (0.14/96.062)*(Sal/1.80655); # in mol/kg-SW

    return TS
end

"""
    default_Ca(Sal, WhichKs) -> Ca

Modern default Ca concentration (mol/kg-SW)
"""
function default_Ca(Sal, WhichKs)
    if WhichKs == Val(6) || WhichKs == Val(7)
        #
        # *** CalculateCaforGEOSECS:
        # Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
        # (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
        Ca = 0.01026*Sal/35;
        # Culkin gives Ca = (.0213/40.078)*(Sal/1.80655) in mol/kg-SW
        # which corresponds to Ca = .01030*Sal/35.
    else
        # CalculateCa:
        # '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
        # '       this is .010285*Sali/35
        Ca = 0.02128/40.087*(Sal/1.80655);# ' in mol/kg-SW
    end

    return Ca
end

function get_fH( E, WhichKs)
    # CalculatefH
    # Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
    if WhichKs == Val(8)        
        fH = 1.0 # this shouldn't occur in the program for this case
    elseif WhichKs == Val(7)
            
        fH = (1.29 - 0.00204*  E.TempK + (0.00046 -
            0.00000148*E.TempK)*E.Sal*E.Sal)
        # Peng et al, Tellus 39B:439-458, 1987:
        # They reference the GEOSECS report, but round the value
        # given there off so that it is about .008 (1#) lower. It
        # doesn't agree with the check value they give on p. 456.
    else
            
        fH = (1.2948 - 0.002036*E.TempK + (0.0004607 -
            0.000001475*E.TempK)*E.Sal^2)
        # Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
        # v. 3, 1982 (p. 80);
    end

    return fH
end

"""
    Kdata

Equilibrium constant at 1 atm and pressure correction data.

K = eqb const at 1 atm 
Pressure correction is calculated as 
lnKfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT
K_atpressure  = K * exp(lnKfac)
"""
struct Kdata{T1, T2, T3}
    K::T1
    deltaV::T2    # cm3/mole
    Kappa::T3     # cm3/mole/bar  
end

"apply pressure correction and return pressure-corrected K"
function pressure_correct(Kin::Kdata, Pbar, RT)
    lnKpfac = (-Kin.deltaV + 0.5*Kin.Kappa*Pbar)*Pbar/RT
    return Kin.K*exp(lnKpfac) # ~ Kin.K * (1 + lnKfac) ?
end

"""
    lnKdata

Equilibrium constant at 1 atm and pressure correction data, stored as lnK and facK.

exp(lnK)*facK = eqb const at 1 atm 
Pressure correction is calculated as 
lnKfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT
K_atpressure  = exp(lnK) * facK * exp(lnKfac)
"""
struct lnKdata{T1, T2, T3, T4}
    lnK::T1      # K = exp(lnK)*facK
    facK::T2  
    deltaV::T3    # cm3/mole
    Kappa::T4     # cm3/mole/bar  
end

"apply pressure correction and return corrected K"
function pressure_correct(Kin::lnKdata, Pbar, RT)
    lnKpfac = (-Kin.deltaV + 0.5*Kin.Kappa*Pbar)*Pbar/RT
    return Kin.facK*exp(Kin.lnK + lnKpfac) # ~ Kin.K * (1 + lnKfac) ?
end

#***************************************************************************
    #CorrectKsForPressureNow:
    # Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
    #       the free scale) are on the SWS scale.
    #       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
    #       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
    #       (the pH scales are the same in this case); the other Ks don't matter.
    #
    #
    # No salinity dependence is given for the pressure coefficients here.
    # It is assumed that the salinity is at or very near Sali = 35.
    # These are valid for the SWS pH scale, but the difference between this and
    # the total only yields a difference of .004 pH units at 1000 bars, much
    # less than the uncertainties in the values.
    #****************************************************************************
    # The sources used are:
    # Millero, 1995:
    #       Millero, F. J., Thermodynamics of the carbon dioxide system in the
    #       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    #       See table 9 and eqs. 90-92, p. 675.
    #       TYPO: a factor of 10^3 was left out of the definition of Kappa
    #       TYPO: the value of R given is incorrect with the wrong units
    #       TYPO: the values of the a's for H2S and H2O are from the 1983
    #                values for fresh water
    #       TYPO: the value of a1 for B(OH)3 should be +.1622
    #        Table 9 on p. 675 has no values for Si.
    #       There are a variety of other typos in Table 9 on p. 675.
    #       There are other typos in the paper, and most of the check values
    #       given don't check.
    # Millero, 1992:
    #       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
    #       CRC Press, 1992. See chapter 6.
    #       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
    #               79, and 96 have typos).
    # Millero, 1983:
    #       Millero, Frank J., Influence of pressure on chemical processes in
    #       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
    #       Chester, R., Academic Press, 1983.
    #       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
    #       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
    #       these two are necessary to match the values given in Table 43.24
    # Millero, 1979:
    #       Millero, F. J., The thermodynamics of the carbon dioxide system
    #       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
    #       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
    # Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
    #       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
    #       This matches the GEOSECS results and is in Edmond and Gieskes.
    # Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
    #       boric acid, and the pH of seawater, Limnology and Oceanography
    #       13:403-417, 1968.
    # Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
    #       seawater with respect to calcium carbonate under in situ conditions,
    #       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
    #****************************************************************************
    # These references often disagree and give different fits for the same thing.
    # They are not always just an update either; that is, Millero, 1995 may agree
    #       with Millero, 1979, but differ from Millero, 1983.
    # For WhichKs# = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
    #       KP3, and KSi as for the other cases. Peng et al didn't consider the
    #       case of P different from 0. GEOSECS did consider pressure, but didn't
    #       include Phos, Si, or OH, so including the factors here won't matter.
    # For WhichKs# = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
    #       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
    #       including the factors won't matter.
    #****************************************************************************
    #       deltaVs are in cm3/mole
    #       Kappas are in cm3/mole/bar
    #****************************************************************************


    

function getK_K0(E)
# CalculateK0:
    # Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    TempK100  = E.TempK/100;
    logTempK100 = E.logTempK - log(100.0)
    lnK0 = (-60.2409 + 93.4517 / TempK100 + 23.3585 * logTempK100 +
                E.Sal * (0.023517 - 0.023656 * TempK100 + 0.0047036 * TempK100 ^2))
    K0   = exp(lnK0);                  # this is in mol/kg-SW/atm
                 # K, deltaV, Kappa, TempCfac
    return (K0 = Kdata(K0, 0.0, 0.0),)  
end

function getK_KS(E, WhoseKSO4)
      # CalculateKS:

      if WhoseKSO4 == Val(1) || WhoseKSO4 == Val(3)
        # Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
        # The goodness of fit is .021.
        # It was given in mol/kg-H2O. I convert it to mol/kg-SW.
        # TYPO on p. 121: the constant e9 should be e8.
        # This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
        lnKS = (-4276.1/E.TempK + 141.328 - 23.093*E.logTempK +
        (-13856 /E.TempK + 324.57 - 47.986*E.logTempK)*sqrt(E.IonS) +
            (35474 /E.TempK - 771.54 + 114.723*E.logTempK)*E.IonS +
            (-2698 /E.TempK)*sqrt(E.IonS)*E.IonS + (1776 /E.TempK)*E.IonS^2)
        
        KS = (exp(lnKS)            # this is on the free pH scale in mol/kg-H2O
            * E.volkg);   # convert to mol/kg-SW       
    elseif WhoseKSO4 == Val(2) || WhoseKSO4 == Val(4)
        # Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
        # KS was found by titrations with a hydrogen electrode
        # of artificial seawater containing sulfate (but without F)
        # at 3 salinities from 20 to 45 and artificial seawater NOT
        # containing sulfate (nor F) at 16 salinities from 15 to 45,
        # both at temperatures from 5 to 40 deg C.
        # KS is on the Free pH scale (inherently so).
        # It was given in mol/kg-H2O. I convert it to mol/kg-SW.
        # He finds log(beta) which = my pKS;
        # his beta is an association constant.
        # The rms error is .0021 in pKS, or about .5# in KS.
        # This is equation 20 on p. 33:
        pKS = 647.59 / E.TempK - 6.3451 + 0.019085*E.TempK - 0.5208*sqrt(E.IonS);
        KS = (10 ^(-pKS)          # this is on the free pH scale in mol/kg-H2O
            * E.volkg);    # convert to mol/kg-SW       
    else
        error("invalid WhoseKSO4 ", WhoseKSO4)
    end

    # PressureEffectsOnKS:
    #       This is from Millero, 1995, which is the same as Millero, 1983.
    #       It is assumed that KS is on the free pH scale.
    deltaV = -18.03 + 0.0466*E.TempC + 0.000316*E.TempC^2;
    Kappa = (-4.53 + 0.09*E.TempC)/1000;
    #lnKSfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;

    return  Kdata(KS, deltaV, Kappa)  
end

function getK_KF(E)
# CalculateKF:
    # Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    lnKF = 1590.2/E.TempK - 12.641 + 1.525*sqrt(E.IonS)
    KF   = (exp(lnKF)                 # this is on the free pH scale in mol/kg-H2O
        *E.volkg);          # convert to mol/kg-SW
    # Another expression exists for KF: Perez and Fraga 1987. Not used here since ill defined for low salinity. (to be used for S: 10-40, T: 9-33)
    # Nonetheless, P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
    # lnKF = 874/TempK - 9.68 + 0.111*Sal^0.5; 
    # KF   = exp(lnKF);                   # this is on the free pH scale in mol/kg-SW

    # PressureEffectsOnKF:
    #       This is from Millero, 1995, which is the same as Millero, 1983.
    #       It is assumed that KF is on the free pH scale.
    deltaV = -9.78 - 0.009*E.TempC - 0.000942*E.TempC^2;
    Kappa = (-3.91 + 0.054*E.TempC)/1000;
    #lnKFfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    
    return Kdata(KF, deltaV, Kappa)
end

function getK_KB(E, WhichKs, pHfacs)
     # CalculateKB:
     if WhichKs == Val(8)
        # Pure water case
        KB = 0.0
        lnKB = -Inf
        Kfac = 0.0
        deltaV = 0.0
        Kappa = 0.0 #; this doesn't matter since TB = 0 for this case
    elseif WhichKs == Val(6) || WhichKs == Val(7)
        # This is for GEOSECS and Peng et al.
        # Lyman, John, UCLA Thesis, 1957
        # fit by Li et al, JGR 74:5507-5525, 1969:
        logKB = -9.26 + 0.00886*E.Sal + 0.01*E.TempC;
        # KB = (10 ^(logKB)  # this is on the NBS scale
        #    /pHfacs.fH)               # convert to the SWS scale
        lnKB = logKB*log(10.0)
        Kfac = 1.0/pHfacs.fH
        # lnKBfac = (27.5 - 0.095*TempC)*Pbar/RT;
        deltaV = -(27.5 - 0.095*E.TempC)
        Kappa = 0.0
    else
        # Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
        lnKBtop = (-8966.9 - 2890.53*E.sqrSal - 77.942*E.Sal +
            1.728*E.sqrSal*E.Sal - 0.0996*E.Sal^2)
        lnKB = (lnKBtop/E.TempK + 148.0248 + 137.1942*E.sqrSal +
            1.62142*E.Sal + (-24.4344 - 25.085*E.sqrSal - 0.2474 *
            E.Sal)*E.logTempK + 0.053105*E.sqrSal*E.TempK)
        #KB = (exp(lnKB)    # this is on the total pH scale in mol/kg-SW
        #    /pHfacs.SWStoTOT)         # convert to SWS pH scale
        Kfac = 1.0/pHfacs.SWStoTOT
        #***PressureEffectsOnKB:
        #               This is from Millero, 1979.
        #               It is from data of Culberson and Pytkowicz, 1968.
        deltaV  = -29.48 + 0.1622*E.TempC - 0.002608*E.TempC^2;
        #               Millero, 1983 has:
        #                 'deltaV = -28.56 + .1211*TempCi - .000321*TempCi*TempCi
        #               Millero, 1992 has:
        #                 'deltaV = -29.48 + .1622*TempCi + .295*(Sali - 34.8)
        #               Millero, 1995 has:
        #                 'deltaV = -29.48 - .1622*TempCi - .002608*TempCi*TempCi
        #                 'deltaV = deltaV + .295*(Sali - 34.8); # Millero, 1979
        Kappa   = -2.84/1000; # Millero, 1979
        #               Millero, 1992 and Millero, 1995 also have this.
        #                 'Kappa = Kappa + .354*(Sali - 34.8)/1000: # Millero,1979
        #               Millero, 1983 has:
        #                 'Kappa = (-3 + .0427*TempCi)/1000
        # lnKBfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    end

    
    # return (KB = Kdata(KB, deltaV, Kappa),)
    return (KB = lnKdata(lnKB, Kfac, deltaV, Kappa),)
end


function getK_KW(E, WhichKs)
    # CalculateKW:

    if WhichKs == Val(7)
            # Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            lnKW = (148.9802 - 13847.26/E.TempK - 23.6521*E.logTempK +
                (-79.2447 + 3298.72/E.TempK + 12.0408*E.logTempK)*
                E.sqrSal - 0.019813*E.Sal)
            # KW = exp(lnKW); # this is on the SWS pH scale in (mol/kg-SW)^2
    elseif WhichKs == Val(8)
            # Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
            # refit data of Harned and Owen, The Physical Chemistry of
            # Electrolyte Solutions, 1958
            lnKW = 148.9802 - 13847.26/E.TempK - 23.6521*E.logTempK
            # KW = exp(lnKW); # this is on the SWS pH scale in (mol/kg-SW)^2
    elseif WhichKs == Val(6)
            # KW = 0.0 # GEOSECS doesn't include OH effects
            lnKW = -Inf
    else
        # Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
        # his check value of 1.6 umol/kg-SW should be 6.2
        lnKW = (148.9802 - 13847.26/E.TempK - 23.6521*E.logTempK +
            (-5.977 + 118.67/E.TempK + 1.0495*E.logTempK)*
            E.sqrSal - 0.01615*E.Sal)
        # KW = exp(lnKW); # this is on the SWS pH scale in (mol/kg-SW)^2        
    end

    # CorrectKWForPressure:

    if WhichKs == Val(8)
        # PressureEffectsOnKWinFreshWater:
        #               This is from Millero, 1983.
        deltaV  =  -25.6 + 0.2324*E.TempC - 0.0036246*E.TempC^2;
        Kappa   = (-7.33 + 0.1368*E.TempC - 0.001233 *E.TempC^2)/1000;
        # lnKWfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
        #               NOTE the temperature dependence of KappaK1 and KappaKW
        #               for fresh water in Millero, 1983 are the same.
    else
        # GEOSECS doesn't include OH term, so this won't matter.
        # Peng et al didn't include pressure, but here I assume that the KW correction
        #       is the same as for the other seawater cases.
        # PressureEffectsOnKW:
        #               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
        deltaV  = -20.02 + 0.1119*E.TempC - 0.001409*E.TempC^2;
        #               Millero, 1992 and Millero, 1995 have:
        Kappa   = (-5.13 + 0.0794*E.TempC)/1000; # Millero, 1983
        #               Millero, 1995 has this too, but Millero, 1992 is different.
        # lnKWfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
        #               Millero, 1979 does not list values for these.
    end

    # return (KW = Kdata(KW, deltaV, Kappa),)
    return (KW = lnKdata(lnKW, 1.0, deltaV, Kappa),)
end

function getK_KP(E, WhichKs, pHfacs)
    # CalculateKP1KP2KP3:

    if WhichKs == Val(7)
        # KP1 = 0.02*1.0
        lnKP1 = log(0.02*1.0)
        KP1fac = 1.0
        # Peng et al don't include the contribution from this term,
        # but it is so small it doesn't contribute. It needs to be
        # kept so that the routines work ok.
        # KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
        # Limnology and Oceanography 12:243-252, 1967:
        # these are only for sals 33 to 36 and are on the NBS scale
        lnKP2 = -9.039 - 1450 /E.TempK
        KP2fac = 1.0/pHfacs.fH
        #KP2 = (exp(-9.039 - 1450 /E.TempK) # this is on the NBS scale
        #    /pHfacs.fH)                          # convert to SWS scale
        lnKP3 = 4.466 - 7276 /E.TempK
        KP3fac = 1.0/pHfacs.fH
        # KP3 = (exp(4.466 - 7276 /E.TempK)  # this is on the NBS scale
        #    /pHfacs.fH)                          # convert to SWS scale
        
    elseif WhichKs == Val(6) || WhichKs == Val(8)
        # KP1 = 0.0; KP2 = 0.0; KP3 =0.0;
        lnKP1 = -Inf; lnKP2 = -Inf; lnKP3 = -Inf
        KP1fac = 0.0; KP2fac = 0.0; KP3fac = 0.0
        # Neither the GEOSECS choice nor the freshwater choice
        # include contributions from phosphate or silicate.
    else
        # Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
        # KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
        lnKP1 = -4576.752/E.TempK + 115.54 - 18.453*E.logTempK + (-106.736/E.TempK +
            0.69171)*E.sqrSal + (-0.65643/E.TempK - 0.01844)*E.Sal;
        # KP1 = exp(lnKP1);
        lnKP2 = -8814.715/E.TempK + 172.1033 - 27.927*E.logTempK + (-160.34/E.TempK +
            1.3566)*E.sqrSal + (0.37335/E.TempK - 0.05778)*E.Sal;
        #KP2 = exp(lnKP2);
        lnKP3 = -3070.75/E.TempK - 18.126 + (17.27039/E.TempK + 2.81197)*E.sqrSal +
            (-44.99486/E.TempK - 0.09984)*E.Sal
        #KP3 = exp(lnKP3);
        KP1fac = 1.0; KP2fac = 1.0; KP3fac = 1.0
    end

    # CorrectKP1KP2KP3ForPressure:
    # These corrections don't matter for the GEOSECS choice (WhichKs# = 6) and
    #       the freshwater choice (WhichKs# = 8). For the Peng choice I assume
    #       that they are the same as for the other choices (WhichKs# = 1 to 5).
    # The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
    #       same as Millero, 1983.
    # PressureEffectsOnKP1:
    deltaV = -14.51 + 0.1211*E.TempC - 0.000321*E.TempC^2;
    Kappa  = (-2.67 + 0.0427*E.TempC)/1000;
    # lnKP1fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    # KP1d = Kdata(KP1, deltaV, Kappa)
    KP1d = lnKdata(lnKP1, KP1fac, deltaV, Kappa)

    # PressureEffectsOnKP2:
    deltaV = -23.12 + 0.1758*E.TempC - 0.002647*E.TempC^2;
    Kappa  = (-5.15 + 0.09  *E.TempC)/1000;
    # lnKP2fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    # KP2d = Kdata(KP2, deltaV, Kappa)
    KP2d = lnKdata(lnKP2, KP2fac, deltaV, Kappa)
    # PressureEffectsOnKP3:
    deltaV = -26.57 + 0.202 *E.TempC - 0.003042*E.TempC^2;
    Kappa  = (-4.08 + 0.0714*E.TempC)/1000;
    # lnKP3fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    # KP3d = Kdata(KP3, deltaV, Kappa)
    KP3d = lnKdata(lnKP3, KP3fac, deltaV, Kappa)

    return (KP1=KP1d, KP2=KP2d, KP3=KP3d)
end

function getK_KSi(E, WhichKs, pHfacs)
    
    # CalculateKSi:

    if WhichKs == Val(7)
        # Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
        # The Chemical Society (London), Special Publ. 17:751, 1964:
        # KSi = (0.0000000004              # this is on the NBS scale
        #    /pHfacs.fH);                          # convert to SWS scale
        lnKSi = log(0.0000000004)
        KSifac = 1.0/pHfacs.fH
    elseif WhichKs == Val(6) || WhichKs == Val(8)
        # KSi = 0.0;
        lnKSi = -Inf
        KSifac = 0.0
        # Neither the GEOSECS choice nor the freshwater choice
        # include contributions from phosphate or silicate.
    else
        # Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
        # KSi was given on the SWS pH scale in molal units. 
        lnKSi = -8904.2/E.TempK + 117.4 - 19.334*E.logTempK + (-458.79/E.TempK +
            3.5913)*sqrt(E.IonS) + (188.74/E.TempK - 1.5998)*E.IonS +
            (-12.1652/E.TempK + 0.07871)*E.IonS^2
        # KSi = (exp(lnKSi)                # this is on the SWS pH scale in mol/kg-H2O
        #    *E.volkg)        # convert to mol/kg-SW
        KSifac = E.volkg
    end

    # PressureEffectsOnKSi:
    #  The only mention of this is Millero, 1995 where it is stated that the
    #    values have been estimated from the values of boric acid. HOWEVER,
    #    there is no listing of the values in the table.
    #    I used the values for boric acid from above.
    deltaV = -29.48 + 0.1622*E.TempC - 0.002608*E.TempC^2;
    Kappa  = -2.84/1000;
    # lnKSifac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;

    # return (KSi = Kdata(KSi, deltaV, Kappa),)
    return (KSi = lnKdata(lnKSi, KSifac, deltaV, Kappa),)
end    

"CalculateK1K2 - split into separate methods that dispatch on type of WhichKs"
function getK_K1K2(E, WhichKs, pHfacs)
    error("invalid WhichKs=", WhichKs)
end

function getK_K1K2(E, ::Val{1}, pHfacs)
    # CalculateK1K2:
            
    # ROY et al, Marine Chemistry, 44:249-267, 1993
    # (see also: Erratum, Marine Chemistry 45:337, 1994
    # and Erratum, Marine Chemistry 52:183, 1996)
    # Typo: in the abstract on p. 249: in the eq. for lnK1* the
    # last term should have S raised to the power 1.5.
    # They claim standard deviations (p. 254) of the fits as
    # .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
    # They also claim (p. 258) 2s precisions of .004 in pK1 and
    # .006 in pK2. These are consistent, but Andrew Dickson
    # (personal communication) obtained an rms deviation of about
    # .004 in pK1 and .003 in pK2. This would be a 2s precision
    # of about 2% in K1 and 1.5# in K2.
    # T:  0-45  S:  5-45. Total Scale. Artificial sewater.
    # This is eq. 29 on p. 254 and what they use in their abstract:
    lnK1 = (2.83655 - 2307.1266/E.TempK - 1.5529413*E.logTempK +
        (-0.20760841 - 4.0484/E.TempK)*E.sqrSal + 0.08468345*E.Sal -
        0.00654208*E.sqrSal*E.Sal)
    K1 = (exp(lnK1)               # this is on the total pH scale in mol/kg-H2O
        *E.volkg   # convert to mol/kg-SW
        /pHfacs.SWStoTOT)                 # convert to SWS pH scale
    # This is eq. 30 on p. 254 and what they use in their abstract:
    lnK2 = -9.226508 - 3351.6106/E.TempK - 0.2005743*E.logTempK +
        (-0.106901773 - 23.9722/E.TempK)*E.sqrSal + 0.1130822*E.Sal -
        0.00846934*E.sqrSal*E.Sal;
    K2 = (exp(lnK2)               # this is on the total pH scale in mol/kg-H2O
        *E.volkg   # convert to mol/kg-SW
        /pHfacs.SWStoTOT);                 # convert to SWS pH scale

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{2}, pHfacs)

    # GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
    # The 2s precision in pK1 is .011, or 2.5% in K1.
    # The 2s precision in pK2 is .02, or 4.5% in K2.
    # This is in Table 5 on p. 1652 and what they use in the abstract:
    pK1 = (812.27/E.TempK + 3.356 - 0.00171*E.Sal*E.logTempK
        + 0.000091*E.Sal^2)
    K1 = 10  ^(-pK1); # this is on the SWS pH scale in mol/kg-SW
    #
    # This is in Table 5 on p. 1652 and what they use in the abstract:
    pK2 = (1450.87/E.TempK + 4.604 - 0.00385*E.Sal*E.logTempK
        + 0.000182*E.Sal^2)
    K2 = 10 ^(-pK2); # this is on the SWS pH scale in mol/kg-SW
    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{3}, pHfacs)
    # HANSSON refit BY DICKSON AND MILLERO
    # Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    # (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    # refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
    # and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
    # on the SWS pH scale in mol/kg-SW.
    # Hansson gave his results on the Total scale (he called it
    # the seawater scale) and in mol/kg-SW.
    # Typo in DM on p. 1739 in Table 4: the equation for pK2*
    # for Hansson should have a .000132 *S^2
    # instead of a .000116 *S^2.
    # The 2s precision in pK1 is .013, or 3% in K1.
    # The 2s precision in pK2 is .017, or 4.1% in K2.
    # This is from Table 4 on p. 1739.
    pK1 = 851.4/E.TempK + 3.237 - 0.0106*E.Sal + 0.000105*E.Sal^2;
    K1 = 10 ^(-pK1); # this is on the SWS pH scale in mol/kg-SW
    #
    # This is from Table 4 on p. 1739.
    pK2 = (-3885.4/E.TempK + 125.844 - 18.141*E.logTempK
        - 0.0192*E.Sal + 0.000132*E.Sal^2)
    K2 = 10 ^(-pK2); # this is on the SWS pH scale in mol/kg-SW
    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{4}, pHfacs)
    # MEHRBACH refit BY DICKSON AND MILLERO
    # Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    # (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    # refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
    # on the SWS pH scale in mol/kg-SW.
    # Mehrbach et al gave results on the NBS scale.
    # The 2s precision in pK1 is .011, or 2.6% in K1.
    # The 2s precision in pK2 is .020, or 4.6% in K2.
    # Valid for salinity 20-40.
    # This is in Table 4 on p. 1739.
    pK1 = (3670.7/E.TempK - 62.008 + 9.7944*E.logTempK
        - 0.0118*E.Sal + 0.000116*E.Sal^2)
    K1 = 10 ^(-pK1); # this is on the SWS pH scale in mol/kg-SW
    #
    # This is in Table 4 on p. 1739.
    pK2 = 1394.7/E.TempK + 4.777 - 0.0184*E.Sal + 0.000118*E.Sal^2;
    K2 = 10 ^(-pK2); # this is on the SWS pH scale in mol/kg-SW

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{5}, pHfacs)
    # HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
    # Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
    # (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    # refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
    # Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
    # and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
    # on the SWS pH scale in mol/kg-SW.
    # Typo in DM on p. 1740 in Table 5: the second equation
    # should be pK2* =, not pK1* =.
    # The 2s precision in pK1 is .017, or 4% in K1.
    # The 2s precision in pK2 is .026, or 6% in K2.
    # Valid for salinity 20-40.
    # This is in Table 5 on p. 1740.
    pK1 = 845 /E.TempK + 3.248 - 0.0098*E.Sal + 0.000087*E.Sal^2
    K1 = 10 ^(-pK1); # this is on the SWS pH scale in mol/kg-SW
    #
    # This is in Table 5 on p. 1740.
    pK2 = 1377.3/E.TempK + 4.824 - 0.0185*E.Sal + 0.000122*E.Sal^2;
    K2 = 10 ^(-pK2); # this is on the SWS pH scale in mol/kg-SW

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Union{Val{6}, Val{7}}, pHfacs)
    # GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
    # Limnology and Oceanography, 18(6):897-907, 1973.
    # I.e., these are the original Mehrbach dissociation constants.
    # The 2s precision in pK1 is .005, or 1.2% in K1.
    # The 2s precision in pK2 is .008, or 2% in K2.
    pK1 = (- 13.7201 + 0.031334*E.TempK + 3235.76/E.TempK
        + 1.3e-5*E.Sal*E.TempK - 0.1032*E.Sal^0.5)
    K1 = (10 ^(-pK1)         # this is on the NBS scale
        /pHfacs.fH)                     # convert to SWS scale
    pK2 = (5371.9645 + 1.671221*E.TempK + 0.22913*E.Sal + 18.3802*log10(E.Sal)
        - 128375.28/E.TempK - 2194.3055*log10(E.TempK) - 8.0944e-4*E.Sal*E.TempK
        - 5617.11*log10(E.Sal)/E.TempK + 2.136*E.Sal/E.TempK) # pK2 is not defined for Sal=0, since log10(0)=-inf
    K2 = (10 ^(-pK2)         # this is on the NBS scale
        /pHfacs.fH);                     # convert to SWS scale

    #               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
    #               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
    #               Culberson and Pytkowicz, L and O 13:403-417, 1968:
    #               but the fits are the same as those in
    #               Edmond and Gieskes, GCA, 34:1261-1291, 1970
    #               who in turn quote Li, personal communication
    # lnK1fac = (24.2 - 0.085*E.TempC)*Pbar/RT;
    deltaV = -(24.2 - 0.085*E.TempC)
    K1d = Kdata(K1, deltaV, 0.0)
    # lnK2fac = (16.4 - 0.04 *E.TempC)*Pbar/RT;
    deltaV = -(16.4 - 0.04 *E.TempC)
    K2d = Kdata(K2, deltaV, 0.0)

    return(K1=K1d, K2=K2d)
end

function getK_K1K2(E, ::Val{8}, pHfacs)
    # PURE WATER CASE
    # Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
    # K1 from refit data from Harned and Davis,
    # J American Chemical Society, 65:2030-2037, 1943.
    # K2 from refit data from Harned and Scholes,
    # J American Chemical Society, 43:1706-1709, 1941.
    # This is only to be used for Sal=0 water (note the absence of S in the below formulations)
    # These are the thermodynamic Constants:
    lnK1 = 290.9097 - 14554.21/E.TempK - 45.0575*E.logTempK;
    K1 = exp(lnK1);
    lnK2 = 207.6548 - 11843.79/E.TempK - 33.6485*E.logTempK;
    K2 = exp(lnK2);

    #***PressureEffectsOnK1inFreshWater:
    #               This is from Millero, 1983.
    deltaV  = -30.54 + 0.1849 *E.TempC - 0.0023366*E.TempC^2;
    Kappa   = (-6.22 + 0.1368 *E.TempC - 0.001233 *E.TempC^2)/1000;
    K1d = Kdata(K1, deltaV, Kappa)
    # lnK1fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    #***PressureEffectsOnK2inFreshWater:
    #               This is from Millero, 1983.
    deltaV  = -29.81 + 0.115*E.TempC - 0.001816*E.TempC^2;
    Kappa   = (-5.74 + 0.093*E.TempC - 0.001896*E.TempC^2)/1000;
    # lnK2fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    K2d = Kdata(K2, deltaV, Kappa)

    return(K1=K1d, K2=K2d)
end

function getK_K1K2(E, ::Val{9}, pHfacs)    
    # From Cai and Wang 1998, for estuarine use.
    # Data used in this work is from:
    # K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
    # K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
    # Sigma of residuals between fits and above data: �0.015, +0.040 for K1 and K2, respectively.
    # Sal 0-40, Temp 0.2-30
    # Limnol. Oceanogr. 43(4) (1998) 657-668
    # On the NBS scale
    # Their check values for F1 don't work out, not sure if this was correctly published...
    F1 = 200.1/E.TempK + 0.3220;
    pK1 = 3404.71/E.TempK + 0.032786*E.TempK - 14.8435 - 0.071692*F1*E.Sal^0.5 + 0.0021487*E.Sal;
    K1  = (10 ^-pK1         # this is on the NBS scale
        /pHfacs.fH)                    # convert to SWS scale (uncertain at low Sal due to junction potential);
    F2 = -129.24/E.TempK + 1.4381;
    pK2 = 2902.39/E.TempK + 0.02379*E.TempK - 6.4980 - 0.3191*F2*E.Sal^0.5 + 0.0198*E.Sal;
    K2  = (10 ^-pK2         # this is on the NBS scale
        /pHfacs.fH)                    # convert to SWS scale (uncertain at low Sal due to junction potential);

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{10}, pHfacs)    
        
    # From Lueker, Dickson, Keeling, 2000
    # This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work.
    # Mar. Chem. 70 (2000) 105-119
    # Total scale and kg-sw
    pK1 = 3633.86/E.TempK-61.2172+9.6777*E.logTempK-0.011555*E.Sal+0.0001152*E.Sal^2;
    #K1  = 10^-pK1...           # this is on the total pH scale in mol/kg-SW
    # K1 = (exp(-pK1*log(10))   # Matlab 10,^ is slow
    #    /pHfacs.SWStoTOT);                # convert to SWS pH scale
    lnK1 = -pK1*log(10)
    facK = 1.0/pHfacs.SWStoTOT
    
    pK2 = 471.78/E.TempK+25.929 -3.16967*E.logTempK-0.01781 *E.Sal+0.0001122*E.Sal^2;
    # K2  = 10^-pK2...           # this is on the total pH scale in mol/kg-SW
    #K2 = (exp(-pK2*log(10))    # Matlab 10^ is slow
    #    /pHfacs.SWStoTOT)                # convert to SWS pH scale
    lnK2 = -pK2*log(10)
   
    deltaV1, Kappa1, deltaV2, Kappa2 = deltaVKappaK1K2_default(E)
    return (K1 = lnKdata(lnK1, facK, deltaV1, Kappa1),
            K2 = lnKdata(lnK2, facK, deltaV2, Kappa2),)
end

function getK_K1K2(E, ::Val{11}, pHfacs)    
    # Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
    # sigma for pK1 is reported to be 0.0056
    # sigma for pK2 is reported to be 0.010
    # This is from the abstract and pages 2536-2537
    pK1 =  -43.6977 - 0.0129037*E.Sal + 1.364e-4*E.Sal^2 + 2885.378/E.TempK +  7.045159*E.logTempK;
    pK2 = (-452.0940 + 13.142162*E.Sal - 8.101e-4*E.Sal^2 + 21263.61/E.TempK + 68.483143*E.logTempK
        + (-581.4428*E.Sal + 0.259601*E.Sal^2)/E.TempK - 1.967035*E.Sal*E.logTempK)
    K1 = 10 ^-pK1; # this is on the SWS pH scale in mol/kg-SW
    K2 = 10 ^-pK2; # this is on the SWS pH scale in mol/kg-SW

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{12}, pHfacs)
    # Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
    # Calculated from overdetermined WOCE-era field measurements
    # sigma for pK1 is reported to be 0.005
    # sigma for pK2 is reported to be 0.008
    # This is from page 1715
    pK1 =  6.359 - 0.00664*E.Sal - 0.01322*E.TempC + 4.989e-5*E.TempC^2;
    pK2 =  9.867 - 0.01314*E.Sal - 0.01904*E.TempC + 2.448e-5*E.TempC^2;
    K1 = 10 ^-pK1; # this is on the SWS pH scale in mol/kg-SW
    K2 = 10 ^-pK2; # this is on the SWS pH scale in mol/kg-SW

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{13}, pHfacs)
    # From Millero 2006 work on pK1 and pK2 from titrations
    # Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
    # S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
    pK1_0 = -126.34048 + 6320.813/E.TempK + 19.568224*E.logTempK;
    A_1   = 13.4191*E.Sal^0.5 + 0.0331*E.Sal - 5.33e-5*E.Sal^2;
    B_1   = -530.123*E.Sal^0.5 - 6.103*E.Sal;
    C_1   = -2.06950*E.Sal^0.5;
    pK1= A_1 + B_1/E.TempK + C_1*E.logTempK + pK1_0; # pK1 sigma = 0.0054
    K1 = 10 ^-(pK1);
    pK2_0= -90.18333 + 5143.692/E.TempK + 14.613358*E.logTempK;
    A_2   = 21.0894*E.Sal^0.5 + 0.1248*E.Sal - 3.687e-4*E.Sal^2;
    B_2   = -772.483*E.Sal^0.5 - 20.051*E.Sal;
    C_2   = -3.3336*E.Sal^0.5;
    pK2= A_2 + B_2/E.TempK + C_2*E.logTempK + pK2_0; #pK2 sigma = 0.011
    K2 = 10 ^-(pK2);

    return pressureK1K2_default(E, K1, K2)
end

function getK_K1K2(E, ::Val{14}, pHfacs)
    # From Millero, 2010, also for estuarine use.
    # Marine and Freshwater Research, v. 61, p. 139�142.
    # Fits through compilation of real seawater titration results:
    # Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
    # Constants for K's on the SWS;
    # This is from page 141
    pK10 = -126.34048 + 6320.813/E.TempK + 19.568224*E.logTempK;
    # This is from their table 2, page 140.
    A1 = 13.4038*E.Sal^0.5 + 0.03206*E.Sal - 5.242e-5*E.Sal^2;
    B1 = -530.659*E.Sal^0.5 - 5.8210*E.Sal;
    C1 = -2.0664*E.Sal^0.5;
    pK1 = pK10 + A1 + B1/E.TempK + C1*E.logTempK;
    K1 = 10 ^-pK1;
    # This is from page 141
    pK20 =  -90.18333 + 5143.692/E.TempK + 14.613358*E.logTempK;
    # This is from their table 3, page 140.
    A2 = 21.3728*E.Sal^0.5 + 0.1218*E.Sal - 3.688e-4*E.Sal^2;
    B2 = -788.289*E.Sal^0.5 - 19.189*E.Sal;
    C2 = -3.374*E.Sal^0.5;
    pK2 = pK20 + A2 + B2/E.TempK + C2*E.logTempK;
    K2 = 10 ^-pK2;
    return pressureK1K2_default(E, K1, K2)
end

function deltaVKappaK1K2_default(E)
    #***PressureEffectsOnK1:
    #               These are from Millero, 1995.
    #               They are the same as Millero, 1979 and Millero, 1992.
    #               They are from data of Culberson and Pytkowicz, 1968.
    deltaV1  = -25.5 + 0.1271*E.TempC;
    #                 'deltaV = deltaV - .151*(Sali - 34.8); # Millero, 1979
    Kappa1   = (-3.08 + 0.0877*E.TempC)/1000;
    #                 'Kappa = Kappa  - .578*(Sali - 34.8)/1000.; # Millero, 1979
    #lnK1fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    #               The fits given in Millero, 1983 are somewhat different.
    
    #***PressureEffectsOnK2:
    #               These are from Millero, 1995.
    #               They are the same as Millero, 1979 and Millero, 1992.
    #               They are from data of Culberson and Pytkowicz, 1968.
    deltaV2  = -15.82 - 0.0219*E.TempC;
    #                  'deltaV = deltaV + .321*(Sali - 34.8); # Millero, 1979
    Kappa2   = (1.13 - 0.1475*E.TempC)/1000;
    #                 'Kappa = Kappa - .314*(Sali - 34.8)/1000: # Millero, 1979
    # lnK2fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    #               The fit given in Millero, 1983 is different.
    #               Not by a lot for deltaV, but by much for Kappa. #
    

    return(deltaV1, Kappa1, deltaV2, Kappa2)
end

function pressureK1K2_default(E, K1, K2)
    deltaV1, Kappa1, deltaV2, Kappa2 = deltaVKappaK1K2_default(E)
    
    K1d = Kdata(K1, deltaV1, Kappa1)    
    K2d = Kdata(K2, deltaV2, Kappa2)

    return(K1=K1d, K2=K2d)
end



function getK_KH2S(E, pHfacs)
    # SD 2014-7-28 Calculate KH2S
    # Original data from
    # Millero, F. J., Plese, T., & Fernandez, M. (1988). Limnol. Oceanogr, 33(2), 269�274.
    # number here from
    # Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507�546, doi:10.1007/s10498-009-9084-1
    lnKH2S= (225.838 + 0.3449*E.sqrSal - 0.0274*E.Sal) - 13275.3/E.TempK - 34.6435*E.logTempK;
    # KH2S = exp(lnKH2S)/pHfacs.SWStoTOT;         # convert to SW pH scale
    facKH2S = 1.0/pHfacs.SWStoTOT

    # SD 2014-7-28 Correct KH2S for pressure
    # From Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507�546, doi:10.1007/s10498-009-9084-1
    #         a0      a1               a2
    deltaV = -14.80 + 0.0020*E.TempC - 0.4000e-3*E.TempC^2;
    #         b0      b1               (b2)
    Kappa  =  (2.89   + 0.0540*E.TempC)/1000;
    # lnKH2Sfac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;

    return (KH2S=lnKdata(lnKH2S, facKH2S, deltaV, Kappa),)
end

function getK_KNH3(E, pHfacs)
     # SD 2015-01-03 Calculate KNH3
    # numbers from 
    # Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507�546, doi:10.1007/s10498-009-9084-1
    lnKNH4= (-0.25444 + 0.46532*E.sqrSal - 0.01992*E.Sal) + (-6285.33-123.7184*E.sqrSal+3.17556*E.Sal)/E.TempK + 0.0001635*E.TempK;
    # KNH4 = exp(lnKNH4);          # parameters are for SW pH scale

    
    # SD 2015-01-03 Correct KNH4 for pressure
    # From Hofmann AF etal  (2010)  Aquatic Geochemistry 16, 507�546, doi:10.1007/s10498-009-9084-1
    #         a0      a1               a2
    deltaV = -26.43 + 0.0889*E.TempC - 0.9050e-3*E.TempC^2;
    #         b0      b1               (b2)
    Kappa  =  (-5.03   + 0.0814*E.TempC)/1000;
    # lnKNH4fac = (-deltaV + 0.5*Kappa*Pbar)*Pbar/RT;
    return (KNH3 = lnKdata(lnKNH4, 1.0, deltaV, Kappa),)
end


function getK_KCaAr(E, WhichKs, pHfacs)
    # add Ca solubility constants, for conditions defined by C
    # modified version of CO2SYS CaSolubility to replace global variables

    #function varargout=CaSolubility(Sal, TempC, Pdbar, TC, pH)
    #global K1 K2 TempK logTempK sqrSal Pbar RT WhichKs
    # '***********************************************************************
    # ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
    # ' Inputs: WhichKs#, Sal, TempCi, Pdbari, TCi, pHi, K1, K2
    # ' Outputs: OmegaCa, OmegaAr
    # ' This calculates omega, the solubility ratio, for calcite and aragonite.
    # ' This is defined by: Omega = [CO3--]*[Ca++]/Ksp,
    # '       where Ksp is the solubility product (either KCa or KAr).
    # '***********************************************************************
    # ' These are from:
    # ' Mucci, Alphonso, The solubility of calcite and aragonite in seawater
    # '       at various salinities, temperatures, and one atmosphere total
    # '       pressure, American Journal of Science 283:781-799, 1983.
    # ' Ingle, S. E., Solubility of calcite in the ocean,
    # '       Marine Chemistry 3:301-319, 1975,
    # ' Millero, Frank, The thermodynamics of the carbonate system in seawater,
    # '       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
    # ' Ingle et al, The solubility of calcite in seawater at atmospheric pressure
    # '       and 35#o salinity, Marine Chemistry 1:295-307, 1973.
    # ' Berner, R. A., The solubility of calcite and aragonite in seawater in
    # '       atmospheric pressure and 34.5#o salinity, American Journal of
    # '       Science 276:713-730, 1976.
    # ' Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
    # ' Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
    # '       boric acid, and the pHi of seawater, Limnology and Oceanography
    # '       13:403-417, 1968.
    # '***********************************************************************

   

    if WhichKs == Val(6) || WhichKs == Val(7)
        #
        # *** CalculateKCaforGEOSECS:
        # Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
        # (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        # but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
        KCa = (0.0000001*(-34.452 - 39.866*E.Sal^(1. /3) +
            110.21*log(E.Sal)/log(10) - 0.0000075752*E.TempK^2))
        # this is in (mol/kg-SW)^2
        #
        # *** CalculateKArforGEOSECS:
        # Berner, R. A., American Journal of Science 276:713-730, 1976:
        # (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
        KAr = 1.45*KCa;# ' this is in (mol/kg-SW)^2
        # Berner (p. 722) states that he uses 1.48.
        # It appears that 1.45 was used in the GEOSECS calculations
        #
        # *** CalculatePressureEffectsOnKCaKArGEOSECS:
        # Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
        # (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
        # but their paper is not even on this topic).
        # The fits appears to be new in the GEOSECS report.
        # I can't find them anywhere else.
        # KCa = KCa*exp((36   - 0.2 *C.TempC)*C.Pbar/C.RT);
        deltaV = -(36   - 0.2 *E.TempC)
        KCad = Kdata(KCa, deltaV, 0.0)
        # KAr = KAr*exp((33.3 - 0.22*C.TempC)*C.Pbar/C.RT);
        deltaV = -(33.3 - 0.22*E.TempC)
        KArd = Kdata(KAr, deltaV, 0.0)         
    else
        # CalciteSolubility:
        # '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKCa = -171.9065 - 0.077993*E.TempK + 2839.319/E.TempK;
        logKCa = logKCa + 71.595*E.logTempK/log(10);
        logKCa = logKCa + (-0.77712 + 0.0028426*E.TempK + 178.34/E.TempK)*E.sqrSal;
        logKCa = logKCa - 0.07711*E.Sal + 0.0041249*E.sqrSal*E.Sal;
        # '       sd fit = .01 (for Sal part, not part independent of Sal)
        #KCa = 10^(logKCa);# ' this is in (mol/kg-SW)^2
        lnKCa = logKCa*log(10) 
        
        # AragoniteSolubility:
        # '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
        logKAr = -171.945 - 0.077993*E.TempK + 2903.293/E.TempK;
        logKAr = logKAr + 71.595*E.logTempK/log(10);
        logKAr = logKAr + (-0.068393 + 0.0017276*E.TempK + 88.135/E.TempK)*E.sqrSal;
        logKAr = logKAr - 0.10018*E.Sal + 0.0059415*E.sqrSal*E.Sal;
        # '       sd fit = .009 (for Sal part, not part independent of Sal)
        # KAr    = exp(logKAr*log(10)); # 10^(logKAr);# ' this is in (mol/kg-SW)^2
        lnKAr = logKAr*log(10)

        # PressureCorrectionForCalcite:
        # '       Ingle, Marine Chemistry 3:301-319, 1975
        # '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
        # '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
        deltaVKCa = -48.76 + 0.5304*E.TempC;
        KappaKCa  = (-11.76 + 0.3692*E.TempC)/1000;
        # lnKCafac  = (-deltaVKCa + 0.5*KappaKCa*C.Pbar)*C.Pbar/C.RT;
        # KCa       = KCa*exp(lnKCafac);
        KCad = lnKdata(lnKCa, 1.0, deltaVKCa, KappaKCa)
        # PressureCorrectionForAragonite:
        # '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
        # '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
        # '       and 10^3 for Kappa factor)
        deltaVKAr = deltaVKCa + 2.8;
        KappaKAr  = KappaKCa;
        # lnKArfac  = (-deltaVKAr + 0.5*KappaKAr*E.Pbar)*E.Pbar/E.RT;
        # KAr       = KAr*exp(lnKArfac);
        KArd = lnKdata(lnKAr, 1.0, deltaVKAr, KappaKAr)
    end

    # Return as NamedTuple
    return (KCa = KCad, KAr = KArd)

end

function get_FugFac(E, WhichKs)
    # CalculateFugacityConstants:
    if WhichKs == Val(6) || WhichKs == Val(7)
        # GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
        FugFac = 1.0;
    else
        # This assumes that the pressure is at one atmosphere, or close to it.
        # Otherwise, the Pres term in the exponent affects the results.
        #       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
        #       Delta and B in cm3/mol  
        Delta = (57.7 - 0.118*E.TempK);
        b = -1636.75 + 12.0408*E.TempK - 0.0327957*E.TempK^2 + 3.16528*0.00001*E.TempK^3; 
        # For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
        P1atm = 1.01325; # in bar
        FugFac = exp((b + 2 *Delta)*P1atm/E.RT);
    end
    
    return FugFac
end

function get_VPFac(E)
    # CalculateVPFac:
    # Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
    #       seawater, Marine Chemistry 8:347-359, 1980.
    # They fit the data of Goff and Gratch (1946) with the vapor pressure
    #       lowering by sea salt as given by Robinson (1954).
    # This fits the more complicated Goff and Gratch, and Robinson equations
    #       from 273 to 313 deg K and 0 to 40 Sali with a standard error
    #       of .015%, about 5 uatm over this range.
    # This may be on IPTS-29 since they didn't mention the temperature scale,
    #       and the data of Goff and Gratch came before IPTS-48.
    # The references are:
    # Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
    #       to 212 deg F, Transactions of the American Society of Heating and
    #       Ventilating Engineers 52:95-122, 1946.
    # Robinson, Journal of the Marine Biological Association of the U. K.
    #       33:449-455, 1954.
    #       This is eq. 10 on p. 350.
    #       This is in atmospheres.
    # VPWP = exp(24.4543 - 67.4509*(100 /E.TempK) - 4.8489*log(E.TempK/100));
    lnVPWP = 24.4543 - 67.4509*(100 /E.TempK) - 4.8489*(E.logTempK -log(100))
    # VPCorrWP = exp(-0.000544*E.Sal);
    lnVPCorrWP = -0.000544*E.Sal
    # VPSWWP = VPWP*VPCorrWP;
    VPSWWP = exp(lnVPWP + lnVPCorrWP)
    VPFac = 1 - VPSWWP; # this assumes 1 atmosphere

    return VPFac
end
