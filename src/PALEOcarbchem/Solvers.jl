
"""
    ResultNames::Tuple

All available CO2SYS result names
"""
const ResultNames = (
    :pHfree, :H, :OH, :TA, :dTAdpH,
    :TS, :HSO4, 
    :TF, :HF, 
    :TCi, :CO2, :HCO3, :CO3, :CAlk, :fCO2, :pCO2, :xCO2dryinp,
    :TB, :BAlk,
    :TP, :H3PO4, :H2PO4, :HPO4, :PO4, :PengCorrection, :PAlk,
    :TSi, :SiAlk,
    :TH2S, :H2S, :HSAlk,
    :TNH3, :NH4, :NH3Alk,
    :Ca, :OmegaCA, :OmegaAR
)

""" 
    RNameIdx::NamedTuple

Indices of results in a Vector of length([`ResultNames`](@ref))
"""
const RNameIdx = NamedTuple{ResultNames}((1:length(ResultNames)))

"""
    getResultDescription(resultName::AbstractString) -> (resultName, units, description)

Get human-readable description of `resultName`
"""
function getResultDescription(resultName::AbstractString)
    ResultDescriptions = (
        # water
        ("H",           "mol kg-1", "Hydrogen ion (free) concentration"),
        ("OH",          "mol kg-1", "Hydroxyl ion concentration"),
        # sulphate
        ("TS",          "mol kg-1", "Sulphate total concentration"),
        ("HSO4",          "mol kg-1", "HSO4- concentration"),
        # fluoride
        ("TF",          "mol kg-1", "Fluoride total concentration"),
        ("HF",          "mol kg-1", "Fluoride HF concentration"),
        # carbon
        ("TCi",         "mol kg-1", "DIC total concentration"),
        ("CO2",         "mol kg-1", "CO2 concentration"),
        ("HCO3",        "mol kg-1", "bicarbonate ion HCO3- concentration"),
        ("CO3",         "mol kg-1", "carbonate ion CO3-- concentration"),
        ("CAlk",        "mol kg-1", "carbon contribution to total alkalinity (HCO3 + 2*CO3)"),
        ("fCO2",        "atm",      "CO2 fugacity"),
        ("pCO2",        "atm",      "CO2 partial pressure (fugacity corrected)"),
        ("xCO2dryinp",  "",         "mixing ratio of CO2 in dry air at 1 atm (always > pCO2 due to H2O vapour pressure)"),
        # boron
        ("TB",          "mol kg-1", "Boron total concentration"),
        ("BAlk",        "mol kg-1", "Boron B(OH)_4^- concentration (= boron contribution to total alkalinity)"),
        # phosphorus
        ("TP",          "mol kg-1", "Phosphate total concentration"),
        ("H3PO4",       "mol kg-1", "Phosphate H3PO4 concentration"),
        ("H2PO4",       "mol kg-1", "Phosphate H2PO4- concentration"),
        ("HPO4",        "mol kg-1", "Phosphate HPO4-- concentration"),
        ("PO4",         "mol kg-1", "Phosphate PO4--- concentration"),
        ("PAlk",        "mol kg-1", "Phosphate contribution to total alkalinity (-H3PO4 + HPO4 + 2*PO4)"),
        # silica
        ("TSi",         "mol kg-1", "Silicate total concentration"),
        ("SiAlk",       "mol kg-1", "Silicate H3SO4- concentration (= silicate contribution to total alkalinity)"),
        # sulphide
        ("TH2S",        "mol kg-1", "Sulphide total concentration"),
        ("H2S",         "mol kg-1", "Sulphide H2S species concentration"),
        ("HSAlk",       "mol kg-1", "Sulphide HS- concentration (= sulphide contribution to total alkalinity)"),
        # ammonia
        ("TNH3",        "mol kg-1", "Ammonia + ammonium total concentration (NH3 + NH_4^+)"),
        ("NH4",         "mol kg-1", "Ammonium NH4- concentration"),
        ("NH3Alk",      "mol kg-1", "Ammonia NH3 concentration (= ammonia contribution to total alkalinity)"),
        # calcium and carbonate saturation
        ("Ca",          "mol kg-1", "Ca concentration"),
        ("OmegaCA",     "",         "calcite saturation"),
        ("OmegaAR",     "",         "aragonite saturation"),
    )

    for resultdesc in ResultDescriptions
        if resultName == resultdesc[1]
            return resultdesc
        end
    end
    return nothing
end
    


"""
    calculateTAfromTCpHfree!(res, C, Options, concs, pHfree; do_dTAdpH=Val(false)) -> (TA, dTAdpH)

Calculate TAlk, and speciation, given pH and conserved concentrations (total DIC, Si, P, SO4, B, F, ...)

# Returns:
- `TA`:  Total Alk, mol/kg-sw      
- `dTAdpH`: derivative, or NaN if do_dTAdpH=Val(false)
 
# Arguments:
- `res`:  (output) Vector `res` of length length([`ResultNames`](@ref)) with details of TA contributions etc
- `C`: constants from [`calc_constants!`](@ref). NB: must be on Free pH scale.
- `concs::NamedTuple`: (mol kg-sw)  total concentrations for sulphate, fluoride, and each optional component of alkalinity enabled in `C`
- `pHfree`: pH on free scale

# Implementation
modified from Matlab CO2SYS CalculateTAfromTCpH(TCi, pHi) to use Free pH scale
"""
function calculateTAfromTCpHfree!(
    r, C, Options, concs, pHfree; 
    do_dTAdpH=Val(false),
)
    CI = CNameIdx
    RI = RNameIdx # short name

    Options.pHScale == Val(3) || error("unsupported pHScale " , Options.pHScale, " != Val(3)")
         
    r[RI.pHfree] = pHfree

    H         = exp(-pHfree*log(10));  # 10^(-pHfree);  # Matlab 10^ is slow
    r[RI.H]   = H

    # Accumlate components to TA (and dTAdH)

    # Water (always included)
    OH          = C[CI.KW]/H;
    r[RI.OH]    = OH
    TA          = OH - H   # set TA as first component, later components add to TA
    if do_dTAdpH == Val(true)
        dOHdH   = - OH/H;
        dTAdH   = dOHdH - 1.0 # set dTAdH as first component, later components add to TA
    end

    # sulphate 
    if Options.Components.S == Val(true)
        TS      = concs.TS
        HSO4    = TS/(1 + C[CI.KS]/H);# ' since KS is on the free scale
        TA      += -HSO4
        r[RI.TS]= TS; r[RI.HSO4]=HSO4       
        if do_dTAdpH == Val(true)
            dHSO4dH = (C[CI.KS]*HSO4)/(H^2 *(C[CI.KS]/H + 1));
            dTAdH   += -dHSO4dH 
        end
    end

    # fluoride
    if Options.Components.F == Val(true)
        TF      = concs.TF
        HF      = TF/(1 + C[CI.KF]/H);# ' since KF is on the free scale
        TA      += -HF
        r[RI.TF] = TF; r[RI.HF] = HF
        if do_dTAdpH == Val(true)           
            dHFdH   = (C[CI.KF]*HF)/(H^2 *(C[CI.KF]/H + 1))        
            dTAdH   += -dHFdH
        end
    end

    # Inorganic carbon
    if Options.Components.Ci == Val(true)
        TCi = concs.TCi
        # Inorganic carbon
        CDenom      = (C[CI.K1]*H + H*H + C[CI.K1]*C[CI.K2]);
        HCO3      = TCi*C[CI.K1]*H  / CDenom;
        CO3       = TCi*C[CI.K1]*C[CI.K2] / CDenom;
        #CO2       = TCi*H*H /CDenom;
        CO2        =TCi-HCO3-CO3;           # CO2   zero level of protons
        CAlk        = HCO3 + 2*CO3;
                
        # CO2 fugacity CO2SYS CalculatefCO2fromTCpH
        fCO2        = CO2/C[CI.K0]
        # pCO2 CO2SYS L351
        pCO2        = fCO2 / C[CI.FugFac]
        xCO2dryinp  = pCO2/ C[CI.VPFac] #  this assumes pTot = 1 atm
        # NamedTuple to include in output
        r[RI.TCi]=TCi; r[RI.CO2]=CO2; r[RI.HCO3]=HCO3; r[RI.CO3]=CO3; r[RI.CAlk]=CAlk;
        r[RI.fCO2]=fCO2;r[RI.pCO2]=pCO2; r[RI.xCO2dryinp]=xCO2dryinp

        # Total alkalinity
        TA          += CAlk
        if do_dTAdpH==Val(true)
            # derivative
            dCAlkdH = -(C[CI.K1]*TCi*(H^2 + 4*C[CI.K2]*H + C[CI.K1]*C[CI.K2]))/CDenom^2;
            dTAdH   += dCAlkdH 
        end
    end

    # Boron
    if Options.Components.B == Val(true)
        TB        = concs.TB
        BAlk      = TB*C[CI.KB]/(C[CI.KB] + H)
        TA        += BAlk
        r[RI.TB]=TB; r[RI.BAlk]=BAlk
        if do_dTAdpH==Val(true)
            dBAlkdH = - BAlk/(C[CI.KB]+H);
            dTAdH += dBAlkdH
        end
    end

    # Phosphate
    if Options.Components.P == Val(true)
        TP = concs.TP
        PhosBot     = H*H*H + C[CI.KP1]*H*H + C[CI.KP1]*C[CI.KP2]*H + C[CI.KP1]*C[CI.KP2]*C[CI.KP3]
        H3PO4     = TP*H*H*H/PhosBot;        # neutral    Alk -1
        H2PO4     = TP*C[CI.KP1]*H*H/PhosBot;       # 1-         zero level of protons
        HPO4      = TP*C[CI.KP1]*C[CI.KP2]*H/PhosBot;    # 2-         Alk +1
        PO4       = TP*C[CI.KP1]*C[CI.KP2]*C[CI.KP3]/PhosBot;  # 3-         Alk +2
        PAlk      = - H3PO4 + HPO4 + 2*PO4;     
        # The vector 'PengCorrection' is used to modify the value of TA, for those
        # cases where WhichKs==7, since PAlk(Peng) = PAlk(Dickson) + TP.
        # Thus, PengCorrection is 0 for all cases where WhichKs is not 7
        if Options.WhichKs == Val(7)
            PengCorrection=TP
        else
            PengCorrection=zero(TP)
        end
        PAlk    += PengCorrection;   # CO2SYS L336
        TA      += PAlk

        r[RI.TP]=TP; r[RI.H3PO4]=H3PO4; r[RI.H2PO4]=H2PO4; r[RI.HPO4]=HPO4; r[RI.PO4]=PO4; 
        r[RI.PengCorrection]=PengCorrection; r[RI.PAlk]=PAlk
        if do_dTAdpH == Val(true)
            dPAlkdH = (-(C[CI.KP1]*TP*(4*C[CI.KP2]*H*H*H + (H^2)^2 + 
                C[CI.KP1]*C[CI.KP2]^2 *C[CI.KP3] + C[CI.KP1]*C[CI.KP2]*H^2 +
                9*C[CI.KP2]*C[CI.KP3]*H^2 + 4*C[CI.KP1]*C[CI.KP2]*C[CI.KP3]*H))
                /(H*H*H + C[CI.KP1]*H^2 + C[CI.KP1]*C[CI.KP2]*H + C[CI.KP1]*C[CI.KP2]*C[CI.KP3])^2)
            dTAdH   += dPAlkdH
        end
    end

    # Silicate
    if Options.Components.Si == Val(true)
        TSi = concs.TSi
        SiAlk     = TSi*C[CI.KSi]/(C[CI.KSi] + H)
        TA      += SiAlk
        if do_dTAdpH == Val(true)
            dSiAlkdH = - SiAlk/(C[CI.KSi] + H)
            dTAdH   += dSiAlkdH
        end
        r[RI.TSi]=TSi; r[RI.SiAlk]=SiAlk
    end

    # Sulphide
    if Options.Components.H2S == Val(true)
        TH2S = concs.TH2S
        HSAlk     = TH2S*C[CI.KH2S]/(C[CI.KH2S] + H)
        H2S       = TH2S - HSAlk
        TA       += HSAlk
        if do_dTAdpH == Val(true)
            dHSAlkdH = - HSAlk/(C[CI.KH2S] + H);
            dTAdH   += dHSAlkdH
        end
        r[RI.TH2S]=TH2S; r[RI.HSAlk]=HSAlk; r[RI.H2S]=H2S
    end 

    # Ammonia
    if Options.Components.NH3 == Val(true)
        TNH3      = concs.TNH3
        NH3Alk    = TNH3*C[CI.KNH3]/(C[CI.KNH3] + H)
        NH4       = TNH3 - NH3Alk
        TA        += NH3Alk
        if do_dTAdpH == Val(true)
            dNH3AlkdH = - NH3Alk/(C[CI.KNH3] + H)
            dTAdH   += dNH3AlkdH
        end
        r[RI.TNH3]=TNH3; r[RI.NH3Alk]=NH3Alk; r[RI.NH4]=NH4
    end

    # Carbonate saturation (not a contribution to alkalinity !)
    if Options.Components.Omega == Val(true)
        Ca      = concs.Ca
        (OmegaCA, OmegaAR) = calculateOmega(C, CO3, Ca)
        r[RI.Ca]=Ca; r[RI.OmegaCA]=OmegaCA; r[RI.OmegaAR]=OmegaAR
    end

    r[RI.TA] = TA

    if do_dTAdpH == Val(true)        
        r[RI.dTAdpH] = -log(10)*H*dTAdH
    else
        r[RI.dTAdpH] = zero(H)*NaN
    end

    return (TA, r[RI.dTAdpH])
end 

"""
    calculatepHfromTATC!(
        res, C, Options, conc_TAx, concs; 
        pHstart=8.0, pHTol=100*eps()
    ) -> (pHfree, steps)

Call [`calculateTAfromTCpHfree!`](@ref) to iteratively solve for pH given (i) a starting value and pH tolerance, 
and (ii) Alk and conserved concentrations (total DIC, Si, P, SO4, B, F, ...)

Intended for use eg in an ocean model to enable a single Newton-Raphson step each model timestep
(ie pHstart from previous value in that grid cell, pHtol set to some large number)

# Returns:
- `pHfree`: pH on free scale
- `steps`: number of Newton iterations

# Arguments:
- `res`:  (output) Vector `res` of length length([`ResultNames`](@ref)) with details of TA contributions etc
- `C`: constants from [`calc_constants!`](@ref). NB: must be on Free pH scale.
- `conc_TAx`: total Alk, mol kg-sw
- `concs`: other input total concentrations, see [`calculateTAfromTCpHfree!`](@ref)
- `pHstart`: starting value (free pH scale)
- `pHtol`: tolerance (accuracy required)

# Implementation
Modified from CO2SYS  SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis.
NB: recoded here to use [`calculateTAfromTCpHfree!`](@ref), hence works on free pH scale.

Units: mol / kg-sw 
"""
function calculatepHfromTATC!(
    res, C, Options, conc_TAx, concs;
    pHstart=8.0,
    pHTol=100*eps()
)

    pHfree      = pHstart
    deltapH     = pHTol+1
    maxsteps    = 1000
    
    step = 0 # count calls to CalculateTA...

    TAx = conc_TAx
    # first step
    (TA, dTAdpH) = calculateTAfromTCpHfree!(
        res, C, Options, concs, pHfree,
        do_dTAdpH=Val(true),
    )
    step += 1
    Residual  = TAx - TA;           
    deltapH   = Residual/dTAdpH # this is Newton's method      
    # println("step=", step, "  deltapH=", deltapH)
    

    # iterate to convergence
    # use any() to allow SIMD datatypes 
    while any(abs(deltapH) > pHTol) && (step < maxsteps)
       # to keep the jump from being too big;
        deltapH = max(deltapH, -1.0)
        deltapH = min(deltapH, 1.0)
        pHfree       = pHfree + deltapH; # Is on the same scale as K1 and K2 were calculated...
       
        (TA, dTAdpH) = calculateTAfromTCpHfree!(
            res, C, Options, concs, pHfree,
            do_dTAdpH=Val(true),
        )
        step += 1
        Residual  = TAx - TA;           
        deltapH   = Residual/dTAdpH # this is Newton's method              
        # println("step=", step, "  deltapH=", deltapH, "  pHTol=", pHTol, "  bodgeany=", bodgeany(abs(deltapH) > pHTol))
    end

    step < maxsteps || @warn "UpdatepHfromTATC maxsteps $maxsteps exceeded" pHstart pHfree deltapH TAx concs
    
    all(isfinite(pHfree)) || @warn "UpdatepHfromTATC failed" pHstart pHfree deltapH

    return (pHfree, step)
end

"""
    calculateOmega(C, CO3,Ca) -> (OmegaCA, OmegaAR)

Calculate carbonate saturation.

Returns omega, the solubility ratio, for calcite and aragonite.
This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
where Ksp is the solubility product (either KCa or KAr).

# Arguments:
- C    - constants from [`calc_constants!`](@ref)
- CO3  - carbonate ion concentration, mol/kg-sw
- Ca   - calcium concentration, mol/kg-sw
"""
function calculateOmega(C, CO3,Ca)
    
    OmegaCA = CO3*Ca/C[CNameIdx.KCa]
    OmegaAR = CO3*Ca/C[CNameIdx.KAr]

    return (OmegaCA, OmegaAR)
end

"""
    mappHscale(C, pHin, scalein, scaleout, TS, TF) -> pHout

Map pH scale at pressure, temperature, salinity defined by constants C

NB: Total, SW scale are not well defined unless using default contemporary values for `TS` (sulphate) and `TF` (fluorine)

# Arguments:
- `C`:  constants Vector from [`calc_constants!`](@ref)
- `pHin`: input pH on scale `scalein`
- `scalein`, `scaleout` :
  -  Val(1) = Total scale
  -  Val(2) = Seawater scale
  -  Val(3) = Free scale
  -  not implemented: Val(4) = NBS scale
- `TS:` total sulphate (mol kg-1)
- `TF:` total fluoride (mol kg-1)
"""
function mappHscale(C, pHin, scalein, scaleout, TS, TF)
         
    SWStoTOT  = (1 + TS/C[CNameIdx.KS])/(1 + TS/C[CNameIdx.KS] + TF/C[CNameIdx.KF])
    FREEtoTOT =  1 + TS/C[CNameIdx.KS]
    pHfreetopHtot = log(FREEtoTOT)/log(0.1) #  pHfree = pHtot - log(FREEtoTOT)/log(0.1);
    pHfreetopHsws = log(FREEtoTOT/SWStoTOT)/log(0.1)  #  pHsws  = pHtot - log(SWStoTOT) /log(0.1)

    # Map pHin to free scale
    if scalein == Val(1)        
        toFree = -pHfreetopHtot
    elseif scalein == Val(2)        
        toFree = -pHfreetopHsws
    elseif scalein == Val(3)
        toFree = 0.0
    else
        error("invalid pH scalein=",scalein)
    end

    # Map free scale to pHout
    if scaleout == Val(1)
        fromFree = pHfreetopHtot
    elseif scaleout == Val(2)
        fromFree = pHfreetopHsws
    elseif scaleout == Val(3)
        fromFree = 0;                
    else
        error("invalid pH scaleout ",scaleout)
    end
    
    pHout = pHin + toFree + fromFree
    return pHout
end
        
