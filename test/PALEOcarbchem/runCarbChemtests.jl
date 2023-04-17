
using Test
using BenchmarkTools
using Logging

import PALEOaqchem
import PALEOaqchem.PALEOcarbchem

@testset "CO2SYS constants" begin

    # test check value against Matlab PALEO version
    # Matlab >> C = paleo_CO2SYSv.calc_constants(5.0, 1000.0, 30.0, 3, 10, 1)
    Carray = Vector{Float64}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = 5.0, 1000.0, 30.0
    PALEOcarbchem.calc_constants!(Carray, TempC, Pdbar, Sal, Options=Options)
    # convert to NamedTuple
    C = NamedTuple{PALEOcarbchem.ConstNames}(Carray)

    C_Matlab = (
       # TempC= 5,
       # Pbar= 100,
       #  Sal= 30,
    # pHScale= 3
    # WhichKs= 10
   # WhoseKSO4= 1
# RGasConstant= 83.1451,
 #      TempK= 278.15,
 #   logTempK= 5.6281605365213,
 #         RT= 23126.809565,
 #         TB= 0.000356314285714286,
          TF= 5.85650054471719e-05,
          TS= 0.0242018006853087,
          fH= 0.773872475,
      FugFac= 0.995889920302133,
       VPFac= 0.991538485813587,
          K0= 0.0536895248802576,
          K1= 8.48716384625353e-07,
          K2= 4.28286066341259e-10,
          KW= 7.74208946778645e-15,
          KB= 1.36576373371171e-09,
          KF= 0.00329431377216681,
          KS= 0.216365651814509,
         KP1= 0.023049183305505,
         KP2= 6.7309137521641e-07,
         KP3= 5.70778009179582e-10,
         KSi= 1.71579301949622e-10,
        KH2S= 1.31794700245033e-07,
  #      KNH4= 1.0691821328762e-10,
        KNH3= 1.0691821328762e-10,
#   FREEtoTOT= 1.1118560200399,
#    SWStoTOT= 0.984262504631478,
#pHfreetopHtot= -0.0460485518587636,
#pHfreetopHsws= -0.0529376108369824,
        #  Ca= 0.00881534545787119,
         KCa= 4.38467838091862e-07,
         KAr= 6.94158466390701e-07,
    )

    # sigdigits = 14
    sigdigits = 13
    for k in keys(C_Matlab)
        @test check_round(getfield(C,k),  getfield(C_Matlab,k), sigdigits, message=k)
    end
end

@testset "CO2SYS modern default concs" begin

    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1))

    Sal = 30.0
    modernC = PALEOcarbchem.calc_modern_default_concs(30.0, Options=Options)

    modernC_Matlab = (
          TB= 0.000356314285714286,
          TF= 5.85650054471719e-05,
          TS= 0.0242018006853087,
          Ca= 0.00881534545787119,
    )

    sigdigits = 13
    for k in keys(modernC_Matlab)
        @test check_round(getfield(modernC, k),  getfield(modernC_Matlab, k), sigdigits, message=k)
    end

end

@testset "CO2SYS Calculate" begin

    # test check value against Matlab PALEO version
    # Matlab C = paleo_CO2SYSv.calc_constants(25.0, 0.0, 35.0)
    C = Vector{Float64}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = 25.0, 0.0, 35.0
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar, Sal, Options=Options)

    modernC = PALEOcarbchem.calc_modern_default_concs(Sal, Options=Options)

    # Matlab
    # [TA, r] = paleo_CO2SYSv.CalculateTAfromTCpHfree(C, 8.0, 2000e-6, 1e-3, 1e-6, C.TS, C.TB, C.TF, 1e-6, 1e-6, true)

    r_Matlab =    (  pHfree= 8,
                    TCi= 0.002,
                    TSi= 0.001,
                TH2S= 1e-06,
                TNH3= 1e-06,
                    TP= 1e-06,
                    TS= 0.0282354341328601,
                    TB= 0.0004157,
                    TF= 6.83258396883673e-05,
                    H= 9.99999999999998e-09,
                HCO3= 0.00182913941036246,
                    CO3= 0.000154374423636126,
                    CO2= 1.6486166001411e-05,
                CAlk= 0.00213788825763472,
                BAlk= 6.84606236173229e-05,
                    OH= 4.69746810439883e-06,
                H3PO4= 5.48753087051033e-15,
                H2PO4= 1.03905993343373e-08,
                HPO4= 8.79005163290632e-07,
                    PO4= 1.106042318875e-07,
                PAlk= 1.1002136215781e-06,
                SiAlk= 3.10201867146702e-05,
                HSAlk= 9.60144986057519e-07,
                NH3Alk= 4.05866587336644e-08,
                HSO4= 2.81503971113345e-09,
                    HF= 2.88841800321286e-10,
                dTAdpH = 0.000583697672714085, # Matlab dTAdPH
                PengCorrection= 0,
                    TA= 0.00224415437745597,
                fCO2= 0.000580664787040614,
                pCO2= 0.000582522767470251,
                xCO2dryinp= 0.000600944914072491)

    CI = PALEOcarbchem.CNameIdx
    Concs = (TCi=[2000e-6], TS=modernC.TS, TF=modernC.TF, TSi=[1e-3], TP=[1e-6],  TB=modernC.TB, TH2S=[1e-6], TNH3=[1e-6], Ca=modernC.Ca)
  
    resarray = Vector{Float64}(undef, length(PALEOcarbchem.RNameIdx))
    (TA, dTAdpH) = PALEOcarbchem.calculateTAfromTCpHfree!(resarray, C, Options, Concs, 8.0, do_dTAdpH=Val(true))
    r = NamedTuple{PALEOcarbchem.ResultNames}(resarray)
    
    # sigdigits = 14
    # sigdigits = 13
    sigdigits = 11
    for k in keys(r_Matlab)
        @test check_round(getfield(r,k),  getfield(r_Matlab,k), sigdigits, message=k)
    end

    # check Newton solve for pHfree
    # Matlab [pHfree, r] = paleo_CO2SYSv.UpdatepHfromTATC(C, 2000e-6, 2000e-6, 8.0, eps, 1e-3, 1e-6, C.TS, C.TB, C.TF, 1e-6, 1e-6)

    pHfree_Matlab = 7.39104603058702
    (pHfree, steps) = PALEOcarbchem.calculatepHfromTATC!(resarray, C, Options, 2000e-6, Concs, pHTol=eps())
    @test check_round(pHfree, pHfree_Matlab, 14, message="pHfree")

end