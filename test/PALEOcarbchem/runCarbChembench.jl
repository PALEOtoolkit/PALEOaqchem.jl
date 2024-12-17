
using Test
using BenchmarkTools
using Logging

import PALEOaqchem
import PALEOaqchem.PALEOcarbchem

import SIMD

function bench_all()
    println("bench_Constants():")
    bench_Constants()
    println()

    println("bench_Constants_simd():")
    bench_Constants_simd()
    println()

    println("bench_TAfromTCpHfree()")
    bench_TAfromTCpHfree()
    println()

    println("bench_TAfromTCpHfree_simd()")
    bench_TAfromTCpHfree_simd()
    println()

    println("bench_pHfromTATC()")
    bench_pHfromTATC()
    println()

    println("bench_pHfromTATC_simd()")
    bench_pHfromTATC()
    println()

    return nothing
end

function bench_Constants()    
    C = Vector{Float64}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = 25.0, 0.0, 35.0
    @btime PALEOcarbchem.calc_constants!($C, $TempC, $Pdbar, $Sal, Options=$Options)
end

function bench_TAfromTCpHfree()
    C = Vector{Float64}(undef, length(PALEOcarbchem.CNameIdx))

    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = 25.0, 0.0, 35.0
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar, Sal, Options=Options)
    modernC = PALEOcarbchem.calc_modern_default_concs(Sal, Options=Options)

    Concs = (TCi=2000e-6, TS=modernC.TS, TF=modernC.TF, TSi=1e-3, TP=1e-6,  TB=modernC.TB, TH2S=1e-6, TNH3=1e-6, Ca=modernC.Ca)

    resarray = Vector{Float64}(undef, length(PALEOcarbchem.RNameIdx))
    println("do_dTAdpH=false")
    @btime (TA, dTAdpH) = PALEOcarbchem.calculateTAfromTCpHfree!($resarray, $C, $Options, $Concs, 8.0, do_dTAdpH=Val(false))
    println("do_dTAdpH=true")
    @btime (TA, dTAdpH) = PALEOcarbchem.calculateTAfromTCpHfree!($resarray, $C, $Options, $Concs, 8.0, do_dTAdpH=Val(true))

end

function bench_pHfromTATC()
    C = Vector{Float64}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = 25.0, 0.0, 35.0
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar, Sal, Options=Options)
    modernC = PALEOcarbchem.calc_modern_default_concs(Sal, Options=Options)

    Concs = (TCi=2000e-6, TS=modernC.TS, TF=modernC.TF, TSi=1e-3, TP=1e-6,  TB=modernC.TB, TH2S=1e-6, TNH3=1e-6, Ca=modernC.Ca)

    resarray = Vector{Float64}(undef, length(PALEOcarbchem.RNameIdx))
   
    @btime (pHfree, steps) = PALEOcarbchem.calculatepHfromTATC!($resarray, $C, $Options, 2000e-6, $Concs, pHTol=2*eps(Float64))

end

###########################################
# Benchmark SIMD packed data types
#############################################

# Short names for packed vector Types
const FP64P4 = SIMD.Vec{4, Float64}
const FP64P2 = SIMD.Vec{2, Float64}
const FP32P8 = SIMD.Vec{8, Float32}
const FP32P4 = SIMD.Vec{4, Float32}

# const ES = FP32P8 # SIMD width to test
const ES = FP32P4 # SIMD width to test


function check_Constants_simd()
    C = Vector{ES}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    
    TempC, Pdbar, Sal = ES(25.0), ES(0.0), ES(35.0)
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar , Sal, Options=Options)

    return C
end

function bench_Constants_simd()    
    C = Vector{ES}(undef, length(PALEOcarbchem.CNameIdx))
    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = ES(25.0), ES(0.0), ES(35.0)
    @btime PALEOcarbchem.calc_constants!($C,$TempC, $Pdbar , $Sal, Options=$Options)
end

function bench_TAfromTCpHfree_simd()

    C = Vector{ES}(undef, length(PALEOcarbchem.CNameIdx))

    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = ES(25.0), ES(0.0), ES(35.0)
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar, Sal, Options=Options)
    modernC = PALEOcarbchem.calc_modern_default_concs(Sal, Options=Options)

    Concs = (TCi=ES(2000e-6), TS=modernC.TS, TF=modernC.TF, TSi=ES(1e-3), TP=ES(1e-6),  TB=modernC.TB, TH2S=ES(1e-6), TNH3=ES(1e-6), Ca=modernC.Ca)

    pHstart = ES(8.0)
    resarray = Vector{ES}(undef, length(PALEOcarbchem.RNameIdx))
    println("do_dTAdpH=false")
    @btime (TA, dTAdpH) = PALEOcarbchem.calculateTAfromTCpHfree!($resarray, $C, $Options, $Concs, $pHstart, do_dTAdpH=Val(false))
    println("do_dTAdpH=true")
    @btime (TA, dTAdpH) = PALEOcarbchem.calculateTAfromTCpHfree!($resarray, $C, $Options, $Concs, $pHstart, do_dTAdpH=Val(true))

end

function bench_pHfromTATC_simd()
    C = Vector{ES}(undef, length(PALEOcarbchem.CNameIdx))

    Options=(; WhichKs=Val(10), WhoseKSO4=Val(1), pHScale=Val(3), Components=PALEOcarbchem.ComponentsAll)
    TempC, Pdbar, Sal = ES(25.0), ES(0.0), ES(35.0)
    PALEOcarbchem.calc_constants!(C, TempC, Pdbar, Sal, Options=Options)
    modernC = PALEOcarbchem.calc_modern_default_concs(Sal, Options=Options)

    Concs = (TCi=ES(2000e-6), TS=modernC.TS, TF=modernC.TF, TSi=ES(1e-3), TP=ES(1e-6),  TB=modernC.TB, TH2S=ES(1e-6), TNH3=ES(1e-6), Ca=modernC.Ca)

    resarray = Vector{ES}(undef, length(PALEOcarbchem.RNameIdx))
   
    TA = ES(2000e-6)
    pHstart = ES(8.0)
    @btime (pHfree, steps) = PALEOcarbchem.calculatepHfromTATC!($resarray, $C, $Options, $TA, $Concs, pHstart=$pHstart, pHTol=2*eps(Float32))

end

##################################################
# optional fast log etc functions
######################################################

const USE_SLEEF = true
@static if USE_SLEEF
    import SLEEF_jll

    @info "Defining Vectorized SIMD functions log, exp, log10 functions from Sleef library $(SLEEF_jll.libsleef)"

    # raw tuple data types held in Vec.data
    const FP64P4_d = NTuple{4,Core.VecElement{Float64}}
    const FP64P2_d = NTuple{2,Core.VecElement{Float64}}
    const FP32P8_d = NTuple{8,Core.VecElement{Float32}}
    const FP32P4_d = NTuple{4,Core.VecElement{Float32}}

    # exp Sleef Vectorized double/single precision base-e exponential functions functions with 1.0 ULP error bound
    sleefexp(v::FP64P4_d)   = ccall((:Sleef_expd4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.exp(v::FP64P4)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP64P2_d)   = ccall((:Sleef_expd2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.exp(v::FP64P2)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP32P8_d)   = ccall((:Sleef_expf8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.exp(v::FP32P8)     = SIMD.Vec(sleefexp(v.data))
    sleefexp(v::FP32P4_d)   = ccall((:Sleef_expf4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.exp(v::FP32P4)     = SIMD.Vec(sleefexp(v.data))

    # log Sleef Vectorized double/single precision natural logarithmic functions with 1.0 ULP error bound
    sleeflog(v::FP64P4_d)   = ccall((:Sleef_logd4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.log(v::FP64P4)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP64P2_d)   = ccall((:Sleef_logd2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.log(v::FP64P2)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP32P8_d)   = ccall((:Sleef_logf8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.log(v::FP32P8)     = SIMD.Vec(sleeflog(v.data))
    sleeflog(v::FP32P4_d)   = ccall((:Sleef_logf4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.log(v::FP32P4)     = SIMD.Vec(sleeflog(v.data))
    # this should work, but isn't specific enough to override the fallbacks in SIMD.jl
    # Base.log(v::SIMD.Vec{N,T}) where {N,T} = SIMD.Vec(sleeflog(v.data))

    # Vectorized double/single precision base-10 logarithmic functions with 1.0 ULP error bound
    sleeflog10(v::FP64P4_d)   = ccall((:Sleef_log10d4_u10, SLEEF_jll.libsleef), FP64P4_d, (FP64P4_d,), v)
    Base.log10(v::FP64P4)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP64P2_d)   = ccall((:Sleef_log10d2_u10, SLEEF_jll.libsleef), FP64P2_d, (FP64P2_d,), v)
    Base.log10(v::FP64P2)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP32P8_d)   = ccall((:Sleef_log10f8_u10, SLEEF_jll.libsleef), FP32P8_d, (FP32P8_d,), v)
    Base.log10(v::FP32P8)     = SIMD.Vec(sleeflog10(v.data))
    sleeflog10(v::FP32P4_d)   = ccall((:Sleef_log10f4_u10, SLEEF_jll.libsleef), FP32P4_d, (FP32P4_d,), v)
    Base.log10(v::FP32P4)     = SIMD.Vec(sleeflog10(v.data))
else
    @info "Using defaults (probably slow scalar fallbacks) for SIMD log, exp, log10 functions"
end

