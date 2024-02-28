module MolecularDiffusion

import PALEOboxes as PB

"""
    create_solute_diffusivity_func(species_or_constant::AbstractString) -> f_diffcoeff

Molecular diffusivity for `species_or_constant` in water

Based on data compiled by [Boudreau1997](@cite) Ch. 4

If `species_or_constant` is a species name, returns a function `f_diffcoeff`, which calculates:
    
    f_diffcoeff(TdegK, Pbar, S) -> Dmol (cm^2/s)

If `species_or_constant` can be parsed into a `Float64`, 
returns a function which provides this constant value (units cm^2/s).

The effect of pressure and salinity is approximated by assuming viscosity 
(calculated by [`MolecularDiffusion.sw_dynamic_viscosity`](@ref)) is the primary consideration, 
[Boudreau1997](@cite) eqn. 4.107 and Table 4.10.
"""
function create_solute_diffusivity_func(species_or_constant::AbstractString)

  
    # helper function to convert T units
    TdegC(TdegK) = TdegK - PB.Constants.k_CtoK

    P1atm = 1.01325 # bar

    # Stokes-Einstein correction for P and S
    # converts Pfrom,0.0 -> Pto, S 
    function corfacPS(TdegC, PFrombar, PTobar, S)
        # viscosity of pure water and seawater at specified T, P, S
        mu0 = sw_dynamic_viscosity(TdegC,  PFrombar, 0.0)
        muS = sw_dynamic_viscosity(TdegC,  PTobar, S)
        return mu0/muS
    end

    function DmolO2(TdegK, Pbar, S)
        # [Boudreau1997](@cite) Fig 4.5 linear correlation
        mu0 = sw_dynamic_viscosity(TdegC(TdegK),  P1atm, 0.0)
        Dmol0 = 1e-5*(0.2604 + 0.006383*(TdegK/mu0))
        Dmol = Dmol0*corfacPS(TdegC(TdegK), P1atm, Pbar, S)
        return Dmol
    end

    function DmolCO2(TdegK, Pbar, S)
        # [Boudreau1997](@cite) Fig 4.5 linear correlation
        mu0 = sw_dynamic_viscosity(TdegC(TdegK),  P1atm, 0.0)
        Dmol0 = 1e-5*(0.194 + 0.005089*(TdegK/mu0))
        Dmol = Dmol0*corfacPS(TdegC(TdegK), P1atm, Pbar, S)
        return Dmol
    end

    function DmolArrhenius(TdegK, Pbar, S, A, Ea)
        # [Boudreau1997](@cite) Eq 4.60 Arrhenius formulation
        Dmol0 = 1e-5*A*exp(-Ea*1000.0/(PB.Constants.k_Rgas*TdegK))
        # correct assuming Dmol0 is for 1 atm, zero S
        Dmol = Dmol0*corfacPS(TdegC(TdegK), P1atm, Pbar, S)
        return Dmol
    end 

    function DmolWilkeChang(TdegK, Pbar, S, Vb)
        # [Boudreau1997](@cite) Eq 4.57 
        mu0 = sw_dynamic_viscosity(TdegC(TdegK),  P1atm, 0.0)
        Dmol0 = 4.72e-9*TdegK/(mu0*1e-2*Vb^0.6)
        Dmol = Dmol0*corfacPS(TdegC(TdegK), P1atm, Pbar, S)
        return Dmol
    end 

    function DmollinT(TdegK, Pbar, S, m0, m1)
        Dmol0 = 1e-6*(m0+TdegC(TdegK)*m1)
        # correct assuming Dmol0 is for 1 atm, zero S
        Dmol = Dmol0*corfacPS(TdegC(TdegK), P1atm, Pbar, S)
        return Dmol
    end 

    # return appropriate function for species_or_constant

    const_diffusivity = tryparse(Float64, species_or_constant)

    if !isnothing(const_diffusivity)
        Dmol_const = let const_diffusivity=const_diffusivity 
            (TdegK, Pbar, S) -> const_diffusivity
        end
        return Dmol_const
    elseif species_or_constant == "O2"
        return DmolO2
    elseif species_or_constant == "CO2"
        return DmolCO2
    elseif Symbol(species_or_constant) in keys(diffArrhenius)
        DmolArrheniusSp = let (A, Ea) = getfield(diffArrhenius, Symbol(species_or_constant))
            (TdegK, Pbar, S) -> DmolArrhenius(TdegK, Pbar, S, A, Ea)
        end 
        return DmolArrheniusSp
    elseif Symbol(species_or_constant) in keys(diffWilkeChang)
        DmolWilkeChangSp = let Vb = getfield(diffWilkeChang, Symbol(species_or_constant))
            (TdegK, Pbar, S) -> DmolWilkeChang(TdegK, Pbar, S, Vb)
        end
        return DmolWilkeChangSp
    elseif Symbol(species_or_constant) in keys(difflinT)
        DmolLinSp = let (m0, m1) = getfield(difflinT, Symbol(species_or_constant))
            (TdegK, Pbar, S) -> DmollinT(TdegK, Pbar, S, m0, m1)
        end
        return DmolLinSp
    else 
        error("unknown species_or_constant $species_or_constant")       
    end

    error("unreachable reached")
end


"""
    sw_dynamic_viscosity(TdegC, Pbar, S) -> mu

Seawater dynamic viscosity ([Boudreau1997](@cite) p94) 
'claimed to accurate to within 7%  over T (0 - 30C), P (1 - 1000 bar), S (0 36)'

Returns `mu  (10^-2 g/cm/s)  (centipoise)`

# Arguments:
- `T`: deg C (0 - 30C) temperature 
- `P`: bar (1 1000) pressure
- `S`: salinity (0 - 36)

See check values for salinity effect Table 4.10 p.125
"""
function sw_dynamic_viscosity(TdegC, Pbar, S)
   
    mu = (1.7910
        -6.144e-2*TdegC + 1.4510e-3*TdegC^2 - 1.6826e-5*TdegC^3
        -1.5290e-4*Pbar + 8.3885e-8*Pbar^2
        + 2.4727e-3*S
        + TdegC*(6.0574e-6*Pbar-2.6760e-9*Pbar^2)
        + S*(4.8429e-5*TdegC-4.7172e-6*TdegC^2 + 7.5986e-8*TdegC^3)
    )

    return mu
end


# Arrhenius diffusivity
# [Boudreau1997](@cite) Eq 4.50 and Table 4.4 - Arrhenius formulation
# H2, CH4, He, Ne, Kr, Xe, Rn: Jahne et al (1987)
# DMS: Saltzman et al. (1993)
# Ar: Ohsumi and Horibe (1984)
const diffArrhenius = (
    #       A           Ea
    #     1e-5 cm^2/s   kJ / mol
    He  =   (818.0,     11.70),
    Ne  =   (1608.0,    14.84),
    Ar  =   (7238.0,    19.81),
    Kr  =   (6393.0,    20.20),
    Xe  =   (9007.0,    21.61),
    Rn  =   (15877.0,   23.26),
    H2  =   (3338.0,    16.06), 
    CH4 =   (3047.0,    18.36),
    DMS =   (2000.0,    18.10),
)
    
#  from Wilke and Chang (1955) as modified by Hayduk and Laudie (1974)
#  Eqn 4.57 and Table 4.3 in [Boudreau1997](@cite)
const diffWilkeChang = (
    #       molar volume cm^3/mol
    NO  =   23.6,
    N2  =   34.7,
    H2S =   35.2,
    N2O =   36.0,
    NH3 =   24.5,
    CO  =   34.5,
    SO2 =   43.8,
)

# Linear regression against temperature for anions and cations
# [Boudreau1997](@cite) Table 4.7, 4.8 p115-6
const difflinT = (
    #         m0        m1

    # Table 4.7 p115 +1 cations 
    H      = (54.4,   1.555),
    Li     = (4.43,   0.241),
    Na     = (6.06,   0.297),
    K      = (6.06,   0.297),
    Cs     = (10.3,   0.416),
    Ag     = (7.82,   0.359),
    NH4    = (9.50,   0.413),

    # Table 4.7 p115 +2 cations
    Ba     = (4.06,   0.176),
    Be     = (2.57,   0.140),
    Ca     = (3.60,   0.179),
    Cd     = (3.31,   0.152),
    Co     = (3.31,   0.152),
    Cu     = (3.39,   0.158),
    Fe     = (3.31,   0.150),
    Hg     = (3.63,   0.208),
    Mg     = (3.43,   0.144),
    Mn     = (3.18,   0.155),
    Ni     = (3.36,   0.130),
    Sr     = (3.69,   0.169),
    Pb     = (4.46,   0.198),
    Ra     = (3.91,   0.199),
    Zn     = (3.31,   0.151),

    # Table 4.7 p115 +3 cations
    Al     = (2.79,   0.172),
    Ce     = (2.95,   0.131),
    La     = (2.78,   0.136),
    Pu     = (2.71,   0.120),       

   # Table 4.8 p116 -1 anions
    OH     = (25.9,   1.094),
    Br     = (10.0,   0.441),
    Cl     = (9.60,   0.438),
    F      = (6.29,   0.343),
    HCO3   = (5.06,   0.275),
    H2PO4  = (4.02,   0.223),
    HS     = (10.4,   0.273),
    HSO3   = (6.35,   0.280),
    HSO4   = (5.99,   0.307),
    I      = (9.81,   0.432),
    IO3    = (4.66,   0.252),
    NO2    = (10.3,   0.331),
    NO3    = (9.50,   0.388),

    # Table 4.8 p116 -2 anions
    CO3    = (4.33,   0.199),
    HPO4   = (3.26,   0.177),
    SO3    = (4.53,   0.249),
    SO4    = (4.88,   0.232),
        
    # Table 4.8 p116 -3 anions
    PO4    = (2.62,   0.143),

)

end # module