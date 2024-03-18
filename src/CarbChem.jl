module CarbChem

import PALEOboxes as PB
using PALEOboxes.DocStrings

import SIMD

import ..PALEOcarbchem
import ..PALEOaqchem

import Infiltrator  # Julia debugger

"""
    ReactionCO2SYS

Carbonate chemistry using PALEOcarbchem.

Two different pH solver approaches are supported, set by `Parameter` `solve_pH`:
- `solve`: iterative solution for pH using an internal Newton solver, at given `TAlk` (eg `TAlk` is a StateExplicit Variable, using explicit ODE solver)
- `speciation`: calculate speciation given `pHfree`.
- `speciationTAlk`: calculate speciation given `pHfree`, also add alkalinity contributions to `TAlk_calc`
  (eg for use in combination with a `ReactionConstraintReservoir` that provides a total variable `TAlk` and `H` (as `pHfree`) as a primary species, for
  use with a DAE solver to solve the algebraic constraint for `pHfree` at a given `TAlk`).

NB: Two options have been removed from `Parameter` `solve_pH`:
- `constraint`: additional State Variable for pH and algebraic Constraint for TAlk (requires DAE solver).
  [replace with `speciationTAlk` and a `ReactionConstraintReservoir`]
- `implicit`: TAlk is a Total Variable, a function of a State Variable for pH (requires DAE solver)
  [replace with `speciationTAlk` and a `ReactionImplicitReservoir`]

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionCO2SYS{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParInt("WhichKs", 10,
            description="PALEOcarbchem WhichKs choice of eqb constant data"),
        PB.ParInt("WhoseKSO4", 1,
            description="PALEOcarbchem WhoseKSO4 choice of eqb constant data"),
        PB.ParStringVec("components", ["Ci", "S", "F", "B", "Omega"], 
            allowed_values=["Ci", "S", "F", "P", "B", "Si", "H2S", "NH3", "Omega"], # needs to match PALEOcarbchem
            description="PALEOcarbchem optional contributions to TAlk"),
        PB.ParStringVec("defaultconcs", ["TS", "TF", "TB", "Ca"], 
            allowed_values=["TS", "TF", "TB", "Ca"],
            description="PALEOcarbchem default concentrations (modern values calculated from salinity)"),

        # PB.ParStringVec("output_consts", String[], 
        #     description="PALEOcarbchem eqb constants etc to include as output variables"),

        PB.ParStringVec("outputs", ["pCO2", "xCO2dryinp"],
            allowed_values=["H", "OH", "TS", "HSO4", "TF", "HF", "TCi", "CO2", "HCO3", "CO3", "CAlk", "fCO2", "pCO2", "xCO2dryinp",
                "TB", "BAlk", "TP", "H3PO4", "H2PO4", "HPO4", "PO4", "PAlk", "TSi", "SiAlk", "TH2S", "H2S", "HSAlk", "TNH3", "NH4", "NH3Alk",
                "Ca", "OmegaCA", "OmegaAR"], # needs to match PALEOcarbchem.ResultNames
            description="PALEOcarbchem output concentrations etc to include as output variables"),

        PB.ParBool("output_pHtot", true, 
            description="true to output pHtot (requires TS_conc and TF_conc)"),

        PB.ParString("solve_pH", "solve", allowed_values=["solve", "speciation", "speciationTAlk"],
            description="'solve' to solve iteratively for pH, "*
                        "'speciation' to just provide speciation at supplied [pHfree], "*
                        "'speciationTAlk' to provide speciation and add to TAlk_calc"),

        PB.ParDouble("pHtol", 100*eps(), units="",
            description="with parameter solve_pH=solve only (internal Newton solver): pH tolerance for iterative solution"),

        PB.ParString("simd_width", "1", allowed_values=["1", "FP64P2", "FP64P4", "FP32P4", "FP32P8"],
            description="with parameter solve_pH=solve only (internal Newton solver): use SIMD (\"1\" - no SIMD, \"FP64P4\" - 4 x Float64, etc)"),
    )
end

function PB.register_methods!(rj::ReactionCO2SYS)

    io = IOBuffer()
    println(io, "register_methods! $(PB.fullname(rj)) ReactionCO2SYS")

    # variables needed by several methods
    var_volume = PB.VarDep("volume",           "m^3",      "Cell volume")
    base_vars = PB.VariableReaction[
        var_volume,
        PB.VarDep("rho_ref",          "kg m^-3",  "density conversion factor"),
        PB.VarDep("temp",             "Kelvin",   "temperature"),
        PB.VarDep("pressure",         "dbar",     "pressure"),
        PB.VarDep("sal",              "psu",      "salinity"),
    ]

    # add method to calculate default concentrations
    default_conc_names = rj.pars.defaultconcs.v.*"_conc"
    PB.setfrozen!(rj.pars.defaultconcs)
    println(io, "    modern default concentration Variables $(default_conc_names)")
    default_conc_vars =[
        PB.VarProp(conc_name, "mol m^-3", conc_name*" modern default scaled from salinity") for conc_name in default_conc_names
    ]
    DefaultConcOptions=(; WhichKs=Val(rj.pars.WhichKs[]), WhoseKSO4=Val(rj.pars.WhoseKSO4[]))
    PB.setfrozen!(rj.pars.WhichKs, rj.pars.WhoseKSO4)
    PB.add_method_do!(
        rj, 
        do_modern_default_concs,
        (   
            PB.VarList_namedtuple(base_vars),
            PB.VarList_namedtuple(default_conc_vars),
        ),
        p = (DefaultConcOptions, ),
    )

    # add methods to calculate carb chem
 
    if rj.pars.solve_pH[] == "solve"
        var_pHfree = PB.VarProp("pHfree", "", "-log10(hydrogen ion concentration), this is used to keep previous value as a starting value for internal Newton solver";
            attributes=(:initial_value=>8.0,))
        PB.add_method_setup_initialvalue_vars_default!(rj, [var_pHfree]; filterfn=v->true)
        push!(base_vars,
            var_pHfree,
            PB.VarDep("TAlk_conc",      "mol m^-3", "total TAlk concentration"),
        )

    elseif rj.pars.solve_pH[] == "speciation"
        var_pHfree = PB.VarDep("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        push!(base_vars, var_pHfree)

    elseif rj.pars.solve_pH[] == "speciationTAlk"
        var_pHfree = PB.VarDep("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        var_TAlk_calc = PB.VarContrib("TAlk_calc", "mol", "accumulated TAlk contributions")
        push!(base_vars, var_pHfree, var_TAlk_calc)

    else
        error(String(take!(io)), "\n", "unknown solve_pH=$(rj.pars.solve_pH[])")
    end
    PB.setfrozen!(rj.pars.solve_pH)

    if rj.pars.output_pHtot[]
        push!(base_vars, PB.VarProp("pHtot",           "",         "-log10([H] + [HS]), pH on total scale"))
    end
    PB.setfrozen!(rj.pars.output_pHtot)

    # get PALEOcarbchem components struct and required input concentrations
    (components, input_names) = PALEOcarbchem.get_components_inputs(rj.pars.components)
    PB.setfrozen!(rj.pars.components)
    input_conc_names = input_names.*"_conc"
    println(io, "    adding input concentration Variables $(input_conc_names)")
    input_vars = [
        # NB: set :field_data=>ScalarData to take just the concentration, not isotopic composition 
        PB.VarDep(vn, "mol m^-3", vn*"input concentration"; attributes=(:field_data=>PB.ScalarData,)) for vn in input_conc_names
    ]
    # field names without _conc suffix so we can create input names needed for PALEOcarbchem
    input_vars_nt = NamedTuple{Tuple(Symbol.(input_names))}(input_vars) 

    # get PALEOcarbchem output required, and create property variables
    output_vars=PB.VariableReaction[]
    idx_output_in_res = Int[]
    output_isconc = Bool[]
    for outputname in rj.pars.outputs
        outputdesc = PALEOcarbchem.getResultDescription(outputname)
        
        if isnothing(outputdesc)
            error(String(take!(io)), "\n", "invalid output Variable $(outputname) not present in PALEOcarbchem output")
        else
            push!(idx_output_in_res, getfield(PALEOcarbchem.RNameIdx, Symbol(outputname)))
            (outputdescname, outputunits, outputlongname) = outputdesc            
            if outputunits=="mol kg-1"
                push!(output_isconc, true)               
                convunits = "mol m-3"
            else
                push!(output_isconc, false)
                convunits = outputunits
            end
            outputlocalname = "%reaction%"*outputname
            println(io, "    adding additional output Variable $(outputlocalname) ($(convunits)) '$(outputlongname)'")
            push!(output_vars, PB.VarProp(outputlocalname,    convunits, outputlongname; attributes=PALEOaqchem.R_conc_attributes_advectfalse))
        end
    end
    PB.setfrozen!(rj.pars.outputs)
    println(io, "    outputs $(rj.pars.outputs.v) at indices $idx_output_in_res in res $(PALEOcarbchem.RNameIdx)")

    # var_output_consts = NamedTuple{}() # TODO
    CO2SYSOptions=(; WhichKs=Val(rj.pars.WhichKs[]), WhoseKSO4=Val(rj.pars.WhoseKSO4[]), pHScale=Val(3), Components=components)
    PB.add_method_do!(
        rj, 
        do_carbchem,
        (   
            PB.VarList_namedtuple(base_vars),
            PB.VarList_namedtuple_fields(input_vars_nt),  # field names without _conc suffix
            PB.VarList_namedtuple(output_vars),
        ),
        p = (
                Val(rj.pars.solve_pH[] == "solve"), # encode solution method as a Type
                CO2SYSOptions, 
                idx_output_in_res,
                output_isconc
            ),
        preparefn = prepare_do_carbchem,
    )

    @info String(take!(io))

    return nothing
end

function do_modern_default_concs(
    m::PB.ReactionMethod,
    pars,
    (vars, default_conc_vars),
    cellrange::PB.AbstractCellRange,
    deltat
)

    (CO2SYSOptions, ) = m.p

    @inbounds for i in cellrange.indices
        # salinity limits will be applied by calc_modern_default_concs
        modernC = PALEOcarbchem.calc_modern_default_concs(
            vars.sal[i],
            Options=CO2SYSOptions,
        )

        # conversion factor to  model units (mol m-3) to PALEOcarbchem units (mol kg-1)       
        rhofac = vars.rho_ref[i]    # kg m-3          
        #                                                  mol m-3                     = mol kg-1   * kg m-3
        hasfield(typeof(default_conc_vars), :TF_conc) && (default_conc_vars.TF_conc[i] = modernC.TF * rhofac)
        hasfield(typeof(default_conc_vars), :TS_conc) && (default_conc_vars.TS_conc[i] = modernC.TS * rhofac)
        hasfield(typeof(default_conc_vars), :TB_conc) && (default_conc_vars.TB_conc[i] = modernC.TB * rhofac)
        hasfield(typeof(default_conc_vars), :Ca_conc) && (default_conc_vars.Ca_conc[i] = modernC.Ca * rhofac)
    end

    return nothing
end


function prepare_do_carbchem(m::PB.ReactionMethod, (vars, input_concs, outputs))

    rj = m.reaction
  
    # work out BufType for temporary arrays
    PB.setfrozen!(rj.pars.simd_width)
    if rj.pars.simd_width[] == "1"
         # may need to be an AD type etc, so check all Variables to look for something that is not Float64
        BufType = Float64
        for v in merge(vars, input_concs, outputs)
            et = eltype(v)
            if et != Float64
                BufType = et
                break
            end
        end  
    elseif rj.pars.simd_width[] == "FP64P2"
        BufType = PB.SIMDutils.FP64P2
    elseif rj.pars.simd_width[] == "FP64P4"
        BufType = PB.SIMDutils.FP64P4
    elseif rj.pars.simd_width[] == "FP32P4"
        BufType = PB.SIMDutils.FP32P4
    elseif rj.pars.simd_width[] == "FP32P8"
        BufType = PB.SIMDutils.FP32P8
    else
        error("ReactionCO2SYS unknown simd_width ", rj.pars.simd_width[])
    end
   
    @info "prepare_do_carbchem: $(PB.fullname(rj)) simd_width=$(rj.pars.simd_width[]) using BufType $BufType"

    # create per-thread buffers for temporary results
    Constants_tbuf  = [Vector{BufType}(undef, length(PALEOcarbchem.CNameIdx)) for t in 1:Threads.nthreads()]    
    res_tbuf    = [Vector{BufType}(undef, length(PALEOcarbchem.RNameIdx)) for t in 1:Threads.nthreads()]
    input_concs_cell_tbuf = [PB.IteratorUtils.named_tuple_ref(keys(input_concs), BufType) for t in 1:Threads.nthreads()]

    buffers = (Constants_tbuf, res_tbuf, input_concs_cell_tbuf)

    return (vars, input_concs, outputs, buffers)
end

function do_carbchem(
    m::PB.ReactionMethod,
    pars,
    (vars, input_concs, outputs, buffers),
    cellrange::PB.AbstractCellRange,
    deltat
)

    (
        solve_pH,
        CO2SYSOptions, 
        idx_output_in_res,
        output_isconc
    ) = m.p

    Constants_tbuf, res_tbuf, input_concs_cell_tbuf = buffers

    Constants           = Constants_tbuf[Threads.threadid()] # get per thread buffer. threadid()=1 for a single-threaded app
    res                 = res_tbuf[Threads.threadid()]
    input_concs_cell    = input_concs_cell_tbuf[Threads.threadid()]

    BufType = eltype(Constants)

    for i in PB.SIMDutils.SIMDIter(cellrange.indices, BufType)
        # SIMD code with type conversion and scalar fallback:
        #   v = vec[i]  --> PB.SIMDutils.vgatherind(BufType, vec, i) # where v is a packed SIMD vector
        #   vec[i] = v  --> PB.SIMDutils.vscatterind!(v, vec, i)
        #   vec[i] += v --> PB.SIMDutils.vaddind!(v, vec, i)
       
        # conversion factors to/from model units (mol m-3) to PALEOcarbchem units (mol kg-1)       
        rhofac = PB.SIMDutils.vgatherind(BufType, vars.rho_ref, i)    # kg m-3        
        r_rhofac = 1.0/rhofac               # m3 kg-1
        
        # temperature and salinity limits will be applied by CO2SYS_Constants
        PALEOcarbchem.calc_constants!(
            Constants,
            PB.SIMDutils.vgatherind(BufType, vars.temp, i) - PB.Constants.k_CtoK, # degrees C
            PB.SIMDutils.vgatherind(BufType, vars.pressure, i),   # dbar
            PB.SIMDutils.vgatherind(BufType, vars.sal, i), # psu
            Options=CO2SYSOptions,
        )
      
        # Set 'input_concs_cell' with 'input_concs' for this cell
        # NB: a Julia bug (v1.7), map allocates so can't just create a NamedTuple
        # workaround using a 'mutable NamedTuple' buffer 'input_concs_cell' which has same fields as 'input_concs'
        for (ic_cell, ic) in PB.IteratorUtils.zipstrict(values(input_concs_cell), values(input_concs))
           ic_cell[] = r_rhofac*PB.SIMDutils.vgatherind(BufType, ic, i)
        end
        # TODO allocates ? (Julia 1.9.0-rc2)
        # PB.IteratorUtils.foreach_tuple(values(input_concs_cell), values(input_concs)) do ic_cell, ic
        #    ic_cell[] = r_rhofac*PB.SIMDutils.vgatherind(BufType, ic, i)
        # end

        if solve_pH == Val(true)  # this will be a compile-time switch as solve_pH Type is part of method signature
                    #   mol       / m^3         / kg m-3    
            pHstart = PB.SIMDutils.vgatherind(BufType, vars.pHfree, i)
            (pHfree, steps) = PALEOcarbchem.calculatepHfromTATC!(
                res,
                Constants,
                CO2SYSOptions,
                r_rhofac*PB.SIMDutils.vgatherind(BufType, vars.TAlk_conc, i),
                input_concs_cell,
                pHstart=pHstart,
                pHTol=pars.pHtol[]
            )
            
            # save pHfree solution as initial guess for next timestep
            PB.SIMDutils.vscatterind!(pHfree, vars.pHfree, i)
        else
            pHfree = PB.SIMDutils.vgatherind(BufType, vars.pHfree, i)
            (TA, _) = PALEOcarbchem.calculateTAfromTCpHfree!(
                res,
                Constants,
                CO2SYSOptions,
                input_concs_cell,
                pHfree,
                do_dTAdpH=Val(false),
            )

            if hasfield(typeof(vars), :TAlk_calc)
                vvolume = PB.SIMDutils.vgatherind(BufType, vars.volume, i)
                #    mol m-3                    mol kg-1 kg m-3                                           
                PB.SIMDutils.vaddind!(TA*rhofac*vvolume, vars.TAlk_calc, i)
            end
        end

       
        if pars.output_pHtot[]
            # TODO not clear whether this is defined for non-contemporary sulphate (and fluoride) ?
            PB.SIMDutils.vscatterind!(
                PALEOcarbchem.mappHscale(
                    Constants,
                    pHfree, 
                    Val(3),  # scalein - Free
                    Val(1),  # scaleout - Total
                    input_concs_cell.TS[], 
                    input_concs_cell.TF[],
                ),
                vars.pHtot, 
                i,
            )            
        end

        # iterate through outputs (a NamedTuple of Variables with fieldnames matching those of res) and set
        for (j, idx) in enumerate(idx_output_in_res)              
            if output_isconc[j]
                # mol m-3    mol kg-1
                PB.SIMDutils.vscatterind!(res[idx]*rhofac, outputs[j], i)
            else
                PB.SIMDutils.vscatterind!(res[idx], outputs[j], i)
            end            
        end

    end

    return nothing
end

end
