module CarbChem

import PALEOboxes as PB
using PALEOboxes.DocStrings

import SIMD

import ..PALEOcarbchem

# import Infiltrator  # Julia debugger

"""
    ReactionCO2SYS

Carbonate chemistry using PALEOcarbchem.

Three different pH solver approaches are supported, set by `Parameter` `solve_pH`:
- `solve`: iterative solution for pH, TAlk is a StateExplicit Variable (eg using explicit ODE solver)
- `constraint`: additional State Variable for pH and algebraic Constraint for TAlk (requires DAE solver).
- `implicit`: TAlk is a Total Variable, a function of a State Variable for pH (requires DAE solver)

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
            description="PALEOcarbchem output concentrations etc to include as output variables"),

        PB.ParBool("output_pHtot", true, 
            description="true to output pHtot (requires TS_conc and TF_conc)"),

        PB.ParString("solve_pH", "solve", allowed_values=["solve", "constraint", "implicit"],
            description="'solve' to solve iteratively for pH, "*
                        "'constraint' to include [pHfree] as a state variable and TAlk as a constraint, "*
                        "'implicit' to include [pHfree] as a state variable and TAlk an implicit variable"),

        PB.ParDouble("pHtol", 100*eps(), units="",
            description="pH tolerance for iterative solution"),

        PB.ParString("simd_width", "1", allowed_values=["1", "FP64P2", "FP64P4", "FP32P4", "FP32P8"],
            description="use SIMD (\"1\" - no SIMD, \"FP64P4\" - 4 x Float64, etc)"),
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
    setup_vars = PB.VariableReaction[var_volume, ]

    if rj.pars.solve_pH[] == "solve"
        var_pHfree = PB.VarProp("pHfree",         "", "-log10(hydrogen ion concentration)")
        push!(setup_vars, PB.VarDep(var_pHfree))
        push!(base_vars,
            var_pHfree,
            PB.VarDep("TAlk_conc",      "mol m^-3", "total TAlk concentration"),
        )

    elseif rj.pars.solve_pH[] == "constraint"
        var_pHfree = PB.VarState("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        var_pHfree_constraint = PB.VarConstraint("pHfree_constraint", "mol", "pH constraint variable for DAE solver")
        var_TAlk_calc = PB.VarTarget("%reaction%TAlk_calc", "mol", "accumulated TAlk contributions")
        push!(setup_vars, var_pHfree, var_pHfree_constraint)
        push!(base_vars, var_pHfree, PB.VarContrib(var_TAlk_calc))
        constraint_vars = [
            PB.VarDep("TAlk",      "mol m^-3", "total TAlk"),            
            var_TAlk_calc,
            var_pHfree_constraint,
        ]
        PB.add_method_do!(
            rj, 
            do_carbchem_pH_constraint,
            (PB.VarList_namedtuple(constraint_vars), ), 
        )

    elseif rj.pars.solve_pH[] == "implicit"
        var_pHfree = PB.VarState("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        var_TAlk = PB.VarTotal(      "TAlk",        "mol",      "TAlk total variable")
        push!(setup_vars, var_pHfree, var_TAlk)
        push!(base_vars,
            var_pHfree,
            PB.VarProp("%reaction%TAlk_calc", "mol", "ReactionCO2SYS TAlk contributions",
                attributes=(:initialize_to_zero=>true, )),
            var_TAlk,
            # PB.VarContrib("TAlk",          "mol", "TAlk Total Variable"),
        )
    else
        error(String(take!(io)), "\n", "unknown solve_pH=$(rj.pars.solve_pH[])")
    end
    PB.setfrozen!(rj.pars.solve_pH)

    if rj.pars.output_pHtot[]
        push!(base_vars, PB.VarProp("pHtot",           "",         "-log10([H] + [HS]), pH on total scale"))
    end
    PB.setfrozen!(rj.pars.output_pHtot)

    PB.add_method_setup!(
        rj, 
        setup_carbchem,
        (PB.VarList_namedtuple(setup_vars), ),
    )

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
            push!(output_vars, PB.VarProp(outputlocalname,    convunits, outputlongname))
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

    PB.add_method_initialize_zero_vars_default!(rj) # constraint and implicit variables

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

function setup_carbchem(
    m::PB.ReactionMethod,
    pars,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    attribute_name in (:initial_value, :norm_value) || return

    initial_pHfree = Dict(:initial_value=>8.0, :norm_value=>1.0)[attribute_name]

    if pars.solve_pH[] == "solve"
        # pHfree is a VarProp used to keep previous value as a starting value
        if attribute_name == :initial_value
            vars.pHfree[cellrange.indices] .= initial_pHfree
            @info "PALEOcarbchem initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        end

    elseif pars.solve_pH[] == "constraint"
        # pHfree is a State Variable
        vars.pHfree[cellrange.indices] .= initial_pHfree
        @info "PALEOcarbchem initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        if attribute_name == :norm_value           
            vars.pHfree_constraint[cellrange.indices] .= 1.0
            @info "PALEOcarbchem initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree_constraint", 30)) = 1.0"
        end

    elseif pars.solve_pH[] == "implicit"
        # pHfree is a State Variable
        vars.pHfree[cellrange.indices] .= initial_pHfree
        @info "PALEOcarbchem initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        initial_TAlk_conc = PB.get_domvar_attribute(PB.get_variable(m, "TAlk"), attribute_name)
        vars.TAlk[cellrange.indices] .= initial_TAlk_conc*vars.volume[cellrange.indices]  
        @info "PALEOcarbchem initialize_state:$(rpad(attribute_name,20)) $(rpad("TAlk", 30)) = $initial_TAlk_conc * volume"
    else
        error("unknown solve_pH $(pars.solve_pH[])")
    end   

    return nothing
end

function prepare_do_carbchem(m::PB.ReactionMethod, (vars, input_concs, outputs))

    rj = m.reaction
  
    # work out BufType for temporary arrays
    PB.setfrozen!(rj.pars.simd_width)
    if rj.pars.simd_width[] == "1"
         # may be an AD type etc, so get type from a Variable that is likely to be a state Variable
        BufType = nothing
        vnames_ad = (:TAlk, :TAlk_calc, :TAlk_conc)
        for vname in vnames_ad
            if vname in propertynames(vars)
                BufType = eltype(getfield(vars, vname))
                break;
            end
        end
        !isnothing(BufType) || error("ReactionCO2SYS simd_width=\"1\" couldnt determine BufType, tried $vnames_ad")       
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

            vvolume = PB.SIMDutils.vgatherind(BufType, vars.volume, i)
            #    mol m-3                    mol kg-1 kg m-3                                           
            PB.SIMDutils.vaddind!(TA*rhofac*vvolume, vars.TAlk_calc, i)

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

    if pars.solve_pH[] == "implicit"
        @inbounds for i in cellrange.indices
            vars.TAlk[i] += vars.TAlk_calc[i]
        end
    end

    return nothing
end

function do_carbchem_pH_constraint(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)

    @inbounds for i in cellrange.indices
        vars.pHfree_constraint[i] =   vars.TAlk_calc[i] - vars.TAlk[i]
    end

    return nothing
end


end
