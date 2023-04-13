module CarbChem

import PALEOboxes as PB
using PALEOboxes.DocStrings

import SIMD

import PALEOco2sys

# import Infiltrator  # Julia debugger

"""
    ReactionCO2SYS

Carbonate chemistry using paleo_CO2SYS.

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
            description="PALEOco2sys WhichKs choice of eqb constant data"),
        PB.ParInt("WhoseKSO4", 1,
            description="PALEOco2sys WhoseKSO4 choice of eqb constant data"),
        PB.ParStringVec("components", ["Ci", "S", "F", "B", "Omega"], 
            allowed_values=["Ci", "S", "F", "P", "B", "Si", "H2S", "NH3", "Omega"], # needs to match paleo_CO2SYS
            description="PALEOco2sys optional contributions to TAlk"),
        PB.ParStringVec("defaultconcs", ["TS", "TF", "TB", "Ca"], 
            allowed_values=["TS", "TF", "TB", "Ca"],
            description="PALEOco2sys default concentrations (modern values calculated from salinity)"),

        PB.ParStringVec("output_consts", String[], 
            description="PALEOco2sys eqb constants etc to include as output variables"),

        PB.ParStringVec("outputs", ["pCO2", "xCO2dryinp"], 
            description="PALEOco2sys output concentrations etc to include as output variables"),

        PB.ParBool("output_pHtot", true, 
            description="true to output pHtot (requires TS and TF)"),

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

    var_volume = PB.VarDep("volume",           "m^3",      "Cell volume")
    basevariables = PB.VariableReaction[
        var_volume,
        PB.VarDep("rho_ref",          "kg m^-3",  "density conversion factor"),
        PB.VarDep("temp",             "Kelvin",   "temperature"),
        PB.VarDep("pressure",         "dbar",     "pressure"),
        PB.VarDep("sal",              "psu",      "salinity"),
    ]

    setup_vars = PB.VariableReaction[var_volume]

    if rj.pars.solve_pH[] == "solve"
        var_pHfree = PB.VarProp("pHfree",         "", "-log10(hydrogen ion concentration)")

        push!(setup_vars, PB.VarDep(var_pHfree))

        push!(basevariables,
            var_pHfree,
            PB.VarDep("TAlk_conc",      "mol m^-3", "total TAlk concentration"),
        )
    elseif rj.pars.solve_pH[] == "constraint"
        var_pHfree = PB.VarState("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        var_pHfree_constraint = PB.VarConstraint("pHfree_constraint", "mol", "pH constraint variable for DAE solver")
        var_TAlk_calc = PB.VarTarget("%reaction%TAlk_calc", "mol", "accumulated TAlk contributions")

        push!(setup_vars, var_pHfree, var_pHfree_constraint)
        
        push!(basevariables, var_pHfree, PB.VarContrib(var_TAlk_calc))

        constraint_vars = [
            PB.VarDep("TAlk",      "mol m^-3", "total TAlk"),            
            var_TAlk_calc,
            var_pHfree_constraint,
        ]

    elseif rj.pars.solve_pH[] == "implicit"
        var_pHfree = PB.VarState("pHfree",         "", "-log10(hydrogen ion concentration) state variable")
        var_TAlk = PB.VarTotal(      "TAlk",        "mol",      "TAlk total variable")

        push!(setup_vars, var_pHfree, var_TAlk)

        push!(basevariables,
            var_pHfree,
            PB.VarProp("%reaction%TAlk_calc", "mol", "CarbChem TAlk contributions",
                attributes=(:initialize_to_zero=>true, )),
            var_TAlk,
            # PB.VarContrib("TAlk",          "mol", "TAlk Total Variable"),
        )
    else
        error("unknown solve_pH=$(rj.pars.solve_pH[])")
    end
    PB.setfrozen!(rj.pars.solve_pH)

    PB.add_method_setup!(
        rj, 
        setup_carbchem,
        (PB.VarList_namedtuple(setup_vars), ),
    )


    if rj.pars.output_pHtot[]
        push!(basevariables, PB.VarProp("pHtot",           "",         "-log10([H] + [HS]), pH on total scale"))
    end
    PB.setfrozen!(rj.pars.output_pHtot)

    # Vector of requested default concentrations followed by non-default input concentrations 
    definput_conc_vars = PB.VariableReaction[]
    definput_conc_names = String[]

    # add any requested default concentrations
    idx_default_conc_in_C = Tuple(getfield(PALEOco2sys.CNameIdx, Symbol(f)) for f in rj.pars.defaultconcs)
    @info "CarbChem default_concs $(rj.pars.defaultconcs.v) at indices $idx_default_conc_in_C in C $(PALEOco2sys.CNameIdx)"
    for conc in rj.pars.defaultconcs
        @info "CarbChem adding additional modern default Variable $(conc)_conc"
        push!(definput_conc_names, conc)
        push!(definput_conc_vars, PB.VarProp("%reaction%"*conc*"_conc", "mol m^-3", conc*" modern default scaled from salinity"))
    end
    PB.setfrozen!(rj.pars.defaultconcs)
   
    # get additional PALEOco2sys components struct, and create required additional input concentrations
    (components, additional_inputs) = PALEOco2sys.get_components_inputs(rj.pars.components)
    PB.setfrozen!(rj.pars.components)
    CO2SYSOptions=(; WhichKs=Val(rj.pars.WhichKs[]), WhoseKSO4=Val(rj.pars.WhoseKSO4[]), pHScale=Val(3), Components=components)
    
    default_conc_vars_dict = Dict(v.localname => v for v in definput_conc_vars)
    for conc in additional_inputs
        varname = conc*"_conc"
        if haskey(default_conc_vars_dict, varname)
            # will get this input from default            
        else
            @info "CarbChem adding additional dependency Variable $(varname)"
            # NB: set :field_data=>ScalarData to take just the concentration, not isotopic composition 
            push!(definput_conc_names, conc)          
            push!(definput_conc_vars, PB.VarDep(varname,        "mol m^-3", conc*" concentration"; attributes=(:field_data=>PB.ScalarData,)))
        end
    end
    definput_conc_nt = NamedTuple{Tuple(Symbol.(definput_conc_names))}(definput_conc_vars)
  
    # get PALEOco2sys output required, and create property variables
    output_vars=PB.VariableReaction[]
    output_isconc = Bool[]
    for outputname in rj.pars.outputs
        outputdesc = PALEOco2sys.getResultDescription(outputname)
        
        if isnothing(outputdesc)
            error("CarbChem invalid output Variable $(outputname) not present in paleo_CO2SYS output")
        else
            (outputdescname, outputunits, outputlongname) = outputdesc            
            if outputunits=="mol kg-1"
                push!(output_isconc, true)               
                convunits = "mol m-3"
            else
                push!(output_isconc, false)
                convunits = outputunits
            end
            outputlocalname = "%reaction%"*outputname
            @info "CarbChem adding additional output Variable $(outputlocalname) ($(convunits)) '$(outputlongname)'"
            push!(output_vars, PB.VarProp(outputlocalname,    convunits, outputlongname))
        end
    end
    PB.setfrozen!(rj.pars.outputs)

    idx_output_in_res = Tuple(getfield(PALEOco2sys.RNameIdx, Symbol(f) ) for f in rj.pars.outputs)
    @info "CarbChem outputs $(rj.pars.outputs.v) at indices $idx_output_in_res in res $(PALEOco2sys.RNameIdx)"

    # var_output_consts = NamedTuple{}() # TODO

    PB.add_method_do!(
        rj, 
        # rj.pars.simd_width[] == "1" ? do_carbchem_scalar : do_carbchem_simd,
        do_carbchem_simd,
        (   
            PB.VarList_namedtuple(basevariables),
            PB.VarList_namedtuple_fields(definput_conc_nt), # fieldnames without _conc suffix
            PB.VarList_namedtuple(output_vars),
        ),
        p = (
                Val(rj.pars.solve_pH[] == "solve"), # encode solution method as a Type
                CO2SYSOptions, 
                idx_default_conc_in_C,
                idx_output_in_res,
                output_isconc
            ),
        preparefn = prepare_do_carbchem,
    )

    if rj.pars.solve_pH[] == "constraint"
        PB.add_method_do!(
            rj, 
            do_carbchem_pH_constraint,
            (PB.VarList_namedtuple(constraint_vars), ), 
        )
    end

    PB.add_method_initialize_zero_vars_default!(rj) # constraint and implicit

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
            @info "PALEOco2sys initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        end

    elseif pars.solve_pH[] == "constraint"
        # pHfree is a State Variable
        vars.pHfree[cellrange.indices] .= initial_pHfree
        @info "PALEOco2sys initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        if attribute_name == :norm_value           
            vars.pHfree_constraint[cellrange.indices] .= 1.0
            @info "PALEOco2sys initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree_constraint", 30)) = 1.0"
        end

    elseif pars.solve_pH[] == "implicit"
        # pHfree is a State Variable
        vars.pHfree[cellrange.indices] .= initial_pHfree
        @info "PALEOco2sys initialize_state:$(rpad(attribute_name,20)) $(rpad("pHfree", 30)) = $initial_pHfree"
        initial_TAlk_conc = PB.get_domvar_attribute(PB.get_variable(m, "TAlk"), attribute_name)
        vars.TAlk[cellrange.indices] .= initial_TAlk_conc*vars.volume[cellrange.indices]  
        @info "PALEOco2sys initialize_state:$(rpad(attribute_name,20)) $(rpad("TAlk", 30)) = $initial_TAlk_conc * volume"
    else
        error("unknown solve_pH $(pars.solve_pH[])")
    end   

    return nothing
end

function prepare_do_carbchem(m::PB.ReactionMethod, (vars, definput_concs, outputs))

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
        !isnothing(BufType) || error("CarbChem simd_width=\"1\" couldnt determine BufType, tried $vnames_ad")       
    elseif rj.pars.simd_width[] == "FP64P2"
        BufType = PB.SIMDutils.FP64P2
    elseif rj.pars.simd_width[] == "FP64P4"
        BufType = PB.SIMDutils.FP64P4
    elseif rj.pars.simd_width[] == "FP32P4"
        BufType = PB.SIMDutils.FP32P4
    elseif rj.pars.simd_width[] == "FP32P8"
        BufType = PB.SIMDutils.FP32P8
    else
        error("CarbChem unknown simd_width ", rj.pars.simd_width[])
    end
   
    @info "prepare_do_carbchem: $(PB.fullname(rj)) simd_width=$(rj.pars.simd_width[]) using BufType $BufType"

    # create per-thread buffers for temporary results
    Carray_tbuf  = [Vector{BufType}(undef, length(PALEOco2sys.CNameIdx)) for t in 1:Threads.nthreads()]    
    res_tbuf    = [Vector{BufType}(undef, length(PALEOco2sys.RNameIdx)) for t in 1:Threads.nthreads()]
    input_concs_cell_tbuf = [PB.IteratorUtils.named_tuple_ref(keys(definput_concs), BufType) for t in 1:Threads.nthreads()]

    buffers = (Carray_tbuf, res_tbuf, input_concs_cell_tbuf)

    return (vars, definput_concs, outputs, buffers)
end


function do_carbchem_simd(
    m::PB.ReactionMethod,
    pars,
    (vars, definput_concs, outputs, buffers),
    cellrange::PB.AbstractCellRange,
    deltat
)

    (
        solve_pH,
        CO2SYSOptions, 
        idx_default_conc_in_C,   
        idx_output_in_res,
        output_isconc
    ) = m.p

    Carray_tbuf, res_tbuf, input_concs_cell_tbuf = buffers

    Carray          = Carray_tbuf[Threads.threadid()] # get per thread buffer. threadid()=1 for a single-threaded app
    res             = res_tbuf[Threads.threadid()]
    input_concs_cell= input_concs_cell_tbuf[Threads.threadid()]

    BufType = eltype(Carray)

    for indvec in PB.SIMDutils.SIMDIter(cellrange.indices, BufType)
        # SIMD code with type conversion and scalar fallback:
        #   v = vec[indvec]  --> PB.SIMDutils.vgatherind(BufType, vec, indvec)
        #   vec[indvec] = v  --> PB.SIMDutils.vscatterind!(v, vec, indvec)
        #   vec[indvec] += v --> PB.SIMDutils.vaddind!(v, vec, indvec) 
       
        # conversion factors to/from model units (mol m-3) to PALEOco2sys units (mol kg-1)       
        rhofac = PB.SIMDutils.vgatherind(BufType, vars.rho_ref, indvec)    # kg m-3        
        r_rhofac = 1.0/rhofac               # m3 kg-1

        vvolume = PB.SIMDutils.vgatherind(BufType, vars.volume, indvec)
        
        # temperature and salinity limits will be applied by CO2SYS_Constants
        PALEOco2sys.calc_constants!(
            Carray,
            PB.SIMDutils.vgatherind(BufType, vars.temp, indvec) - PB.Constants.k_CtoK, # degrees C
            PB.SIMDutils.vgatherind(BufType, vars.pressure, indvec),   # dbar
            PB.SIMDutils.vgatherind(BufType, vars.sal, indvec), # psu
            Options=CO2SYSOptions,
        )
      
        # set defaults from C     
        for (j, idx) in enumerate(idx_default_conc_in_C)
            # mol m-3   mol kg-1        kg m-3
            PB.SIMDutils.vscatterind!(Carray[idx]*rhofac, definput_concs[j], indvec)
        end     
     
       
        # NamedTuple with concs for this cell
        # NB: a Julia bug (v1.7), map allocates so can't easily create a NamedTuple
        # workaround using a 'mutable NamedTuple' buffer
        for j in eachindex(values(input_concs_cell))
            input_concs_cell[j][] = r_rhofac*PB.SIMDutils.vgatherind(BufType, definput_concs[j], indvec)
        end

        if solve_pH == Val(true)  # this will be a compile-time switch as solve_pH Type is part of method signature
                    #   mol       / m^3         / kg m-3    
            pHstart = PB.SIMDutils.vgatherind(BufType, vars.pHfree, indvec)
            (pHfree, steps) = PALEOco2sys.calculatepHfromTATC!(
                res,
                Carray,
                CO2SYSOptions,
                r_rhofac*PB.SIMDutils.vgatherind(BufType, vars.TAlk_conc, indvec),
                input_concs_cell,
                pHstart=pHstart,
                pHTol=pars.pHtol[]
            )
            
            # save pHfree solution as initial guess for next timestep
            PB.SIMDutils.vscatterind!(pHfree, vars.pHfree, indvec)
        else
            pHfree = PB.SIMDutils.vgatherind(BufType, vars.pHfree, indvec)
            (TA, dTAdpH) = PALEOco2sys.calculateTAfromTCpHfree!(
                res,
                Carray,
                CO2SYSOptions,
                input_concs_cell,
                pHfree,
                do_dTAdpH=Val(false),
            )

            #    mol m-3                    mol kg-1 kg m-3                                           
            PB.SIMDutils.vaddind!(TA*rhofac*vvolume, vars.TAlk_calc, indvec)

        end

       
        if pars.output_pHtot[]
            # TODO not clear whether this is defined for non-contemporary sulphate (and fluoride) ?
            PB.SIMDutils.vscatterind!(
                PALEOco2sys.mappHscale(
                    Carray,
                    pHfree, 
                    Val(3),  # scalein - Free
                    Val(1),  # scaleout - Total
                    input_concs_cell.TS[], 
                    input_concs_cell.TF[],
                ),
                vars.pHtot, 
                indvec,
            )            
        end

        # iterate through outputs (a NamedTuple of Variables with fieldnames matching those of res) and set
        for (j, idx) in enumerate(idx_output_in_res)              
            if output_isconc[j]
                # mol m-3    mol kg-1
                PB.SIMDutils.vscatterind!(res[idx]*rhofac, outputs[j], indvec)
            else
                PB.SIMDutils.vscatterind!(res[idx], outputs[j], indvec)
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
