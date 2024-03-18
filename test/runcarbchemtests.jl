
using Test
using Logging
using Sundials
using OrdinaryDiffEq # for QNDF
using BenchmarkTools

import PALEOboxes as PB
import PALEOocean
import PALEOaqchem
import PALEOmodel

@testset "carbchem" begin
skipped_testsets = [
    # "Ocean 3 box carb chem",
    # "Ocean 3 box carb chem Hfree DAE", 
    # "Ocean 3 box carb chem Hfree implicit",   
]

configfile = joinpath(@__DIR__, "configcarbchem.yaml")


!("Ocean 3 box carb chem" in skipped_testsets) && @testset "Ocean 3 box carb chem" begin

    model = PB.create_model_from_config(configfile, "test_carbchem")

    initial_state, modeldata = PALEOmodel.initialize!(model)
   
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # Check CO2SYS 
    
    PB.do_deriv(modeldata.dispatchlists_all)

    PB.TestUtils.bench_model(
        model, modeldata; 
        bench_whole=false, domainname="ocean", reactionname="carbchem", methodname="do_carbchem"
    )
   
    println("ocean  model created variables after do_stateandeqb:")
    ocean_domain = PB.get_domain(model, "ocean")
    vars = PB.get_variables(ocean_domain, hostdep=false)
    for var in vars
        println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
    end

    var_pHtot = PB.get_variable(ocean_domain, "pHtot")    
    @test round.(PB.get_data(var_pHtot, modeldata), sigdigits=12) == round.([7.3540128696428155, 7.6020215057634175, 7.525099896712386], sigdigits=12)

end

!("Ocean 3 box carb chem Hfree DAE" in skipped_testsets) && @testset "Ocean 3 box carb chem Hfree DAE" begin

    model = PB.create_model_from_config(
        configfile, 
        "test_carbchem_constraint"; 
    )

    # add pH, TAlk to state variables and apply constraint on TAlk      

    initial_state, modeldata = PALEOmodel.initialize!(model)
    
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # Test ODE solvers with mass matrix, and DAE solvers, with/without Jacobian
    
    for solver in ["integrateDAE", "integrate", "integrateDAEForwardDiff", "integrateDAEForwardDiffsparse"]    
   
        println("===========================================================")
        println("integrate to approx steady state with solver ", solver)
        println("===========================================================")

        if solver=="integrateDAE"  # DAE with no Jacobian
            PALEOmodel.ODE.integrateDAE(run, initial_state, modeldata, (0, 1e5), alg=IDA())  # nevals 157     nsteps 12
        elseif solver=="integrate"  # ODE with mass matrix, no Jacobian
            # PALEOmodel.ODE.integrate(run, initial_state, modeldata, (0, 1e5), alg=IDA())  # fails, mass matrix not supported
            PALEOmodel.ODE.integrate(run, initial_state, modeldata, (0, 1e5), alg=QNDF(autodiff=false)) # nevals 7069     nsteps 24

            # PALEOmodel.ODE.integrate(run, initial_state, modeldata, (0, 1e5), alg=Rodas5(autodiff=false)) # inefficient ? nevals=1035, nsteps=44
        elseif solver=="integrateDAEForwardDiff"  # DAE with ForwardDiff dense Jacobian
            PALEOmodel.ODE.integrateDAE(run, initial_state, modeldata, (0, 1e5), alg=IDA(),
                                       jac_ad=:ForwardDiff) # nevals 24 njacs 11 nsteps 12
        elseif solver=="integrateDAEForwardDiffsparse"  # DAE with ForwardDiff sparse Jacobian
            PALEOmodel.ODE.integrateDAEForwardDiff(run, initial_state, modeldata, (0, 1e5)) # nevals 24, njacs 11, nsteps 12
        else
            error("unrecognized solver='$(solver)'")
        end

        println("ocean  model created variables after integrate to steady-state with $(solver):")
        ocean_domain = PB.get_domain(model, "ocean")
        vars = PB.get_variables(ocean_domain, hostdep=false)
        for var in vars
            println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
        end

        var_pHtot = PB.get_variable(ocean_domain, "pHtot")
        sigdigits = 5
        @test round.(PB.get_data(var_pHtot, modeldata), sigdigits=sigdigits) == round.([7.3540128696428155, 7.6020215057634175, 7.525099896712386], sigdigits=sigdigits)
    end
end

!("Ocean 3 box carb chem Hfree implicit" in skipped_testsets) && @testset "Ocean 3 box carb chem Hfree implicit" begin

    model = PB.create_model_from_config(
        configfile, 
        "test_carbchem_implicit";
    )

    initial_state, modeldata = PALEOmodel.initialize!(model)
    
    run = PALEOmodel.Run(model=model, output = PALEOmodel.OutputWriters.OutputMemory())

    # Test DAE solvers, with Jacobian
    
    for solver in ["integrateDAEForwardDiff", "integrateDAEForwardDiffsparse"]    
   
        println("===========================================================")
        println("integrate to approx steady state with solver ", solver)
        println("===========================================================")

        if solver=="integrateDAEForwardDiff"  # DAE with ForwardDiff dense Jacobian
            PALEOmodel.ODE.integrateDAE(run, initial_state, modeldata, (0, 1e5), alg=IDA(),
                                       jac_ad=:ForwardDiff) # nevals 24 njacs 11 nsteps 12
        elseif solver=="integrateDAEForwardDiffsparse"  # DAE with ForwardDiff sparse Jacobian
            PALEOmodel.ODE.integrateDAEForwardDiff(run, initial_state, modeldata, (0, 1e5)) # nevals 24, njacs 11, nsteps 12
        else
            error("unrecognized solver='$(solver)'")
        end

        println("ocean  model created variables after integrate to steady-state with $(solver):")
        ocean_domain = PB.get_domain(model, "ocean")
        vars = PB.get_variables(ocean_domain, hostdep=false)
        for var in vars
            println("\t", PB.fullname(var), " = ", PB.get_data(var, modeldata))
        end

        # implicit variable TAlk is initialised from pHfree, not a specified TAlk_conc
        var_TAlk_conc = PB.get_variable(ocean_domain, "TAlk_conc")
        var_pHtot = PB.get_variable(ocean_domain, "pHtot")
        sigdigits = 5
        @test round.(PB.get_data(var_TAlk_conc, modeldata); sigdigits) == round.([2.1939990876584, 2.099727277183955, 2.1263041720076443]; sigdigits)
        @test round.(PB.get_data(var_pHtot, modeldata); sigdigits) == [7.9042, 7.9547, 7.9602]
    end
end

end

