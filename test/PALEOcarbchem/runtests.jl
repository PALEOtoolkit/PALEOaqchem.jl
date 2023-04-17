using Test
using Documenter
import PALEOaqchem
import PALEOaqchem.PALEOcarbchem

function check_round(a, b, sigdigits; message="")
    lhs = round(a, sigdigits=sigdigits)
    rhs = round(b, sigdigits=sigdigits)
    result = lhs == rhs
    if !result
        println("check_round (", message, ") ",a, " != ",b, " sigdigits=", sigdigits, " round(a) ", lhs, " != round(b) ", rhs)
    end
    return result
end

@testset "PALEOcarbchem all" begin

include("runCarbChemtests.jl")

include("runCarbChembench.jl")

@test bench_all() === nothing

# doctest(PALEOcarbchem; manual=false)  

end # @testset
